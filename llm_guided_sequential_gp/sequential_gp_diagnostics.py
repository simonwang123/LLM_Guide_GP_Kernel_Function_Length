"""Diagnose one revealed sequential point using its sealed M_t prediction.

Run without arguments for sample 21. For later steps, provide the revealed
training prefix (through t), the M_t model directory, and the next sample
index. This script only prepares evidence; it neither refits the GP nor calls
the LLM. The truth CSV is read only up to the requested row. The diagnostic
context contains only aggregate statistics and current-point diagnostics.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import pickle
from pathlib import Path

import numpy as np

try:
    from . import sequential_gp_data as data_io
except ImportError:
    import sequential_gp_data as data_io


HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "outputs" / "semiparametric_guarded_20_20261109"
DEFAULT_MODEL_DIR = HERE / "outputs" / "models" / "samples_0020"
VARIABLES = (
    "upstream_slope", "downstream_slope", "crest_width", "modulus_zone_height",
    "secondary_top_width", "secondary_height", "secondary_bottom_width",
    "secondary_right_slope_distance",
)
RESPONSES = (
    "crest_max_settlement_m", "dam_max_settlement_m",
    "dam_max_horizontal_displacement_m",
)
UNIT_FIELDS = tuple(f"{name}_unit" for name in VARIABLES)
BOUNDS = {
    "upstream_slope": (1.3, 1.7),
    "downstream_slope": (1.3, 1.7),
    "crest_width": (6.0, 12.0),
    "modulus_zone_height": (5.0, 20.0),
    "secondary_top_width": (10.0, 30.0),
    "secondary_height": (20.0, 40.0),
    "secondary_bottom_width": (25.0, 70.0),
    "secondary_right_slope_distance": (2.0, 15.0),
}
PRIOR_CATEGORIES = frozenset({"short", "base", "long"})


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_hash(value: object) -> str:
    raw = json.dumps(value, ensure_ascii=False, sort_keys=True,
                     separators=(",", ":"), allow_nan=False).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def unique_pairs(pairs: list[tuple[str, object]]) -> dict:
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON number: {value}")


def read_json(path: Path) -> dict:
    try:
        result = json.loads(path.read_text(encoding="utf-8-sig"),
                            object_pairs_hook=unique_pairs,
                            parse_constant=reject_constant)
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(result, dict):
        raise ValueError(f"JSON root must be an object: {path}")
    return result


def csv_rows(path: Path, required: tuple[str, ...], *, target: int | None = None) -> list[dict[str, str]]:
    """When target is given, stop at that row; never read later FE truth."""
    if not path.is_file():
        raise ValueError(f"missing CSV file: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames
        if fields is None or len(fields) != len(set(fields)) or set(required) - set(fields):
            raise ValueError(f"missing/duplicate CSV columns in {path}; required: {required}")
        rows = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"malformed CSV row in {path}")
            if target is None:
                rows.append(row)
            elif int(row["sample_index"]) == target:
                return [row]
    if not rows:
        raise ValueError(f"no requested rows in {path}" + (f" for sample {target}" if target else ""))
    return rows


def numeric(value: str | int | float, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} is not numeric") from error
    if not math.isfinite(result):
        raise ValueError(f"{label} is not finite")
    return result


def unit_design(row: dict[str, str], label: str) -> np.ndarray:
    x = np.asarray([numeric(row[name], f"{label}/{name}") for name in UNIT_FIELDS])
    if np.any(x < -1e-10) or np.any(x > 1 + 1e-10):
        raise ValueError(f"{label} normalized design is outside [0, 1]")
    physical = set(VARIABLES) & set(row)
    if physical:
        if physical != set(VARIABLES):
            raise ValueError(f"{label} has incomplete physical columns")
        bounds = np.asarray(list(BOUNDS.values()))
        values = np.asarray([numeric(row[name], f"{label}/{name}") for name in VARIABLES])
        expected = (values - bounds[:, 0]) / (bounds[:, 1] - bounds[:, 0])
        if not np.allclose(x, expected, rtol=0, atol=1e-9):
            raise ValueError(f"{label} physical and normalized designs disagree")
    return x


def identity(row: dict[str, str], label: str) -> tuple[int, str]:
    try:
        index = int(row["sample_index"])
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} has invalid sample_index") from error
    design_id = row["design_id"].strip()
    if index < 1 or not design_id:
        raise ValueError(f"{label} has invalid sample_index/design_id")
    return index, design_id


def lengthscale_rule_mode(rules: dict) -> str:
    """Keep legacy lower bounds distinct from the new soft prior categories."""
    version = rules.get("schema_version")
    if type(version) is not int:
        raise ValueError("unsupported length-scale rules schema")
    if version == 1:
        key = "lengthscale_tiers"
        allowed = {"Tier-1", "Tier-2", "Tier-3"}
        mode = "legacy_hard_lower_bound"
    elif version == 2:
        key = "lengthscale_prior_categories"
        allowed = PRIOR_CATEGORIES
        mode = "soft_log_normal"
    else:
        raise ValueError("unsupported length-scale rules schema")
    mapping = rules.get(key)
    if not isinstance(mapping, dict) or set(mapping) != set(RESPONSES):
        raise ValueError(f"{key} must contain all three responses")
    for response in RESPONSES:
        values = mapping[response]
        if (not isinstance(values, dict) or set(values) != set(VARIABLES)
                or any(type(category) is not str or category not in allowed
                       for category in values.values())):
            raise ValueError(f"{key}/{response} must contain eight valid variable values")
    return mode


def model_prediction(fitted: tuple, x: np.ndarray) -> tuple[float, float, float]:
    model, mean, scale, _ = fitted
    prediction, std = model.predict(x.reshape(1, -1), return_std=True)
    nugget = float(model.kernel_.k2.noise_level)
    observed_std = float(scale * std[0])
    latent_std = float(scale * math.sqrt(max(float(std[0] ** 2) - nugget, 0.0)))
    return float(mean + scale * prediction[0]), observed_std, latent_std


def model_mean(fitted: tuple, x: np.ndarray) -> float:
    model, mean, scale, _ = fitted
    result = float(mean + scale * model.predict(x.reshape(1, -1))[0])
    if not math.isfinite(result):
        raise ValueError("GP returned a non-finite prediction mean")
    return result


def nearest_context(x: np.ndarray, train_x: np.ndarray,
                    train_rows: list[dict[str, str]]) -> tuple[dict, list[dict]]:
    distances = np.linalg.norm(train_x - x, axis=1)
    order = np.argsort(distances, kind="stable")[:3]
    train_distances = np.linalg.norm(train_x[:, None, :] - train_x[None, :, :], axis=2)
    np.fill_diagonal(train_distances, np.inf)
    sorted_train = np.sort(train_distances, axis=1)
    q90_first = float(np.quantile(sorted_train[:, 0], 0.9))
    q90_third = float(np.quantile(sorted_train[:, 2], 0.9))
    first_far = float(distances[order[0]]) > q90_first
    third_far = float(distances[order[2]]) > q90_third
    coverage = "sparse" if first_far and third_far else "mixed" if first_far or third_far else "typical"
    neighbors = []
    for rank, i in enumerate(order, 1):
        row = train_rows[int(i)]
        neighbors.append({
            "rank": rank,
            "sample_index": int(row["sample_index"]),
            "design_id": row["design_id"].strip(),
            "distance_euclidean_unit_space": float(distances[i]),
            "design_unit": train_x[i].tolist(),
            "delta_unit_new_minus_neighbor": (x - train_x[i]).tolist(),
        })
    return {
        "method": "Euclidean distance in fixed [0,1]^8 space",
        "k": 3,
        "nearest_distance": float(distances[order[0]]),
        "very_close_screening_distance": 0.08,
        "has_very_close_historical_point": bool(distances[order[0]] < 0.08),
        "distance_q90_reference": "leave-one-out first/third-neighbor distances among current training designs",
        "first_neighbor_q90": q90_first,
        "third_neighbor_q90": q90_third,
        "first_neighbor_far": first_far,
        "third_neighbor_far": third_far,
        "coverage_class": coverage,
        "interpretation": (
            "The 0.08 screen means very close, while the q90 flags measure relative "
            "sparsity for this training set. Typical relative spacing does not imply local "
            "comparability or accuracy, and sparse does not mean outside the design bounds."
        ),
    }, neighbors


def taylor_attribution(fitted: tuple, x: np.ndarray, neighbors: list[dict]) -> tuple[dict, dict]:
    """Rank model-prediction changes; never attribute observed residual causally."""
    if x.shape != (len(VARIABLES),) or len(neighbors) != 3:
        raise ValueError("Taylor attribution requires one 8D point and three neighbors")
    step = 1e-4
    zero_tol = 1e-12
    gradient = np.empty(len(VARIABLES))
    for j in range(len(VARIABLES)):
        low, high = x.copy(), x.copy()
        low[j] = max(0.0, x[j] - step)
        high[j] = min(1.0, x[j] + step)
        gradient[j] = (model_mean(fitted, high) - model_mean(fitted, low)) / (high[j] - low[j])
    current_mean = model_mean(fitted, x)
    contributions = []
    comparisons = []
    detailed_neighbors = []
    for neighbor in neighbors:
        design = np.asarray(neighbor["design_unit"], dtype=float)
        delta = x - design
        signed = gradient * delta
        actual = current_mean - model_mean(fitted, design)
        first_order = float(np.sum(signed))
        remainder = actual - first_order
        absolute_total = float(np.sum(np.abs(signed)))
        poor = abs(remainder) > 0.25 * max(absolute_total, zero_tol)
        contributions.append(signed)
        comparisons.append({
            "neighbor_rank": neighbor["rank"],
            "predicted_increment_m": actual,
            "first_order_increment_m": first_order,
            "remainder_m": remainder,
            "sum_absolute_contributions_m": absolute_total,
            "poor_approximation": bool(poor),
        })
        detailed_neighbors.append({
            **neighbor,
            "signed_contributions_m": dict(zip(VARIABLES, signed.tolist(), strict=True)),
            "prediction_comparison": comparisons[-1],
        })
    by_neighbor = np.vstack(contributions)
    mean_absolute = np.mean(np.abs(by_neighbor), axis=0)
    total = float(np.sum(mean_absolute))
    no_dominant = total <= zero_tol
    top3 = []
    if not no_dominant:
        for j in sorted(range(len(VARIABLES)), key=lambda index: (-mean_absolute[index], index))[:3]:
            signed = by_neighbor[:, j]
            top3.append({
                "variable": VARIABLES[j],
                "mean_absolute_contribution_m": float(mean_absolute[j]),
                "share_of_all_absolute_contributions": float(mean_absolute[j] / total),
                "mean_signed_contribution_m": float(np.mean(signed)),
                "direction_counts": {
                    "positive": int(np.sum(signed > zero_tol)),
                    "negative": int(np.sum(signed < -zero_tol)),
                    "zero": int(np.sum(np.abs(signed) <= zero_tol)),
                },
            })
    public = {
        "method": "First-order Taylor terms at the new point; three Euclidean nearest neighbors, equal weight",
        "neighbor_distances_unit_space": [item["distance_euclidean_unit_space"] for item in neighbors],
        "mean_absolute_contribution_total_m": total,
        "no_dominant_contribution": bool(no_dominant),
        "top3": top3,
        "approximation_checks": comparisons,
    }
    detailed = {
        "new_design_unit": dict(zip(VARIABLES, x.tolist(), strict=True)),
        "gradient_m_per_unit": dict(zip(VARIABLES, gradient.tolist(), strict=True)),
        "mean_absolute_contributions_m": dict(zip(VARIABLES, mean_absolute.tolist(), strict=True)),
        "neighbors": detailed_neighbors,
    }
    return public, detailed


def run(train: Path, designs: Path, truth_path: Path, model_dir: Path,
        output: Path, sample_index: int | None, ledger_path: Path | None,
        prior_context_path: Path | None) -> dict:
    train_rows = csv_rows(train, ("sample_index", "design_id", *UNIT_FIELDS, *RESPONSES))
    if len(train_rows) < 20:
        raise ValueError("training CSV must contain at least 20 revealed samples")
    expected = len(train_rows) + 1
    target = expected if sample_index is None else sample_index
    if target != expected:
        raise ValueError(f"target sample must be {expected} for this training prefix")
    identities = [identity(row, "training") for row in train_rows]
    if [index for index, _ in identities] != list(range(1, expected)):
        raise ValueError("training CSV must be a contiguous revealed prefix")
    if len({design_id for _, design_id in identities}) != len(identities):
        raise ValueError("duplicate training design_id")
    train_x = np.vstack([unit_design(row, "training") for row in train_rows])
    if np.unique(train_x, axis=0).shape[0] != len(train_x):
        raise ValueError("duplicate training design vectors")

    design = csv_rows(designs, ("sample_index", "design_id", *UNIT_FIELDS), target=target)[0]
    prediction_path = model_dir / "blind_prediction.csv"
    prediction_rows = csv_rows(prediction_path, ("sample_index", "design_id"))
    if len(prediction_rows) != 1:
        raise ValueError("blind_prediction.csv must contain exactly one sealed prediction")
    prediction = prediction_rows[0]
    revealed = csv_rows(truth_path, ("sample_index", "design_id", *RESPONSES), target=target)[0]
    ids = [identity(row, label) for label, row in (
        ("design", design), ("prediction", prediction), ("truth", revealed)
    )]
    if len(set(ids)) != 1 or ids[0][0] != target or ids[0][1] in {item[1] for item in identities}:
        raise ValueError("design, prediction, and revealed truth do not identify the same new point")
    x = unit_design(design, "new design")
    if np.any(np.all(np.isclose(train_x, x, rtol=0, atol=1e-12), axis=1)):
        raise ValueError("new design duplicates a training design")

    report_path = model_dir / "fit_report.json"
    rules_path = model_dir / "rules_used.json"
    model_path = model_dir / "gp_models.pkl"
    report, rules = read_json(report_path), read_json(rules_path)
    if report.get("training_sample_count") != len(train_rows) or report.get("train_sha256") != file_hash(train):
        raise ValueError("model report does not match the supplied revealed training prefix")
    if report.get("blind_prediction_sample_index") != target:
        raise ValueError("model report does not match the sealed prediction sample")
    if report.get("blind_prediction_source") != str(designs.resolve()):
        raise ValueError("design CSV is not the source used for the sealed blind prediction")
    if report.get("rules_sha256") != canonical_hash(rules):
        raise ValueError("model report and rules_used.json disagree")
    if rules.get("responses") != list(RESPONSES) or rules.get("design_variables") != list(VARIABLES):
        raise ValueError("model rules have unexpected response or variable names")
    prior_mode = lengthscale_rule_mode(rules)
    if set(report.get("responses", {})) != set(RESPONSES):
        raise ValueError("model report does not contain all three responses")
    if report.get("fit_settings", {}).get("kernel") != "Constant * ARD Matern(nu=2.5) + WhiteKernel":
        raise ValueError("this diagnostic only supports the current pure-data ARD GP model")

    # This pickle must come from the locally trusted sequential_gp_model.py run.
    try:
        with model_path.open("rb") as stream:
            saved_model = pickle.load(stream)
    except (ModuleNotFoundError, ImportError) as error:
        raise ValueError(
            "cannot load the GP snapshot with this Python environment; run with the "
            "same environment used for training (normally the project .venv)"
        ) from error
    if (saved_model.get("training_sample_count") != len(train_rows)
            or saved_model.get("train_sha256") != file_hash(train)
            or tuple(saved_model.get("responses", ())) != RESPONSES
            or tuple(saved_model.get("variables", ())) != VARIABLES
            or saved_model.get("rules") != rules):
        raise ValueError("saved GP model does not match report, rules, or training data")
    fitted = saved_model.get("fitted")
    if not isinstance(fitted, list) or len(fitted) != len(RESPONSES):
        raise ValueError("saved GP model does not contain three fitted responses")

    truth = {name: numeric(revealed[name], f"truth/{name}") for name in RESPONSES}
    train_y = np.asarray([[numeric(row[name], f"training/{name}") for name in RESPONSES]
                          for row in train_rows])
    statistics = data_io.screening_statistics(
        np.vstack((train_x, x)),
        np.vstack((train_y, np.asarray([truth[name] for name in RESPONSES]))),
    )
    predicted, observed_std, latent_std = {}, {}, {}
    for j, name in enumerate(RESPONSES):
        mu, std, latent = model_prediction(fitted[j], x)
        for suffix, calculated in (("mean_m", mu), ("observed_std_m", std), ("latent_std_m", latent)):
            saved = numeric(prediction.get(f"{name}_{suffix}"), f"prediction/{name}_{suffix}")
            if not math.isclose(saved, calculated, rel_tol=1e-8, abs_tol=1e-10):
                raise ValueError(f"sealed prediction differs from saved M_t model: {name}/{suffix}")
        if std <= 0:
            raise ValueError(f"non-positive observed prediction std for {name}")
        predicted[name], observed_std[name], latent_std[name] = mu, std, latent

    coverage, neighbors = nearest_context(x, train_x, train_rows)
    diagnostics = {}
    taylor_details = {}
    for response_index, name in enumerate(RESPONSES):
        info = report["responses"][name]
        if set(info.get("length_scales_unit_space", {})) != set(VARIABLES):
            raise ValueError(f"invalid fitted length scales for {name}")
        residual = truth[name] - predicted[name]
        scale = numeric(info["training_std_m"], f"report/{name}/training_std_m")
        # Protect z against machine-scale sigma without hiding the original sigma.
        sigma_floor = max(1e-12, 1e-6 * scale)
        sigma_used = max(observed_std[name], sigma_floor)
        z = abs(residual) / sigma_used
        sparse = coverage["coverage_class"] == "sparse"
        issue = (
            "large_standardized_error_sparse_region" if sparse and z > 2.0 else
            "small_standardized_error_sparse_region" if sparse and z <= 1.8 else
            "moderate_standardized_error_sparse_region" if sparse else
            "large_standardized_error" if z > 2.0 else
            "no_strong_standardized_anomaly"
        )
        lengths = info["length_scales_unit_space"]
        fitted_lengths = np.asarray(
            [numeric(lengths[v], f"length_scale/{name}/{v}") for v in VARIABLES]
        )
        if np.any(fitted_lengths <= 0):
            raise ValueError(f"fitted length scales must be positive for {name}")
        weighted_distances = np.linalg.norm((train_x - x) / fitted_lengths, axis=1)
        attribution, taylor_details[name] = taylor_attribution(fitted[response_index], x, neighbors)
        diagnostics[name] = {
            "truth_m": truth[name],
            "blind_mean_m": predicted[name],
            "blind_observed_std_m": observed_std[name],
            "blind_latent_std_m": latent_std[name],
            "blind_pi95_low_m": predicted[name] - 1.96 * observed_std[name],
            "blind_pi95_high_m": predicted[name] + 1.96 * observed_std[name],
            "signed_error_m": residual,
            "absolute_error_m": abs(residual),
            "relative_error_percent": (100 * abs(residual) / abs(truth[name])
                                       if abs(truth[name]) > 1e-12 else None),
            "z_absolute_error_over_observed_std": z,
            "sigma_floor_m": sigma_floor,
            "sigma_floor_applied": observed_std[name] < sigma_floor,
            "within_nominal_95_prediction_interval": abs(residual) <= 1.96 * observed_std[name],
            "issue_class": issue,
            "kernel_scaled_nearest_distance": float(np.min(weighted_distances)),
            "kernel_scaled_nearest_sample_index": int(train_rows[int(np.argmin(weighted_distances))]["sample_index"]),
            "fitted_length_scales_unit_space": lengths,
            "taylor_attribution": attribution,
        }
        if prior_mode == "soft_log_normal":
            diagnostics[name]["lengthscale_prior_categories"] = rules[
                "lengthscale_prior_categories"
            ][name]
        else:
            diagnostics[name].update({
                "lengthscale_tiers": rules["lengthscale_tiers"][name],
                "length_scale_at_lower_bound": info["length_scale_at_lower_bound"],
                "length_scale_at_upper_bound": info["length_scale_at_upper_bound"],
                "nugget_tier": rules["nugget_tiers"][name],
                "nugget_variance_standardized": info["nugget_variance_standardized"],
                "variable_attribution": (
                    "insufficient: all three neighbors are relatively far and joint increments "
                    "cannot identify a responsible variable" if sparse else
                    "not established by this diagnostic; even typical relative spacing is not "
                    "a controlled single-variable comparison"
                ),
            })

    prior = None
    if prior_context_path is not None:
        prior_context = read_json(prior_context_path)
        if (prior_context.get("initial_count") != 20
                or prior_context.get("design_variables") != list(VARIABLES)
                or prior_context.get("responses") != list(RESPONSES)):
            raise ValueError("prior context does not match the 20-point experiment")
        initial_fields = ("design_id", *VARIABLES, *RESPONSES)
        initial_payload = [{field: row[field] for field in initial_fields}
                           for row in train_rows[:20]]
        if prior_context.get("initial_observations_sha256") != canonical_hash(initial_payload):
            raise ValueError("prior context belongs to different first-20 observations")
        prior = {"prior_sign": prior_context["prior_sign"],
                 "prior_tier": prior_context["prior_tier"],
                 "note": "engineering direction/strength, not hard GP monotonicity"}

    ledger = None
    ledger_entries = None
    if ledger_path is not None:
        saved_ledger = read_json(ledger_path)
        entries = saved_ledger.get("entries")
        if not isinstance(entries, list):
            raise ValueError("ledger JSON must contain an entries array")
        if any(type(item.get("sample_index")) is not int or item["sample_index"] >= target
               for item in entries if isinstance(item, dict)) or any(not isinstance(item, dict) for item in entries):
            raise ValueError("ledger must contain only prior-step object entries")
        ledger_entries = entries[-5:]
        # Previous free-text LLM reasons can quote historical coordinates.
        # Keep those in the local ledger, but expose only structured decisions.
        ledger = [{
            "sample_index": entry["sample_index"],
            "outcome": entry["outcome"],
            "proposed_candidates": [{key: item[key] for key in (
                "response", "variable", "from_category", "to_category")}
                for item in entry["proposed_candidates"]],
            "adopted_candidates": [{key: item[key] for key in (
                "response", "variable", "from_category", "to_category")}
                for item in entry["adopted_candidates"]],
        } for entry in ledger_entries]

    context = {
        "schema_version": 3,
        "stage": "sequential",
        "sample_index": target,
        "design_id": ids[0][1],
        "training_sample_count": len(train_rows),
        "statistics": statistics,
        "coverage": coverage,
        "responses": diagnostics,
        "current_model": {
            "kernel": report["fit_settings"]["kernel"],
            "lengthscale_rule_mode": prior_mode,
        },
        "engineering_prior": prior,
        "recent_interventions": ledger,
        "recent_interventions_status": "provided" if ledger_path is not None else "not_provided",
        "diagnostic_limit": (
            "No future truths or independent test data are used. The three-neighbor Taylor "
            "ranking describes GP prediction increments, not causes of the observed residual. "
            "A proposed prior change still requires independent model selection."
        ),
    }
    if prior_mode == "legacy_hard_lower_bound":
        context["current_model"].update({
            "directional_prior_enabled": False,
            "local_knowledge_veto_available": False,
        })
    signature = canonical_hash({
        "diagnostic_format_version": 3,
        "training_sha256": file_hash(train),
        "model_sha256": file_hash(model_path),
        "report_sha256": file_hash(report_path),
        "prediction_sha256": file_hash(prediction_path),
        "rules_sha256": file_hash(rules_path),
        "design_row": design,
        "revealed_truth_row": revealed,
        "ledger_entries": ledger_entries,
        "prior_context_sha256": file_hash(prior_context_path) if prior_context_path else None,
    })
    manifest_path = output / "diagnostic_manifest.json"
    if output.exists() and any(output.iterdir()):
        old = read_json(manifest_path) if manifest_path.is_file() else {}
        if old.get("input_signature_sha256") != signature:
            raise ValueError(f"output contains a different diagnosis; choose another --output: {output}")
    output.mkdir(parents=True, exist_ok=True)
    (output / "diagnostic_context.json").write_text(
        json.dumps(context, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    (output / "taylor_contributions.json").write_text(
        json.dumps({"schema_version": 1, "sample_index": target, "responses": taylor_details},
                   ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    manifest_path.write_text(json.dumps({
        "schema_version": 3,
        "input_signature_sha256": signature,
        "sample_index": target,
        "train_source": str(train.resolve()),
        "design_source": str(designs.resolve()),
        "truth_source": str(truth_path.resolve()),
        "model_dir": str(model_dir.resolve()),
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Sample {target} ({ids[0][1]}), training prefix: {len(train_rows)}")
    print(f"Relative sparsity: {coverage['coverage_class']}; "
          f"very-close neighbor: {'yes' if coverage['has_very_close_historical_point'] else 'no'}; "
          "nearest distances: "
          + ", ".join(f"{item['distance_euclidean_unit_space']:.4f}" for item in neighbors))
    for name, item in diagnostics.items():
        print(f"  {name}: error={item['signed_error_m']:+.6g} m, "
              f"sigma={item['blind_observed_std_m']:.6g} m, "
              f"z={item['z_absolute_error_over_observed_std']:.3f}, "
              f"issue={item['issue_class']}")
    print(f"Context: {output / 'diagnostic_context.json'}")
    print(f"Taylor details: {output / 'taylor_contributions.json'}")
    return context


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train", type=Path, default=DATA_DIR / "initial20.csv")
    parser.add_argument("--designs", type=Path, default=DATA_DIR / "sequential_designs.csv")
    parser.add_argument("--truth", type=Path, default=DATA_DIR / "sequential_truth.csv")
    parser.add_argument("--model-dir", type=Path, default=DEFAULT_MODEL_DIR)
    parser.add_argument("--sample-index", type=int,
                        help="must equal the revealed training count plus one")
    parser.add_argument("--ledger", type=Path,
                        help="optional JSON object with an entries array; only earlier steps")
    parser.add_argument("--prior-context", type=Path,
                        default=DATA_DIR / "initial20_context.json",
                        help="initial engineering-prior context")
    parser.add_argument("--output", type=Path,
                        help="diagnostic directory; defaults to outputs/diagnostics/sample_NNNN")
    args = parser.parse_args()
    sample_index = args.sample_index
    if sample_index is None:
        sample_index = len(csv_rows(args.train, ("sample_index", "design_id"))) + 1
    output = args.output or HERE / "outputs" / "diagnostics_summary_taylor_v4" / f"sample_{sample_index:04d}"
    try:
        run(args.train.resolve(), args.designs.resolve(), args.truth.resolve(),
            args.model_dir.resolve(), output.resolve(), sample_index,
            args.ledger.resolve() if args.ledger else None,
            args.prior_context.resolve())
    except (ValueError, OSError, KeyError, IndexError, TypeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
