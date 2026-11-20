"""Fit three independent ARD Gaussian processes from revealed observations.

No arguments: fit M_20 from the prepared initial data and validated LLM priors,
then make a blind prediction for sample 21. Later steps may pass --train and
--rules. The caller must provide only responses already revealed at that step.

Use --loo when a proposed rule change needs a same-data LOO comparison with
the old rule. This script reports metrics but never adopts or rolls back rules.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import pickle
import warnings
from dataclasses import asdict
from pathlib import Path

import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from threadpoolctl import threadpool_limits

try:
    from llm_guided_sequential_gp.sequential_gp_prior import (
        CATEGORY_LOG_SHIFT, LengthscaleMAPOptimizer, PriorConfig,
        lengthscale_theta_indices, log_centers,
    )
except ModuleNotFoundError:
    from sequential_gp_prior import (
        CATEGORY_LOG_SHIFT, LengthscaleMAPOptimizer, PriorConfig,
        lengthscale_theta_indices, log_centers,
    )


HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "outputs" / "semiparametric_guarded_20_20261109"
DEFAULT_DECISION = HERE / "outputs" / "prompt_requests" / "initial" / "validation_result.json"
RESPONSES = (
    "crest_max_settlement_m",
    "dam_max_settlement_m",
    "dam_max_horizontal_displacement_m",
)
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
VARIABLES = tuple(BOUNDS)
UNIT_FIELDS = tuple(f"{name}_unit" for name in VARIABLES)
NUMERICAL_JITTER = 1e-8
INITIAL_COUNT = 20
DEFAULT_PRIOR_CONFIG = PriorConfig()
# The still-unmigrated policy imports these names. They are never used by the
# version-2 fitter and can be removed when that policy is converted.
LENGTH_FLOORS = {"Tier-1": 0.15, "Tier-2": 0.60, "Tier-3": 2.00}
NUGGET_FLOORS = {"Nugget-Low": 1e-4, "Nugget-Mid": 1e-3, "Nugget-High": 1e-2}


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_sha256(value: object) -> str:
    encoded = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def unique_json_pairs(pairs: list[tuple[str, object]]) -> dict:
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_json_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON number: {value}")


def read_json(path: Path) -> dict:
    try:
        value = json.loads(
            path.read_text(encoding="utf-8-sig"),
            object_pairs_hook=unique_json_pairs,
            parse_constant=reject_json_constant,
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(f"cannot read JSON file {path}: {error}") from error
    if not isinstance(value, dict):
        raise ValueError(f"JSON root must be an object: {path}")
    return value


def read_csv(path: Path, required: tuple[str, ...]) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"missing CSV file: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames
        if fields is None or len(fields) != len(set(fields)):
            raise ValueError(f"missing or duplicate CSV headers: {path}")
        missing = set(required) - set(fields)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        rows = list(reader)
    if not rows or any(None in row for row in rows):
        raise ValueError(f"empty CSV or extra row fields: {path}")
    return rows


def parse_rows(rows: list[dict[str, str]], *, training: bool) -> tuple[np.ndarray, np.ndarray | None]:
    indices = []
    ids = []
    try:
        for row in rows:
            indices.append(int(row["sample_index"]))
            ids.append(row["design_id"].strip())
        x = np.asarray([[float(row[name]) for name in UNIT_FIELDS] for row in rows])
        y = (
            np.asarray([[float(row[name]) for name in RESPONSES] for row in rows])
            if training else None
        )
    except (ValueError, TypeError, KeyError) as error:
        raise ValueError("sample index, normalized design, or response is not numeric") from error
    if any(not value for value in ids) or len(set(ids)) != len(ids):
        raise ValueError("design_id must be nonempty and unique")
    if len(set(indices)) != len(indices):
        raise ValueError("sample_index values must be unique")
    if not np.isfinite(x).all() or np.any(x < -1e-10) or np.any(x > 1 + 1e-10):
        raise ValueError("normalized designs must be finite and inside [0, 1]")
    if np.unique(x, axis=0).shape[0] != len(x):
        raise ValueError("duplicate normalized design vectors")
    if y is not None and not np.isfinite(y).all():
        raise ValueError("training responses must be finite")
    physical_fields = set(VARIABLES) & set(rows[0])
    if physical_fields and physical_fields != set(VARIABLES):
        raise ValueError("physical design columns must be complete when supplied")
    if physical_fields:
        physical = np.asarray([[float(row[name]) for name in VARIABLES] for row in rows])
        bounds = np.asarray(list(BOUNDS.values()))
        recomputed = (physical - bounds[:, 0]) / (bounds[:, 1] - bounds[:, 0])
        if not np.isfinite(physical).all() or not np.allclose(x, recomputed, atol=1e-9, rtol=0):
            raise ValueError("physical and normalized design columns disagree")
    if training:
        if len(rows) < INITIAL_COUNT:
            raise ValueError(f"training data must contain at least {INITIAL_COUNT} rows")
        if indices != list(range(1, len(rows) + 1)):
            raise ValueError("training data must contain a contiguous revealed prefix")
    else:
        if set(RESPONSES) & set(rows[0]):
            raise ValueError("blind prediction input must not contain response columns")
    return x, y


def validate_rules(rules: dict) -> dict:
    expected = {"schema_version", "responses", "design_variables", "lengthscale_prior_categories"}
    if (not isinstance(rules, dict) or set(rules) != expected
            or type(rules.get("schema_version")) is not int or rules["schema_version"] != 2):
        raise ValueError("rules must use the version-2 soft length-scale prior schema")
    if rules["responses"] != list(RESPONSES) or rules["design_variables"] != list(VARIABLES):
        raise ValueError("rules response or design-variable order does not match this experiment")
    categories = rules["lengthscale_prior_categories"]
    if not isinstance(categories, dict) or set(categories) != set(RESPONSES):
        raise ValueError("rules must have one length-scale prior mapping per response")
    for response in RESPONSES:
        values = categories[response]
        if not isinstance(values, dict) or set(values) != set(VARIABLES):
            raise ValueError(f"{response} must have eight named prior categories")
        if any(type(value) is not str or value not in CATEGORY_LOG_SHIFT
               for value in values.values()):
            raise ValueError(f"{response} has an illegal length-scale prior category")
    return rules


def rules_from_initial_decision(path: Path, initial_rows: list[dict[str, str]]) -> dict:
    result = read_json(path)
    answer = result.get("answer")
    if result.get("valid") is not True or result.get("stage") != "initial":
        raise ValueError("initial LLM output must be a validated initial-stage decision")
    if not isinstance(answer, dict) or answer.get("status") != "OK":
        raise ValueError("initial LLM decision is not OK")
    context = read_json(path.parent / "context.json")
    if result.get("context_sha256") != canonical_sha256(context):
        raise ValueError("initial decision does not match its saved prompt context")
    relevant = ("design_id", *VARIABLES, *RESPONSES)
    initial_payload = [{name: row[name] for name in relevant} for row in initial_rows]
    if context.get("initial_observations_sha256") != canonical_sha256(initial_payload):
        raise ValueError("initial decision belongs to different first-20 observations")
    if context.get("schema_version") != 2 or "initial_observations" in context:
        raise ValueError("initial context must be the statistics-only version")
    if answer.get("schema_version") != 2:
        raise ValueError("initial answer must use version-2 prior categories")
    assignments = answer.get("assignments")
    expected = [(response, variable) for response in RESPONSES for variable in VARIABLES]
    if not isinstance(assignments, list) or len(assignments) != len(expected):
        raise ValueError("initial LLM decision must contain 24 assignments")
    categories = {response: {} for response in RESPONSES}
    for row, (response, variable) in zip(assignments, expected, strict=True):
        if (not isinstance(row, dict) or set(row) != {"response", "variable", "prior_category", "reason"}
                or (row.get("response"), row.get("variable")) != (response, variable)):
            raise ValueError(f"initial assignment order differs at {response}/{variable}")
        if type(row.get("reason")) is not str or not row["reason"].strip():
            raise ValueError(f"initial assignment {response}/{variable} has no reason")
        categories[response][variable] = row.get("prior_category")
    return validate_rules({
        "schema_version": 2,
        "responses": list(RESPONSES),
        "design_variables": list(VARIABLES),
        "lengthscale_prior_categories": categories,
    })


def fit_one(
    x: np.ndarray, y: np.ndarray, rules: dict, response: str,
    *, restarts: int, seed: int, prior_config: PriorConfig = DEFAULT_PRIOR_CONFIG,
) -> tuple[GaussianProcessRegressor, float, float, list[str]]:
    mean = float(np.mean(y))
    scale = float(np.std(y, ddof=1))
    if not np.isfinite(scale) or scale <= 1e-12 * max(1.0, float(np.max(np.abs(y)))):
        raise ValueError(f"{response} has nearly constant training values")
    centers = log_centers(rules["lengthscale_prior_categories"][response],
                          VARIABLES, prior_config)
    common_start = prior_config.base_factor * math.sqrt(len(VARIABLES))
    kernel = ConstantKernel(1.0, (1e-2, 1e2)) * Matern(
        length_scale=np.full(len(VARIABLES), common_start),
        length_scale_bounds=(prior_config.length_scale_lower,
                             prior_config.length_scale_upper),
        nu=2.5,
    ) + WhiteKernel(
        noise_level=prior_config.noise_initial,
        noise_level_bounds=(prior_config.noise_lower, prior_config.noise_upper),
    )
    optimizer = LengthscaleMAPOptimizer(
        centers, prior_config.log_std, lengthscale_theta_indices(kernel, len(VARIABLES))
    )
    model = GaussianProcessRegressor(
        kernel=kernel,
        alpha=NUMERICAL_JITTER,
        normalize_y=False,
        optimizer=optimizer,
        n_restarts_optimizer=restarts,
        random_state=seed,
    )
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        model.fit(x, (y - mean) / scale)
    # The snapshot is for prediction. Keep it independent of how the fitting
    # script was invoked, so sklearn does not pickle a custom optimizer class.
    model.optimizer = None
    return model, mean, scale, [str(item.message) for item in captured]


def predict_one(
    fitted: tuple[GaussianProcessRegressor, float, float, list[str]],
    x: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    model, mean, scale, _ = fitted
    prediction, observed_std = model.predict(x, return_std=True)
    noise_variance = float(model.kernel_.k2.noise_level)
    latent_variance = np.maximum(observed_std**2 - noise_variance, 0.0)
    return (
        mean + scale * prediction,
        scale * observed_std,
        scale * np.sqrt(latent_variance),
    )


def fit_all(x: np.ndarray, y: np.ndarray, rules: dict, *, restarts: int, seed: int,
            prior_config: PriorConfig = DEFAULT_PRIOR_CONFIG) -> list:
    with threadpool_limits(limits=1):
        return [
            fit_one(x, y[:, index], rules, response, restarts=restarts,
                    seed=seed + index, prior_config=prior_config)
            for index, response in enumerate(RESPONSES)
        ]


def loo_validation(
    x: np.ndarray, y: np.ndarray, rules: dict, reference_scale: np.ndarray,
    *, restarts: int, seed: int, prior_config: PriorConfig = DEFAULT_PRIOR_CONFIG,
) -> tuple[dict, list[dict]]:
    rows = []
    metrics = {}
    with threadpool_limits(limits=1):
        for response_index, response in enumerate(RESPONSES):
            predicted = np.empty(len(x))
            predicted_std = np.empty(len(x))
            fold_warnings = 0
            for held in range(len(x)):
                train_mask = np.arange(len(x)) != held
                model = fit_one(
                    x[train_mask], y[train_mask, response_index], rules, response,
                    restarts=restarts, seed=seed + response_index,
                    prior_config=prior_config,
                )
                fold_warnings += len(model[3])
                mu, std, _ = predict_one(model, x[held : held + 1])
                predicted[held], predicted_std[held] = mu[0], std[0]
                rows.append({
                    "sample_index": held + 1,
                    "response": response,
                    "observed_m": float(y[held, response_index]),
                    "loo_mean_m": float(mu[0]),
                    "loo_observed_std_m": float(std[0]),
                })
            rmse = float(np.sqrt(np.mean((y[:, response_index] - predicted) ** 2)))
            metrics[response] = {
                "loo_rmse_m": rmse,
                "loo_nrmse": rmse / float(reference_scale[response_index]),
                "nrmse_denominator_initial20_std_m": float(reference_scale[response_index]),
                "folds": len(x),
                "fold_warning_count": fold_warnings,
            }
    return metrics, rows


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: dict) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def run(
    train_path: Path, reference_path: Path, output: Path, *,
    rules_path: Path | None, initial_decision: Path,
    predict_path: Path | None, predict_index: int | None,
    restarts: int, seed: int, loo: bool,
    prior_config: PriorConfig = DEFAULT_PRIOR_CONFIG,
) -> dict:
    if restarts < 0:
        raise ValueError("--restarts must be nonnegative")
    prior_config.validate(len(VARIABLES))
    train_rows = read_csv(train_path, ("sample_index", "design_id", *UNIT_FIELDS, *RESPONSES))
    x, y = parse_rows(train_rows, training=True)
    reference_rows = read_csv(reference_path, ("sample_index", "design_id", *UNIT_FIELDS, *RESPONSES))
    if len(reference_rows) != INITIAL_COUNT:
        raise ValueError("reference initial dataset must have exactly 20 rows")
    x_ref, y_ref = parse_rows(reference_rows, training=True)
    if (
        [row["design_id"] for row in train_rows[:INITIAL_COUNT]]
        != [row["design_id"] for row in reference_rows]
        or not np.allclose(x[:INITIAL_COUNT], x_ref, rtol=0, atol=1e-12)
        or not np.allclose(y[:INITIAL_COUNT], y_ref, rtol=0, atol=1e-12)
    ):
        raise ValueError("training prefix does not match the fixed initial 20 points")
    reference_scale = np.std(y_ref, axis=0, ddof=1)
    if np.any(reference_scale <= 0) or not np.isfinite(reference_scale).all():
        raise ValueError("initial 20 responses must have positive standard deviations")

    if rules_path is None:
        if len(x) != INITIAL_COUNT:
            raise ValueError("use --rules for training beyond the initial 20 points")
        rules = rules_from_initial_decision(initial_decision, train_rows)
        rule_source = initial_decision
    else:
        rules = validate_rules(read_json(rules_path))
        rule_source = rules_path
    rules_hash = canonical_sha256(rules)

    prediction_rows = []
    prediction_source = None
    if predict_path is not None:
        blind_rows = read_csv(predict_path, ("sample_index", "design_id", *UNIT_FIELDS))
        blind_x, _ = parse_rows(blind_rows, training=False)
        target = predict_index if predict_index is not None else len(x) + 1
        if target != len(x) + 1:
            raise ValueError("blind prediction target must be the next sample after training")
        selected = [i for i, row in enumerate(blind_rows) if int(row["sample_index"]) == target]
        if len(selected) != 1:
            raise ValueError(f"prediction input needs exactly one sample_index={target}")
        index = selected[0]
        prediction_rows = [(blind_rows[index], blind_x[index : index + 1])]
        prediction_source = predict_path

    signature = canonical_sha256({
        "train_sha256": file_sha256(train_path),
        "reference_sha256": file_sha256(reference_path),
        "rules_sha256": rules_hash,
        "prediction_sha256": file_sha256(prediction_source) if prediction_source else None,
        "prediction_index": predict_index if predict_index is not None else len(x) + 1,
        "restarts": restarts, "seed": seed, "loo": loo,
        "prior_config": asdict(prior_config), "model_schema_version": 2,
    })
    report_path = output / "fit_report.json"
    if report_path.exists():
        prior_report = read_json(report_path)
        if prior_report.get("input_signature_sha256") != signature:
            raise ValueError(f"output contains a different fit; choose another --output: {output}")
    elif output.exists() and any(output.iterdir()):
        raise ValueError(f"output is not an empty model directory: {output}")

    fitted = fit_all(x, y, rules, restarts=restarts, seed=seed,
                     prior_config=prior_config)
    fit_info = {}
    for response, (model, mean, scale, warning_list) in zip(RESPONSES, fitted, strict=True):
        centers = log_centers(rules["lengthscale_prior_categories"][response],
                              VARIABLES, prior_config)
        indices = lengthscale_theta_indices(model.kernel_, len(VARIABLES))
        optimizer = LengthscaleMAPOptimizer(centers, prior_config.log_std, indices)
        lengths = np.asarray(model.kernel_.k1.k2.length_scale)
        noise = float(model.kernel_.k2.noise_level)
        raw_lml = float(model.log_marginal_likelihood(model.kernel_.theta))
        penalty = optimizer.penalty(model.kernel_.theta)
        map_score = raw_lml - penalty
        if (not math.isfinite(map_score)
                or not math.isclose(map_score, model.log_marginal_likelihood_value_,
                                    rel_tol=1e-8, abs_tol=1e-8)):
            raise ValueError(f"MAP score disagrees with fitted optimizer objective for {response}")
        fit_info[response] = {
            "training_mean_m": mean,
            "training_std_m": scale,
            "fitted_kernel": str(model.kernel_),
            "lengthscale_prior_categories": rules["lengthscale_prior_categories"][response],
            "lengthscale_prior_centers_unit_space": dict(zip(
                VARIABLES, np.exp(centers).tolist(), strict=True
            )),
            "length_scales_unit_space": dict(zip(VARIABLES, lengths.tolist(), strict=True)),
            "length_scale_at_optimizer_lower_bound": dict(zip(
                VARIABLES, (lengths <= prior_config.length_scale_lower * 1.01).tolist(), strict=True
            )),
            "length_scale_at_optimizer_upper_bound": dict(zip(
                VARIABLES, (lengths >= prior_config.length_scale_upper * 0.99).tolist(), strict=True
            )),
            "raw_log_marginal_likelihood": raw_lml,
            "lengthscale_prior_penalty": penalty,
            "map_log_objective_without_prior_constant": map_score,
            "nugget_variance_standardized": noise,
            "nugget_variance_m2": noise * scale**2,
            "warnings": warning_list,
        }

    predictions = []
    for source_row, point in prediction_rows:
        item = {"sample_index": int(source_row["sample_index"]), "design_id": source_row["design_id"]}
        for response, model in zip(RESPONSES, fitted, strict=True):
            mu, observed_std, latent_std = predict_one(model, point)
            item.update({
                f"{response}_mean_m": float(mu[0]),
                f"{response}_observed_std_m": float(observed_std[0]),
                f"{response}_latent_std_m": float(latent_std[0]),
                f"{response}_pi95_low_m": float(mu[0] - 1.96 * observed_std[0]),
                f"{response}_pi95_high_m": float(mu[0] + 1.96 * observed_std[0]),
            })
        predictions.append(item)

    loo_metrics, loo_rows = (
        loo_validation(x, y, rules, reference_scale, restarts=restarts,
                       seed=seed, prior_config=prior_config)
        if loo else (None, [])
    )
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / "rules_used.json", rules)
    with (output / "gp_models.pkl").open("wb") as stream:
        pickle.dump({
            "schema_version": 2,
            "responses": RESPONSES,
            "variables": VARIABLES,
            "fitted": fitted,
            "rules": rules,
            "training_sample_count": len(x),
            "train_sha256": file_sha256(train_path),
            "prior_config": asdict(prior_config),
        }, stream, protocol=pickle.HIGHEST_PROTOCOL)
    if predictions:
        write_csv(output / "blind_prediction.csv", list(predictions[0]), predictions)
    if loo:
        write_csv(output / "loo_predictions.csv", list(loo_rows[0]), loo_rows)
    report = {
        "schema_version": 2,
        "training_sample_count": len(x),
        "train_file": str(train_path),
        "train_sha256": file_sha256(train_path),
        "reference_initial_file": str(reference_path),
        "reference_initial_sha256": file_sha256(reference_path),
        "rule_source": str(rule_source),
        "rule_source_sha256": file_sha256(rule_source),
        "rules_sha256": rules_hash,
        "input_signature_sha256": signature,
        "fit_settings": {
            "kernel": "Constant * ARD Matern(nu=2.5) + WhiteKernel",
            "lengthscale_prior_mode": "soft_log_normal_map",
            "lengthscale_prior_base_factor": prior_config.base_factor,
            "lengthscale_prior_log_std": prior_config.log_std,
            "lengthscale_prior_category_log_shifts": CATEGORY_LOG_SHIFT,
            "length_scale_optimizer_bounds": [prior_config.length_scale_lower,
                                              prior_config.length_scale_upper],
            "nugget_optimizer_bounds_standardized_variance": [prior_config.noise_lower,
                                                               prior_config.noise_upper],
            "nugget_optimizer_initial_standardized_variance": prior_config.noise_initial,
            "numerical_jitter_standardized_variance": NUMERICAL_JITTER,
            "optimizer_restarts": restarts,
            "random_seed": seed,
            "pi95_definition": "nominal Gaussian prediction interval: mean +/- 1.96 * observed_std",
        },
        "responses": fit_info,
        "blind_prediction_source": str(prediction_source) if prediction_source else None,
        "blind_prediction_sample_index": predictions[0]["sample_index"] if predictions else None,
        "loo_metrics": loo_metrics,
        "scope": "Only supplied training rows are fitted. LOO is diagnostic; this script does not adopt or roll back rules.",
    }
    write_json(report_path, report)
    print(f"Training: {len(x)} revealed samples from {train_path}")
    print(f"Rules: {rule_source}")
    for response in RESPONSES:
        info = fit_info[response]
        bound_hits = [v for v, hit in info["length_scale_at_optimizer_upper_bound"].items() if hit]
        print(f"  {response}: nugget={info['nugget_variance_standardized']:.6g}, "
              f"MAP penalty={info['lengthscale_prior_penalty']:.6g}, "
              f"numeric upper-bound hits={bound_hits}")
    if predictions:
        print(f"Blind prediction: sample {predictions[0]['sample_index']} ({predictions[0]['design_id']})")
        for response in RESPONSES:
            print(f"  {response}: mean={predictions[0][f'{response}_mean_m']:.6g} m, "
                  f"observed std={predictions[0][f'{response}_observed_std_m']:.6g} m")
    if loo_metrics:
        for response in RESPONSES:
            print(f"  LOO {response}: NRMSE={loo_metrics[response]['loo_nrmse']:.6g}")
    print(f"Output: {output}")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fit three independent ARD GPs from revealed observations.",
        epilog="No arguments fit M_20 and predict sample 21. Use --loo for same-data rule comparison.",
    )
    parser.add_argument("--train", type=Path, default=DATA_DIR / "initial20.csv",
                        help="CSV containing only the currently revealed training prefix")
    parser.add_argument("--reference-initial", type=Path, default=DATA_DIR / "initial20.csv",
                        help="fixed first-20 CSV defining the NRMSE denominator")
    parser.add_argument("--initial-decision", type=Path, default=DEFAULT_DECISION,
                        help="validated cold-start LLM result, used when --rules is omitted")
    parser.add_argument("--rules", type=Path,
                        help="previously adopted rules_used.json or an equivalent edited rules file")
    parser.add_argument("--predict-designs", type=Path, default=DATA_DIR / "sequential_designs.csv",
                        help="design-only CSV containing the next sample")
    parser.add_argument("--predict-sample-index", type=int,
                        help="must equal the current training count plus one")
    parser.add_argument("--no-predict", action="store_true",
                        help="fit without producing a blind prediction")
    parser.add_argument("--loo", action="store_true",
                        help="refit every held-out fold for rule-comparison metrics")
    parser.add_argument("--restarts", type=int, default=3)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--prior-base-factor", type=float,
                        default=DEFAULT_PRIOR_CONFIG.base_factor,
                        help="fixed c in the log-length-scale prior center c*sqrt(8)")
    parser.add_argument("--prior-log-std", type=float,
                        default=DEFAULT_PRIOR_CONFIG.log_std,
                        help="fixed standard deviation of the Gaussian log-length-scale prior")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    train_path = args.train.resolve()
    try:
        # The output name follows the supplied training prefix. It is computed
        # before fitting so different shadow rules can use separate --output dirs.
        count = len(read_csv(train_path, ("sample_index", "design_id", *UNIT_FIELDS, *RESPONSES)))
        output = (args.output or HERE / "outputs" / "models" / f"samples_{count:04d}").resolve()
        predict_path = None if args.no_predict else args.predict_designs.resolve()
        if predict_path is not None:
            design_rows = read_csv(predict_path, ("sample_index", "design_id", *UNIT_FIELDS))
            parse_rows(design_rows, training=False)
            available_indices = [int(row["sample_index"]) for row in design_rows]
            if count >= max(available_indices):
                predict_path = None
        run(
            train_path,
            args.reference_initial.resolve(),
            output,
            rules_path=args.rules.resolve() if args.rules else None,
            initial_decision=args.initial_decision.resolve(),
            predict_path=predict_path,
            predict_index=args.predict_sample_index,
            restarts=args.restarts,
            seed=args.seed,
            loo=args.loo,
            prior_config=PriorConfig(base_factor=args.prior_base_factor,
                                     log_std=args.prior_log_std),
        )
    except (ValueError, OSError, np.linalg.LinAlgError) as error:
        parser.exit(2, f"GP fitting failed: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
