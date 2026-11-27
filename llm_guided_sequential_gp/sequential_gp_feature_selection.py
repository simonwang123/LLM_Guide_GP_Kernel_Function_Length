"""Screen a validated physical feature pool and replay sequential GP predictions."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

if __package__:
    from . import sequential_gp_feature_candidates as candidates
    from . import sequential_gp_model as gp
    from . import sequential_gp_primitives as primitives
else:
    import sequential_gp_feature_candidates as candidates
    import sequential_gp_model as gp
    import sequential_gp_primitives as primitives


ABS_STD_FLOOR = 1e-8
REL_STD_FLOOR = 1e-8
REDUNDANCY_LIMIT = 0.95
FEATURE_CAP = 2
COLD_START_MIN_DELTA_MAP = (0.50, 1.00)
COLD_START_MIN_PHYSICAL_SHARE = 0.20
COLD_START_MAX_PENALTY_INCREASE = 1.0
COLD_START_MIN_ANA_LOO_GAIN = 0.0
COLD_START_MIN_ANA_NLPD_GAIN = 0.05
IDENTITY = ("sample_index", "design_id")
DESIGN_FIELDS = (*IDENTITY, *gp.VARIABLES, *gp.UNIT_FIELDS)
FORMULA_FIELDS = ("name", "expression", "physical_rationale")
FORMULA_ERRORS = (ValueError, SyntaxError, TypeError, OverflowError)


def design_namespace(rows, config):
    fields = ["design_id", *gp.VARIABLES]
    designs = pd.DataFrame(rows, columns=fields).astype(dict.fromkeys(gp.VARIABLES, float))
    return primitives.build_formula_namespace(designs, config=config)


def feature_values(feature, namespace, count):
    formula = {field: feature[field] for field in FORMULA_FIELDS}
    return candidates.evaluate_candidate(formula, namespace, count, min_std=0)[1]


def selected_formulas(report: dict) -> dict[str, list[dict]]:
    responses = report.get("responses")
    if report.get("schema_version") not in (1, 2) or not isinstance(responses, dict) or set(responses) != set(gp.RESPONSES):
        raise ValueError("expected a feature_selection.json with the three constraint responses")
    selected = {}
    for response in gp.RESPONSES:
        result = responses[response]
        names = result.get("selected") if isinstance(result, dict) else None
        screening = result.get("screening") if isinstance(result, dict) else None
        if not isinstance(names, list) or not isinstance(screening, list):
            raise ValueError(f"invalid selected features or screening records for {response}")
        records = {}
        for item in screening:
            if not isinstance(item, dict) or any(not isinstance(item.get(key), str) for key in FORMULA_FIELDS):
                raise ValueError(f"invalid screening record for {response}")
            name = item["name"]
            if name in records:
                raise ValueError(f"duplicate screening record: {response}/{name}")
            records[name] = item
        if (any(not isinstance(name, str) or name not in records for name in names)
                or len(set(names)) != len(names)):
            raise ValueError(f"invalid selected features for {response}")
        selected[response] = [{key: records[name][key] for key in FORMULA_FIELDS} for name in names]
    return selected


def load_pool(path, initial_rows, namespace):
    pool = gp.read_json(path)
    initial_ids = [row["design_id"].strip() for row in initial_rows]
    responses = pool.get("responses")
    if (pool.get("schema_version") != 1 or pool.get("design_ids") != initial_ids
            or not isinstance(responses, dict) or set(responses) != set(gp.RESPONSES)):
        raise ValueError("candidate_pool.json must match the initial 20 designs and response schema")
    initial_namespace = {key: values[:len(initial_ids)] for key, values in namespace.items()}
    valid = {response: [] for response in responses}
    rejected = {response: [] for response in responses}
    for response, items in responses.items():
        if not isinstance(items, list) or not candidates.MIN_CANDIDATES <= len(items) <= candidates.MAX_CANDIDATES:
            raise ValueError(
                f"{response} needs {candidates.MIN_CANDIDATES}–{candidates.MAX_CANDIDATES} saved candidate records"
            )
        names = set()
        for feature in items:
            if not isinstance(feature, dict) or set(feature) != {*FORMULA_FIELDS, "unit", "values"}:
                raise ValueError(f"invalid candidate record for {response}")
            if any(not isinstance(feature[key], str) for key in (*FORMULA_FIELDS, "unit")):
                raise ValueError(f"candidate formula fields and unit must be strings for {response}")
            name = feature["name"]
            if name in names:
                raise ValueError(f"duplicate candidate name: {response}/{name}")
            names.add(name)
            try:
                saved = np.asarray(feature["values"], dtype=float)
            except (TypeError, ValueError, OverflowError) as exc:
                raise ValueError(f"invalid saved values: {response}/{name}") from exc
            if saved.shape != (len(initial_ids),) or not np.isfinite(saved).all():
                raise ValueError(f"invalid saved values: {response}/{name}")

            formula = {key: feature[key] for key in FORMULA_FIELDS}
            try:
                evaluated, measured = candidates.evaluate_candidate(formula, initial_namespace, len(initial_ids), min_std=0)
            except FORMULA_ERRORS as exc:
                rejected[response].append({
                    **formula, "unit": feature["unit"], "status": "invalid_pool_formula", "reason": str(exc),
                })
                continue
            if evaluated["unit"] != feature["unit"] or not np.allclose(measured, saved, atol=1e-12, rtol=1e-10):
                raise ValueError(f"initial geometry or unit mismatch: {response}/{name}")
            valid[response].append(feature)
    return valid, rejected


def migration_outcomes(previous, result):
    screening = {item["name"]: item for item in result["screening"]}
    outcomes = []
    for formula in previous:
        name = formula["name"]
        item = screening.get(name)
        if item is None:
            outcome = {"status": "absent_from_pool"}
        elif any(formula[key] != item[key] for key in FORMULA_FIELDS):
            outcome = {"status": "formula_changed"}
        elif name in result["selected"]:
            outcome = {"status": "retained"}
        elif item["status"] == "eligible":
            outcome = {"status": "not_selected_by_validation", "stop_reason": result["stop_reason"]}
        else:
            outcome = {key: item[key] for key in ("status", "reason", "retained_instead") if key in item}
        outcomes.append({**formula, **outcome})
    return {"previous_features": outcomes, "selected": result["selected"]}


def std_floor(values):
    return np.maximum(ABS_STD_FLOOR, REL_STD_FLOOR * np.max(np.abs(values), axis=0))


def screen_candidates(x, y, pool, namespace, *, scaling):
    values, scaled_values = {}, {}
    audit = {}
    ref_namespace = scaling.namespace if pool else None
    for feature in pool:
        name = feature["name"]
        item = audit[name] = {field: feature[field] for field in (*FORMULA_FIELDS, "unit")}
        try:
            column = feature_values(feature, namespace, len(x))
        except FORMULA_ERRORS as exc:
            item.update(status="invalid_formula_values", reason=str(exc))
            continue
        formula = {field: feature[field] for field in FORMULA_FIELDS}
        try:
            candidates.check_feasible_domain_sanity(formula, column, ref_namespace)
        except FORMULA_ERRORS as exc:
            item.update(status="invalid_feasible_domain", reason=str(exc))
            continue
        with np.errstate(over="ignore", invalid="ignore"):
            variance = float(np.var(column))
        floor = float(std_floor(column))
        item.update(variance=variance if math.isfinite(variance) else None, std_floor=floor)
        if not math.isfinite(variance):
            item.update(status="invalid_variance")
            continue
        elif math.sqrt(variance) <= floor:
            item.update(status="near_constant")
            continue
        try:
            scaled = scaling.transform(formula, column)
        except FORMULA_ERRORS as exc:
            item.update(status="invalid_reference_scaling", reason=str(exc))
            continue
        item["reference_bounds"] = scaling.bounds(formula)
        if scaled.std() <= std_floor(scaled):
            item.update(status="near_constant_scaled")
            continue
        values[name] = column
        scaled_values[name] = scaled
        item.update(status="eligible")

    # The design variables remain in all models. Their correlations are diagnostic.
    eligible = list(values)
    if not eligible:
        return values, list(audit.values())
    width = x.shape[1]
    matrix = pd.DataFrame(np.column_stack([x, *values.values(), y]))
    correlations = matrix.corr(method="spearman").to_numpy()
    candidate_rhos = correlations[width:-1, width:-1]
    response_rhos = correlations[width:-1, -1]
    if not np.isfinite(response_rhos).all():
        raise ValueError("response Spearman correlations are undefined")
    for i, name in enumerate(eligible):
        item = audit[name]
        item["rho_response"] = float(response_rhos[i])
        item["rho_original"] = {
            key: float(rho) if math.isfinite(rho) else None
            for key, rho in zip(gp.VARIABLES, correlations[width + i, :width], strict=True)
        }
        item["rho_candidates"] = dict(zip(eligible, candidate_rhos[i].tolist(), strict=True))

    retained = []
    # Stronger response association wins. Exact ties preserve candidate input order.
    for i in sorted(range(len(eligible)), key=lambda i: -abs(response_rhos[i])):
        name = eligible[i]
        conflict = next((j for j in retained if abs(candidate_rhos[i, j]) > REDUNDANCY_LIMIT), None)
        if conflict is None:
            retained.append(i)
        else:
            audit[name].update(status="correlated_candidate", retained_instead=eligible[conflict],
                               rho_with_retained=float(candidate_rhos[i, conflict]))
    return {eligible[i]: scaled_values[eligible[i]] for i in retained}, list(audit.values())


def score(y, prediction, reference_scale):
    residual = y - prediction
    rmse = float(np.sqrt(np.mean(residual**2)))
    return {"rmse_m": rmse, "nrmse": rmse / reference_scale,
            "mae_m": float(np.mean(np.abs(residual)))}


def relative_gain(base, new):
    return 1.0 - new / base if base > 0 else None


def stage_a_screen(x, y, extra, response, reference_scale, settings):
    """Run Stage-1 full-fit Bayesian evidence and analytical GP LOO."""
    if __package__:
        from . import sequential_gp_additive as additive
    else:
        import sequential_gp_additive as additive

    full_restarts = max(int(settings["restarts"]), 8)
    ev = additive.fit_full_evidence(
        x, y, extra, response, reference_scale, restarts=full_restarts, seed=settings["seed"],
    )
    ana = ev["ana_loo"]
    summary = {
        "raw_lml": ev["raw_lml"],
        "map_score": ev["map_score"],
        "prior_penalty": float(ev["raw_lml"] - ev["map_score"]),
        "physical_signal_share_mean": ev["physical_signal_share_mean"],
        "physical_length_scales": list(ev["length_scales"]["physical"]),
        "ana_rmse_m": ana["rmse_m"],
        "ana_nrmse": ana["nrmse"],
        "ana_mae_m": ana["mae_m"],
        "ana_nlpd": ana["nlpd"],
        "warning_count": ev["warning_count"],
        "optimizer_failure_count": ev["optimizer_failure_count"],
    }
    return summary, np.asarray(ana["predicted_m"], dtype=float)


def stage_a_calibration(x, y, extra, response, reference_scale, settings):
    """Run Stage-2 refitted Repeated Stratified 5x3 CV NLPD calibration guard."""
    if __package__:
        from . import sequential_gp_additive as additive
    else:
        import sequential_gp_additive as additive

    cv_restarts = max(1, min(int(settings["restarts"]), 5))
    rep = additive.evaluate_stage_a_calibration(
        x, y, extra, response, reference_scale, restarts=cv_restarts, seed=settings["seed"],
    )
    return {
        "rep5f_rmse_m": rep["rmse_m"],
        "rep5f_nrmse": rep["nrmse"],
        "rep5f_mae_m": rep["mae_m"],
        "rep5f_nlpd": rep["nlpd"],
        "rep5f_repeat_nlpds": rep["repeat_nlpds"],
        "rep5f_warning_count": rep["warning_count"],
    }, np.asarray(rep["predicted_m"], dtype=float)


def select_response(x, y, pool, namespace, response, reference_scale, settings, *, scaling):
    if __package__:
        from . import sequential_gp_additive as additive
    else:
        import sequential_gp_additive as additive

    eligible, audit = screen_candidates(x, y, pool, namespace, scaling=scaling)
    extra = np.empty((len(x), 0))
    base_screen, base_ana_pred = stage_a_screen(x, y, extra, response, reference_scale, settings)
    base_cal, base_rep_pred = stage_a_calibration(x, y, extra, response, reference_scale, settings)
    baseline = {**base_screen, **base_cal}
    current, current_pred = baseline, base_rep_pred
    selected, rounds = [], []
    stop = "feature_limit"
    for step_idx, min_delta_map in enumerate(COLD_START_MIN_DELTA_MAP, start=1):
        if not eligible:
            stop = "no_more_candidates"
            break
        trials = {}
        for name, column in eligible.items():
            trial_extra = np.column_stack([extra, column])
            screen, ana_pred = stage_a_screen(x, y, trial_extra, response, reference_scale, settings)
            delta_map = float(screen["map_score"] - current["map_score"])
            delta_lml = float(screen["raw_lml"] - current["raw_lml"])
            penalty_increase = float(screen["prior_penalty"] - current["prior_penalty"])
            ana_gain = relative_gain(current["ana_nrmse"], screen["ana_nrmse"]) or 0.0
            ana_nlpd_gain = float(current["ana_nlpd"] - screen["ana_nlpd"])
            phys_share = float(screen["physical_signal_share_mean"])
            phys_lengths = screen.get("physical_length_scales", ())
            ard_active = (
                len(phys_lengths) == trial_extra.shape[1]
                and all(isinstance(v, (int, float)) and math.isfinite(v) and v <= additive.MAX_ACTIVE_PHYSICAL_LENGTH_SCALE
                        for v in phys_lengths)
            )
            passes_gate = (
                delta_map >= min_delta_map
                and phys_share >= COLD_START_MIN_PHYSICAL_SHARE
                and ard_active
                and penalty_increase <= COLD_START_MAX_PENALTY_INCREASE
                and ana_gain >= COLD_START_MIN_ANA_LOO_GAIN
                and ana_nlpd_gain >= COLD_START_MIN_ANA_NLPD_GAIN
                and screen["optimizer_failure_count"] <= current["optimizer_failure_count"]
            )
            if passes_gate:
                cal, rep_pred = stage_a_calibration(x, y, trial_extra, response, reference_scale, settings)
                rep_nlpd_gain = float(current["rep5f_nlpd"] - cal["rep5f_nlpd"])
                rep_wins = int(sum(c < b for c, b in zip(cal["rep5f_repeat_nlpds"], current["rep5f_repeat_nlpds"])))
                metrics = {**screen, **cal}
                prediction = ana_pred
            else:
                rep_nlpd_gain = None
                rep_wins = 0
                metrics = {**screen, "rep5f_nlpd": None, "rep5f_nrmse": None}
                prediction = ana_pred

            gate = {
                "delta_map": delta_map,
                "delta_lml": delta_lml,
                "prior_penalty_increase": penalty_increase,
                "ana_loo_gain": ana_gain,
                "ana_nlpd_gain": ana_nlpd_gain,
                "rep5f_nlpd_gain": rep_nlpd_gain,
                "rep5f_nlpd_repeat_wins": rep_wins,
                "ard_active": ard_active,
                "passes_gate": passes_gate,
            }
            trials[name] = metrics, prediction, gate

        highest_map = max(trials, key=lambda name: trials[name][2]["delta_map"])
        passing = [name for name, (_, _, gate) in trials.items() if gate["passes_gate"]]
        winner = max(
            passing,
            key=lambda name: (
                trials[name][2]["delta_map"],
                trials[name][2]["ana_loo_gain"],
                trials[name][2]["rep5f_nlpd_gain"] or 0.0,
            ),
            default=None,
        )
        best_candidate = highest_map if winner is None else winner
        best, prediction, gate = trials[best_candidate]
        rounds.append({
            "step": step_idx,
            "best_candidate": best_candidate,
            "highest_map_candidate": highest_map,
            "accepted": winner is not None,
            "previous_map_score": current["map_score"],
            "previous_ana_nrmse": current["ana_nrmse"],
            "previous_rep5f_nlpd": current["rep5f_nlpd"],
            "trials": {name: {**metrics, **checks} for name, (metrics, _, checks) in trials.items()},
        })
        if winner is None:
            stop = "insufficient_joint_evidence"
            break
        current, current_pred = best, prediction
        extra = np.column_stack([extra, eligible.pop(winner)])
        selected.append(winner)

    result = {
        "screening": audit, "baseline_validation": baseline, "selected_validation": current,
        "delta_map": float(current["map_score"] - baseline["map_score"]),
        "ana_loo_gain": relative_gain(baseline["ana_nrmse"], current["ana_nrmse"]),
        "ana_nlpd_gain": float(baseline["ana_nlpd"] - current["ana_nlpd"]),
        "rep5f_nlpd_gain": float(baseline["rep5f_nlpd"] - current["rep5f_nlpd"]),
        "rounds": rounds, "selected": selected, "stop_reason": stop,
    }
    return result, base_ana_pred, current_pred, extra


def load_sequential_rows(data_dir):
    designs = gp.read_csv(data_dir / "sequential_designs.csv", DESIGN_FIELDS)
    truth = gp.read_csv(data_dir / "sequential_truth.csv", (*IDENTITY, *gp.RESPONSES))
    if len(designs) != len(truth):
        raise ValueError("sequential designs and truth must have the same length")
    if [int(row["sample_index"]) for row in designs] != list(
        range(gp.INITIAL_COUNT + 1, gp.INITIAL_COUNT + len(designs) + 1)
    ) or any(
        tuple(design[key] for key in IDENTITY) != tuple(observed[key] for key in IDENTITY)
        for design, observed in zip(designs, truth, strict=True)
    ):
        raise ValueError("sequential designs and truth must align by sample index and design ID")
    return designs, truth


def replay(x, y, selections, designs, truth, config, namespace, *, scaling):
    """Freeze the formulas, refit both arms with each newly revealed observation."""
    if __package__:
        from . import sequential_gp_additive as additive
    else:
        import sequential_gp_additive as additive

    records = []
    has_physical_features = any(formulas for formulas, _ in selections)
    blind_x, _ = gp.parse_rows(designs, training=False)
    if np.unique(np.vstack([x, blind_x]), axis=0).shape[0] != len(x) + len(blind_x):
        raise ValueError("revealed prefix and future designs share input coordinates")
    for i, design in enumerate(designs):
        target = blind_x[i:i + 1]
        target_namespace = design_namespace([design], config) if has_physical_features else None
        target_records = []
        for j, (response, (formulas, settings)) in enumerate(zip(gp.RESPONSES, selections, strict=True)):
            fitted, params = additive.fit_response(x, y[:, j], response, **settings)
            mean, _, _ = gp.predict_one(fitted, target)
            base, base_warnings = float(mean[0]), len(params["fit_warnings"])
            augmented, augmented_warnings, fallback = base, base_warnings, ""
            try:
                x_use, point_use, active = additive.model_inputs(
                    x, target, formulas, namespace, target_namespace, scaling=scaling,
                )
            except FORMULA_ERRORS as exc:
                fallback = str(exc)
            else:
                if active:
                    fitted, params = additive.fit_response(x_use, y[:, j], response, **settings)
                    mean, _, _ = gp.predict_one(fitted, point_use)
                    augmented, augmented_warnings = float(mean[0]), len(params["fit_warnings"])
            target_records.append({
                "sample_index": int(design["sample_index"]), "design_id": design["design_id"],
                "response": response, "baseline_mean_m": base, "selected_mean_m": augmented,
                "baseline_warning_count": base_warnings, "selected_warning_count": augmented_warnings,
                "fallback_reason": fallback,
            })
        # Read this response only after both models have predicted it.
        actual = np.asarray([float(truth[i][response]) for response in gp.RESPONSES])
        if not np.isfinite(actual).all():
            raise ValueError(f"non-finite truth for sample {design['sample_index']}")
        for row, value in zip(target_records, actual.tolist(), strict=True):
            records.append({**row, "observed_m": value})
        x = np.vstack([x, target])
        y = np.vstack([y, actual])
        if target_namespace is not None:
            namespace = {name: np.concatenate((values, target_namespace[name])) for name, values in namespace.items()}
    return records


def replay_metrics(records, response, reference_scale):
    rows = [row for row in records if row["response"] == response]
    actual, baseline, selected = np.array([
        [row[field] for row in rows] for field in ("observed_m", "baseline_mean_m", "selected_mean_m")
    ])
    base = score(actual, baseline, reference_scale)
    augmented = score(actual, selected, reference_scale)
    return {
        "count": len(rows), "first_target": rows[0]["sample_index"], "last_target": rows[-1]["sample_index"],
        "baseline": base, "selected": augmented,
        "relative_gain": relative_gain(base["nrmse"], augmented["nrmse"]),
        "targets_with_lower_absolute_error": int(np.sum(np.abs(actual - selected) < np.abs(actual - baseline))),
        "fallback_count": sum(bool(row["fallback_reason"]) for row in rows),
    }


def run(
    data_dir: Path, pool_path: Path, output: Path, *, train_count: int, replay_after: bool,
    restarts: int, seed: int, geometry_template: Path | None, previous_selection: Path | None = None,
) -> None:
    if __package__:
        from . import sequential_gp_additive as additive
    else:
        import sequential_gp_additive as additive

    if train_count < gp.INITIAL_COUNT or restarts < 0 or not 0 <= seed <= 2**32 - len(gp.RESPONSES):
        raise ValueError("need t ≥ 20, nonnegative restarts and a valid 32-bit seed for all responses")
    if output.exists():
        raise ValueError(f"choose a new output directory: {output}")
    initial_path = data_dir / "initial20.csv"
    input_hashes = {path: gp.file_sha256(path) for path in (initial_path, pool_path)}
    previous_features = None
    if previous_selection is not None:
        if train_count != gp.INITIAL_COUNT:
            raise ValueError(f"--previous-selection requires t={gp.INITIAL_COUNT}")
        input_hashes[previous_selection] = gp.file_sha256(previous_selection)
        previous_report = gp.read_json(previous_selection)
        if previous_report.get("train_count") != gp.INITIAL_COUNT:
            raise ValueError(f"--previous-selection must be a report from the initial {gp.INITIAL_COUNT} samples")
        previous_features = selected_formulas(previous_report)
        recorded = previous_report.get("inputs_sha256", {})
        if not isinstance(recorded, dict):
            raise ValueError("previous feature selection has invalid input hashes")
        for path, filename in ((initial_path, "initial20.csv"), (pool_path, "candidate_pool.json")):
            expected = [digest for source, digest in recorded.items()
                        if source.replace("\\", "/").rsplit("/", 1)[-1] == filename]
            if expected and input_hashes[path] not in expected:
                raise ValueError(f"{path} differs from the recorded source of the previous feature selection")
        if geometry_template is None and previous_report.get("geometry_template"):
            geometry_template = Path(previous_report["geometry_template"])
    initial = gp.read_csv(initial_path, (*DESIGN_FIELDS, *gp.RESPONSES))
    _, initial_y = gp.parse_rows(initial, training=True)
    if len(initial) != gp.INITIAL_COUNT:
        raise ValueError("initial20.csv must contain exactly 20 samples")
    designs, truth = [], []
    if train_count > gp.INITIAL_COUNT or replay_after:
        input_hashes.update({data_dir / name: gp.file_sha256(data_dir / name)
                             for name in ("sequential_designs.csv", "sequential_truth.csv")})
        designs, truth = load_sequential_rows(data_dir)
    if train_count > len(initial) + len(designs):
        raise ValueError(f"requested t={train_count} exceeds the available sample count")
    extra_count = train_count - len(initial)
    training_rows = initial + [
        {**design, **{key: observed[key] for key in gp.RESPONSES}}
        for design, observed in zip(designs[:extra_count], truth[:extra_count], strict=True)
    ]
    if len({row["design_id"].strip() for row in [*initial, *designs]}) != len(initial) + len(designs):
        raise ValueError("design IDs must be unique across initial and sequential rows")
    x, y = gp.parse_rows(training_rows, training=True)
    reference_scales = initial_y.std(axis=0, ddof=1)
    if not np.isfinite(reference_scales).all() or np.any(reference_scales <= 0):
        raise ValueError("NRMSE requires positive, finite initial response standard deviations")
    config = primitives.DamGeometryConfig(geometry_template) if geometry_template else primitives.DamGeometryConfig()
    namespace = design_namespace(training_rows, config)
    pool, rejected = load_pool(pool_path, initial, namespace)
    scaling = candidates.ReferenceFeatureScaling(config)
    results, val_rows, selections = {}, [], []
    with threadpool_limits(limits=1):
        for j, response in enumerate(gp.RESPONSES):
            settings = dict(restarts=restarts, seed=seed + j)
            scale = float(reference_scales[j])
            print(f"t={train_count} | {response}: screening and Stage-A 3D validation (MAP + AnaLOO + Rep5x3 NLPD)", flush=True)
            result, baseline, prediction, _ = select_response(
                x, y[:, j], pool[response], namespace, response, scale, settings, scaling=scaling,
            )
            result["screening"] = rejected[response] + result["screening"]
            results[response] = {"initial20_std_m": scale, **result}
            by_name = {feature["name"]: feature for feature in pool[response]}
            selections.append(([by_name[name] for name in result["selected"]], settings))
            for row, actual, base, selected in zip(training_rows, y[:, j], baseline, prediction, strict=True):
                val_rows.append({
                    "sample_index": int(row["sample_index"]), "design_id": row["design_id"],
                    "response": response, "observed_m": float(actual),
                    "baseline_mean_m": float(base), "selected_mean_m": float(selected),
                })
        online = []
        if replay_after and extra_count < len(designs):
            print(f"Predicting samples {train_count + 1}–{len(initial) + len(designs)} in sequence", flush=True)
            online = replay(x, y, selections, designs[extra_count:], truth[extra_count:], config, namespace, scaling=scaling)

    for j, response in enumerate(gp.RESPONSES):
        results[response]["sequential_replay"] = (
            replay_metrics(online, response, float(reference_scales[j])) if online else None
        )
    report = {
        "schema_version": 2, "train_count": train_count, "replay_requested": replay_after,
        "inputs_sha256": {str(path.resolve()): digest for path, digest in input_hashes.items()},
        "geometry_template": str(config.template_path.resolve()),
        "geometry_template_sha256": candidates.geometry_template_sha256(config.template_path),
        "settings": {
            "model": "matched prior 8D baseline GP for empty feature sets; asymmetric additive physical and residual GP otherwise",
            "optimizer_restarts": restarts, "base_seed": seed,
            "original_scaling": "fixed design bounds mapped to [0, 1]",
            "feature_scaling": candidates.reference_scaling_settings(),
            "nrmse_denominator": "sample standard deviation of the first 20 responses",
            "absolute_std_floor": ABS_STD_FLOOR, "relative_std_floor": REL_STD_FLOOR,
            "candidate_redundancy_abs_rho_strictly_greater_than": REDUNDANCY_LIMIT,
            "feature_cap": FEATURE_CAP,
            "minimum_delta_map": COLD_START_MIN_DELTA_MAP,
            "minimum_physical_signal_share": COLD_START_MIN_PHYSICAL_SHARE,
            "max_active_physical_length_scale": additive.MAX_ACTIVE_PHYSICAL_LENGTH_SCALE,
            "max_prior_penalty_increase": COLD_START_MAX_PENALTY_INCREASE,
            "minimum_analytical_loo_gain": COLD_START_MIN_ANA_LOO_GAIN,
            "minimum_analytical_nlpd_gain": COLD_START_MIN_ANA_NLPD_GAIN,
            "selection_rule": "highest delta_map among candidates passing MAP evidence, physical share, penalty stability, and Analytical GP LOO (RMSE & NLPD)",
            "feasible_reference": candidates.FEASIBLE_REFERENCE_SETTINGS,
        },
        "responses": results,
    }
    if previous_features is not None:
        report["migration"] = {
            "previous_selection": str(previous_selection.resolve()),
            "responses": {response: migration_outcomes(previous_features[response], result)
                          for response, result in results.items()},
        }
    for path, digest in input_hashes.items():
        if gp.file_sha256(path) != digest:
            raise ValueError(f"input changed during feature selection: {path}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".feature-selection-", dir=output.parent) as staging:
        ready = Path(staging)
        gp.write_csv(ready / "validation_predictions.csv", list(val_rows[0]), val_rows)
        if online:
            gp.write_csv(ready / "sequential_predictions.csv", list(online[0]), online)
        gp.write_json(ready / "feature_selection.json", report)
        ready.rename(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True,
                        help="prepared experiment directory with initial20.csv and sequential CSVs")
    parser.add_argument("--candidate-pool", type=Path, required=True)
    parser.add_argument("--previous-selection", type=Path,
                        help="previous t=20 feature_selection.json to audit when reselecting from its candidate pool")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--train-count", type=int, default=gp.INITIAL_COUNT)
    parser.add_argument("--replay", action="store_true", help="evaluate on samples t+1 onward in their original order")
    parser.add_argument("--restarts", type=int, default=15)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--geometry-template", type=Path)
    args = parser.parse_args()
    try:
        run(args.data_dir, args.candidate_pool, args.output, train_count=args.train_count,
            replay_after=args.replay, restarts=args.restarts, seed=args.seed,
            geometry_template=args.geometry_template, previous_selection=args.previous_selection)
    except (OSError, TypeError, ValueError) as exc:
        parser.exit(2, f"{exc}\n")
    print(args.output / "feature_selection.json")


if __name__ == "__main__":
    main()
