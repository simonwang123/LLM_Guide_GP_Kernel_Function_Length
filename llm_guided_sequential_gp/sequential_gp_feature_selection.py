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
MIN_GAINS = (0.01, 0.02, 0.02)
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


def load_pool(path, initial_rows, namespace):
    pool = gp.read_json(path)
    initial_ids = [row["design_id"].strip() for row in initial_rows]
    responses = pool.get("responses")
    if (pool.get("schema_version") != 1 or pool.get("design_ids") != initial_ids
            or not isinstance(responses, dict) or set(responses) != set(gp.RESPONSES)):
        raise ValueError("candidate_pool.json must match the initial 20 designs and response schema")
    initial_namespace = {key: values[:len(initial_ids)] for key, values in namespace.items()}
    for response, items in responses.items():
        if not isinstance(items, list) or not candidates.MIN_CANDIDATES <= len(items) <= candidates.MAX_CANDIDATES:
            raise ValueError(f"{response} needs 5–6 validated candidates")
        names = set()
        for feature in items:
            if not isinstance(feature, dict) or set(feature) != {*FORMULA_FIELDS, "unit", "values"}:
                raise ValueError(f"invalid candidate record for {response}")
            measured = feature_values(feature, initial_namespace, len(initial_ids))
            name = feature["name"]
            saved = np.asarray(feature["values"], dtype=float)
            if name in names or saved.shape != measured.shape or not np.allclose(
                measured, saved, atol=1e-12, rtol=1e-10
            ):
                raise ValueError(f"duplicate candidate name or initial geometry mismatch: {response}/{name}")
            names.add(name)
    return responses


def std_floor(values):
    return np.maximum(ABS_STD_FLOOR, REL_STD_FLOOR * np.max(np.abs(values), axis=0))


def screen_candidates(x, y, pool, namespace):
    values = {}
    audit = {}
    for feature in pool:
        name = feature["name"]
        item = audit[name] = {field: feature[field] for field in (*FORMULA_FIELDS, "unit")}
        try:
            column = feature_values(feature, namespace, len(x))
        except FORMULA_ERRORS as exc:
            item.update(status="invalid_formula_values", reason=str(exc))
            continue
        with np.errstate(over="ignore", invalid="ignore"):
            variance = float(np.var(column))
        floor = float(std_floor(column))
        item.update(variance=variance if math.isfinite(variance) else None, std_floor=floor)
        if not math.isfinite(variance):
            item.update(status="invalid_variance")
        elif math.sqrt(variance) <= floor:
            item.update(status="near_constant")
        else:
            values[name] = column
            item.update(status="eligible")

    # The design variables remain in all models. Their correlations are diagnostic.
    eligible = list(values)
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
    return {eligible[i]: values[eligible[i]] for i in retained}, list(audit.values())


def predict_fold(x, y, extra, target_x, target_extra, response, settings):
    if extra.shape[1]:
        lo = extra.min(axis=0)
        span = extra.max(axis=0) - lo
        variable = extra.std(axis=0) > std_floor(extra)
        # A feature may vary in D_t but be near constant in one LOO training fold.
        # A constant coordinate makes no contribution to the GP kernel distance.
        x = np.column_stack([x, (extra[:, variable] - lo[variable]) / span[variable]])
        target_x = np.column_stack([
            target_x, (target_extra[:, variable] - lo[variable]) / span[variable],
        ])
    fitted = gp.fit_standard(x, y, response, **settings)
    mean, _, _ = gp.predict_one(fitted, target_x)
    return float(mean[0]), len(fitted[3])


def score(y, prediction, reference_scale):
    residual = y - prediction
    rmse = float(np.sqrt(np.mean(residual**2)))
    return {"rmse_m": rmse, "nrmse": rmse / reference_scale,
            "mae_m": float(np.mean(np.abs(residual)))}


def relative_gain(base, new):
    return 1.0 - new / base if base > 0 else None


def loo_score(x, y, extra, response, reference_scale, settings):
    predictions = np.empty(len(y))
    warnings = 0
    for held in range(len(y)):
        train = np.arange(len(y)) != held
        predictions[held], count = predict_fold(
            x[train], y[train], extra[train], x[held:held + 1], extra[held:held + 1],
            response, settings,
        )
        warnings += count
    return {**score(y, predictions, reference_scale), "warning_count": warnings}, predictions


def select_response(x, y, pool, namespace, response, reference_scale, settings):
    eligible, audit = screen_candidates(x, y, pool, namespace)
    extra = np.empty((len(x), 0))
    baseline, baseline_pred = loo_score(x, y, extra, response, reference_scale, settings)
    current, current_pred = baseline, baseline_pred
    selected, rounds = [], []
    stop = "feature_limit"
    for threshold in MIN_GAINS:
        if not eligible or current["nrmse"] == 0:
            stop = "no_more_candidates" if not eligible else "zero_current_error"
            break
        trials = {}
        for name, column in eligible.items():
            trial_extra = np.column_stack([extra, column])
            trials[name] = loo_score(x, y, trial_extra, response, reference_scale, settings)
        winner = min(trials, key=lambda name: trials[name][0]["nrmse"])
        best, prediction = trials[winner]
        previous = current["nrmse"]
        accepted = best["nrmse"] <= previous * (1 - threshold)
        rounds.append({
            "step": len(selected) + 1, "threshold": threshold, "best_candidate": winner,
            "accepted": accepted, "previous_loo_nrmse": previous,
            "trials": {name: {**metrics, "relative_gain": relative_gain(previous, metrics["nrmse"])}
                       for name, (metrics, _) in trials.items()},
        })
        if not accepted:
            stop = "insufficient_gain"
            break
        current, current_pred = best, prediction
        extra = np.column_stack([extra, eligible.pop(winner)])
        selected.append(winner)

    result = {
        "screening": audit, "baseline_loo": baseline, "selected_loo": current,
        "relative_loo_gain": relative_gain(baseline["nrmse"], current["nrmse"]),
        "rounds": rounds, "selected": selected, "stop_reason": stop,
    }
    return result, baseline_pred, current_pred, extra


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


def replay(x, y, selections, designs, truth, config):
    """Freeze the formulas, refit both arms with each newly revealed observation."""
    records = []
    extra_history = [extra for _, extra, _ in selections]
    blind_x, _ = gp.parse_rows(designs, training=False)
    if np.unique(np.vstack([x, blind_x]), axis=0).shape[0] != len(x) + len(blind_x):
        raise ValueError("revealed prefix and future designs share input coordinates")
    for i, design in enumerate(designs):
        target = blind_x[i:i + 1]
        namespace = design_namespace([design], config) if any(items for items, _, _ in selections) else None
        target_records = []
        for j, response in enumerate(gp.RESPONSES):
            formulas, _, settings = selections[j]
            empty = np.empty((len(x), 0))
            base, base_warnings = predict_fold(x, y[:, j], empty, target, empty[:1], response, settings)
            augmented, augmented_warnings, fallback = base, base_warnings, ""
            if formulas:
                try:
                    values = np.array([feature_values(f, namespace, 1)[0] for f in formulas])
                except FORMULA_ERRORS as exc:
                    values = np.full(len(formulas), np.nan)
                    fallback = str(exc)
                extra_history[j] = np.vstack([extra_history[j], values])
                if not fallback and not np.isfinite(extra_history[j]).all():
                    fallback = "a selected formula is undefined for an earlier revealed design"
                if not fallback:
                    augmented, augmented_warnings = predict_fold(
                        x, y[:, j], extra_history[j][:-1], target, values[None, :], response, settings,
                    )
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
    restarts: int, seed: int, geometry_template: Path | None,
) -> None:
    if train_count < gp.INITIAL_COUNT or restarts < 0 or not 0 <= seed <= 2**32 - len(gp.RESPONSES):
        raise ValueError("need t ≥ 20, nonnegative restarts and a valid 32-bit seed for all responses")
    if output.exists():
        raise ValueError(f"choose a new output directory: {output}")
    initial_path = data_dir / "initial20.csv"
    initial = gp.read_csv(initial_path, (*DESIGN_FIELDS, *gp.RESPONSES))
    _, initial_y = gp.parse_rows(initial, training=True)
    if len(initial) != gp.INITIAL_COUNT:
        raise ValueError("initial20.csv must contain exactly 20 samples")
    designs, truth = load_sequential_rows(data_dir) if train_count > gp.INITIAL_COUNT or replay_after else ([], [])
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
    pool = load_pool(pool_path, initial, namespace)
    results, loo_rows, selections = {}, [], []
    with threadpool_limits(limits=1):
        for j, response in enumerate(gp.RESPONSES):
            settings = dict(restarts=restarts, seed=seed + j)
            scale = float(reference_scales[j])
            print(f"t={train_count} | {response}: screening and refitted LOO", flush=True)
            result, baseline, prediction, extra = select_response(
                x, y[:, j], pool[response], namespace, response, scale, settings,
            )
            results[response] = {"initial20_std_m": scale, **result}
            by_name = {feature["name"]: feature for feature in pool[response]}
            selections.append(([by_name[name] for name in result["selected"]], extra, settings))
            for row, actual, base, selected in zip(training_rows, y[:, j], baseline, prediction, strict=True):
                loo_rows.append({
                    "sample_index": int(row["sample_index"]), "design_id": row["design_id"],
                    "response": response, "observed_m": float(actual),
                    "baseline_mean_m": float(base), "selected_mean_m": float(selected),
                })
        online = []
        if replay_after and extra_count < len(designs):
            print(f"Predicting samples {train_count + 1}–{len(initial) + len(designs)} in sequence", flush=True)
            online = replay(x, y, selections, designs[extra_count:], truth[extra_count:], config)

    for j, response in enumerate(gp.RESPONSES):
        results[response]["sequential_replay"] = (
            replay_metrics(online, response, float(reference_scales[j])) if online else None
        )
    source_paths = [initial_path, pool_path]
    if designs:
        source_paths += [data_dir / "sequential_designs.csv", data_dir / "sequential_truth.csv"]
    report = {
        "schema_version": 1, "train_count": train_count, "replay_requested": replay_after,
        "inputs_sha256": {str(path.resolve()): gp.file_sha256(path) for path in source_paths},
        "geometry_template": str(config.template_path.resolve()),
        "settings": {
            "model": "no-prior Constant * ARD Matérn 5/2 + WhiteKernel",
            "optimizer_restarts": restarts, "base_seed": seed,
            "original_scaling": "fixed design bounds mapped to [0, 1]",
            "feature_scaling": "min-max fitted within each training fold, without clipping",
            "nrmse_denominator": "sample standard deviation of the first 20 responses",
            "absolute_std_floor": ABS_STD_FLOOR, "relative_std_floor": REL_STD_FLOOR,
            "candidate_redundancy_abs_rho_strictly_greater_than": REDUNDANCY_LIMIT,
            "minimum_relative_gains": MIN_GAINS,
        },
        "responses": results,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".feature-selection-", dir=output.parent) as staging:
        ready = Path(staging)
        gp.write_csv(ready / "loo_predictions.csv", list(loo_rows[0]), loo_rows)
        if online:
            gp.write_csv(ready / "sequential_predictions.csv", list(online[0]), online)
        gp.write_json(ready / "feature_selection.json", report)
        ready.rename(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True,
                        help="prepared experiment directory with initial20.csv and sequential CSVs")
    parser.add_argument("--candidate-pool", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--train-count", type=int, default=gp.INITIAL_COUNT)
    parser.add_argument("--replay", action="store_true", help="evaluate on samples t+1 onward in their original order")
    parser.add_argument("--restarts", type=int, default=3)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--geometry-template", type=Path)
    args = parser.parse_args()
    try:
        run(args.data_dir, args.candidate_pool, args.output, train_count=args.train_count,
            replay_after=args.replay, restarts=args.restarts, seed=args.seed,
            geometry_template=args.geometry_template)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"{exc}\n")
    print(args.output / "feature_selection.json")


if __name__ == "__main__":
    main()
