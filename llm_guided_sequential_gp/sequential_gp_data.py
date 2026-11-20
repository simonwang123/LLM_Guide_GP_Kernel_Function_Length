"""Prepare a leakage-aware 20+N dataset for the sequential GP experiment.

Run without arguments to prepare the existing historical 20+40 experiment, or
pass --sequential-count to select a shorter prefix (for example, 10 points).
The source run must contain at least 20+N completed rows. This script does not
fit a GP or call an LLM.

Outputs intentionally separate future designs from their finite-element truth:
initial20.csv, sequential_designs.csv, sequential_truth.csv,
initial20_context.json, and data_summary.json.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import rankdata


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN = ROOT / "runs" / "semiparametric_guarded_20_20261109"
INITIAL_COUNT = 20
SEQUENTIAL_COUNT = 40
TEST_COUNT = 40

# Kept here so this experimental script runs directly without importing another
# experimental script or changing geostudio_automation. These values match the
# current DesignVariables/knowledge calibration schema in the project.
DESIGN_BOUNDS = {
    "upstream_slope": (1.3, 1.7),
    "downstream_slope": (1.3, 1.7),
    "crest_width": (6.0, 12.0),
    "modulus_zone_height": (5.0, 20.0),
    "secondary_top_width": (10.0, 30.0),
    "secondary_height": (20.0, 40.0),
    "secondary_bottom_width": (25.0, 70.0),
    "secondary_right_slope_distance": (2.0, 15.0),
}
DESIGN_VARS = tuple(DESIGN_BOUNDS)
RESPONSES = (
    "crest_max_settlement_m",
    "dam_max_settlement_m",
    "dam_max_horizontal_displacement_m",
)
PRIOR_SIGN = (
    (-1, -1, -1),
    (-1, -1, -1),
    (+1, +1, +1),
    (-1, -1, +1),
    (+1, +1, +1),
    (+1, +1, +1),
    (+1, +1, +1),
    (+1, +1, -1),
)
PRIOR_TIER = (
    (0.5, 0.1, 0.1),
    (1.0, 1.0, 1.0),
    (0.1, 0.1, 0.1),
    (1.0, 0.5, 0.5),
    (1.0, 1.0, 0.5),
    (0.1, 1.0, 1.0),
    (1.0, 1.0, 0.1),
    (0.5, 1.0, 0.5),
)


def read_csv(path: Path, required: tuple[str, ...]) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"missing input file: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        headers = reader.fieldnames
        if headers is None or len(headers) != len(set(headers)):
            raise ValueError(f"missing or duplicate CSV headers: {path}")
        missing = set(required) - set(headers)
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    if not rows or any(None in row for row in rows):
        raise ValueError(f"empty CSV or row with extra fields: {path}")
    return rows


def parse_designs(rows: list[dict[str, str]], label: str) -> tuple[np.ndarray, np.ndarray]:
    ids = [row["design_id"].strip() for row in rows]
    if any(not value for value in ids) or len(ids) != len(set(ids)):
        raise ValueError(f"{label} has empty or duplicate design_id values")
    try:
        x = np.array([[float(row[name]) for name in DESIGN_VARS] for row in rows])
        y = np.array([[float(row[name]) for name in RESPONSES] for row in rows])
    except (ValueError, TypeError) as error:
        raise ValueError(f"{label} has a non-numeric design or response") from error
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError(f"{label} has missing or non-finite design/response values")
    bounds = np.array(list(DESIGN_BOUNDS.values()))
    if np.any(x < bounds[:, 0] - 1e-12) or np.any(x > bounds[:, 1] + 1e-12):
        raise ValueError(f"{label} has a design outside declared physical bounds")
    if np.unique(x, axis=0).shape[0] != len(x):
        raise ValueError(f"{label} has duplicate design vectors")
    normalized = (x - bounds[:, 0]) / (bounds[:, 1] - bounds[:, 0])
    return normalized, y


def ordered_history(
    rows: list[dict[str, str]], history: list[dict[str, str]],
    sequential_count: int = SEQUENTIAL_COUNT,
) -> list[dict[str, str]]:
    if type(sequential_count) is not int or sequential_count < 0:
        raise ValueError("sequential_count must be a nonnegative integer")
    total = INITIAL_COUNT + sequential_count
    source_count = len(history)
    if len(rows) != source_count or source_count < total:
        raise ValueError(
            f"source needs at least {total} matching dataset/history rows; "
            f"got {len(rows)} and {source_count}"
        )
    by_id = {row["design_id"].strip(): row for row in rows}
    if len(by_id) != source_count:
        raise ValueError("dataset contains duplicate design_id values")
    history_ids = [row["design_id"].strip() for row in history]
    if len(set(history_ids)) != source_count or set(history_ids) != set(by_id):
        raise ValueError("dataset and history have different or duplicate design IDs")
    if any(row["status"].strip() != "succeeded" for row in history[:total]):
        raise ValueError("all selected history rows must be successful finite-element evaluations")
    attempts = [int(row["attempt"]) for row in history]
    if attempts != sorted(attempts) or len(set(attempts)) != source_count:
        raise ValueError("optimization_history.csv must be in unique attempt order")
    if any(row["phase"].strip() not in {"initial", "imported_initial"} for row in history[:INITIAL_COUNT]):
        raise ValueError("the first 20 history rows must be initial designs")
    if any(row["phase"].strip() in {"initial", "imported_initial"} for row in history[INITIAL_COUNT:]):
        raise ValueError("history has an initial design after the first 20 rows")
    ordered = [by_id[design_id] for design_id in history_ids[:total]]
    for data_row, history_row in zip(ordered, history[:total], strict=True):
        for name in RESPONSES:
            if not math.isclose(
                float(data_row[name]), float(history_row[name]), rel_tol=0, abs_tol=1e-10
            ):
                raise ValueError(f"dataset/history response mismatch: {data_row['design_id']}/{name}")
    return ordered


def screening_statistics(x: np.ndarray, y: np.ndarray) -> dict:
    """Two noncausal screens for one revealed training prefix."""
    if (x.ndim != 2 or y.ndim != 2 or x.shape[1] != len(DESIGN_VARS)
            or y.shape != (len(x), len(RESPONSES)) or len(x) < INITIAL_COUNT
            or not np.isfinite(x).all() or not np.isfinite(y).all()):
        raise ValueError("screening statistics require at least 20 finite 8D designs and 3 responses")
    n = len(x)
    ranks = rankdata(np.column_stack((x, y)), axis=0)
    partial = np.empty((8, 3))
    quad_gain = np.empty((8, 3))
    for i in range(8):
        other = [j for j in range(8) if j != i]
        controls = np.column_stack((np.ones(n), ranks[:, other]))
        x_residual = ranks[:, i] - controls @ np.linalg.lstsq(
            controls, ranks[:, i], rcond=None
        )[0]
        for m in range(3):
            response_rank = ranks[:, 8 + m]
            y_residual = response_rank - controls @ np.linalg.lstsq(
                controls, response_rank, rcond=None
            )[0]
            partial[i, m] = (
                np.corrcoef(x_residual, y_residual)[0, 1]
                if np.std(x_residual) > 1e-10 and np.std(y_residual) > 1e-10
                else np.nan
            )

            # Compare one-variable linear and quadratic fits by leave-one-out
            # prediction. This is a screening statistic, not a main-effect claim.
            predictions = np.empty((2, n))
            for held in range(n):
                train = np.arange(n) != held
                for degree in (1, 2):
                    features = np.vander(x[train, i], degree + 1, increasing=True)
                    coefficients = np.linalg.lstsq(features, y[train, m], rcond=None)[0]
                    predictions[degree - 1, held] = np.polynomial.polynomial.polyval(
                        x[held, i], coefficients
                    )
            denominator = np.sum((y[:, m] - np.mean(y[:, m])) ** 2)
            quad_gain[i, m] = (
                (np.sum((y[:, m] - predictions[0]) ** 2)
                 - np.sum((y[:, m] - predictions[1]) ** 2)) / denominator
                if denominator > 0 else np.nan
            )
    return {
        "partial_spearman_design_response": finite_or_none(partial),
        "delta_loo_r2_quad_univariate": finite_or_none(quad_gain),
        "interpretation": (
            f"Statistics use {n} revealed points. Partial rank correlation controls "
            "the other seven input ranks linearly. Quadratic gain compares "
            "univariate leave-one-out predictive R²; neither identifies a causal "
            "effect in coupled geometry. Null means the statistic is undefined."
        ),
    }


def finite_or_none(values: np.ndarray) -> list[list[float | None]]:
    return [[float(v) if np.isfinite(v) else None for v in row] for row in values]


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(run: Path, output: Path, test_dataset: Path | None,
            sequential_count: int = SEQUENTIAL_COUNT) -> dict:
    dataset_path = run / "dataset.csv"
    history_path = run / "optimization_history.csv"
    dataset = read_csv(dataset_path, ("design_id", *DESIGN_VARS, *RESPONSES))
    history = read_csv(history_path, ("attempt", "design_id", "phase", "status", *RESPONSES))
    ordered = ordered_history(dataset, history, sequential_count)
    x, y = parse_designs(ordered, "historical dataset")
    strata = np.minimum(np.floor(np.clip(x[:INITIAL_COUNT], 0, 1) * INITIAL_COUNT).astype(int), INITIAL_COUNT - 1)
    lhs_counts = [len(np.unique(strata[:, i])) for i in range(8)]
    if any(count != INITIAL_COUNT for count in lhs_counts):
        raise ValueError(f"initial 20 designs do not occupy all LHS strata: {lhs_counts}")

    metadata = {}
    for name in ("analysis_name", "step", "cvar_alpha"):
        if name in ordered[0]:
            values = {row[name] for row in ordered}
            if len(values) != 1:
                raise ValueError(f"historical dataset mixes {name} values")
            metadata[name] = values.pop()

    test_rows = None
    test_x = None
    if test_dataset is not None:
        test_rows = read_csv(test_dataset, ("design_id", *DESIGN_VARS, *RESPONSES))
        if len(test_rows) != TEST_COUNT:
            raise ValueError(f"independent test dataset must contain {TEST_COUNT} rows")
        test_x, _ = parse_designs(test_rows, "independent test dataset")
        if {row["design_id"].strip() for row in ordered} & {
            row["design_id"].strip() for row in test_rows
        }:
            raise ValueError("historical and test datasets share design IDs")
        all_x = np.vstack((x, test_x))
        if np.unique(all_x, axis=0).shape[0] != len(all_x):
            raise ValueError("historical and test datasets share design vectors")
        for name, expected in metadata.items():
            if name in test_rows[0] and any(row[name] != expected for row in test_rows):
                raise ValueError(f"test dataset has inconsistent {name}")

    # The prompt context has no future responses, holdout metrics, or full-run hash.
    initial_payload = [
        {name: row[name] for name in ("design_id", *DESIGN_VARS, *RESPONSES)}
        for row in ordered[:INITIAL_COUNT]
    ]
    initial_hash = hashlib.sha256(
        json.dumps(initial_payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    context = {
        "schema_version": 2,
        "initial_count": INITIAL_COUNT,
        "initial_observations_sha256": initial_hash,
        "design_variables": list(DESIGN_VARS),
        "design_bounds_physical": DESIGN_BOUNDS,
        "input_scale": "fixed physical bounds mapped to [0, 1]",
        "responses": list(RESPONSES),
        "prior_sign": PRIOR_SIGN,
        "prior_tier": PRIOR_TIER,
        "prior_tier_meaning": "static knowledge strength (0.1, 0.5, 1.0), not GP length-scale tiers",
        "statistics": screening_statistics(x[:INITIAL_COUNT], y[:INITIAL_COUNT]),
    }

    output.mkdir(parents=True, exist_ok=True)
    normalized_fields = [f"{name}_unit" for name in DESIGN_VARS]
    initial_fields = ["sample_index", "design_id", *DESIGN_VARS, *normalized_fields, *RESPONSES]
    design_fields = ["sample_index", "design_id", *DESIGN_VARS, *normalized_fields]
    truth_fields = ["sample_index", "design_id", *RESPONSES]
    initial_output, design_output, truth_output = [], [], []
    for index, row in enumerate(ordered):
        design = {"sample_index": index + 1, "design_id": row["design_id"]}
        design.update({name: row[name] for name in DESIGN_VARS})
        design.update(dict(zip(normalized_fields, x[index].tolist(), strict=True)))
        truth = {"sample_index": index + 1, "design_id": row["design_id"]}
        truth.update({name: row[name] for name in RESPONSES})
        if index < INITIAL_COUNT:
            initial_output.append({**design, **{name: row[name] for name in RESPONSES}})
        else:
            design_output.append(design)
            truth_output.append(truth)

    write_csv(output / "initial20.csv", initial_fields, initial_output)
    write_csv(output / "sequential_designs.csv", design_fields, design_output)
    write_csv(output / "sequential_truth.csv", truth_fields, truth_output)
    if test_rows is not None:
        holdout_dir = output / "holdout"
        holdout_dir.mkdir(exist_ok=True)
        holdout_output = []
        for index, row in enumerate(test_rows):
            item = {"test_index": index + 1, "design_id": row["design_id"]}
            item.update({name: row[name] for name in DESIGN_VARS})
            item.update(dict(zip(normalized_fields, test_x[index].tolist(), strict=True)))
            item.update({name: row[name] for name in RESPONSES})
            holdout_output.append(item)
        write_csv(
            holdout_dir / "test_holdout.csv",
            ["test_index", "design_id", *DESIGN_VARS, *normalized_fields, *RESPONSES],
            holdout_output,
        )
    (output / "initial20_context.json").write_text(
        json.dumps(context, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8"
    )
    summary = {
        "schema_version": 1,
        "source_run": str(run),
        "dataset": str(dataset_path),
        "dataset_sha256": sha256(dataset_path),
        "history": str(history_path),
        "history_sha256": sha256(history_path),
        "sample_counts": {"initial": INITIAL_COUNT, "sequential": sequential_count,
                          "independent_test": len(test_rows) if test_rows else 0},
        "source_sample_count": len(history),
        "initial_lhs_occupied_strata_per_variable": dict(zip(DESIGN_VARS, lhs_counts, strict=True)),
        "response_unit": "m",
        "result_selection": metadata,
        "test_dataset": str(test_dataset) if test_dataset else None,
        "test_dataset_sha256": sha256(test_dataset) if test_dataset else None,
        "outputs": {name: str(output / name) for name in (
            "initial20.csv", "sequential_designs.csv", "sequential_truth.csv", "initial20_context.json"
        )},
        "causal_boundary": "Only initial20_context.json may inform cold-start LLM decisions. "
                           "Sequential truth must be revealed one row at a time; "
                           "the independent test set must never enter LLM or model-selection gates.",
    }
    if test_rows is not None:
        summary["outputs"]["test_holdout.csv"] = str(output / "holdout" / "test_holdout.csv")
    (output / "data_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8"
    )
    print(f"Input dataset: {dataset_path}")
    print(f"Input order:   {history_path}")
    if test_dataset is not None:
        print(f"Input test:    {test_dataset}")
    print(f"Selected: {INITIAL_COUNT} initial + {sequential_count} sequential + "
          f"{summary['sample_counts']['independent_test']} independent test samples")
    print(f"Initial LHS occupied strata: {lhs_counts}")
    print(f"First initial sample: {initial_output[0]['design_id']}, "
          f"x_unit={[round(float(v), 4) for v in x[0]]}, "
          f"y_m={[round(float(v), 6) for v in y[0]]}")
    if design_output:
        print(f"Next sequential design: {design_output[0]['design_id']}, "
              f"x_unit={[round(float(v), 4) for v in x[INITIAL_COUNT]]}")
    print("Future truth is stored separately in sequential_truth.csv; "
          "the initial LLM context contains only samples 1-20.")
    print(f"Output directory: {output}")
    for path in summary["outputs"].values():
        print(f"  {path}")
    print(f"  {output / 'data_summary.json'}")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=DEFAULT_RUN,
                        help="directory containing dataset.csv and optimization_history.csv")
    parser.add_argument("--test-dataset", type=Path,
                        help="optional, separate 40-row solved finite-element test dataset")
    parser.add_argument("--sequential-count", type=int, default=SEQUENTIAL_COUNT,
                        help="number of sequential rows to select after the initial 20 (default: 40)")
    parser.add_argument("--output", type=Path,
                        help="output directory; defaults to a separate summary_taylor_v4 directory")
    args = parser.parse_args()
    run = args.run.resolve()
    default_name = f"{run.name}_20_plus_{args.sequential_count}_summary_taylor_v4"
    output = (args.output or Path(__file__).resolve().parent / "outputs" / default_name).resolve()
    try:
        prepare(run, output, args.test_dataset.resolve() if args.test_dataset else None,
                args.sequential_count)
    except ValueError as error:
        parser.exit(2, f"Data validation failed: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
