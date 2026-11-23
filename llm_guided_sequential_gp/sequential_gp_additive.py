"""Fit one physical-plus-residual ARD Matérn GP per constraint response."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from threadpoolctl import threadpool_limits

if __package__:
    from . import sequential_gp_feature_selection as features
    from . import sequential_gp_model as gp
    from . import sequential_gp_primitives as primitives
else:
    import sequential_gp_feature_selection as features
    import sequential_gp_model as gp
    import sequential_gp_primitives as primitives


PHYSICAL_INITIAL_VARIANCE = 0.8
RESIDUAL_INITIAL_VARIANCE = 0.2
SIGNAL_VARIANCE_BOUNDS = (1e-2, 1e2)


class BlockMatern(Matern):
    """Apply a Matérn kernel to one block of the combined input matrix."""

    def __init__(self, length_scale=1.0, length_scale_bounds=(1e-5, 1e5), nu=1.5, columns=()):
        self.columns = columns
        super().__init__(length_scale=length_scale, length_scale_bounds=length_scale_bounds, nu=nu)

    def __call__(self, X, Y=None, eval_gradient=False):
        x_block = np.atleast_2d(X)[:, self.columns]
        y_block = None if Y is None else np.atleast_2d(Y)[:, self.columns]
        return super().__call__(x_block, y_block, eval_gradient=eval_gradient)


def additive_kernel(physical_count: int):
    """Build the covariance for inputs ordered as [original design, physical features]."""
    if physical_count < 1:
        raise ValueError("an additive GP needs at least one physical coordinate")
    design_count = len(gp.VARIABLES)
    config = gp.DEFAULT_PRIOR_CONFIG
    start = config.base_factor * math.sqrt(design_count)
    length_bounds = (config.length_scale_lower, config.length_scale_upper)
    physical = BlockMatern(
        length_scale=np.full(physical_count, start), length_scale_bounds=length_bounds,
        nu=2.5, columns=tuple(range(design_count, design_count + physical_count)),
    )
    residual = BlockMatern(
        length_scale=np.full(design_count, start), length_scale_bounds=length_bounds,
        nu=2.5, columns=tuple(range(design_count)),
    )
    physical_signal = ConstantKernel(PHYSICAL_INITIAL_VARIANCE, SIGNAL_VARIANCE_BOUNDS) * physical
    residual_signal = ConstantKernel(RESIDUAL_INITIAL_VARIANCE, SIGNAL_VARIANCE_BOUNDS) * residual
    noise = WhiteKernel(config.noise_initial, (config.noise_lower, config.noise_upper))
    return (physical_signal + residual_signal) + noise


def fit_response(x: np.ndarray, y: np.ndarray, response: str, *, restarts: int, seed: int):
    """Fit the GP from normalized input coordinates and return its signal parameters."""
    physical_count = x.shape[1] - len(gp.VARIABLES)
    if physical_count == 0:
        fitted = gp.fit_standard(x, y, response, restarts=restarts, seed=seed)
        signal = fitted[0].kernel_.k1
        physical_variance, residual_variance = 0.0, float(signal.k1.constant_value)
        physical_lengths, residual_lengths = [], np.atleast_1d(signal.k2.length_scale).tolist()
    else:
        kernel = additive_kernel(physical_count)
        fitted = gp._fit_gp(x, y, response, kernel, "fmin_l_bfgs_b", restarts=restarts, seed=seed)
        signal = fitted[0].kernel_.k1
        physical_variance = float(signal.k1.k1.constant_value)
        residual_variance = float(signal.k2.k1.constant_value)
        physical_lengths = np.atleast_1d(signal.k1.k2.length_scale).tolist()
        residual_lengths = np.atleast_1d(signal.k2.k2.length_scale).tolist()
    total = physical_variance + residual_variance
    return fitted, {
        "model": "additive" if physical_count else "standard_8d",
        "fitted_signal_variance_standardized": {"physical": physical_variance, "residual": residual_variance},
        "fitted_signal_share": {"physical": physical_variance / total, "residual": residual_variance / total},
        "length_scales": {"physical": physical_lengths, "residual": dict(zip(gp.VARIABLES, residual_lengths, strict=True))},
        "nugget_variance_standardized": float(fitted[0].kernel_.k2.noise_level),
        "fit_warnings": fitted[3],
    }


def selected_formulas(report: dict) -> dict[str, list[dict]]:
    responses = report.get("responses")
    if report.get("schema_version") != 1 or not isinstance(responses, dict) or set(responses) != set(gp.RESPONSES):
        raise ValueError("expected a feature_selection.json with the three constraint responses")
    selected = {}
    for response in gp.RESPONSES:
        names = responses[response]["selected"]
        records = {item["name"]: item for item in responses[response]["screening"]}
        if (not isinstance(names, list)
                or any(not isinstance(name, str) or name not in records for name in names)
                or len(set(names)) != len(names)):
            raise ValueError(f"invalid selected features for {response}")
        selected[response] = [
            {key: records[name][key] for key in features.FORMULA_FIELDS} for name in names
        ]
    return selected


def training_prefix(data_dir: Path, count: int) -> list[dict]:
    initial = gp.read_csv(data_dir / "initial20.csv", (*features.DESIGN_FIELDS, *gp.RESPONSES))
    if len(initial) != gp.INITIAL_COUNT:
        raise ValueError("initial20.csv must contain exactly 20 samples")
    additional = count - gp.INITIAL_COUNT
    if not additional:
        return initial
    designs, truth = features.load_sequential_rows(data_dir)
    if additional > len(designs):
        raise ValueError(f"requested t={count} exceeds the available sample count")
    return initial + [
        {**design, **{key: observed[key] for key in gp.RESPONSES}}
        for design, observed in zip(designs[:additional], truth[:additional], strict=True)
    ]


def physical_values(formulas, namespace, count):
    if not formulas:
        return np.empty((count, 0))
    return np.column_stack([features.feature_values(item, namespace, count) for item in formulas])


def scale_physical(z):
    """Fit the physical-coordinate scaling using training rows only."""
    keep = z.std(axis=0) > features.std_floor(z)
    z = z[:, keep]
    low, span = z.min(axis=0), np.ptp(z, axis=0)
    return keep, low, span


def model_inputs(x, point, formulas, namespace, target_namespace):
    """Scale each varying physical feature using the revealed training prefix."""
    if not formulas:
        return x, point, []
    z = physical_values(formulas, namespace, len(x))
    keep, low, span = scale_physical(z)
    active = [item for item, varying in zip(formulas, keep, strict=True) if varying]
    x = np.column_stack((x, (z[:, keep] - low) / span))
    if point is not None and active:
        target = physical_values(active, target_namespace, len(point))
        point = np.column_stack((point, (target - low) / span))
    return x, point, active


def loo_response(x, y, z, response, reference_scale, *, restarts, seed):
    """Refit scaling and both GP kernels within each leave-one-out fold."""
    prediction = np.empty(len(y))
    warning_count = 0
    indices = np.arange(len(y))
    for held in indices:
        train = indices != held
        x_train, point = x[train], x[held:held + 1]
        if z.shape[1]:
            keep, low, span = scale_physical(z[train])
            x_train = np.column_stack((x_train, (z[train][:, keep] - low) / span))
            point = np.column_stack((point, (z[held:held + 1, keep] - low) / span))
        fitted, params = fit_response(x_train, y[train], response, restarts=restarts, seed=seed)
        prediction[held] = gp.predict_one(fitted, point)[0][0]
        warning_count += len(params["fit_warnings"])
    return {
        **features.score(y, prediction, reference_scale), "folds": len(y),
        "warning_count": warning_count, "predicted_m": prediction.tolist(),
        "residual_m": (y - prediction).tolist(),
    }


def run(
    data_dir: Path, selection_path: Path, output: Path, *, train_count: int | None,
    predict_path: Path | None, geometry_template: Path | None, restarts: int, seed: int,
) -> Path:
    if output.exists():
        raise ValueError(f"choose a new output directory: {output}")
    if restarts < 0 or not 0 <= seed <= 2**32 - len(gp.RESPONSES):
        raise ValueError("need nonnegative restarts and a valid 32-bit seed for all responses")
    selection = gp.read_json(selection_path)
    selected_at = selection.get("train_count")
    if type(selected_at) is not int or selected_at < gp.INITIAL_COUNT:
        raise ValueError("feature_selection.json has an invalid selection sample count")
    count = selected_at if train_count is None else train_count
    if count < selected_at:
        raise ValueError(f"training count must be at least the feature selection count ({selected_at})")
    formulas = selected_formulas(selection)
    train = training_prefix(data_dir, count)
    x, y = gp.parse_rows(train, training=True)
    template = geometry_template or Path(selection["geometry_template"])
    config = primitives.DamGeometryConfig(template)
    namespace = features.design_namespace(train, config) if any(formulas.values()) else None

    point = None
    target = None
    target_namespace = None
    if predict_path is not None:
        designs = gp.read_csv(predict_path, features.DESIGN_FIELDS)
        blind, _ = gp.parse_rows(designs, training=False)
        indices = [i for i, row in enumerate(designs) if int(row["sample_index"]) == count + 1]
        if len(indices) != 1:
            raise ValueError(f"prediction CSV needs exactly one sample_index={count + 1}")
        index = indices[0]
        point = blind[index:index + 1]
        target = designs[index]
        if target["design_id"].strip() in {row["design_id"].strip() for row in train} or np.any(np.all(x == point, axis=1)):
            raise ValueError("the next prediction design already appears in the training prefix")
        if namespace is not None:
            target_namespace = features.design_namespace([target], config)

    result = {}
    predictions = []
    with threadpool_limits(limits=1):
        for j, response in enumerate(gp.RESPONSES):
            try:
                x_use, point_use, active = model_inputs(x, point, formulas[response], namespace, target_namespace)
            except features.FORMULA_ERRORS as exc:
                raise ValueError(f"{response}: {exc}") from exc
            fitted, params = fit_response(x_use, y[:, j], response, restarts=restarts, seed=seed + j)
            result[response] = {
                "selected_features": formulas[response], "active_features": [item["name"] for item in active],
                **params,
            }
            if point_use is not None:
                mean, observed_std, latent_std = gp.predict_one(fitted, point_use)
                predictions.append({
                    "sample_index": count + 1, "design_id": target["design_id"], "response": response,
                    "mean_m": float(mean[0]), "observed_std_m": float(observed_std[0]),
                    "latent_std_m": float(latent_std[0]),
                })

    report = {
        "schema_version": 1, "training_sample_count": count,
        "feature_selection_sample_count": selected_at, "feature_selection": str(selection_path.resolve()),
        "feature_selection_sha256": gp.file_sha256(selection_path), "data_dir": str(data_dir.resolve()),
        "training_rows_sha256": gp.canonical_sha256(train),
        "geometry_template": str(template.resolve()),
        "initial_signal_variance_standardized": {"physical": PHYSICAL_INITIAL_VARIANCE, "residual": RESIDUAL_INITIAL_VARIANCE},
        "signal_variance_bounds": SIGNAL_VARIANCE_BOUNDS, "optimizer_restarts": restarts, "base_seed": seed,
        "responses": result,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".additive-gp-", dir=output.parent) as staging:
        ready = Path(staging)
        gp.write_json(ready / "additive_fit_report.json", report)
        if predictions:
            gp.write_csv(ready / "blind_prediction.csv", list(predictions[0]), predictions)
        ready.rename(output)
    return output / "additive_fit_report.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True, help="directory with the prepared experiment CSVs")
    parser.add_argument("--feature-selection", type=Path, required=True, help="existing feature_selection.json")
    parser.add_argument("--output", type=Path, required=True, help="new output directory")
    parser.add_argument("--train-count", type=int, help="revealed sample count (default: feature selection count)")
    parser.add_argument("--predict", type=Path, help="response-free CSV containing the next design")
    parser.add_argument("--geometry-template", type=Path, help="override the template path saved in feature_selection.json")
    parser.add_argument("--restarts", type=int, default=3)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    try:
        report = run(args.data_dir, args.feature_selection, args.output, train_count=args.train_count,
                     predict_path=args.predict, geometry_template=args.geometry_template,
                     restarts=args.restarts, seed=args.seed)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(2, f"{exc}\n")
    print(report)


if __name__ == "__main__":
    main()
