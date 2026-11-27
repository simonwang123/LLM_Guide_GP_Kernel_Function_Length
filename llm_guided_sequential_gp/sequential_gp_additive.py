"""Fit one physical-plus-residual ARD Matérn GP per constraint response."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from threadpoolctl import threadpool_limits
from scipy.optimize import minimize
from sklearn.exceptions import ConvergenceWarning
import warnings

if __package__:
    from . import sequential_gp_feature_candidates as candidates
    from . import sequential_gp_feature_selection as features
    from . import sequential_gp_model as gp
    from . import sequential_gp_primitives as primitives
    from . import sequential_gp_validation as validation
else:
    import sequential_gp_feature_candidates as candidates
    import sequential_gp_feature_selection as features
    import sequential_gp_model as gp
    import sequential_gp_primitives as primitives
    import sequential_gp_validation as validation


# === 1. 拆分物理/残差方差与长度尺度安全边界（允许物理核自由退化，残差核采用非对称防坍缩先验） ===
PHYSICAL_INITIAL_VARIANCE = 0.5
RESIDUAL_INITIAL_VARIANCE = 0.5
PHYSICAL_VARIANCE_BOUNDS = (1e-4, 100.0)   # 允许无效物理特征自然坍缩至 ~0%，交由门控淘汰
RESIDUAL_VARIANCE_BOUNDS = (1e-3, 100.0)   # 与标准 GP 同量级上限，避免高方差工况触顶
SIGNAL_VARIANCE_BOUNDS = PHYSICAL_VARIANCE_BOUNDS  # 保留导出名兼容报告字段
MIN_PHYSICAL_SIGNAL_SHARE = 0.10
MAX_ACTIVE_PHYSICAL_LENGTH_SCALE = 8.0

PHYSICAL_LENGTH_SCALE_LOWER = 0.15
RESIDUAL_LENGTH_SCALE_LOWER = 0.15         # 放宽残差核下界至 0.15，兼顾灵活性与数值稳定性
PHYSICAL_PRIOR_LOG_STD = math.log(3.0)     # 物理核：宽松双边正则，保留灵活性
RESIDUAL_PRIOR_LOG_STD = math.log(2.0)     # 残差核短尺度防坍缩先验（< center 时生效）
RESIDUAL_PRIOR_UPPER_LOG_STD = math.log(4.0)  # 残差核长尺度宽松先验（>= center 时允许 ARD 屏蔽无关维度）
OPTIMIZER_METHOD = "L-BFGS-B"

class BlockMatern(Matern):
    """Apply a Matérn kernel to one block of the combined input matrix."""

    def __init__(self, length_scale=1.0, length_scale_bounds=(1e-5, 1e5), nu=1.5, columns=()):
        self.columns = columns
        super().__init__(length_scale=length_scale, length_scale_bounds=length_scale_bounds, nu=nu)

    def __call__(self, X, Y=None, eval_gradient=False):
        x_block = np.atleast_2d(X)[:, self.columns]
        y_block = None if Y is None else np.atleast_2d(Y)[:, self.columns]
        return super().__call__(x_block, y_block, eval_gradient=eval_gradient)


def baseline_8d_kernel():
    """Build the matched 8D baseline kernel sharing residual bounds and prior."""
    design_count = len(gp.VARIABLES)
    config = gp.DEFAULT_PRIOR_CONFIG
    start = config.base_factor * math.sqrt(design_count)
    signal = ConstantKernel(1.0, (1e-2, 1e2)) * Matern(
        length_scale=np.full(design_count, start),
        length_scale_bounds=(RESIDUAL_LENGTH_SCALE_LOWER, config.length_scale_upper),
        nu=2.5,
    )
    noise = WhiteKernel(config.noise_initial, (config.noise_lower, config.noise_upper))
    return signal + noise


def _baseline_lengthscale_indices(kernel, design_count: int) -> np.ndarray:
    """Locate residual length-scale indices in the matched 8D baseline kernel."""
    offset = 0
    for param in kernel.hyperparameters:
        if param.fixed:
            continue
        if param.name == "k1__k2__length_scale":
            if param.n_elements != design_count:
                raise ValueError("无法匹配 8D 基线核的长度尺度超参数索引")
            return np.arange(offset, offset + design_count, dtype=int)
        offset += param.n_elements
    raise ValueError("8D 基线核缺少可训练的长度尺度超参数")


def additive_kernel(physical_count: int):
    """Build the covariance for inputs ordered as [original design, physical features]."""
    if physical_count < 1:
        raise ValueError("an additive GP needs at least one physical coordinate")
    design_count = len(gp.VARIABLES)
    config = gp.DEFAULT_PRIOR_CONFIG
    start = config.base_factor * math.sqrt(design_count)
    physical = BlockMatern(
        length_scale=np.full(physical_count, start),
        length_scale_bounds=(PHYSICAL_LENGTH_SCALE_LOWER, config.length_scale_upper),
        nu=2.5,
        columns=tuple(range(design_count, design_count + physical_count)),
    )
    residual = BlockMatern(
        length_scale=np.full(design_count, start),
        length_scale_bounds=(RESIDUAL_LENGTH_SCALE_LOWER, config.length_scale_upper),
        nu=2.5,
        columns=tuple(range(design_count)),
    )
    physical_signal = ConstantKernel(PHYSICAL_INITIAL_VARIANCE, PHYSICAL_VARIANCE_BOUNDS) * physical
    residual_signal = ConstantKernel(RESIDUAL_INITIAL_VARIANCE, RESIDUAL_VARIANCE_BOUNDS) * residual
    noise = WhiteKernel(config.noise_initial, (config.noise_lower, config.noise_upper))
    return (physical_signal + residual_signal) + noise

def _additive_parameter_indices(kernel, physical_count: int, design_count: int):
    """Locate the length scales and signal variances in sklearn's log parameters."""
    indices = {}
    offset = 0
    for param in kernel.hyperparameters:
        if param.fixed:
            continue
        indices[param.name] = np.arange(offset, offset + param.n_elements, dtype=int)
        offset += param.n_elements
    physical_lengths = indices.get("k1__k1__k2__length_scale", ())
    residual_lengths = indices.get("k1__k2__k2__length_scale", ())
    physical_variance = indices.get("k1__k1__k1__constant_value", ())
    residual_variance = indices.get("k1__k2__k1__constant_value", ())
    if (len(physical_lengths) != physical_count or len(residual_lengths) != design_count
            or len(physical_variance) != 1 or len(residual_variance) != 1):
        raise ValueError("无法匹配加性核的长度尺度与信号方差超参数索引")
    return physical_lengths, residual_lengths, (int(physical_variance[0]), int(residual_variance[0]))

class AdditiveLengthscaleMAPOptimizer:
    """Fit GP length scales with asymmetric anti-collapse residual regularization."""

    def __init__(
        self,
        phys_idx: np.ndarray,
        res_idx: np.ndarray,
        center: float,
        *,
        variance_indices: tuple[int, int] | None = None,
        phys_log_std: float = PHYSICAL_PRIOR_LOG_STD,
        res_log_std: float = RESIDUAL_PRIOR_LOG_STD,
        res_upper_log_std: float = RESIDUAL_PRIOR_UPPER_LOG_STD,
    ):
        self.phys_idx = np.asarray(phys_idx, dtype=int)
        self.res_idx = np.asarray(res_idx, dtype=int)
        self.variance_indices = variance_indices
        self.center = float(center)
        self.phys_log_std = float(phys_log_std)
        self.res_log_std = float(res_log_std)
        self.res_upper_log_std = float(res_upper_log_std)

    def __call__(self, obj_func, initial_theta: np.ndarray, bounds: np.ndarray) -> tuple[np.ndarray, float]:
        def objective(theta: np.ndarray) -> tuple[float, np.ndarray]:
            value, grad = obj_func(theta, eval_gradient=True)
            derivative = np.asarray(grad, dtype=float).copy()
            pen = 0.0
            if len(self.res_idx):
                # 8D 残差核非对称先验：短尺度严格防坍缩，长尺度宽松放行 ARD 屏蔽无关维度
                res_dev = theta[self.res_idx] - self.center
                res_std = np.where(res_dev < 0.0, self.res_log_std, self.res_upper_log_std)
                inv_var = np.where(np.isfinite(res_std), 1.0 / (res_std ** 2), 0.0)
                pen += 0.5 * float(np.sum((res_dev ** 2) * inv_var))
                derivative[self.res_idx] += res_dev * inv_var
            if len(self.phys_idx):
                # 物理主核：宽松双边正则，保留物理特征自适应不同尺度的灵活性
                phys_dev = theta[self.phys_idx] - self.center
                pen += 0.5 * float(phys_dev @ phys_dev) / (self.phys_log_std ** 2)
                derivative[self.phys_idx] += phys_dev / (self.phys_log_std ** 2)
            return float(value + pen), derivative

        result = minimize(objective, initial_theta, method=OPTIMIZER_METHOD, jac=True, bounds=bounds)
        if not np.isfinite(result.fun) or not np.isfinite(result.x).all():
            raise ValueError("Additive MAP optimizer produced a non-finite result")
        if not result.success:
            warnings.warn(f"Additive MAP optimizer did not converge: {result.message}", ConvergenceWarning, stacklevel=2)
        return np.asarray(result.x, dtype=float), float(result.fun)


def fit_response(x: np.ndarray, y: np.ndarray, response: str, *, restarts: int, seed: int):
    """Fit the GP from normalized input coordinates and return its signal parameters."""
    design_count = len(gp.VARIABLES)
    physical_count = x.shape[1] - design_count
    center = math.log(gp.DEFAULT_PRIOR_CONFIG.base_factor * math.sqrt(design_count))
    if physical_count == 0:
        kernel = baseline_8d_kernel()
        res_idx = _baseline_lengthscale_indices(kernel, design_count)
        optimizer = AdditiveLengthscaleMAPOptimizer(
            np.empty(0, dtype=int), res_idx, center,
            res_log_std=RESIDUAL_PRIOR_LOG_STD,
            res_upper_log_std=RESIDUAL_PRIOR_UPPER_LOG_STD,
        )
        fitted = gp._fit_gp(x, y, response, kernel, optimizer, restarts=restarts, seed=seed)
        signal = fitted[0].kernel_.k1
        physical_variance, residual_variance = 0.0, float(signal.k1.constant_value)
        physical_lengths, residual_lengths = [], np.atleast_1d(signal.k2.length_scale).tolist()
    else:
        kernel = additive_kernel(physical_count)
        phys_idx, res_idx, variance_indices = _additive_parameter_indices(kernel, physical_count, design_count)
        optimizer = AdditiveLengthscaleMAPOptimizer(
            phys_idx, res_idx, center, variance_indices=variance_indices,
            phys_log_std=PHYSICAL_PRIOR_LOG_STD,
            res_log_std=RESIDUAL_PRIOR_LOG_STD,
            res_upper_log_std=RESIDUAL_PRIOR_UPPER_LOG_STD,
        )
        fitted = gp._fit_gp(x, y, response, kernel, optimizer, restarts=restarts, seed=seed)
        signal = fitted[0].kernel_.k1
        physical_variance = float(signal.k1.k1.constant_value)
        residual_variance = float(signal.k2.k1.constant_value)
        physical_lengths = np.atleast_1d(signal.k1.k2.length_scale).tolist()
        residual_lengths = np.atleast_1d(signal.k2.k2.length_scale).tolist()
    total = physical_variance + residual_variance
    physical_share = physical_variance / total if total > 0 else 0.0
    residual_share = residual_variance / total if total > 0 else 0.0
    return fitted, {
        "model": "additive" if physical_count else "standard_8d",
        "fitted_signal_variance_standardized": {"physical": physical_variance, "residual": residual_variance},
        "fitted_signal_share": {"physical": physical_share, "residual": residual_share},
        "length_scales": {"physical": physical_lengths, "residual": dict(zip(gp.VARIABLES, residual_lengths, strict=True))},
        "nugget_variance_standardized": float(fitted[0].kernel_.k2.noise_level),
        "fit_warnings": fitted[3],
    }


selected_formulas = features.selected_formulas


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


def scaled_physical_values(formulas, namespace, count, *, scaling):
    raw = physical_values(formulas, namespace, count)
    for index, formula in enumerate(formulas):
        raw[:, index] = scaling.transform(formula, raw[:, index])
    return raw


def model_inputs(x, point, formulas, namespace, target_namespace, *, scaling):
    """Append fixed reference coordinates that vary on the revealed prefix."""
    if not formulas:
        return x, point, []
    z = scaled_physical_values(formulas, namespace, len(x), scaling=scaling)
    keep = z.std(axis=0) > features.std_floor(z)
    active = [item for item, varying in zip(formulas, keep, strict=True) if varying]
    x = np.column_stack((x, z[:, keep]))
    if point is not None and active:
        target = scaled_physical_values(active, target_namespace, len(point), scaling=scaling)
        point = np.column_stack((point, target))
    return x, point, active


def summarize_fit_warnings(messages: list[str]) -> dict[str, int]:
    """Separate optimizer failures from kernel-bound and other fit warnings."""
    counts = {"optimizer_failure_count": 0, "bound_warning_count": 0, "other_warning_count": 0}
    for message in messages:
        text = message.strip().lower()
        if text.startswith(("additive map optimizer did not converge:", "lbfgs failed to converge",
                            "l-bfgs-b failed to converge")):
            counts["optimizer_failure_count"] += 1
        elif text.startswith("the optimal value found for dimension") and (
            "close to the specified lower bound" in text or "close to the specified upper bound" in text
        ):
            counts["bound_warning_count"] += 1
        else:
            counts["other_warning_count"] += 1
    return {"warning_count": len(messages), **counts}


def _stack_inputs(x: np.ndarray, scaled_z: np.ndarray, train_idx: np.ndarray, test_idx: np.ndarray | None = None):
    """Append physical columns that vary on train_idx without leaking test_idx statistics."""
    x_train = x[train_idx]
    x_test = None if test_idx is None else x[test_idx]
    if scaled_z.shape[1]:
        train_z = scaled_z[train_idx]
        keep = train_z.std(axis=0) > features.std_floor(train_z)
        x_train = np.column_stack((x_train, train_z[:, keep]))
        if x_test is not None:
            x_test = np.column_stack((x_test, scaled_z[test_idx][:, keep]))
    return x_train, x_test


def fit_predict_subset(
    x: np.ndarray,
    scaled_z: np.ndarray,
    y: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    response: str,
    *,
    restarts: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Fit additive/baseline GP on train_idx and predict mean and observed std on test_idx."""
    x_train, x_test = _stack_inputs(x, scaled_z, train_idx, test_idx)
    fitted, params = fit_response(x_train, y[train_idx], response, restarts=restarts, seed=seed)
    mu, obs_std, _ = gp.predict_one(fitted, x_test)
    return mu, np.maximum(obs_std, validation.STD_EPS), params["fit_warnings"]


def fit_full_evidence(
    x: np.ndarray,
    y: np.ndarray,
    scaled_z: np.ndarray,
    response: str,
    reference_scale: float,
    *,
    restarts: int,
    seed: int,
) -> dict:
    """Fit once on the full revealed prefix to extract MAP/LML, share, and analytical GP LOO."""
    full_idx = np.arange(len(y), dtype=int)
    x_full, _ = _stack_inputs(x, scaled_z, full_idx)
    fitted, params = fit_response(x_full, y, response, restarts=restarts, seed=seed)
    evidence = validation.full_fit_evidence(
        fitted, y, reference_scale,
        physical_share=float(params["fitted_signal_share"]["physical"]),
    )
    return {
        **evidence,
        "physical_signal_share_mean": evidence["physical_signal_share"],
        "length_scales": params["length_scales"],
        **summarize_fit_warnings(evidence["fit_warnings"]),
    }


def evaluate_stage_a_calibration(
    x: np.ndarray,
    y: np.ndarray,
    scaled_z: np.ndarray,
    response: str,
    reference_scale: float,
    *,
    restarts: int,
    seed: int,
    n_splits: int = validation.DEFAULT_REP_SPLITS,
    n_repeats: int = validation.DEFAULT_REP_REPEATS,
) -> dict:
    """Run Stage-A refitted Repeated Stratified K-Fold CV as an NLPD calibration guard."""
    rep = validation.repeated_stratified_cv_evaluate(
        lambda tr, te, fold_seed: fit_predict_subset(
            x, scaled_z, y, tr, te, response, restarts=restarts, seed=seed + fold_seed,
        ),
        x, y, reference_scale, n_splits=n_splits, n_repeats=n_repeats, seed=seed,
    )
    return {
        **rep,
        **summarize_fit_warnings(rep["fit_warnings"]),
    }


def evaluate_stage_b(
    x: np.ndarray,
    y: np.ndarray,
    scaled_z: np.ndarray,
    response: str,
    reference_scale: float,
    *,
    restarts: int,
    seed: int,
    cv_restarts: int | None = None,
    n_clusters: int = validation.DEFAULT_CLUSTER_COUNT,
    preq_max_steps: int = validation.DEFAULT_PREQ_MAX_STEPS,
) -> dict:
    """Evaluate all four Stage-B dimensions: MAP/LML, AnaLOO, Cluster-CV, and Causal Prequential."""
    evidence = fit_full_evidence(x, y, scaled_z, response, reference_scale, restarts=restarts, seed=seed)
    sub_restarts = max(1, min(restarts, 4)) if cv_restarts is None else cv_restarts

    cluster = validation.cluster_cv_evaluate(
        lambda tr, te, fold_idx: fit_predict_subset(
            x, scaled_z, y, tr, te, response, restarts=sub_restarts, seed=seed + fold_idx,
        ),
        x, y, reference_scale, n_clusters=n_clusters,
    )
    preq_restarts = max(1, min(sub_restarts, 3))
    prequential = validation.prequential_evaluate(
        lambda tr, te, fold_idx: fit_predict_subset(
            x, scaled_z, y, tr, te, response, restarts=preq_restarts, seed=seed + fold_idx,
        ),
        x, y, reference_scale, start_index=gp.INITIAL_COUNT, max_steps=preq_max_steps,
    )
    all_warnings = [*evidence["fit_warnings"], *cluster["fit_warnings"]]
    if prequential is not None:
        all_warnings.extend(prequential["fit_warnings"])

    return {
        "raw_lml": evidence["raw_lml"],
        "map_score": evidence["map_score"],
        "physical_signal_share_mean": evidence["physical_signal_share_mean"],
        "length_scales": evidence["length_scales"],
        "ana_loo": evidence["ana_loo"],
        "cluster": {k: v for k, v in cluster.items() if k != "fit_warnings"},
        "prequential": None if prequential is None else {k: v for k, v in prequential.items() if k != "fit_warnings"},
        # Expose out-of-cluster metrics as top-level spatial summary for diagnosis compatibility.
        "rmse_m": cluster["rmse_m"],
        "nrmse": cluster["nrmse"],
        "mae_m": cluster["mae_m"],
        "nlpd": cluster["nlpd"],
        "predicted_m": cluster["predicted_m"],
        "predicted_std_m": cluster["predicted_std_m"],
        "residual_m": cluster["residual_m"],
        "folds": cluster["folds"],
        **summarize_fit_warnings(all_warnings),
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
    scaling = candidates.ReferenceFeatureScaling(config)
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
                x_use, point_use, active = model_inputs(
                    x, point, formulas[response], namespace, target_namespace, scaling=scaling,
                )
            except features.FORMULA_ERRORS as exc:
                raise ValueError(f"{response}: {exc}") from exc
            fitted, params = fit_response(x_use, y[:, j], response, restarts=restarts, seed=seed + j)
            result[response] = {
                "selected_features": formulas[response], "active_features": [item["name"] for item in active],
                "physical_feature_bounds": {item["name"]: scaling.bounds(item) for item in active},
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
        "geometry_template": str(config.template_path.resolve()),
        "geometry_template_sha256": candidates.geometry_template_sha256(config.template_path),
        "feature_scaling": candidates.reference_scaling_settings(),
        "feasible_reference": candidates.FEASIBLE_REFERENCE_SETTINGS,
        "initial_signal_variance_standardized": {"physical": PHYSICAL_INITIAL_VARIANCE, "residual": RESIDUAL_INITIAL_VARIANCE},
        "min_physical_signal_share": MIN_PHYSICAL_SIGNAL_SHARE,
        "max_active_physical_length_scale": MAX_ACTIVE_PHYSICAL_LENGTH_SCALE,
        "residual_signal_variance_bounds": RESIDUAL_VARIANCE_BOUNDS,
        "signal_variance_bounds": SIGNAL_VARIANCE_BOUNDS,
        "residual_length_scale_lower": RESIDUAL_LENGTH_SCALE_LOWER,
        "residual_length_scale_prior_log_std": {
            "lower": RESIDUAL_PRIOR_LOG_STD,
            "upper": RESIDUAL_PRIOR_UPPER_LOG_STD,
        },
        "optimizer_restarts": restarts, "base_seed": seed,
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
    parser.add_argument("--restarts", type=int, default=15)
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
