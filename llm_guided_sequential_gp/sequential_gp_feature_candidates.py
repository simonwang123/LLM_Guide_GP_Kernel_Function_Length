"""独立的物理特征冷启动，支持渲染、LLM 调用和离线回复校验。"""

from __future__ import annotations

import argparse
import ast
import json
import math
import os
import re
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import yeojohnson, yeojohnson_normmax

if __package__:
    from . import sequential_gp_data as data_io
    from . import sequential_gp_llm as llm_io
    from . import sequential_gp_primitives as primitives
    from . import sequential_gp_prompt as prompt_io
else:
    import sequential_gp_data as data_io
    import sequential_gp_llm as llm_io
    import sequential_gp_primitives as primitives
    import sequential_gp_prompt as prompt_io


TEMPLATE = Path(__file__).resolve().parent / "prompts" / "gp_feature_primitives.md"
PLACEHOLDER = "{primitive_summary_table}"
UNSAFE_TOKEN = "{unsafe_denominators}"
MIN_CANDIDATES, MAX_CANDIDATES = 5, 6
PROMPT_RESPONSES = dict(zip(
    ("crest_max_settlement", "body_max_settlement", "body_max_horizontal_displacement"),
    data_io.RESPONSES,
    strict=True,
))
DESIGN_VARIABLES = primitives.CANONICAL_DESIGN_VARS
ALIAS_TO_CANONICAL = {f"x{i}": name for i, name in enumerate(DESIGN_VARIABLES, 1)}
ALIAS_TO_CANONICAL["d_ds_perp"] = "secondary_right_slope_distance"
UNIT_POWERS = {"1": 0, "m": 1, "m^2": 2}
SYMBOL_UNITS = dict(zip(DESIGN_VARIABLES, (0, 0, 1, 1, 1, 1, 1, 1), strict=True))
SYMBOL_UNITS.update({name: UNIT_POWERS[meta["unit"]] for name, meta in primitives.PRIMITIVE_METADATA.items()})
# 局部小量在倒数中必须有独立的加法下界，例如总断面比例或分区之和。
UNSAFE_STANDALONE_DENOMINATORS = frozenset({
    "X_c_sec", "d_us_min", "d_ds_perp",
    "A_sec", "A_mod", "A_ds_shell", "A_drain",
    "h_drain", "W_sec_mean", "crest_width",
    "x3", "x4", "x5", "x6", "x7", "x8",
    "modulus_zone_height", "secondary_top_width", "secondary_height",
    "secondary_bottom_width", "secondary_right_slope_distance",
})

FEASIBLE_REFERENCE_COUNT = 2000
FEASIBLE_REFERENCE_SEED = 2026
FEASIBLE_REFERENCE_MAX_DRAWS = 20000
FEASIBLE_MAX_SPAN_RATIO = 3.5
FEASIBLE_MAX_ABS_SKEW = 4.0
FEASIBLE_REFERENCE_SETTINGS = {
    "sample_count": FEASIBLE_REFERENCE_COUNT, "seed": FEASIBLE_REFERENCE_SEED,
    "max_draws": FEASIBLE_REFERENCE_MAX_DRAWS, "max_span_ratio": FEASIBLE_MAX_SPAN_RATIO,
    "max_abs_skew": FEASIBLE_MAX_ABS_SKEW,
}
FEASIBLE_GRID_CACHE: dict[tuple[Path, str, bool, int, int, int], dict[str, np.ndarray]] = {}
REFERENCE_ABS_SPAN_FLOOR = 1e-8
REFERENCE_REL_SPAN_FLOOR = 1e-8
REFERENCE_TAIL_MARGIN = 0.25
REFERENCE_SKEW_THRESHOLD = 0.50
REFERENCE_YEOJOHNSON_BOUNDS = (-2.0, 2.0)
REFERENCE_TARGET_MEAN = 0.5
REFERENCE_TARGET_STD = 1.0 / math.sqrt(12.0)
REFERENCE_TAIL_QUANTILES = (0.01, 0.99)


def reference_scaling_settings() -> dict:
    return {
        "method": "feasible_reference_yeojohnson_isovariance_with_rational_tails",
        "skew_threshold": REFERENCE_SKEW_THRESHOLD,
        "yeojohnson_lambda_bounds": list(REFERENCE_YEOJOHNSON_BOUNDS),
        "target_mean": REFERENCE_TARGET_MEAN,
        "target_std": REFERENCE_TARGET_STD,
        "tail_quantiles": list(REFERENCE_TAIL_QUANTILES),
        "tail_margin": REFERENCE_TAIL_MARGIN,
        "absolute_span_floor": REFERENCE_ABS_SPAN_FLOOR,
        "relative_span_floor": REFERENCE_REL_SPAN_FLOOR,
    }


@dataclass(frozen=True)
class _ReferenceScaleParams:
    lower: float
    upper: float
    pre_center: float
    pre_scale: float
    lmbda: float | None
    post_center: float
    post_std: float
    tail_low: float
    tail_high: float
    tail_margin: float


def _fit_yeojohnson_lambda(values: np.ndarray, ref_skew: float) -> float | None:
    """Fit a bounded Yeo-Johnson exponent only when reference skewness exceeds threshold."""
    if not math.isfinite(ref_skew) or abs(ref_skew) <= REFERENCE_SKEW_THRESHOLD:
        return None
    try:
        with np.errstate(all="ignore"):
            raw_lmbda = float(yeojohnson_normmax(values, brack=REFERENCE_YEOJOHNSON_BOUNDS))
    except (ValueError, RuntimeError, FloatingPointError):
        return None
    if not math.isfinite(raw_lmbda):
        return None
    lo_lam, hi_lam = REFERENCE_YEOJOHNSON_BOUNDS
    return float(np.clip(raw_lmbda, lo_lam, hi_lam))


def _apply_yeojohnson(values: np.ndarray, lmbda: float | None) -> np.ndarray:
    if lmbda is None:
        return values
    with np.errstate(over="ignore", invalid="ignore"):
        return np.asarray(yeojohnson(values, lmbda), dtype=float)


class ReferenceFeatureScaling:
    """Keep formula de-skewing and iso-variance scaling fixed on the same geometry reference."""

    def __init__(self, geometry: primitives.DamGeometryConfig):
        self.geometry = geometry
        self._namespace = None
        self._params: dict[str, _ReferenceScaleParams] = {}

    @property
    def namespace(self) -> dict[str, np.ndarray]:
        if self._namespace is None:
            self._namespace = build_feasible_reference_namespace(self.geometry)
        return self._namespace

    def _fit_params(self, formula: dict) -> _ReferenceScaleParams:
        expression = formula["expression"].strip()
        if expression in self._params:
            return self._params[expression]

        namespace = self.namespace
        count = len(next(iter(namespace.values())))
        item = {key: formula[key] for key in ("name", "expression", "physical_rationale")}
        try:
            _, values = evaluate_candidate(item, namespace, count, min_std=0)
        except ValueError as exc:
            raise ValueError(f"formula is invalid on the feasible reference: {exc}") from exc

        arr = np.asarray(values, dtype=float)
        low, high = float(np.min(arr)), float(np.max(arr))
        raw_span = high - low
        span_floor = max(REFERENCE_ABS_SPAN_FLOOR, REFERENCE_REL_SPAN_FLOOR * max(abs(low), abs(high)))
        if not math.isfinite(raw_span) or raw_span <= span_floor:
            raise ValueError("formula has a zero, near-zero or non-finite span on the feasible reference")

        # 1. Pre-normalize by robust quantiles so Yeo-Johnson operates near O(1) regardless of units.
        q01, q25, q50, q75, q99 = (float(v) for v in np.quantile(arr, [0.01, 0.25, 0.50, 0.75, 0.99]))
        pre_scale = next((span for span in (q75 - q25, q99 - q01, raw_span) if span > span_floor), 0.0)
        if not math.isfinite(pre_scale) or pre_scale <= span_floor:
            raise ValueError("formula has a zero, near-zero or non-finite span on the feasible reference")
        pre_norm = (arr - q50) / pre_scale

        # 2. Conditionally de-skew using bounded Yeo-Johnson when reference skew is substantial.
        ref_std = float(np.std(arr))
        ref_skew = float(np.mean(((arr - float(np.mean(arr))) / ref_std) ** 3)) if ref_std > span_floor else 0.0
        lmbda = _fit_yeojohnson_lambda(pre_norm, ref_skew)
        warped = _apply_yeojohnson(pre_norm, lmbda)
        if not np.isfinite(warped).all():
            lmbda = None
            warped = pre_norm

        # 3. Rescale to match the variance of Uniform[0, 1] design coordinates.
        post_center = float(np.mean(warped))
        post_std = float(np.std(warped))
        post_floor = max(REFERENCE_ABS_SPAN_FLOOR, REFERENCE_REL_SPAN_FLOOR * max(1.0, abs(post_center)))
        if not math.isfinite(post_std) or post_std <= post_floor:
            raise ValueError("formula has a zero, near-zero or non-finite span after de-skewing")

        z_ref = REFERENCE_TARGET_MEAN + ((warped - post_center) / post_std) * REFERENCE_TARGET_STD
        q_lo_p, q_hi_p = REFERENCE_TAIL_QUANTILES
        tail_low, tail_high = (float(v) for v in np.quantile(z_ref, [q_lo_p, q_hi_p]))
        tail_span = tail_high - tail_low
        if not math.isfinite(tail_span) or tail_span <= REFERENCE_ABS_SPAN_FLOOR:
            tail_low, tail_high = float(np.min(z_ref)), float(np.max(z_ref))
            tail_span = tail_high - tail_low
        if not math.isfinite(tail_span) or tail_span <= REFERENCE_ABS_SPAN_FLOOR:
            raise ValueError("formula has a zero, near-zero or non-finite quantile span on the feasible reference")

        params = _ReferenceScaleParams(
            lower=low,
            upper=high,
            pre_center=q50,
            pre_scale=pre_scale,
            lmbda=lmbda,
            post_center=post_center,
            post_std=post_std,
            tail_low=tail_low,
            tail_high=tail_high,
            tail_margin=REFERENCE_TAIL_MARGIN * tail_span,
        )
        self._params[expression] = params
        return params

    def bounds(self, formula: dict) -> dict[str, float]:
        params = self._fit_params(formula)
        return {"lower": params.lower, "upper": params.upper}

    def transform(self, formula: dict, values: np.ndarray) -> np.ndarray:
        params = self._fit_params(formula)
        arr = np.asarray(values, dtype=float)
        if arr.size == 0:
            return arr.copy()
        if not np.isfinite(arr).all():
            raise ValueError("formula values must be finite before reference scaling")

        pre_norm = (arr - params.pre_center) / params.pre_scale
        warped = _apply_yeojohnson(pre_norm, params.lmbda)
        with np.errstate(over="ignore", invalid="ignore"):
            scaled = REFERENCE_TARGET_MEAN + ((warped - params.post_center) / params.post_std) * REFERENCE_TARGET_STD
            low = scaled < params.tail_low
            high = scaled > params.tail_high
            margin = params.tail_margin
            # Keep unit slope at the quantile endpoints and approach finite asymptotes
            # monotonically, so distant outliers remain distinguishable without exploding.
            scaled[low] = params.tail_low - margin + margin / (1.0 + (params.tail_low - scaled[low]) / margin)
            scaled[high] = params.tail_high + margin - margin / (1.0 + (scaled[high] - params.tail_high) / margin)
        if not np.isfinite(scaled).all():
            raise ValueError("formula produced non-finite values after reference scaling")
        return scaled


def geometry_template_sha256(path: Path) -> str:
    """Hash template bytes and relative filenames, including reference parameters."""
    path = path.resolve()
    if path.is_dir():
        files = {item.relative_to(path).as_posix(): data_io.sha256(item)
                 for item in path.rglob("*") if item.is_file()}
        snapshot = {"directory": files}
    else:
        snapshot = {"archive": data_io.sha256(path)}
        parameters = path.with_name("parameters.json")
        if parameters.is_file():
            snapshot["parameters.json"] = data_io.sha256(parameters)
    return prompt_io.digest(prompt_io.canonical_json(snapshot))

def build_feasible_reference_namespace(
    config: primitives.DamGeometryConfig,
    *,
    n_target: int = FEASIBLE_REFERENCE_COUNT,
    seed: int = FEASIBLE_REFERENCE_SEED,
    max_draws: int = FEASIBLE_REFERENCE_MAX_DRAWS,
) -> dict[str, np.ndarray]:
    """在设计边界内采样几何可行断面并缓存命名空间（仅纯几何计算，零有限元开销）。"""
    if not 1 <= n_target <= max_draws:
        raise ValueError("参考样本数必须为正且不能超过抽样上限")
    template = Path(config.template_path).resolve()
    key = (template, geometry_template_sha256(template), config.include_x_aliases, n_target, seed, max_draws)
    if key in FEASIBLE_GRID_CACHE:
        return FEASIBLE_GRID_CACHE[key]

    rng = np.random.default_rng(seed)
    bounds = data_io.DESIGN_BOUNDS
    lows = np.array([bounds[name][0] for name in DESIGN_VARIABLES], dtype=float)
    highs = np.array([bounds[name][1] for name in DESIGN_VARIABLES], dtype=float)

    draws = rng.uniform(lows, highs, size=(max_draws, len(DESIGN_VARIABLES)))
    ref_ns = primitives.build_feasible_formula_namespace(draws, config, n_target)
    FEASIBLE_GRID_CACHE[key] = ref_ns
    return ref_ns


def check_feasible_domain_sanity(
    formula: object,
    train_values: np.ndarray,
    ref_namespace: dict[str, np.ndarray] | None,
    *,
    max_span_ratio: float = FEASIBLE_MAX_SPAN_RATIO,
    max_abs_skew: float = FEASIBLE_MAX_ABS_SKEW,
) -> None:
    """全可行域外推安全阀：拦截长尾/极端偏态/越界倍增的特征公式。"""
    if ref_namespace is None:
        return
    n_ref = len(next(iter(ref_namespace.values())))
    try:
        _, ref_vals = evaluate_candidate(formula, ref_namespace, n_ref, min_std=0.0)
    except ValueError as exc:
        raise ValueError(f"候选特征在可行域参考断面上求值失败: {exc}") from exc

    train_span = float(np.ptp(train_values))
    ref_span = float(np.ptp(ref_vals))
    if train_span > 0 and ref_span > max_span_ratio * train_span:
        raise ValueError(
            f"候选特征全可行域跨度 ({ref_span:.3g}) 超过训练集跨度 ({train_span:.3g}) 的 {max_span_ratio} 倍，存在外推爆雷风险"
        )

    std = float(np.std(ref_vals))
    if std > 0:
        z = (ref_vals - float(np.mean(ref_vals))) / std
        skew = float(np.mean(z ** 3))
        if abs(skew) > max_abs_skew:
            raise ValueError(f"候选特征全可行域偏度过高 (skew={skew:.2f})，建议改用占总和比例或 log1p 变换")

def render(initial_csv: Path, output: Path, geometry_template: Path | None = None) -> dict:
    initial_hash = data_io.sha256(initial_csv)
    rows = data_io.read_csv(
        initial_csv, ("sample_index", "design_id", *DESIGN_VARIABLES, *data_io.RESPONSES)
    )
    indices = list(range(1, data_io.INITIAL_COUNT + 1))
    if len(rows) != len(indices) or [int(row["sample_index"]) for row in rows] != indices:
        raise ValueError("需要按样本编号 1–20 排列的恰好 20 行初始观测")
    _, y = data_io.parse_designs(rows, "initial20.csv")
    design_ids = [row["design_id"].strip() for row in rows]
    x = pd.DataFrame([{name: float(row[name]) for name in DESIGN_VARIABLES} for row in rows])
    x["design_id"] = design_ids
    config = primitives.DamGeometryConfig(geometry_template) if geometry_template else primitives.DamGeometryConfig()
    namespace = primitives.build_formula_namespace(x, config=config)
    namespace = {
        name: values for name, values in namespace.items()
        if name in SYMBOL_UNITS or name in ALIAS_TO_CANONICAL
    }

    summary = [
        "统计仅来自初始样本 1–20（样本量较小，Spearman 相关系数仅作方向参考，切勿机械拟合偶然高相关）。",
        "",
        "| 响应 | 均值 (mm) | 样本标准差 (mm) | 观测范围 (mm) |",
        "|---|---:|---:|---:|",
    ]
    for key, values in zip(PROMPT_RESPONSES, 1000 * y.T, strict=True):
        summary.append(f"| `{key}` | {values.mean():.1f} | {values.std(ddof=1):.1f} | [{values.min():.1f}, {values.max():.1f}] |")
    summary.extend([
        "",
        "| 原语 | 单位 | 观测范围（量级参考） | 定义 | ρ(C1) | ρ(C2) | ρ(C3) |",
        "|---|---|---:|---|---:|---:|---:|",
    ])
    response_ranks = pd.DataFrame(y).rank()
    for name, meta in primitives.PRIMITIVE_METADATA.items():
        values = namespace[name]
        ranks = pd.Series(values).rank()
        correlations = response_ranks.corrwith(ranks) if ranks.nunique() > 1 else [math.nan] * 3
        rho = " | ".join(f"{v:+.2f}" if math.isfinite(v) else "N/A" for v in correlations)
        lo, hi = float(values.min()), float(values.max())
        bounds_str = f"[{lo:.0f}, {hi:.0f}]" if meta["unit"] == "m^2" else f"[{lo:.1f}, {hi:.1f}]"
        summary.append(f"| `{name}` | {meta['unit']} | {bounds_str} | {meta['desc']} | {rho} |")

    summary.extend([
        "",
        "### 常见无量纲构型方向参考（未按相关性排序，仅供物理机理联想）",
        "| 机理方向 | 典型表达式示例 | 物理意义 |",
        "|---|---|---|",
        "| 分区面积占比 | `A_sec / A_tot`、`A_mod / A_tot`、`A_ds_shell / A_tot` | 软弱压缩区、顶帽加强区与下游保护壳相对规模 |",
        "| 刚柔比值（防长尾求和分母） | `A_sec / (A_sec + A_drain)`、`A_mod / (A_mod + A_main)` | 软弱区与高刚度支承骨架相对份额 |",
        "| 软区空间位置与相对净距 | `Y_c_sec / H_dam`、`X_c_sec / B_bot`、`d_us_min / B_bot`、`d_ds_perp / B_bot` | 软弱区高程、偏心及上下游约束净距 |",
        "| 坡比不对称与细长比 | `(x2 - x1) / (x1 + x2)`、`W_sec_mean / B_bot`、`h_drain / H_dam` | 上下游边坡非对称性与局部几何尺度比 |",
    ])

    template = TEMPLATE.read_text(encoding="utf-8")
    if template.count(PLACEHOLDER) != 1:
        raise ValueError("物理特征模板应包含一个原语统计占位符")
    rendered = template.replace(PLACEHOLDER, "\n".join(summary))
    if UNSAFE_TOKEN in rendered:
        rendered = rendered.replace(UNSAFE_TOKEN, ", ".join(sorted(UNSAFE_STANDALONE_DENOMINATORS)))
    prompt = rendered.encode("utf-8")
    # 保存本地求值所需的几何快照。API 请求只使用渲染后的统计摘要。
    context = {"design_ids": design_ids, "geometry_template": str(config.template_path.resolve()),
               "namespace": {name: values.tolist() for name, values in namespace.items()}}
    if data_io.sha256(initial_csv) != initial_hash:
        raise ValueError("初始 CSV 在渲染期间发生变化")
    manifest = {
        "schema_version": 1,
        "stage": "feature_initial",
        "sample_count": len(indices),
        "initial_csv_sha256": initial_hash,
        "template_sha256": prompt_io.digest(template.encode("utf-8")),
        "context_sha256": prompt_io.digest(prompt_io.canonical_json(context)),
        "prompt_sha256": prompt_io.digest(prompt),
    }
    output.mkdir(parents=True, exist_ok=False)
    llm_io.write_json(output / "context.json", context)
    (output / "prompt.md").write_bytes(prompt)
    llm_io.write_json(output / "request_manifest.json", manifest)
    return manifest


def load_request(request: Path) -> tuple[dict, dict, bytes]:
    manifest = prompt_io.exact_keys(prompt_io.read_json_file(request / "request_manifest.json"), {
        "schema_version", "stage", "sample_count", "initial_csv_sha256", "template_sha256",
        "context_sha256", "prompt_sha256",
    }, "物理特征请求清单")
    context = prompt_io.read_json_file(request / "context.json")
    fields = {"design_ids", "namespace"}
    if isinstance(context, dict) and "geometry_template" in context:
        fields.add("geometry_template")
        prompt_io.nonempty_string(context["geometry_template"], "geometry_template")
    context = prompt_io.exact_keys(context, fields, "物理特征上下文")
    prompt = (request / "prompt.md").read_bytes()
    if (manifest.get("schema_version") != 1 or manifest.get("stage") != "feature_initial"
            or manifest.get("sample_count") != data_io.INITIAL_COUNT
            or manifest.get("context_sha256") != prompt_io.digest(prompt_io.canonical_json(context))
            or manifest.get("prompt_sha256") != prompt_io.digest(prompt)):
        raise ValueError("物理特征请求清单与上下文或提示词不一致")
    return manifest, context, prompt


@dataclass(frozen=True)
class FormulaShape:
    unit_power: Fraction
    sources: frozenset[str] = frozenset()
    nonlinear: bool = False
    normalized_primitive: bool = False


def inspect_formula(node: ast.expr) -> FormulaShape:
    """核对原提示词白名单，并计算单位 m 的指数与变量依赖。"""
    if isinstance(node, ast.Name):
        symbol = ALIAS_TO_CANONICAL.get(node.id, node.id)
        if symbol not in SYMBOL_UNITS:
            raise ValueError(f"未知公式变量 {node.id}")
        sources = frozenset() if symbol == "H_dam" else frozenset({symbol})
        return FormulaShape(Fraction(SYMBOL_UNITS[symbol]), sources)
    if isinstance(node, ast.Constant) and type(node.value) in (int, float) and math.isfinite(node.value):
        return FormulaShape(Fraction(0))
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        return inspect_formula(node.operand)
    if isinstance(node, ast.BinOp):
        a, b = inspect_formula(node.left), inspect_formula(node.right)
        sources = a.sources | b.sources
        nonlinear = a.nonlinear or b.nonlinear
        if isinstance(node.op, (ast.Add, ast.Sub)):
            if a.unit_power != b.unit_power:
                raise ValueError("加减项量纲须一致")
            return FormulaShape(a.unit_power, sources, nonlinear)
        if isinstance(node.op, ast.Mult):
            return FormulaShape(a.unit_power + b.unit_power, sources, nonlinear or bool(a.sources and b.sources))
        if isinstance(node.op, ast.Div):
            power = a.unit_power - b.unit_power
            normalized = (power == 0 and isinstance(node.left, ast.Name)
                          and node.left.id in primitives.PRIMITIVE_METADATA
                          and len(a.sources) == 1 and not b.sources)
            return FormulaShape(power, sources, nonlinear or bool(b.sources), normalized)
        if isinstance(node.op, ast.Pow):
            exponent = ast.literal_eval(node.right)
            if type(exponent) not in (int, float) or abs(exponent) > 3:
                raise ValueError("幂指数须为绝对值不大于 3 的数值常量")
            return FormulaShape(a.unit_power * Fraction(str(exponent)), sources, nonlinear or exponent not in (0, 1))
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and not node.keywords:
        function = node.func.id
        parts = [inspect_formula(arg) for arg in node.args]
        if function in {"sqrt", "abs", "log", "log1p"} and len(parts) == 1:
            a = parts[0]
            if function in {"log", "log1p"} and a.unit_power != 0:
                raise ValueError("对数函数需要无量纲输入")
            power = a.unit_power / 2 if function == "sqrt" else a.unit_power
            return FormulaShape(power, a.sources, True)
        if function in {"min", "max"} and len(parts) == 2:
            a, b = parts
            # 将零用于净厚等量的下界，零继承另一项的单位。
            if a.unit_power != b.unit_power and not any(
                isinstance(arg, ast.Constant) and arg.value == 0 for arg in node.args
            ):
                raise ValueError("min/max 参数量纲须一致")
            return FormulaShape(a.unit_power or b.unit_power, a.sources | b.sources, True)
    raise ValueError("公式含有提示词白名单以外的语法或函数")


def _small_denominator_factors(node: ast.expr) -> frozenset[str] | None:
    """Find exposed small factors, with None denoting a structurally zero term."""
    if isinstance(node, ast.Name):
        symbol = ALIAS_TO_CANONICAL.get(node.id, node.id)
        return frozenset({symbol}) if symbol in UNSAFE_STANDALONE_DENOMINATORS else frozenset()
    if isinstance(node, ast.Constant):
        return None if node.value == 0 else frozenset()
    if isinstance(node, ast.UnaryOp):
        return _small_denominator_factors(node.operand)
    if isinstance(node, ast.Call):
        if node.func.id in {"abs", "sqrt"}:
            return _small_denominator_factors(node.args[0])
        if node.func.id in {"min", "max"}:
            left, right = (_small_denominator_factors(arg) for arg in node.args)
            return right if left is None else left if right is None else left & right
        return frozenset()
    if not isinstance(node, ast.BinOp):
        return frozenset()

    if isinstance(node.op, ast.Div):
        return _small_denominator_factors(node.left)
    if isinstance(node.op, ast.Pow):
        return _small_denominator_factors(node.left) if ast.literal_eval(node.right) > 0 else frozenset()

    left, right = _small_denominator_factors(node.left), _small_denominator_factors(node.right)
    if isinstance(node.op, ast.Mult):
        return None if left is None or right is None else left | right
    if isinstance(node.op, ast.Sub) and ast.dump(node.left) == ast.dump(node.right):
        return None
    # Only a factor common to both nonzero terms survives a sum or difference.
    if left is None:
        return right
    if right is None:
        return left
    return left & right


def _check_denominators(tree: ast.Expression, namespace: dict, n: int) -> None:
    evaluator = primitives.SafeFormulaEvaluator(tuple(namespace), min_std=0)
    for node in ast.walk(tree):
        if not isinstance(node, ast.BinOp):
            continue
        if isinstance(node.op, ast.Div):
            base = denominator = node.right
        elif isinstance(node.op, ast.Pow) and (exponent := ast.literal_eval(node.right)) < 0:
            base = node.left
            # Check the same values as division by the corresponding positive power.
            denominator = ast.BinOp(left=base, op=ast.Pow(), right=ast.Constant(value=-exponent))
        else:
            continue

        factors = _small_denominator_factors(base)
        if factors:
            names = sorted({part.id for part in ast.walk(base) if isinstance(part, ast.Name)
                            and ALIAS_TO_CANONICAL.get(part.id, part.id) in factors})
            labels = ", ".join(names)
            raise ValueError(
                f"禁止将局部小量 '{labels}' 作为分母因子或取负幂；"
                "请改用占总断面/之和的比例（如 A_sec / (A_mod + A_sec)）"
            )

        with np.errstate(all="ignore"):
            result = evaluator.evaluate("denominator", ast.unparse(denominator), namespace, n)
        values = result.values
        if not result.is_valid or values is None:
            raise ValueError(f"分母求值失败: {result.error_message}")
        min_abs, max_abs = float(np.min(np.abs(values))), float(np.max(np.abs(values)))
        if np.min(values) <= 0 <= np.max(values) or min_abs < 0.05 * max_abs:
            raise ValueError("分母在样本上变号、触零或动态范围过宽（最小/最大比 < 0.05）")
        if min_abs < evaluator.min_abs_denom:
            raise ValueError("分母存在接近 0 的数值，拒绝该特征以防数值奇异")


def evaluate_candidate(
    candidate: object, namespace: dict, n: int, *, min_std: float = 1e-8,
) -> tuple[dict, np.ndarray]:
    fields = {"name", "expression", "physical_rationale"}
    item = prompt_io.exact_keys(candidate, fields, "候选特征")
    name = prompt_io.nonempty_string(item["name"], "name").strip()
    expression = prompt_io.nonempty_string(item["expression"], "expression").strip()
    rationale = prompt_io.nonempty_string(item["physical_rationale"], "physical_rationale").strip()
    if not re.fullmatch(r"[a-z][a-z0-9]*(?:_[a-z0-9]+)*", name) or len(expression) > 512:
        raise ValueError("候选名须为 snake_case，公式长度最多 512 字符")
    if not any(character.isalpha() for character in rationale):
        raise ValueError("需要提供力学解释")
    tree = ast.parse(expression, mode="eval")
    shape = inspect_formula(tree.body)
    if not (shape.normalized_primitive or (shape.nonlinear and len(shape.sources) >= 2)):
        raise ValueError("公式需要组合多个变化量，或按提示词对衍生几何原语作无量纲化")

    _check_denominators(tree, namespace, n)
    evaluator = primitives.SafeFormulaEvaluator(tuple(namespace), min_std=min_std)
    with np.errstate(all="ignore"):
        result = evaluator.evaluate(name, expression, namespace, n)
    if not result.is_valid or result.values is None:
        raise ValueError(result.error_message)
    values = result.values
    power = shape.unit_power
    unit = "1" if power == 0 else "m" if power == 1 else f"m^{power}"
    return {"name": name, "expression": expression, "physical_rationale": rationale,
            "unit": unit, "values": values.tolist()}, values


def validate_candidate(
    candidate: object, namespace: dict, basis: np.ndarray, ref_namespace: dict[str, np.ndarray] | None = None,
) -> tuple[dict, np.ndarray]:
    feature, values = evaluate_candidate(candidate, namespace, len(basis))
    # 检出别名消去或 x5*x6/x6 等在观测点上退化为原始坐标的线性变换。
    scale = np.max(np.abs(values))
    scaled = values / scale
    std = scaled.std()
    if std == 0:
        raise ValueError("公式在已观测样本上为常量")
    standardized = (scaled - scaled.mean()) / std
    coefficients, *_ = np.linalg.lstsq(basis, standardized, rcond=None)
    if np.max(np.abs(standardized - basis @ coefficients)) < 1e-9:
        raise ValueError(f"公式在 {len(basis)} 个样本上等价于原始坐标的仿射变换")
    if ref_namespace is not None:
        check_feasible_domain_sanity(candidate, values, ref_namespace)
    return feature, standardized


def equivalent_feature(a: np.ndarray, b: np.ndarray, scaled_a: np.ndarray, scaled_b: np.ndarray) -> bool:
    if (np.allclose(scaled_a, scaled_b, atol=1e-9, rtol=0)
            or np.allclose(scaled_a, -scaled_b, atol=1e-9, rtol=0)):
        return True
    # 倒数或经常数缩放的倒数在两个特征均有定义时，逐点乘积为非零常数。
    product = (a / np.max(np.abs(a))) * (b / np.max(np.abs(b)))
    return bool(product[0] != 0 and np.allclose(product, product[0], rtol=1e-9, atol=0))


def validate(request: Path, raw_reply: bytes, *, geometry_template: Path | None = None) -> tuple[dict, Path]:
    manifest, context, _ = load_request(request)
    response_hash = prompt_io.digest(raw_reply)
    output = request / "responses" / response_hash
    output.mkdir(parents=True, exist_ok=True)
    (output / "llm_response_raw.txt").write_bytes(raw_reply)
    report = {"schema_version": 1, "prompt_sha256": manifest["prompt_sha256"],
              "response_sha256": response_hash, "valid": False, "errors": [], "checks": {}}
    candidates = {}
    try:
        parsed = prompt_io.parse_llm_json(raw_reply)
        answer = prompt_io.exact_keys(parsed, set(PROMPT_RESPONSES), "物理特征回复")
        namespace = {name: np.array(values, dtype=float) for name, values in context["namespace"].items()}
        x = np.column_stack([namespace[name] for name in DESIGN_VARIABLES])
        basis = np.column_stack([np.ones(len(x)), x - x.mean(axis=0)])
        path = geometry_template or Path(context.get("geometry_template", primitives.DEFAULT_TEMPLATE_PATH))
        config = primitives.DamGeometryConfig(path)
        regenerated = primitives.compute_geometric_primitives(x, config=config, return_dataframe=False)
        if any(not np.allclose(values, namespace[name], atol=1e-12, rtol=1e-10) for name, values in regenerated.items()):
            raise ValueError("几何模板与请求快照不一致，请用 --geometry-template 指定原模板")
        ref_namespace = build_feasible_reference_namespace(config)
        report["feasible_reference"] = {**FEASIBLE_REFERENCE_SETTINGS, "geometry_template": str(path.resolve())}
        for key, response in PROMPT_RESPONSES.items():
            items = answer[key]
            if not isinstance(items, list) or not MIN_CANDIDATES <= len(items) <= MAX_CANDIDATES:
                report["errors"].append(f"{key} 须提供 {MIN_CANDIDATES}–{MAX_CANDIDATES} 条候选")
                continue
            checks = report["checks"][response] = []
            accepted = candidates[response] = []
            accepted_by_name = {}
            for position, item in enumerate(items, 1):
                check = {"position": position, "name": item.get("name") if isinstance(item, dict) else None}
                try:
                    feature, standardized = validate_candidate(item, namespace, basis, ref_namespace=ref_namespace)
                    name = feature["name"]
                    values = np.asarray(feature["values"])
                    if name in accepted_by_name or any(
                        equivalent_feature(values, previous, standardized, scaled_previous)
                        for previous, scaled_previous in accepted_by_name.values()
                    ):
                        raise ValueError("候选名重复，或特征在初始样本上仿射等价或互为缩放倒数")
                    accepted_by_name[name] = (values, standardized)
                    accepted.append(feature)
                    check["valid"] = True
                except (SyntaxError, TypeError, ValueError, OverflowError) as exc:
                    check.update(valid=False, reason=str(exc))
                checks.append(check)
            if len(accepted) < MIN_CANDIDATES:
                report["errors"].append(f"{key} 仅有 {len(accepted)} 条候选通过校验，至少需要 {MIN_CANDIDATES} 条")
    except (UnicodeError, TypeError, ValueError) as exc:
        report["errors"].append(str(exc))
    report["valid"] = not report["errors"]
    if report["valid"]:
        llm_io.write_json(output / "candidate_pool.json", {
            "schema_version": 1, "prompt_sha256": manifest["prompt_sha256"],
            "response_sha256": response_hash, "design_ids": context["design_ids"],
            "responses": candidates,
        })
    llm_io.write_json(output / "validation_result.json", report)
    return report, output


def call(args: argparse.Namespace) -> tuple[dict, Path]:
    manifest, _, prompt = load_request(args.request)
    default_url, default_style, default_key = llm_io.PROVIDER_DEFAULTS.get(
        args.provider, (None, "chat-completions", None)
    )
    base_url, key_env = args.base_url or default_url, args.api_key_env or default_key
    if not base_url or not key_env or not args.model.strip() or not math.isfinite(args.timeout) or args.timeout <= 0:
        raise ValueError("需要模型、有效超时时间，以及提供商的 API 根地址和密钥环境变量名")
    style = args.api_style or default_style
    url = llm_io.endpoint(base_url, style)
    body = llm_io.request_body(args.model, prompt.decode("utf-8"), style, {})
    audit = {"prompt_sha256": manifest["prompt_sha256"], "provider": args.provider, "model": args.model,
             "endpoint": url, "request_body_sha256": prompt_io.digest(prompt_io.canonical_json(body))}
    if args.dry_run:
        return {**audit, "valid": True, "dry_run": True}, args.request
    key = os.environ.get(key_env)
    if not key:
        raise ValueError(f"请设置 {key_env} 环境变量")
    audit_path = args.request / "api_call.json"
    with audit_path.open("x", encoding="utf-8") as stream:
        json.dump({**audit, "status": "pending"}, stream)
    try:
        provider_reply, status = llm_io.post_json(url, body, key, args.timeout)
        llm_io.write_json(args.request / "provider_response.json", provider_reply)
        reply, provider_id, usage, model = llm_io.extract_answer(provider_reply, style)
        report, output = validate(args.request, reply.encode("utf-8"),
                                  geometry_template=getattr(args, "geometry_template", None))
        audit.update(status="validated" if report["valid"] else "invalid_output", http_status=status,
                     provider_id=provider_id, usage=usage, provider_model=model, validation_dir=str(output))
        return report, output
    except (OSError, ValueError, llm_io.APIError) as exc:
        audit.update(status="error", error=str(exc))
        raise
    finally:
        llm_io.write_json(audit_path, audit)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("render", help="由 initial20.csv 渲染现有物理特征模板")
    prepare.add_argument("--initial-csv", type=Path, required=True)
    prepare.add_argument("--output", type=Path, required=True)
    prepare.add_argument("--geometry-template", type=Path)
    check = commands.add_parser("validate", help="离线校验 LLM 的 JSON 回复")
    check.add_argument("--request", type=Path, required=True)
    check.add_argument("--response", type=Path, required=True)
    check.add_argument("--geometry-template", type=Path, help="旧请求或模板迁移时指定原几何模板")
    send = commands.add_parser("call", help="向提供商发送已渲染的提示词")
    send.add_argument("--request", type=Path, required=True)
    send.add_argument("--provider", default="openai", choices=("openai", "openrouter", "custom"))
    send.add_argument("--model", required=True)
    send.add_argument("--base-url")
    send.add_argument("--api-style", choices=("responses", "chat-completions"))
    send.add_argument("--api-key-env")
    send.add_argument("--geometry-template", type=Path, help="旧请求或模板迁移时指定原几何模板")
    send.add_argument("--timeout", type=float, default=120.0)
    send.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    try:
        if args.command == "render":
            result = render(args.initial_csv, args.output, args.geometry_template)
            print(f"{args.output / 'prompt.md'} | SHA-256 {result['prompt_sha256']}")
            return 0
        if args.command == "validate":
            result, output = validate(args.request, args.response.read_bytes(), geometry_template=args.geometry_template)
        else:
            result, output = call(args)
        if result.get("dry_run"):
            print(json.dumps(result, ensure_ascii=False))
        else:
            print(output / ("candidate_pool.json" if result["valid"] else "validation_result.json"))
        return 0 if result["valid"] else 2
    except (OSError, ValueError, llm_io.APIError) as exc:
        print(str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
