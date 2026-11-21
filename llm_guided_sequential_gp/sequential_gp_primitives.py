"""
确定性断面几何原语库 (Deterministic Cross-Sectional Geometric Primitive Library, DC-GPL)
及 AST 沙箱物理公式求值器。

放置路径: llm_guided_sequential_gp/sequential_gp_primitives.py
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import geostudio_automation
import numpy as np
import pandas as pd
from geostudio_automation.geometry import CORE_MATERIALS, GeometryBuild, build_geometry
from geostudio_automation.model import CONSTANTS, DesignVariables, GeometryError
from geostudio_automation.template import TemplateArchive
from geostudio_automation.xml_model import GeoStudioXmlModel
from shapely import get_parts, union_all
from shapely.geometry import LineString, Polygon
from shapely.geometry.base import BaseGeometry


DEFAULT_TEMPLATE_PATH = (
    Path(geostudio_automation.__file__).resolve().parents[1]
    / "templates"
    / "geo_studio_model_sample"
)


# =============================================================================
# 1. 原始变量与原语元数据定义
# =============================================================================

# 标准 8 维设计变量顺序 (x1 ~ x8)
CANONICAL_DESIGN_VARS: tuple[str, ...] = (
    "upstream_slope",                 # x1: 上游坡比 H/V
    "downstream_slope",               # x2: 下游坡比 H/V
    "crest_width",                    # x3: 物理坝顶宽度 (对称)
    "modulus_zone_height",            # x4: 坝顶增模区高度
    "secondary_top_width",            # x5: 次堆石区顶部宽度
    "secondary_height",               # x6: 次堆石区高度
    "secondary_bottom_width",         # x7: 次堆石区底部宽度
    "secondary_right_slope_distance", # x8: 次堆石区右腰距下游参考线的法向距离
)

# 兼容不同脚本或 CSV 表头的别名映射
DESIGN_VAR_ALIASES: dict[str, tuple[str, ...]] = {
    "upstream_slope": ("upstream_slope", "x1", "m_up"),
    "downstream_slope": ("downstream_slope", "x2", "m_down"),
    "crest_width": ("crest_width", "x3", "B_crest"),
    "modulus_zone_height": ("modulus_zone_height", "x4", "h_mod"),
    "secondary_top_width": ("secondary_top_width", "x5", "w_sec_top"),
    "secondary_height": ("secondary_height", "x6", "h_sec"),
    "secondary_bottom_width": ("secondary_bottom_width", "x7", "w_sec_bot"),
    "secondary_right_slope_distance": (
        "secondary_right_slope_distance",
        "secondary_right_dist",
        "x8",
        "d_right",
    ),
}

# 14 项确定性几何原语 + 1 项常量说明（用于自动生成 Prompt 字典与校验白名单）
PRIMITIVE_METADATA: dict[str, dict[str, str]] = {
    "A_tot": {"unit": "m^2", "desc": "坝体主体域面积（GeoStudio 材料 4、5、6、12 的并集）"},
    "A_mod": {"unit": "m^2", "desc": "坝顶增模区(Material 12)面积（顶部高模量加强帽盖）"},
    "A_sec": {"unit": "m^2", "desc": "次堆石区(Material 5)面积（内部低模量软弱压缩区）"},
    "A_drain": {"unit": "m^2", "desc": "底部排水区(Material 6)面积（底部高刚度支承区）"},
    "A_main": {"unit": "m^2", "desc": "主堆石区(Material 4)面积（高刚度主承重与抗推骨架）"},
    "A_ds_shell": {"unit": "m^2", "desc": "次堆石区顶底高程之间、其右腰下游侧的材料 4 面积"},
    "B_bot": {"unit": "m", "desc": "坝体主体域在坝底高程的水平宽度"},
    "B_mod_bot": {"unit": "m", "desc": "坝顶增模区底面宽度"},
    "h_drain": {"unit": "m", "desc": "底部排水区高度（次堆石区底面至坝基的距离）"},
    "d_ds_perp": {"unit": "m", "desc": "次堆石区右腰到 GeoStudio 下游参考线的法向距离，等于 x8"},
    "d_us_min": {"unit": "m", "desc": "次堆石区左腰距主体域上游边界的最小水平净宽"},
    "X_c_sec": {"unit": "m", "desc": "次堆石区形心相对坝轴线的横坐标，正值偏下游"},
    "Y_c_sec": {"unit": "m", "desc": "次堆石区形心相对坝底的高程"},
    "W_sec_mean": {"unit": "m", "desc": "次堆石区平均水平宽度 0.5*(x5+x7)"},
    "H_dam": {"unit": "m", "desc": "大坝总高度常量"},
}


@dataclass(frozen=True)
class DamGeometryConfig:
    """GeoStudio 建模模板路径及公式命名选项。"""

    template_path: Path = DEFAULT_TEMPLATE_PATH
    include_x_aliases: bool = True


# =============================================================================
# 2. 从 GeoStudio 材料分区计算几何原语
# =============================================================================

def _extract_design_columns(
    X: pd.DataFrame | np.ndarray | Mapping[str, float | np.ndarray],
    feature_names: Sequence[str] | None = None,
) -> dict[str, np.ndarray]:
    """将 DataFrame / ndarray / dict 统一标准化为 8 个 shape=(N,) 的 float64 数组字典。"""
    if isinstance(X, pd.DataFrame):
        available_keys = set(X.columns)
        get_col = lambda k: X[k].to_numpy(dtype=np.float64)
    elif isinstance(X, Mapping):
        available_keys = set(X.keys())
        get_col = lambda k: np.atleast_1d(np.asarray(X[k], dtype=np.float64))
    else:
        arr = np.asarray(X, dtype=np.float64)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        if arr.ndim != 2 or arr.shape[1] != len(CANONICAL_DESIGN_VARS):
            raise ValueError(
                f"期望输入 8 列设计变量，实际收到 shape={arr.shape}。"
            )
        names = tuple(feature_names) if feature_names is not None else CANONICAL_DESIGN_VARS
        if len(names) != len(CANONICAL_DESIGN_VARS) or len(set(names)) != len(names):
            raise ValueError("feature_names 必须提供 8 个不同的设计变量名")
        col_map = {name: arr[:, idx] for idx, name in enumerate(names)}
        available_keys = set(col_map.keys())
        get_col = lambda k: col_map[k]

    resolved: dict[str, np.ndarray] = {}
    for canonical_name in CANONICAL_DESIGN_VARS:
        matched_key = None
        for alias in DESIGN_VAR_ALIASES[canonical_name]:
            if alias in available_keys:
                matched_key = alias
                break
        if matched_key is None:
            raise KeyError(
                f"输入数据缺少设计变量 '{canonical_name}'（支持别名: {DESIGN_VAR_ALIASES[canonical_name]}），"
                f"当前表头为: {sorted(available_keys)}"
            )
        column = np.atleast_1d(get_col(matched_key)).astype(np.float64)
        if column.ndim != 1:
            raise ValueError(f"设计变量 '{canonical_name}' 必须是标量或一维数组")
        resolved[canonical_name] = column

    return dict(zip(resolved, np.broadcast_arrays(*resolved.values()), strict=True))


def _horizontal_limits(shape: BaseGeometry, elevation: float) -> tuple[float, float]:
    min_x, _, max_x, _ = shape.bounds
    section = shape.intersection(
        LineString([(min_x - 1.0, elevation), (max_x + 1.0, elevation)])
    )
    left, _, right, _ = section.bounds
    if section.is_empty:
        raise GeometryError(f"高程 {elevation:g} m 处未找到分区截面")
    return left, right


def _measure_primitives(geometry: GeometryBuild) -> dict[str, float]:
    # 使用实际写入有限元模型的材料分区，合并同一材料的各施工层。
    materials = {
        material: union_all(
            [cell.polygon for cell in geometry.cells if cell.material_id == material]
        )
        for material in CORE_MATERIALS
    }
    main, secondary, drain, modulus = (materials[4], materials[5], materials[6], materials[12])
    core = union_all(list(materials.values()))
    base = CONSTANTS.dam_base_elevation
    _, secondary_bottom, _, secondary_top = secondary.bounds
    sec_left_bottom, sec_right_bottom = _horizontal_limits(secondary, secondary_bottom)
    sec_left_top, sec_right_top = _horizontal_limits(secondary, secondary_top)
    base_left, base_right = _horizontal_limits(core, base)
    mod_left, mod_right = _horizontal_limits(modulus, modulus.bounds[1])

    # 主堆石保护壳限于次堆石区顶底高程之间及右腰下游侧。
    shell_right = core.bounds[2] + 1.0
    shell_window = Polygon([
        (sec_right_bottom, secondary_bottom),
        (shell_right, secondary_bottom),
        (shell_right, secondary_top),
        (sec_right_top, secondary_top),
    ])

    # 真实主体域的坡面可含折点，不能只用次堆石区的顶底两端。
    elevations = {secondary_bottom, secondary_top}
    elevations.update(
        y for part in get_parts(core) for _, y in part.exterior.coords
        if secondary_bottom < y < secondary_top
    )
    min_upstream_width = min(
        sec_left_bottom + (sec_left_top - sec_left_bottom)
        * (y - secondary_bottom) / (secondary_top - secondary_bottom)
        - _horizontal_limits(core, y)[0]
        for y in elevations
    )
    centroid = secondary.centroid
    return {
        "A_tot": core.area,
        "A_mod": modulus.area,
        "A_sec": secondary.area,
        "A_drain": drain.area,
        "A_main": main.area,
        "A_ds_shell": main.intersection(shell_window).area,
        "B_bot": base_right - base_left,
        "B_mod_bot": mod_right - mod_left,
        "h_drain": drain.bounds[3] - base,
        "d_ds_perp": geometry.derived.secondary_right_slope_distance,
        "d_us_min": min_upstream_width,
        "X_c_sec": centroid.x - CONSTANTS.dam_axis_x,
        "Y_c_sec": centroid.y - base,
        "W_sec_mean": 0.5 * (sec_right_bottom - sec_left_bottom + sec_right_top - sec_left_top),
        "H_dam": CONSTANTS.crest_elevation - base,
    }


def compute_geometric_primitives(
    X: pd.DataFrame | np.ndarray | Mapping[str, float | np.ndarray],
    config: DamGeometryConfig | None = None,
    feature_names: Sequence[str] | None = None,
    return_dataframe: bool = True,
) -> pd.DataFrame | dict[str, np.ndarray]:
    """用 GeoStudio 模板重建断面，从材料多边形计算原语。"""
    cfg = config or DamGeometryConfig()
    cols = _extract_design_columns(X, feature_names=feature_names)
    index = X.index if isinstance(X, pd.DataFrame) else None

    rows = []
    if len(next(iter(cols.values()))):
        template = TemplateArchive.load(Path(cfg.template_path))
        model = GeoStudioXmlModel(template.root(), reference_design=template.reference_design)
        for position, row in enumerate(zip(*cols.values(), strict=True)):
            if isinstance(X, pd.DataFrame) and "design_id" in X:
                design_id = str(X["design_id"].iat[position])
            else:
                design_id = f"sample_{position + 1}"
            design = DesignVariables(design_id, *map(float, row))
            try:
                geometry = build_geometry(
                    model.original_regions, design, reference_design=model.reference_design
                )
                rows.append(_measure_primitives(geometry))
            except GeometryError as exc:
                raise GeometryError(f"{design_id}: {exc}") from exc

    result = pd.DataFrame(rows, index=index, columns=list(PRIMITIVE_METADATA), dtype=float)
    if return_dataframe:
        return result
    return {name: result[name].to_numpy(copy=True) for name in result}


def build_formula_namespace(
    X: pd.DataFrame | np.ndarray | Mapping[str, float | np.ndarray],
    config: DamGeometryConfig | None = None,
    feature_names: Sequence[str] | None = None,
) -> dict[str, np.ndarray]:
    """构建供安全公式求值器使用的完整变量命名空间（原始 8 变量 + 14 项原语 + H_dam）。"""
    cfg = config or DamGeometryConfig()
    design_cols = _extract_design_columns(X, feature_names=feature_names)
    prim_cols = compute_geometric_primitives(
        X, config=cfg, feature_names=feature_names, return_dataframe=False
    )
    assert isinstance(prim_cols, dict)

    ns: dict[str, np.ndarray] = {}
    ns.update(design_cols)
    ns.update(prim_cols)

    if cfg.include_x_aliases:
        for idx, name in enumerate(CANONICAL_DESIGN_VARS, start=1):
            ns[f"x{idx}"] = design_cols[name]
        ns["secondary_right_dist"] = design_cols["secondary_right_slope_distance"]

    return ns


# =============================================================================
# 3. AST 安全物理公式解析与校验器（杜绝 eval 注入与数值炸裂）
# =============================================================================

@dataclass
class FormulaEvaluationResult:
    name: str
    expression: str
    is_valid: bool
    values: np.ndarray | None
    used_symbols: tuple[str, ...]
    error_message: str = ""


class SafeFormulaEvaluator:
    """
    基于 Python AST 白名单的向量化公式求值器。
    - 仅允许基础算术 (+, -, *, /, **) 与安全数学函数 (sqrt, abs, log, log1p, min, max, pow)
    - 严格拦截除零、负数开方/对数、NaN/Inf、零方差退化与单变量纯线性变换
    """

    _ALLOWED_FUNCS: dict[str, Any] = {
        "sqrt": np.sqrt,
        "abs": np.abs,
        "log": np.log,
        "log1p": np.log1p,
        "min": np.minimum,
        "max": np.maximum,
        "pow": np.power,
    }

    def __init__(
        self,
        allowed_symbols: Sequence[str],
        min_std: float = 1e-8,
        min_abs_denom: float = 1e-8,
        max_abs_exponent: float = 4.0,
    ) -> None:
        self.allowed_symbols = set(allowed_symbols)
        self.min_std = min_std
        self.min_abs_denom = min_abs_denom
        self.max_abs_exponent = max_abs_exponent

    def _eval_node(
        self,
        node: ast.AST,
        namespace: Mapping[str, np.ndarray],
        used: set[str],
    ) -> np.ndarray | float:
        if isinstance(node, ast.Expression):
            return self._eval_node(node.body, namespace, used)

        if isinstance(node, ast.Constant):
            if isinstance(node.value, (int, float)) and not isinstance(node.value, bool):
                return float(node.value)
            raise ValueError(f"不允许的非数值常量: {node.value!r}")

        if isinstance(node, ast.Name):
            if node.id not in self.allowed_symbols or node.id not in namespace:
                raise ValueError(f"使用了白名单之外的未知符号: '{node.id}'")
            used.add(node.id)
            return namespace[node.id]

        if isinstance(node, ast.UnaryOp):
            operand = self._eval_node(node.operand, namespace, used)
            if isinstance(node.op, ast.UAdd):
                return +operand
            if isinstance(node.op, ast.USub):
                return -operand
            raise ValueError(f"不支持的一元运算符: {type(node.op).__name__}")

        if isinstance(node, ast.BinOp):
            left = self._eval_node(node.left, namespace, used)
            right = self._eval_node(node.right, namespace, used)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                denom = np.asarray(right, dtype=np.float64)
                if np.any(np.abs(denom) < self.min_abs_denom):
                    raise ValueError("分母存在接近 0 的数值，拒绝该特征以防数值奇异")
                return left / right
            if isinstance(node.op, ast.Pow):
                exp_arr = np.asarray(right, dtype=np.float64)
                if np.any(np.abs(exp_arr) > self.max_abs_exponent):
                    raise ValueError(f"幂指数绝对值超过上限 {self.max_abs_exponent}")
                base_arr = np.asarray(left, dtype=np.float64)
                if np.any((base_arr < 0) & (~np.isclose(exp_arr, np.round(exp_arr)))):
                    raise ValueError("出现负数非整数次幂")
                return np.power(left, right)
            raise ValueError(f"不支持的二元运算符: {type(node.op).__name__}")

        if isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name):
                raise ValueError("不允许调用属性方法或复杂可调用对象")
            func_name = node.func.id
            if func_name not in self._ALLOWED_FUNCS:
                raise ValueError(
                    f"不允许调用函数 '{func_name}'，仅支持: {sorted(self._ALLOWED_FUNCS)}"
                )
            if node.keywords:
                raise ValueError("公式函数调用不支持关键字参数")
            args = [self._eval_node(arg, namespace, used) for arg in node.args]

            # 针对定义域敏感函数的严格预检
            if func_name == "sqrt":
                if len(args) != 1 or np.any(np.asarray(args[0]) < 0.0):
                    raise ValueError("sqrt 参数数量错误或存在负数")
            elif func_name == "log":
                if len(args) != 1 or np.any(np.asarray(args[0]) <= 0.0):
                    raise ValueError("log 参数数量错误或存在非正数")
            elif func_name == "log1p":
                if len(args) != 1 or np.any(np.asarray(args[0]) <= -1.0):
                    raise ValueError("log1p 参数数量错误或 <= -1")
            elif func_name in ("min", "max", "pow"):
                if len(args) != 2:
                    raise ValueError(f"{func_name} 需要且仅需要 2 个参数")

            return self._ALLOWED_FUNCS[func_name](*args)

        raise ValueError(f"禁止的语法节点类型: {type(node).__name__}")

    def evaluate(
        self,
        name: str,
        expression: str,
        namespace: Mapping[str, np.ndarray],
        n_samples: int,
    ) -> FormulaEvaluationResult:
        """解析、校验并向量化计算单条 LLM 候选物理特征公式。"""
        used: set[str] = set()
        try:
            tree = ast.parse(expression.strip(), mode="eval")
            raw_val = self._eval_node(tree, namespace, used)
            arr = np.broadcast_to(
                np.asarray(raw_val, dtype=np.float64), (n_samples,)
            ).copy()

            if not np.all(np.isfinite(arr)):
                return FormulaEvaluationResult(
                    name, expression, False, None, tuple(sorted(used)),
                    "求值结果包含 NaN 或 Inf",
                )
            if float(np.std(arr)) < self.min_std:
                return FormulaEvaluationResult(
                    name, expression, False, None, tuple(sorted(used)),
                    f"特征标准差过小 (< {self.min_std})，属于退化常量特征",
                )
            return FormulaEvaluationResult(
                name=name,
                expression=expression.strip(),
                is_valid=True,
                values=arr,
                used_symbols=tuple(sorted(used)),
            )
        except Exception as exc:  # noqa: BLE001
            return FormulaEvaluationResult(
                name=name,
                expression=expression.strip(),
                is_valid=False,
                values=None,
                used_symbols=tuple(sorted(used)),
                error_message=str(exc),
            )


# =============================================================================
# 4. 面向 LLM 的原语统计摘要表生成器（严格仅用前 t 个样本，零泄露）
# =============================================================================

def format_primitive_catalog_markdown(
    X_observed: pd.DataFrame,
    y_observed: pd.DataFrame | None = None,
    config: DamGeometryConfig | None = None,
) -> str:
    """
    基于当前已观测的前 t 个样本，生成喂给 LLM 的原语取值范围与 Spearman 相关性摘要表。
    注意：严禁传入第 t+1..60 个未观测样本！
    """
    prim_df = compute_geometric_primitives(X_observed, config=config, return_dataframe=True)
    assert isinstance(prim_df, pd.DataFrame)

    resp_cols = list(y_observed.columns) if y_observed is not None else []
    header = "| 原语名 | 单位 | 当前观测范围 [min, max] | 物理含义 |"
    sep = "|---|---:|---|---|"
    if resp_cols:
        header += " " + " | ".join(f"ρ({c})" for c in resp_cols) + " |"
        sep += " " + " | ".join("---:" for _ in resp_cols) + " |"

    lines = [header, sep]
    for col in prim_df.columns:
        meta = PRIMITIVE_METADATA[col]
        vmin, vmax = float(prim_df[col].min()), float(prim_df[col].max())
        row = f"| `{col}` | {meta['unit']} | [{vmin:.2f}, {vmax:.2f}] | {meta['desc']} |"
        if resp_cols and y_observed is not None:
            rhos = []
            for rcol in resp_cols:
                if prim_df[col].nunique() <= 1:
                    rhos.append("N/A")
                else:
                    rho = prim_df[col].corr(y_observed[rcol], method="spearman")
                    rhos.append(f"{rho:+.3f}" if np.isfinite(rho) else "N/A")
            row += " " + " | ".join(rhos) + " |"
        lines.append(row)

    return "\n".join(lines)


# =============================================================================
# 5. 端到端自检脚本（自动尝试加载 data/dataset.csv，无文件时使用真实边界回退）
# =============================================================================

if __name__ == "__main__":
    repo_csv = Path(__file__).resolve().parent.parent / "data" / "dataset.csv"
    if repo_csv.exists():
        print(f"[INFO] 加载实际数据集: {repo_csv}")
        df_all = pd.read_csv(repo_csv)
        df_init = df_all.iloc[:20].copy()  # 严格仅取前 20 点做冷启动演示
    else:
        print("[INFO] 未检测到 data/dataset.csv，使用 8 维大坝可行域典型方案演示：")
        df_init = pd.DataFrame(
            [
                {
                    "upstream_slope": 1.45,
                    "downstream_slope": 1.40,
                    "crest_width": 8.0,
                    "modulus_zone_height": 12.0,
                    "secondary_top_width": 16.0,
                    "secondary_height": 32.0,
                    "secondary_bottom_width": 44.0,
                    "secondary_right_slope_distance": 6.0,
                },
                {
                    "upstream_slope": 1.60,
                    "downstream_slope": 1.55,
                    "crest_width": 10.0,
                    "modulus_zone_height": 18.0,
                    "secondary_top_width": 12.0,
                    "secondary_height": 26.0,
                    "secondary_bottom_width": 36.0,
                    "secondary_right_slope_distance": 10.0,
                },
                {
                    "upstream_slope": 1.70,
                    "downstream_slope": 1.50,
                    "crest_width": 10.0,
                    "modulus_zone_height": 15.0,
                    "secondary_top_width": 20.0,
                    "secondary_height": 30.0,
                    "secondary_bottom_width": 60.0,
                    "secondary_right_slope_distance": 5.0,
                },
            ]
        )

    cfg = DamGeometryConfig()
    ns = build_formula_namespace(df_init, config=cfg)
    prim_df = pd.DataFrame({name: ns[name] for name in PRIMITIVE_METADATA}, index=df_init.index)
    print("\n=== 1. 前 3 个样本的确定性几何原语输出 ===")
    print(prim_df.head(3).round(2).to_string())

    # 输出由 GeoStudio 几何得到的关键净距与分区面积。
    print("\n=== 2. 几何拓扑合理性自检 ===")
    print(f"  min(h_drain)  = {prim_df['h_drain'].min():.2f} m (应 > 0)")
    print(f"  min(A_main)   = {prim_df['A_main'].min():.2f} m^2 (应 > 0)")
    print(f"  min(d_us_min) = {prim_df['d_us_min'].min():.2f} m (应 > 0)")

    # 演示 AST 沙箱公式求值器验证 LLM 候选特征
    evaluator = SafeFormulaEvaluator(allowed_symbols=list(ns.keys()))

    demo_formulas = [
        ("soft_area_ratio", "A_sec / A_tot"),
        ("crest_cap_stiffness_ratio", "A_mod / (A_tot * (1.0 + Y_c_sec / H_dam))"),
        ("ds_shell_vs_soft_ratio", "A_ds_shell / (A_sec + 1.0)"),
        ("soft_core_eccentricity", "X_c_sec / B_bot"),
        ("invalid_zero_div_test", "A_sec / (H_dam - 74.0)"),
    ]
    print("\n=== 3. AST 安全公式求值器测试 ===")
    for fname, fexpr in demo_formulas:
        res = evaluator.evaluate(fname, fexpr, ns, n_samples=len(df_init))
        if res.is_valid and res.values is not None:
            print(
                f"  [PASS] {fname:28s} | range=[{res.values.min():.4f}, {res.values.max():.4f}] | symbols={res.used_symbols}"
            )
        else:
            print(f"  [REJECT] {fname:26s} | reason={res.error_message}")
