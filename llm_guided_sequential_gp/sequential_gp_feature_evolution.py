"""Build residual-based physical feature requests and validate their replies."""

from __future__ import annotations

import json
import math
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
from scipy.stats import rankdata

if __package__:
    from . import sequential_gp_additive as additive
    from . import sequential_gp_feature_candidates as candidates
    from . import sequential_gp_feature_selection as features
    from . import sequential_gp_llm as llm
    from . import sequential_gp_model as gp
    from . import sequential_gp_primitives as primitives
    from . import sequential_gp_prompt as prompt_io
else:
    import sequential_gp_additive as additive
    import sequential_gp_feature_candidates as candidates
    import sequential_gp_feature_selection as features
    import sequential_gp_llm as llm
    import sequential_gp_model as gp
    import sequential_gp_primitives as primitives
    import sequential_gp_prompt as prompt_io


PROMPT_TEMPLATE = Path(__file__).resolve().parent / "prompts" / "gp_feature_evolution.md"
CONTEXT_TOKEN = "{{CONTEXT_JSON}}"
BRIEF_TOKEN = "{{DIAGNOSTIC_BRIEF}}"
DENOMINATOR_TOKEN = "{{UNSAFE_DENOMINATORS}}"
REQUEST_SCHEMA_VERSION = 4
TOP_COUNT = 3
PROPOSAL_COUNT = 2
FEATURE_CAP = features.FEATURE_CAP
TRIGGER_MILESTONES = (32, 44)
TRIGGER_MIN_SAMPLE_COUNT = 28
TRIGGER_COOLDOWN = 8
ONLINE_TRIGGER_WINDOW = 4
ONLINE_Z_THRESHOLD = 2.0
ONLINE_RMSE_RATIO = 1.05
LOW_SHARE_PRUNE_THRESHOLD = 0.05
TRIGGER_MODES = ("milestones", "anomaly", "anomaly_or_milestones")
MIN_RESIDUAL_SIGNAL_ABS_RHO = 0.20
SIGNAL_STRONG_ABS_RHO = 0.35
ACTIVE_COLLINEARITY_WARN_RHO = 0.90
ACTIVE_COLLINEARITY_DROP_RHO = 0.98
SHIFT_MIN_EFFECT_Z = 0.50
SHIFT_STRONG_EFFECT_Z = 0.75
SPARSE_SSE_SHARE_PCT = 55
HIGH_SHARE_WEAK_RHO_CUTOFF = 0.25


def should_trigger_evolution(
    sample_count: int,
    *,
    last_trigger_count: int | Mapping[str, int | None] | None,
    recent_online_zscores: dict[str, list[float]] | None = None,
    recent_online_rmse_vs_standard: dict[str, tuple[float, float]] | None = None,
    recent_online_physical_shares: dict[str, list[float]] | None = None,
    has_active_features: Mapping[str, bool] | None = None,
    mode: str = "anomaly_or_milestones",
    milestones: tuple[int, ...] = TRIGGER_MILESTONES,
    min_sample_count: int = TRIGGER_MIN_SAMPLE_COUNT,
    cooldown: int = TRIGGER_COOLDOWN,
    z_threshold: float = ONLINE_Z_THRESHOLD,
    low_share_threshold: float = LOW_SHARE_PRUNE_THRESHOLD,
) -> tuple[bool, list[str]]:
    """
    低频/异常触发控制器（逐响应独立冷却）：
    - 禁止在 N < 28 的极小样本期（如 t=24）以及冷却期内盲目高频改特征；
    - 冷却期按响应独立判定，避免单个响应异常触发锁死其余响应或吞掉里程碑；
    - mode="milestones": 仅在指定节点（默认 t=32, 44）触发未处于冷却期的响应；
    - mode="anomaly": 只在已揭晓的在线预测窗口失配或物理核持续退化时触发；
    - mode="anomaly_or_milestones": 某响应脱离冷却期后，到达里程碑或近窗出现 |z| > 2.0 /
      在线 RMSE 落后同构 8D 基线 GP / 在役物理核方差占比持续低于阈值时触发，并返回需要进化的响应列表。
    """
    if mode not in TRIGGER_MODES:
        raise ValueError(f"unknown evolution trigger mode: {mode}")
    if sample_count < min_sample_count:
        return False, []

    is_milestone = mode != "anomaly" and sample_count in milestones
    flagged = []
    for response in gp.RESPONSES:
        last_t = (
            last_trigger_count.get(response)
            if isinstance(last_trigger_count, Mapping)
            else last_trigger_count
        )
        if last_t is not None and (sample_count - last_t) < cooldown:
            continue
        if is_milestone:
            flagged.append(response)
            continue
        if mode == "milestones":
            continue
        z_vals = (recent_online_zscores or {}).get(response, [])
        z_bad = any(abs(z) > z_threshold for z in z_vals)
        rmse_pair = (recent_online_rmse_vs_standard or {}).get(response)
        # rmse_pair = (current_dynamic_window_rmse, matched_8d_window_rmse)
        rmse_bad = rmse_pair is not None and rmse_pair[0] > rmse_pair[1] * ONLINE_RMSE_RATIO
        share_vals = (recent_online_physical_shares or {}).get(response, [])
        low_share_bad = (
            bool((has_active_features or {}).get(response))
            and len(share_vals) >= ONLINE_TRIGGER_WINDOW
            and all(share < low_share_threshold for share in share_vals)
        )
        if z_bad or rmse_bad or low_share_bad:
            flagged.append(response)
    return bool(flagged), flagged

def write_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        try:
            stream.write(data)
        except BaseException:
            stream.close()
            temporary.unlink(missing_ok=True)
            raise
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def write_json(path: Path, data: dict) -> None:
    content = json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    write_bytes(path, content.encode("utf-8"))


def _spearman_corr(a: np.ndarray, b: np.ndarray) -> float | None:
    """Compute Spearman rank correlation safely, returning None when undefined."""
    a_arr = np.asarray(a, dtype=float)
    b_arr = np.asarray(b, dtype=float)
    if len(a_arr) < 3 or not (np.isfinite(a_arr).all() and np.isfinite(b_arr).all()):
        return None
    ranks = rankdata(np.column_stack((a_arr, b_arr)), axis=0)
    centered = ranks - ranks.mean(axis=0)
    norms = np.linalg.norm(centered, axis=0)
    if norms.min() <= 1e-12:
        return None
    return float(np.clip((centered[:, 0] @ centered[:, 1]) / np.prod(norms), -1.0, 1.0))


def _safe_ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray | None:
    """Divide arrays only when denominator is finite and strictly nonzero."""
    num = np.asarray(numerator, dtype=float)
    den = np.asarray(denominator, dtype=float)
    if not (np.isfinite(num).all() and np.isfinite(den).all()) or np.any(np.abs(den) <= 1e-12):
        return None
    return num / den


def _residual_signal_catalog(namespace: dict, x: np.ndarray | None = None) -> list[tuple[str, str, np.ndarray]]:
    """Pre-compute interpretable dimensionless quantities for residual correlation diagnosis."""
    del x
    specs: list[tuple[str, str, np.ndarray | None]] = [
        ("软区面积占比", "A_sec / A_tot", _safe_ratio(namespace["A_sec"], namespace["A_tot"])),
        ("帽盖面积占比", "A_mod / A_tot", _safe_ratio(namespace["A_mod"], namespace["A_tot"])),
        ("排水区面积占比", "A_drain / A_tot", _safe_ratio(namespace["A_drain"], namespace["A_tot"])),
        ("主堆石面积占比", "A_main / A_tot", _safe_ratio(namespace["A_main"], namespace["A_tot"])),
        ("下游保护壳占比", "A_ds_shell / A_tot", _safe_ratio(namespace["A_ds_shell"], namespace["A_tot"])),
        ("软区/主堆石面积比", "A_sec / A_main", _safe_ratio(namespace["A_sec"], namespace["A_main"])),
        ("帽盖/(帽盖+主堆石)", "A_mod / (A_mod + A_main)", _safe_ratio(namespace["A_mod"], namespace["A_mod"] + namespace["A_main"])),
        ("软区/(软区+排水)", "A_sec / (A_sec + A_drain)", _safe_ratio(namespace["A_sec"], namespace["A_sec"] + namespace["A_drain"])),
        ("软区相对高程", "Y_c_sec / H_dam", _safe_ratio(namespace["Y_c_sec"], namespace["H_dam"])),
        ("软区相对偏心", "X_c_sec / B_bot", _safe_ratio(namespace["X_c_sec"], namespace["B_bot"])),
        ("排水区相对高度", "h_drain / H_dam", _safe_ratio(namespace["h_drain"], namespace["H_dam"])),
        ("上游相对净宽", "d_us_min / B_bot", _safe_ratio(namespace["d_us_min"], namespace["B_bot"])),
        ("下游相对净距", "d_ds_perp / B_bot", _safe_ratio(namespace["d_ds_perp"], namespace["B_bot"])),
        ("软区相对平均宽", "W_sec_mean / B_bot", _safe_ratio(namespace["W_sec_mean"], namespace["B_bot"])),
        ("上下游坡比不对称度", "(x2 - x1) / (x1 + x2)", _safe_ratio(namespace["x2"] - namespace["x1"], namespace["x1"] + namespace["x2"])),
    ]
    return [(label, expr, vals) for label, expr, vals in specs if vals is not None]


def _residual_signal_correlations(
    namespace: dict,
    x: np.ndarray,
    residual: np.ndarray,
    *,
    active_columns: Sequence[tuple[str, np.ndarray]] = (),
    limit: int = 5,
) -> list[dict]:
    """Filter and rank dimensionless signals by Spearman correlation with residual."""
    items = []
    for label, expr, values in _residual_signal_catalog(namespace, x):
        rho = _spearman_corr(values, residual)
        if rho is None or abs(rho) < MIN_RESIDUAL_SIGNAL_ABS_RHO:
            continue
        max_active_name, max_active_rho = None, None
        for act_name, act_vals in active_columns:
            act_rho = _spearman_corr(values, act_vals)
            if act_rho is not None and (max_active_rho is None or abs(act_rho) > abs(max_active_rho)):
                max_active_name, max_active_rho = act_name, act_rho
        if max_active_rho is not None and abs(max_active_rho) > ACTIVE_COLLINEARITY_DROP_RHO:
            continue
        if max_active_rho is not None and abs(max_active_rho) > ACTIVE_COLLINEARITY_WARN_RHO:
            collinearity_note = f"与在役 {max_active_name} 共线 ({max_active_rho:+.2f})，仅可作替换评估"
        elif active_columns:
            collinearity_note = "与在役特征低共线"
        else:
            collinearity_note = "无在役冲突"
        strength = "中强相关" if abs(rho) >= SIGNAL_STRONG_ABS_RHO else "弱相关"
        direction = "正相关(越大越易低估)" if rho > 0 else "负相关(越大越易高估)"
        items.append({
            "label": label,
            "expression": expr,
            "rho_residual": round(rho, 2),
            "strength": strength,
            "direction": direction,
            "max_active_rho": round(max_active_rho, 2) if max_active_rho is not None else None,
            "collinearity_note": collinearity_note,
        })
    items.sort(key=lambda item: -abs(item["rho_residual"]))
    return items[:limit]


def _case_dimensionless_ratios(namespace: dict, idx: int, *, residual_val: float | None = None, observed_val: float | None = None) -> dict[str, float]:
    """Compute compact dimensionless percentages/ratios for a single cross-section."""
    a_tot = float(namespace["A_tot"][idx]) or 1.0
    b_bot = float(namespace["B_bot"][idx]) or 1.0
    h_dam = float(namespace["H_dam"][idx]) or 1.0
    ratios = {
        "soft_pct": round(100.0 * float(namespace["A_sec"][idx]) / a_tot, 1),
        "mod_pct": round(100.0 * float(namespace["A_mod"][idx]) / a_tot, 1),
        "drain_pct": round(100.0 * float(namespace["A_drain"][idx]) / a_tot, 1),
        "ds_shell_pct": round(100.0 * float(namespace["A_ds_shell"][idx]) / a_tot, 1),
        "rel_elev": round(float(namespace["Y_c_sec"][idx]) / h_dam, 2),
        "ecc_ratio": round(float(namespace["X_c_sec"][idx]) / b_bot, 2),
        "us_clear_ratio": round(float(namespace["d_us_min"][idx]) / b_bot, 2),
        "ds_clear_ratio": round(float(namespace["d_ds_perp"][idx]) / b_bot, 2),
    }
    if residual_val is not None and observed_val is not None and abs(observed_val) > 1e-12:
        ratios["rel_err_pct"] = round(100.0 * residual_val / observed_val, 1)
    return ratios


def _shift_strength_label(z: float) -> str:
    az = abs(z)
    if az >= SHIFT_STRONG_EFFECT_Z:
        return "强偏移"
    if az >= SHIFT_MIN_EFFECT_Z:
        return "显著偏移"
    return "弱偏移"


def _residual_ard_summary(length_scales: dict | None) -> dict[str, list[str]] | None:
    """Bucket 8D residual kernel ARD length scales into qualitative sensitivity labels."""
    residual_ls = (length_scales or {}).get("residual")
    if not isinstance(residual_ls, dict) or not residual_ls:
        return None
    center = gp.DEFAULT_PRIOR_CONFIG.base_factor * math.sqrt(len(gp.VARIABLES))
    alias_for = {canonical: alias for alias, canonical in candidates.ALIAS_TO_CANONICAL.items() if alias.startswith("x")}
    focused, relaxed, moderate = [], [], []
    for var_name in gp.VARIABLES:
        ell = residual_ls.get(var_name)
        if not isinstance(ell, (int, float)) or not math.isfinite(ell):
            continue
        tag = alias_for.get(var_name, var_name)
        if ell <= 1.0 * center:
            focused.append(tag)
        elif ell >= 1.80 * center:
            relaxed.append(tag)
        else:
            moderate.append(tag)
    return {"focused": focused, "relaxed": relaxed, "moderate": moderate}


def _tail_primitive_shifts(namespace: dict, residual: np.ndarray) -> tuple[list[dict], dict[str, int]]:
    """Compare only the high-error tails (top 1/3 on each sign) to avoid near-zero jitter noise."""
    pos_idx = np.where(residual > 0)[0]
    neg_idx = np.where(residual < 0)[0]
    if len(pos_idx) < 3 or len(neg_idx) < 3:
        return [], {"underestimated_tail": 0, "overestimated_tail": 0}
    k_pos = max(2, len(pos_idx) // 3)
    k_neg = max(2, len(neg_idx) // 3)
    pos_tail = pos_idx[np.argsort(-residual[pos_idx], kind="stable")[:k_pos]]
    neg_tail = neg_idx[np.argsort(residual[neg_idx], kind="stable")[:k_neg]]

    shifts = []
    for pname in primitives.PRIMITIVE_METADATA:
        vals = np.asarray(namespace[pname], dtype=float)
        std = float(np.std(vals))
        if std <= 1e-12:
            continue
        pos_mean = float(np.mean(vals[pos_tail]))
        neg_mean = float(np.mean(vals[neg_tail]))
        effect_z = (pos_mean - neg_mean) / std
        if abs(effect_z) < SHIFT_MIN_EFFECT_Z:
            continue
        shifts.append({
            "primitive": pname,
            "underestimated_tail_mean": round(pos_mean, 4),
            "overestimated_tail_mean": round(neg_mean, 4),
            "standardized_shift": round(effect_z, 2),
            "strength": _shift_strength_label(effect_z),
        })
    shifts.sort(key=lambda d: -abs(d["standardized_shift"]))
    return shifts[:4], {"underestimated_tail": int(k_pos), "overestimated_tail": int(k_neg)}


def response_context(rows, x, y, namespace, formulas, baseline, *, tried_history=None):
    residual = np.asarray(baseline["residual_m"], dtype=float)
    order = np.argsort(-np.abs(residual), kind="stable")
    hard = order[:TOP_COUNT]

    # 用已揭晓前缀内的归一化 8D 最近邻距离标记采样稀疏程度。
    dists = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=-1)
    np.fill_diagonal(dists, np.inf)
    nn_dist = np.min(dists, axis=1)
    sparse_cutoff = float(np.quantile(nn_dist, 0.85))

    pos_mask = residual > 0
    neg_mask = residual < 0
    directions = np.select([pos_mask, neg_mask], ["underestimated", "overestimated"], default="zero")
    group_counts = {
        "underestimated": int(np.count_nonzero(pos_mask)),
        "overestimated": int(np.count_nonzero(neg_mask)),
        "zero": int(np.count_nonzero(residual == 0)),
    }
    primitive_shifts, tail_counts = _tail_primitive_shifts(namespace, residual)

    active_columns: list[tuple[str, np.ndarray]] = []
    for f in formulas:
        try:
            active_columns.append((f["name"], features.feature_values(f, namespace, len(x))))
        except features.FORMULA_ERRORS:
            continue

    sse = float(np.sum(residual ** 2))
    top3_sse_share_pct = round(100.0 * float(np.sum(residual[hard] ** 2)) / sse) if sse > 0 else 0
    top3_sparse_count = sum(1 for idx in hard if nn_dist[idx] >= sparse_cutoff)

    top_cases = []
    for idx in hard:
        top_cases.append({
            "sample_index": int(rows[idx]["sample_index"]),
            "observed_mm": round(float(y[idx]) * 1000.0, 2),
            "cluster_pred_mm": round(float(baseline["predicted_m"][idx]) * 1000.0, 2),
            "residual_mm": round(float(residual[idx]) * 1000.0, 2),
            "residual_direction": str(directions[idx]),
            "nearest_neighbor_dist_8d": round(float(nn_dist[idx]), 2),
            "is_sparse_in_prefix": bool(nn_dist[idx] >= sparse_cutoff),
            "ratios": _case_dimensionless_ratios(
                namespace, idx, residual_val=float(residual[idx]), observed_val=float(y[idx]),
            ),
            "key_primitives": {
                k: round(float(namespace[k][idx]), 3)
                for k in ("A_sec", "A_mod", "A_ds_shell", "A_drain", "A_tot", "B_bot", "d_us_min", "d_ds_perp", "Y_c_sec", "X_c_sec")
                if k in namespace
            },
        })

    return {
        "current_formulas": [{"name": f["name"], "expression": f["expression"]} for f in formulas],
        "tried_formulas": list(tried_history or ()),
        "cluster_rmse_mm": round(float(baseline["rmse_m"]) * 1000.0, 2),
        "cluster_mae_mm": round(float(baseline["mae_m"]) * 1000.0, 2),
        "cluster_nrmse": round(float(baseline["nrmse"]), 3),
        "cluster_nlpd": round(float(baseline["nlpd"]), 3),
        "prequential_nrmse": round(float(baseline["prequential"]["nrmse"]), 3) if baseline.get("prequential") else None,
        "map_score": round(float(baseline["map_score"]), 2) if "map_score" in baseline else None,
        "physical_signal_share_mean": round(float(baseline.get("physical_signal_share_mean", 0.0)), 2),
        "residual_ard_summary": _residual_ard_summary(baseline.get("length_scales")),
        "residual_group_counts": group_counts,
        "tail_group_counts": tail_counts,
        "top3_sse_share_pct": int(top3_sse_share_pct),
        "top3_sparse_count": int(top3_sparse_count),
        "top_residual_signals": _residual_signal_correlations(
            namespace, x, residual, active_columns=active_columns,
        ),
        "top_group_primitive_shifts": primitive_shifts,
        "top_residual_cases": top_cases,
    }


def make_context(rows, x, y, namespace, active, baselines, *, target_responses=None, tried_history=None):
    targets = list(gp.RESPONSES if target_responses is None else target_responses)
    history = tried_history or {}
    return {
        "sample_count": len(rows), "response_unit": "mm",
        "residual_definition": "observed minus out-of-cluster Ward spatial CV prediction",
        "top_residual_case_limit": TOP_COUNT,
        "target_responses": targets,
        "unsafe_standalone_denominators": sorted(candidates.UNSAFE_STANDALONE_DENOMINATORS),
        "original_variable_units": {name: "m" if candidates.SYMBOL_UNITS[name] == 1 else "1"
                                    for name in gp.VARIABLES},
        "aliases": candidates.ALIAS_TO_CANONICAL,
        "primitive_catalog": primitives.PRIMITIVE_METADATA,
        "responses": {
            response: response_context(
                rows, x, y[:, j], namespace, active[response], baselines[response],
                tried_history=history.get(response),
            )
            for j, response in enumerate(gp.RESPONSES)
        },
    }


def _online_anomaly_reasons(context: dict, resp: str) -> list[str]:
    """Return any online blind-window anomaly reasons for a response."""
    recent = (context.get("trigger") or {}).get("recent_online") or {}
    z_vals = (recent.get("z_scores") or {}).get(resp, [])
    rmse_pair = (recent.get("rmse_dynamic_vs_standard_m") or {}).get(resp)
    share_vals = (recent.get("physical_signal_shares") or {}).get(resp, [])
    has_active = bool((context.get("responses") or {}).get(resp, {}).get("current_formulas"))
    reasons = []
    if z_vals:
        max_abs_z = max(abs(z) for z in z_vals)
        if max_abs_z > ONLINE_Z_THRESHOLD:
            reasons.append(f"近 {len(z_vals)} 步盲测出现校准异常 (max |z| = {max_abs_z:.2f} > {ONLINE_Z_THRESHOLD:.1f})")
    if rmse_pair and rmse_pair[1] > 0 and rmse_pair[0] > rmse_pair[1] * ONLINE_RMSE_RATIO:
        ratio = rmse_pair[0] / rmse_pair[1]
        reasons.append(f"近窗盲测 RMSE 落后同构 8D 基线 ({ratio:.2f}x > {ONLINE_RMSE_RATIO:.2f}x)")
    if has_active and len(share_vals) >= ONLINE_TRIGGER_WINDOW and all(s < LOW_SHARE_PRUNE_THRESHOLD for s in share_vals):
        reasons.append(f"在役物理核近 {len(share_vals)} 步方差占比持续低于 {LOW_SHARE_PRUNE_THRESHOLD:.2f}，触发剪枝或替换评估")
    return reasons


def _format_trigger_note(context: dict, resp: str) -> str:
    """Format trigger provenance (milestone vs recent blind anomaly) for one target response."""
    trigger = context.get("trigger") or {}
    mode = trigger.get("mode", "milestones")
    recent = trigger.get("recent_online") or {}
    resp_indices = (recent.get("sample_indices_by_response") or {}).get(resp) or recent.get("sample_indices") or []
    reasons = []
    if mode != "anomaly" and context.get("sample_count") in TRIGGER_MILESTONES:
        reasons.append(f"到达里程碑节点 t={context['sample_count']}")
    reasons.extend(_online_anomaly_reasons(context, resp))
    win_str = f"（窗口样本 #{resp_indices[0]}–#{resp_indices[-1]}）" if resp_indices else ""
    if not reasons:
        return f"触发来源：`{mode}` 例行检查{win_str}"
    return f"触发来源：`{mode}`（{'；'.join(reasons)}）{win_str}"


def format_diagnostic_brief(context: dict) -> str:
    """Render context into compact Markdown tables so the LLM avoids raw float overload."""
    targets = context.get("target_responses") or context.get("trigger", {}).get("target_responses") or list(gp.RESPONSES)
    trigger_mode = (context.get("trigger") or {}).get("mode", "milestones")
    lines = [
        f"## 当前已观测前缀摘要（样本数 t = {context['sample_count']}，触发模式 `{trigger_mode}`）",
        "",
        f"- **本轮需处理的目标响应**：`{', '.join(targets)}`（仅需在 JSON 中输出这些键）",
        "- **残差定义**：`真值 - 空间聚簇外推预测`（正值 = `低估`，负值 = `高估`；Top-3 画像以相对偏差%展示）",
        "",
        "### 允许使用的变量与几何原语（精简速查表）",
        "| 符号 / 别名 | 单位 | 含义 |",
        "|---|---|---|",
    ]
    var_descs = {
        "upstream_slope": "上游坡比 H/V (`x1`)",
        "downstream_slope": "下游坡比 H/V (`x2`)",
        "crest_width": "坝顶宽度 (`x3`)",
        "modulus_zone_height": "增模区高度 (`x4`)",
        "secondary_top_width": "次堆石区顶宽 (`x5`)",
        "secondary_height": "次堆石区高度 (`x6`)",
        "secondary_bottom_width": "次堆石区底宽 (`x7`)",
        "secondary_right_slope_distance": "次堆石右腰距下游参考线法向距离 (`x8` / `d_ds_perp`)",
    }
    alias_for = {canonical: alias for alias, canonical in context["aliases"].items() if alias.startswith("x")}
    for name, unit in context["original_variable_units"].items():
        alias = alias_for.get(name, "")
        label = f"`{name}` (`{alias}`)" if alias else f"`{name}`"
        lines.append(f"| {label} | {unit} | {var_descs.get(name, name)} |")
    for pname, meta in context["primitive_catalog"].items():
        lines.append(f"| `{pname}` | {meta['unit']} | {meta['desc']} |")

    for resp in targets:
        rdata = context["responses"][resp]
        counts = rdata["residual_group_counts"]
        current_str = ", ".join(f"`{f['name']} = {f['expression']}`" for f in rdata["current_formulas"]) or "无（当前为同先验非对称正则化 8D 基线 GP）"
        signals = rdata.get("top_residual_signals", [])
        shifts = rdata.get("top_group_primitive_shifts", [])
        max_rho = max((abs(s["rho_residual"]) for s in signals), default=0.0)
        phys_share = float(rdata.get("physical_signal_share_mean", 0.0))
        top3_sse_pct = int(rdata.get("top3_sse_share_pct", 0))
        top3_sparse = int(rdata.get("top3_sparse_count", 0))
        sparse_dominated = top3_sse_pct >= SPARSE_SSE_SHARE_PCT and top3_sparse >= 2
        no_cross_signal = max_rho < MIN_RESIDUAL_SIGNAL_ABS_RHO and not shifts
        high_share_quiet = phys_share >= 0.50 and max_rho < HIGH_SHARE_WEAK_RHO_CUTOFF
        anomaly_flagged = bool(_online_anomaly_reasons(context, resp))

        if high_share_quiet and not anomaly_flagged:
            health_note = (
                f"提示：当前在役物理特征已承担主要信号方差（占比 {phys_share:.2f}），且残差中无强全局相关性线索（|ρ| < {HIGH_SHARE_WEAK_RHO_CUTOFF:.2f}），"
                "建议直接返回 `[]` 保持现状，避免强行叠特征或替换破坏已有机制。"
            )
        elif high_share_quiet and anomaly_flagged:
            health_note = (
                f"提示：当前在役特征占比已达 {phys_share:.2f} 且全局相关较弱（|ρ| < {HIGH_SHARE_WEAK_RHO_CUTOFF:.2f}），但近期盲测窗口触发了在线异常；"
                "请优先检查尾部几何偏移是否指向缺失的补充/替代机理，若无可靠机理仍应返回 `[]` 保护在役特征。"
            )
        elif no_cross_signal:
            sparse_desc = (
                f"，且高残差高度集中在个别稀疏边界点（Top-3 占总平方误差约 {top3_sse_pct}%），更像是局部采样稀疏而非全局缺失规律"
                if sparse_dominated else ""
            )
            health_note = (
                f"提示：当前全局残差无显著单调线索（|ρ| < {MIN_RESIDUAL_SIGNAL_ABS_RHO:.2f}）且显著高低估尾部无明显几何偏移{sparse_desc}，"
                "建议直接返回 `[]` 保持现状，切勿硬凑过拟合公式。"
            )
        elif sparse_dominated:
            health_note = (
                f"提示：高残差集中于稀疏边界断面（Top-3 占总平方误差约 {top3_sse_pct}%），但仍存在跨样本残差线索；"
                "请基于全样本单调/偏移方向提出平滑全局无量纲比例，禁止专为孤立点构造特异公式。"
            )
        else:
            health_note = "提示：存在显著的跨断面残差相关性或高低估尾部几何偏移，可提出 1 ~ 2 条机理互补的候选公式（若无把握也可返回 `[]`）。"

        preq_str = f"，因果前向 NRMSE = `{rdata['prequential_nrmse']:.3f}`" if rdata.get("prequential_nrmse") is not None else ""
        lines.extend([
            "",
            f"### 目标响应：`{resp}`",
            f"- **{_format_trigger_note(context, resp)}**",
            f"- **当前在役公式**：{current_str}",
            f"- **模型健康度**：Clust-8 NRMSE = `{rdata['cluster_nrmse']:.3f}`{preq_str}，物理核方差占比 = `{phys_share:.2f}`，Top-3 误差集中度 = `{top3_sse_pct}%`（含 `{top3_sparse}/3` 稀疏边界点）",
            f"- **全量残差方向分布**：低估 `{counts['underestimated']}` 例 / 高估 `{counts['overestimated']}` 例 / 零 `{counts['zero']}` 例",
            f"- **{health_note}**",
        ])

        ard = rdata.get("residual_ard_summary")
        if ard:
            ard_parts = []
            if ard.get("focused"):
                ard_parts.append(f"重点拟合维度：`{', '.join(ard['focused'])}`")
            if ard.get("moderate"):
                ard_parts.append(f"中等敏感维度：`{', '.join(ard['moderate'])}`")
            if ard.get("relaxed"):
                ard_parts.append(f"已放宽维度：`{', '.join(ard['relaxed'])}`")
            lines.append(f"- **8D 残差核 ARD 敏感度概况**：{'；'.join(ard_parts)}")

        tried = rdata.get("tried_formulas", [])
        if tried:
            lines.append("- **历史已试错/淘汰公式（严禁重提等价表达式）**：")
            for item in tried:
                lines.append(f"  - `{item['name']} = {item['expression']}` → {item['outcome']}")

        if signals:
            lines.extend([
                "",
                "**显著全局空间聚簇残差相关性线索（Python 预计算，仅展示 |ρ| ≥ 0.20）**：",
                "| 物理含义 | 预计算无量纲式 | ρ_res | 强度 | 方向提示 | 在役共线检查 |",
                "|---|---|---:|:---:|---|---|",
            ])
            for sig in signals:
                lines.append(
                    f"| {sig['label']} | `{sig['expression']}` | {sig['rho_residual']:+.2f} | "
                    f"{sig.get('strength', '有效')} | {sig['direction']} | {sig.get('collinearity_note', '无')} |"
                )
        else:
            lines.append(f"- **全局残差相关性线索**：全局无显著残差相关性线索（预计算无量纲比值均满足 |ρ| < {MIN_RESIDUAL_SIGNAL_ABS_RHO:.2f}）。")

        if shifts:
            tail_counts = rdata.get("tail_group_counts", {})
            n_u = tail_counts.get("underestimated_tail", "?")
            n_o = tail_counts.get("overestimated_tail", "?")
            shift_parts = [
                f"`{s['primitive']}` ({s['standardized_shift']:+.2f}σ, {s['strength']}, 显著低估组均值{'更大' if s['standardized_shift'] > 0 else '更小'})"
                for s in shifts
            ]
            lines.append(f"- **显著高误差子集原语偏移（低估端前 {n_u} 例 vs 高估端前 {n_o} 例，仅展示 |z| ≥ {SHIFT_MIN_EFFECT_Z:.1f}）**：{'；'.join(shift_parts)}")
        else:
            lines.append(f"- **显著高误差子集原语偏移**：正负显著误差组之间无明显几何原语偏移（|effect_z| < {SHIFT_MIN_EFFECT_Z:.1f}）。")

        cases = rdata.get("top_residual_cases", [])
        if cases:
            lines.extend([
                "",
                "**Top-3 高残差断面归一化画像（无需心算比例）**：",
                "| 样本 | 相对偏差 | 稀疏点 | 软区 | 帽盖 | 排水 | 下游壳 | 软区高程 | 偏心率 | 上游净宽 | 下游净距 |",
                "|---:|---:|:---:|---:|---:|---:|---:|---:|---:|---:|---:|",
            ])
            for c in cases:
                r = c.get("ratios", _case_dimensionless_ratios_fallback(c))
                dir_cn = "低估" if c["residual_direction"] == "underestimated" else "高估" if c["residual_direction"] == "overestimated" else "零"
                sparse_cn = "是" if c["is_sparse_in_prefix"] else "否"
                rel_err = r.get("rel_err_pct")
                if rel_err is None:
                    obs = c.get("observed_mm", 0.0)
                    rel_err = round(100.0 * c.get("residual_mm", 0.0) / obs) if abs(obs) > 1e-9 else 0
                lines.append(
                    f"| #{c['sample_index']} | {rel_err:+.0f}% ({dir_cn}) | {sparse_cn} | "
                    f"{r['soft_pct']:.0f}% | {r['mod_pct']:.0f}% | {r['drain_pct']:.0f}% | {r.get('ds_shell_pct', 0.0):.0f}% | "
                    f"{r['rel_elev']:.2f}H | {r['ecc_ratio']:+.2f}B | {r['us_clear_ratio']:.2f}B | {r.get('ds_clear_ratio', 0.0):.2f}B |"
                )

    return "\n".join(lines)


def _case_dimensionless_ratios_fallback(case: dict) -> dict[str, float]:
    """Compute ratios from key_primitives if a legacy context lacks pre-computed ratios."""
    kp = case.get("key_primitives", {})
    a_tot = kp.get("A_tot", 1.0) or 1.0
    b_bot = kp.get("B_bot", 1.0) or 1.0
    obs = case.get("observed_mm", 0.0)
    res = case.get("residual_mm", 0.0)
    return {
        "soft_pct": round(100.0 * kp.get("A_sec", 0.0) / a_tot, 1),
        "mod_pct": round(100.0 * kp.get("A_mod", 0.0) / a_tot, 1),
        "drain_pct": round(100.0 * kp.get("A_drain", 0.0) / a_tot, 1),
        "ds_shell_pct": round(100.0 * kp.get("A_ds_shell", 0.0) / a_tot, 1),
        "rel_elev": round(kp.get("Y_c_sec", 0.0) / 74.0, 2),
        "ecc_ratio": round(kp.get("X_c_sec", 0.0) / b_bot, 2),
        "us_clear_ratio": round(kp.get("d_us_min", 0.0) / b_bot, 2),
        "ds_clear_ratio": round(kp.get("d_ds_perp", 0.0) / b_bot, 2),
        "rel_err_pct": round(100.0 * res / obs, 1) if abs(obs) > 1e-9 else 0.0,
    }


def render_request(path: Path, context: dict, *, signature: str) -> None:
    template = PROMPT_TEMPLATE.read_text(encoding="utf-8")
    unsafe_list = ", ".join(context.get("unsafe_standalone_denominators") or sorted(candidates.UNSAFE_STANDALONE_DENOMINATORS))
    if BRIEF_TOKEN in template:
        rendered = template.replace(BRIEF_TOKEN, format_diagnostic_brief(context))
        if DENOMINATOR_TOKEN in rendered:
            rendered = rendered.replace(DENOMINATOR_TOKEN, unsafe_list)
    elif CONTEXT_TOKEN in template:
        payload = json.dumps(context, ensure_ascii=False, indent=2, allow_nan=False)
        rendered = template.replace(CONTEXT_TOKEN, payload)
    else:
        raise ValueError("evolution prompt must contain a diagnostic brief or context placeholder")
    prompt = rendered.encode("utf-8")
    write_json(path / "context.json", context)
    write_bytes(path / "prompt.md", prompt)
    write_json(path / "request_manifest.json", {
        "schema_version": REQUEST_SCHEMA_VERSION, "stage": "feature_evolution", "sample_count": context["sample_count"],
        "signature": signature,
        "context_sha256": gp.canonical_sha256(context), "prompt_sha256": gp.file_sha256(path / "prompt.md"),
    })


def validate_reply(raw: bytes, x, namespace, active, *, scaling, target_responses: Sequence[str] | None = None,
                   allow_empty: bool = True):
    targets = tuple(gp.RESPONSES) if target_responses is None else tuple(target_responses)
    options = {response: [] for response in gp.RESPONSES}
    audit = {response: [] for response in gp.RESPONSES}
    try:
        reply = prompt_io.parse_llm_json(raw)
        if not isinstance(reply, dict) or set(reply) != set(targets):
            raise ValueError(f"reply must contain exactly the target response keys: {list(targets)}")
    except (UnicodeError, ValueError) as exc:
        return options, audit, str(exc)

    basis = np.column_stack([np.ones(len(x)), x - x.mean(axis=0)])
    for response in targets:
        rows = reply[response]
        if not isinstance(rows, list) or len(rows) > PROPOSAL_COUNT or (not allow_empty and len(rows) != PROPOSAL_COUNT):
            msg = f"expected 0 to {PROPOSAL_COUNT} candidate formulas" if allow_empty else f"expected exactly {PROPOSAL_COUNT} candidate formulas"
            audit[response].append({"status": "rejected", "reason": msg})
            continue
        if not rows:
            audit[response].append({"status": "abstained", "reason": "llm_returned_empty_list"})
            continue
        ref_namespace = scaling.namespace
        names = {item["name"] for item in active[response]}
        for position, proposed in enumerate(rows, 1):
            check = {
                "position": position,
                "name": proposed.get("name") if isinstance(proposed, dict) else None,
                "expression": proposed.get("expression") if isinstance(proposed, dict) and isinstance(proposed.get("expression"), str) else None,
            }
            try:
                feature, _ = candidates.validate_candidate(proposed, namespace, basis, ref_namespace=ref_namespace)
                column = np.asarray(feature["values"])
                if not np.isfinite(column.std()) or column.std() <= features.std_floor(column):
                    raise ValueError("candidate is near constant on the revealed prefix")
                if feature["name"] in names:
                    raise ValueError("candidate name duplicates an existing or proposed feature")
                formula = {field: feature[field] for field in features.FORMULA_FIELDS}
                column = scaling.transform(formula, column)
                if column.std() <= features.std_floor(column):
                    raise ValueError("scaled candidate is near constant on the revealed prefix")
                names.add(feature["name"])
                options[response].append((formula, column))
                check.update(status="valid", formula=formula, unit=feature["unit"],
                             reference_bounds=scaling.bounds(formula))
            except features.FORMULA_ERRORS as exc:
                check.update(status="rejected", reason=str(exc))
            audit[response].append(check)
    return options, audit, None


@dataclass(frozen=True)
class ApiSettings:
    provider: str
    url: str
    model: str
    style: str
    key_env: str
    timeout: float


def request_reply(path: Path, offline_file: Path | None, api: ApiSettings | None):
    """Use a saved reply, wait for an offline reply, or send one explicitly configured API request."""
    raw_path = path / "reply.json"
    call_path = path / "api_call.json"
    provider_path = path / "provider_response.json"
    if raw_path.is_file():
        return raw_path.read_bytes(), {"status": "reply_received", "source": "saved"}
    if offline_file is not None and offline_file.is_file():
        raw = offline_file.read_bytes()
        write_bytes(raw_path, raw)
        return raw, {"status": "reply_received", "source": "offline"}
    if call_path.is_file():
        call = gp.read_json(call_path)
        if call["status"] == "error":
            return None, {"status": "provider_error", "error": call["error"]}
        if not provider_path.is_file():
            return None, {"status": "awaiting_reply", "reason": "an earlier API call has no saved response"}
    elif api is None:
        return None, {"status": "awaiting_reply"}
    else:
        key = os.environ.get(api.key_env)
        if not key:
            raise ValueError(f"set {api.key_env} to use live feature proposals")
        prompt = (path / "prompt.md").read_text(encoding="utf-8")
        body = llm.request_body(api.model, prompt, api.style, {})
        call = {"status": "pending", "provider": api.provider, "model": api.model, "style": api.style,
                "endpoint": api.url, "request_sha256": gp.canonical_sha256(body)}
        # The marker prevents duplicate API requests during a resume.
        with call_path.open("x", encoding="utf-8") as stream:
            json.dump(call, stream, allow_nan=False)
        try:
            reply, http_status = llm.post_json(api.url, body, key, api.timeout)
            write_json(provider_path, reply)
            call["http_status"] = http_status
        except llm.APIError as exc:
            write_json(call_path, {**call, "status": "error", "error": str(exc)})
            return None, {"status": "provider_error", "error": str(exc)}
    try:
        answer, provider_id, usage, model = llm.extract_answer(gp.read_json(provider_path), call["style"])
    except llm.APIError as exc:
        write_json(call_path, {**call, "status": "error", "error": str(exc)})
        return None, {"status": "provider_error", "error": str(exc)}
    raw = answer.encode("utf-8")
    write_bytes(raw_path, raw)
    write_json(call_path, {**call, "status": "complete", "provider_id": provider_id,
                           "provider_model": model, "usage": usage})
    return raw, {"status": "reply_received", "source": "live"}
