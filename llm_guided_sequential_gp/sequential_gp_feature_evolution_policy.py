"""Compare proposed feature changes with the incumbent using Stage-B 4D GP validation."""

from __future__ import annotations

import numpy as np
from scipy.stats import rankdata

if __package__:
    from . import sequential_gp_additive as additive
    from . import sequential_gp_feature_evolution as evolution
    from . import sequential_gp_feature_selection as features
    from . import sequential_gp_model as gp
    from . import sequential_gp_validation as validation
else:
    import sequential_gp_additive as additive
    import sequential_gp_feature_evolution as evolution
    import sequential_gp_feature_selection as features
    import sequential_gp_model as gp
    import sequential_gp_validation as validation

# === Stage-B 四维联合门控常量 ===
MIN_GAIN_REPLACE = 0.02        # 替换已有特征（维度不变）：Clust+Preq 综合相对增益 >= 2%
MIN_GAIN_ADD = 0.03            # 新增特征（维度 +1）：Clust+Preq 综合相对增益 >= 3%
MIN_GAIN_PRUNE = 0.02          # 剪枝冗余特征：Clust+Preq 综合净改善 >= 2%
MIN_GAIN_CLEAR = 0.04          # 清零回退至 8D：提高迟滞门槛防局部波动误杀
MIN_GAIN_CLEAR_HIGH_SHARE = 0.06
HIGH_PHYSICAL_SHARE_GUARD = 0.50
MIN_DELTA_MAP_ADD = 0.50       # 新增特征必须带来 >= 0.50 nats 的净 MAP 证据增益
MIN_DELTA_MAP_REPLACE = 0.25   # 替换特征需 >= 0.25 nats 的 MAP 证据增益
MAX_PRUNE_MAP_DROP = 0.25      # 常规剪枝允许的最大 MAP 下降；高物理占比清零要求 delta_map >= 0
MAX_ANA_REGRESSION = 0.01      # 解析 GP LOO 允许的微小误差容忍
MAX_EXPLORE_REGRESSION = 0.02  # 稀疏探索点（explore_steps >= 4）不允许显著退化
MAX_WARNING_RATIO = 1.25
MIN_WARNING_ALLOWANCE = 5
MAX_NLPD_INCREASE = 0.02
MIN_GAIN_VS_STANDARD = 0.02


def gate_settings(min_gain: float) -> dict:
    return {
        "validation_protocol": "stage_b_four_dimensional",
        "min_composite_gain": {
            "add": max(MIN_GAIN_ADD, min_gain),
            "replace": max(MIN_GAIN_REPLACE, min_gain),
            "prune": MIN_GAIN_PRUNE,
            "clear": MIN_GAIN_CLEAR,
            "clear_high_share": MIN_GAIN_CLEAR_HIGH_SHARE,
        },
        "min_delta_map": {"add": MIN_DELTA_MAP_ADD, "replace": MIN_DELTA_MAP_REPLACE},
        "max_prune_map_drop": MAX_PRUNE_MAP_DROP,
        "high_physical_share_guard": HIGH_PHYSICAL_SHARE_GUARD,
        "min_physical_share": additive.MIN_PHYSICAL_SIGNAL_SHARE,
        "max_active_physical_length_scale": additive.MAX_ACTIVE_PHYSICAL_LENGTH_SCALE,
        "max_nlpd_increase": MAX_NLPD_INCREASE,
        "max_ana_regression": MAX_ANA_REGRESSION,
        "max_explore_regression": MAX_EXPLORE_REGRESSION,
        "cluster_count": validation.DEFAULT_CLUSTER_COUNT,
        "prequential_max_steps": validation.DEFAULT_PREQ_MAX_STEPS,
        "warning_metric": "optimizer_failure_count",
        "max_warning_ratio": MAX_WARNING_RATIO,
        "min_warning_allowance": MIN_WARNING_ALLOWANCE,
        "min_composite_gain_vs_standard": MIN_GAIN_VS_STANDARD,
        "redundancy_abs_spearman": features.REDUNDANCY_LIMIT,
        "feature_cap": evolution.FEATURE_CAP,
    }


def incumbent_baseline(x, y, namespace, active, scales, *, restarts, seed, scaling):
    results = {}
    for j, response in enumerate(gp.RESPONSES):
        fallback = None
        try:
            values = additive.scaled_physical_values(active[response], namespace, len(x), scaling=scaling)
        except features.FORMULA_ERRORS as exc:
            values = np.empty((len(x), 0))
            fallback = str(exc)
        score = additive.evaluate_stage_b(
            x, y[:, j], values, response, float(scales[j]), restarts=restarts, seed=seed + j,
        )
        results[response] = {**score, "formula_fallback": fallback}
    return results


def spearman(a, b):
    ranks = rankdata(np.column_stack((a, b)), axis=0)
    centered = ranks - ranks.mean(axis=0)
    length = np.linalg.norm(centered, axis=0)
    if length.min() == 0:
        return None
    return float(np.clip((centered[:, 0] @ centered[:, 1]) / np.prod(length), -1, 1))


def remove_redundant_proposals(options, residual):
    if len(options) < 2:
        return options, []
    first, second = options
    rho = spearman(first[1], second[1])
    if rho is None or abs(rho) <= features.REDUNDANCY_LIMIT:
        return options, []
    residual_rhos = [spearman(values, residual) for _, values in options]
    winner = max(range(len(options)), key=lambda i: abs(residual_rhos[i] or 0.0))
    loser = 1 - winner
    return [options[winner]], [{
        "name": options[loser][0]["name"], "status": "redundant_proposal",
        "retained_instead": options[winner][0]["name"], "candidate_pair_rho": rho,
        "rho_residual": residual_rhos[loser], "retained_rho_residual": residual_rhos[winner],
    }]


def _composite_nrmse(score: dict) -> float:
    cluster_n = float(score["cluster"]["nrmse"])
    preq = score.get("prequential")
    if preq is None:
        return cluster_n
    return 0.5 * cluster_n + 0.5 * float(preq["nrmse"])


def _evaluate_stage_b_gate(
    score: dict,
    baseline: dict,
    standard_8d: dict,
    *,
    action_type: str,
    min_gain_override: float | None = None,
    clear_all: bool | None = None,
    expected_physical_count: int | None = None,
) -> tuple[bool, dict]:
    """统一执行 Stage-B 四维联合门控：MAP证据 + 解析LOO + 空间聚簇CV + 因果前向Prequential + NLPD校准。"""
    cluster_gain = validation.relative_gain(baseline["cluster"]["nrmse"], score["cluster"]["nrmse"])
    has_preq = bool(score.get("prequential") and baseline.get("prequential"))
    preq_gain = (
        validation.relative_gain(baseline["prequential"]["nrmse"], score["prequential"]["nrmse"])
        if has_preq else cluster_gain
    )
    composite_gain = 0.5 * cluster_gain + 0.5 * preq_gain if (cluster_gain is not None and preq_gain is not None) else cluster_gain
    ana_gain = validation.relative_gain(baseline["ana_loo"]["nrmse"], score["ana_loo"]["nrmse"])

    explore_gain = None
    if has_preq:
        b_exp = baseline["prequential"].get("explore_nrmse")
        c_exp = score["prequential"].get("explore_nrmse")
        if b_exp is not None and c_exp is not None:
            explore_gain = validation.relative_gain(b_exp, c_exp)

    delta_map = float(score["map_score"] - baseline["map_score"])
    delta_lml = float(score["raw_lml"] - baseline["raw_lml"])

    cluster_nlpd_ok = score["cluster"]["nlpd"] <= baseline["cluster"]["nlpd"] + MAX_NLPD_INCREASE
    preq_nlpd_ok = (
        True if not has_preq
        else score["prequential"]["nlpd"] <= baseline["prequential"]["nlpd"] + MAX_NLPD_INCREASE
    )

    warning_metric = "optimizer_failure_count"
    baseline_warnings = int(baseline[warning_metric])
    candidate_warnings = int(score[warning_metric])
    warn_limit = max(baseline_warnings + MIN_WARNING_ALLOWANCE, int(baseline_warnings * MAX_WARNING_RATIO))
    warn_ok = candidate_warnings <= warn_limit

    common_checks = cluster_nlpd_ok and preq_nlpd_ok and warn_ok

    is_clear_all = False
    if action_type == "prune":
        is_clear_all = (score is standard_8d) if clear_all is None else bool(clear_all)
        base_phys_share = float(baseline.get("physical_signal_share_mean", 0.0))
        high_share_incumbent = base_phys_share >= HIGH_PHYSICAL_SHARE_GUARD
        if is_clear_all:
            req_gain = MIN_GAIN_CLEAR_HIGH_SHARE if high_share_incumbent else MIN_GAIN_CLEAR
            map_ok = delta_map >= (0.0 if high_share_incumbent else -MAX_PRUNE_MAP_DROP)
        else:
            req_gain = MIN_GAIN_PRUNE
            map_ok = delta_map >= -MAX_PRUNE_MAP_DROP
        explore_ok = True if explore_gain is None else explore_gain >= 0.0
        axes_ok = (cluster_gain is not None and cluster_gain >= 0.0
                   and preq_gain is not None and preq_gain >= 0.0)
        ard_active = True
        phys_lengths = (score.get("length_scales") or {}).get("physical")
        eligible = common_checks and map_ok and explore_ok and axes_ok
    else:
        req_gain = MIN_GAIN_ADD if action_type == "add" else MIN_GAIN_REPLACE
        if min_gain_override is not None:
            req_gain = max(req_gain, min_gain_override)
        min_dmap = MIN_DELTA_MAP_ADD if action_type == "add" else MIN_DELTA_MAP_REPLACE
        map_ok = delta_map >= min_dmap
        phys_share = float(score.get("physical_signal_share_mean", 0.0))
        phys_lengths = (score.get("length_scales") or {}).get("physical")
        if phys_lengths is None:
            ard_active = True
        else:
            dim_ok = expected_physical_count is None or len(phys_lengths) == expected_physical_count
            ard_active = (
                dim_ok
                and len(phys_lengths) >= 1
                and all(isinstance(v, (int, float)) and np.isfinite(v) and v <= additive.MAX_ACTIVE_PHYSICAL_LENGTH_SCALE
                        for v in phys_lengths)
            )
        std_cluster_gain = validation.relative_gain(standard_8d["cluster"]["nrmse"], score["cluster"]["nrmse"])
        std_preq_gain = (
            validation.relative_gain(standard_8d["prequential"]["nrmse"], score["prequential"]["nrmse"])
            if (score.get("prequential") and standard_8d.get("prequential")) else std_cluster_gain
        )
        std_comp_gain = (
            0.5 * std_cluster_gain + 0.5 * std_preq_gain
            if (std_cluster_gain is not None and std_preq_gain is not None) else None
        )
        beats_standard = (
            std_comp_gain is not None
            and std_comp_gain >= MIN_GAIN_VS_STANDARD
            and score["map_score"] >= standard_8d["map_score"] + MIN_DELTA_MAP_ADD
        )
        explore_ok = True if explore_gain is None else explore_gain >= -MAX_EXPLORE_REGRESSION
        axes_ok = (
            cluster_gain is not None and cluster_gain >= 0.0
            and preq_gain is not None and preq_gain >= 0.0
            and (ana_gain is None or ana_gain >= -MAX_ANA_REGRESSION)
        )
        eligible = (
            common_checks
            and map_ok
            and phys_share >= additive.MIN_PHYSICAL_SIGNAL_SHARE
            and ard_active
            and beats_standard
            and explore_ok
            and axes_ok
        )

    passed = eligible and composite_gain is not None and composite_gain >= req_gain
    diagnostics = {
        "action_type": action_type,
        "clear_all": is_clear_all,
        "required_composite_gain": req_gain,
        "composite_gain": composite_gain,
        "cluster_gain": cluster_gain,
        "prequential_gain": preq_gain,
        "ana_loo_gain": ana_gain,
        "explore_gain": explore_gain,
        "delta_map": delta_map,
        "delta_lml": delta_lml,
        "map_ok": map_ok,
        "cluster_nlpd_ok": cluster_nlpd_ok,
        "preq_nlpd_ok": preq_nlpd_ok,
        "warn_ok": warn_ok,
        "warning_metric": warning_metric,
        "baseline_warning_count": baseline_warnings,
        "candidate_warning_count": candidate_warnings,
        "warning_limit": warn_limit,
        "physical_signal_share_mean": score.get("physical_signal_share_mean"),
        "physical_length_scales": phys_lengths,
        "ard_active": ard_active,
        "standard_8d_composite_nrmse": _composite_nrmse(standard_8d),
        "gate_passed": passed,
    }
    return passed, diagnostics


def compare_response(x, y, namespace, old, options, baseline, response, scale, *, min_gain, restarts, seed, scaling):
    options, redundant = remove_redundant_proposals(options, baseline["residual_m"])
    old_values = {}
    for item in old:
        try:
            column = features.feature_values(item, namespace, len(x))
            old_values[item["name"]] = scaling.transform(item, column)
        except features.FORMULA_ERRORS as exc:
            old_values[item["name"]] = str(exc)

    if not old:
        standard_8d = baseline
    else:
        standard_8d = additive.evaluate_stage_b(
            x, y, np.empty((len(x), 0)), response, scale, restarts=restarts, seed=seed,
        )

    trials = []
    best_eligible = None

    def _register_trial(
        trial: dict,
        score: dict,
        action_type: str,
        override: float | None = None,
        *,
        clear_all: bool = False,
        expected_physical_count: int | None = None,
    ):
        nonlocal best_eligible
        passed, gate = _evaluate_stage_b_gate(
            score, baseline, standard_8d,
            action_type=action_type, min_gain_override=override, clear_all=clear_all,
            expected_physical_count=expected_physical_count,
        )
        trial.update(status="scored", validation=score, relative_gain=gate["composite_gain"], gate=gate)
        if passed and (best_eligible is None or _composite_nrmse(score) < _composite_nrmse(best_eligible["validation"])):
            best_eligible = trial

    if old:
        clear_trial = {
            "added": None, "removed": [item["name"] for item in old],
            "action": "prune_to_standard_8d", "feature_set": [],
        }
        trials.append(clear_trial)
        _register_trial(clear_trial, standard_8d, action_type="prune", clear_all=True)

        if len(old) >= 2:
            for slot in range(len(old)):
                retained = [item for i, item in enumerate(old) if i != slot]
                invalid = [item["name"] for item in retained if isinstance(old_values[item["name"]], str)]
                prune_trial = {
                    "added": None, "removed": old[slot]["name"],
                    "action": "prune_one", "feature_set": retained,
                }
                trials.append(prune_trial)
                if invalid:
                    prune_trial.update(status="invalid_retained_features", reasons={n: old_values[n] for n in invalid})
                    continue
                matrix = np.column_stack([old_values[f["name"]] for f in retained])
                try:
                    score = additive.evaluate_stage_b(x, y, matrix, response, scale, restarts=restarts, seed=seed)
                except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
                    prune_trial.update(status="fit_failed", reason=str(exc))
                    continue
                _register_trial(prune_trial, score, action_type="prune", clear_all=False)

    for candidate, column in options:
        slots = range(len(old)) if len(old) >= evolution.FEATURE_CAP else (None, *range(len(old)))
        for slot in slots:
            retained = [item for i, item in enumerate(old) if i != slot]
            dropped = old[slot]["name"] if slot is not None else None
            proposed = old.copy()
            if slot is None:
                proposed.append(candidate)
                action_type = "add"
            else:
                proposed[slot] = candidate
                action_type = "replace"

            trial = {
                "added": candidate["name"], "removed": dropped,
                "action": action_type, "feature_set": proposed,
            }
            trials.append(trial)
            invalid = [item["name"] for item in retained if isinstance(old_values[item["name"]], str)]
            if invalid:
                trial.update(status="invalid_retained_features", reasons={n: old_values[n] for n in invalid})
                continue
            correlations = {item["name"]: spearman(column, old_values[item["name"]]) for item in retained}
            trial["rho_retained_features"] = correlations
            if any(rho is not None and abs(rho) > features.REDUNDANCY_LIMIT for rho in correlations.values()):
                trial.update(status="redundant_with_retained_feature")
                continue
            matrix = np.column_stack([
                column if f["name"] == candidate["name"] else old_values[f["name"]] for f in proposed
            ])
            try:
                score = additive.evaluate_stage_b(x, y, matrix, response, scale, restarts=restarts, seed=seed)
            except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
                trial.update(status="fit_failed", reason=str(exc))
                continue
            _register_trial(
                trial, score, action_type=action_type, override=min_gain,
                expected_physical_count=matrix.shape[1],
            )

    adopted = best_eligible is not None
    scored_trials = [t for t in trials if t.get("status") == "scored"]
    best_overall = (
        best_eligible if adopted
        else min(scored_trials, key=lambda t: _composite_nrmse(t["validation"]), default=None)
    )
    decision = "adopted" if adopted else "below_threshold" if scored_trials else "no_eligible_candidate"
    return {
        "decision": decision,
        "shadow_features": old,
        "active_features": best_eligible["feature_set"] if adopted else old,
        "shadow_validation": baseline,
        "standard_8d_validation": standard_8d,
        "best_trial": best_overall,
        "redundant_proposals": redundant,
        "trials": trials,
    }


def decide(x, y, namespace, active, options, baselines, scales, *, min_gain, restarts, seed, scaling,
           target_responses=None):
    targets = set(gp.RESPONSES if target_responses is None else target_responses)
    decisions = {}
    for j, response in enumerate(gp.RESPONSES):
        if response not in targets:
            decisions[response] = {
                "decision": "not_triggered", "shadow_features": active[response],
                "active_features": active[response], "shadow_validation": baselines[response],
                "standard_8d_validation": None, "best_trial": None, "redundant_proposals": [], "trials": [],
            }
            continue
        decisions[response] = compare_response(
            x, y[:, j], namespace, active[response], options[response], baselines[response], response,
            float(scales[j]), min_gain=min_gain, restarts=restarts, seed=seed + j, scaling=scaling,
        )
    return decisions
