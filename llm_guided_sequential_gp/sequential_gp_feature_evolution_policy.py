"""Compare proposed feature changes with the incumbent using matched additive-GP LOO."""

from __future__ import annotations

import numpy as np
from scipy.stats import rankdata

if __package__:
    from . import sequential_gp_additive as additive
    from . import sequential_gp_feature_evolution as evolution
    from . import sequential_gp_feature_selection as features
    from . import sequential_gp_model as gp
else:
    import sequential_gp_additive as additive
    import sequential_gp_feature_evolution as evolution
    import sequential_gp_feature_selection as features
    import sequential_gp_model as gp


def incumbent_loo(x, y, namespace, active, scales, *, restarts, seed):
    results = {}
    for j, response in enumerate(gp.RESPONSES):
        fallback = None
        try:
            values = additive.physical_values(active[response], namespace, len(x))
        except features.FORMULA_ERRORS as exc:
            values = np.empty((len(x), 0))
            fallback = str(exc)
        score = additive.loo_response(x, y[:, j], values, response, float(scales[j]),
                                      restarts=restarts, seed=seed + j)
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


def compare_response(x, y, namespace, old, options, baseline, response, scale, *, min_gain, restarts, seed):
    options, redundant = remove_redundant_proposals(options, baseline["residual_m"])
    old_values = {}
    for item in old if options else []:
        try:
            old_values[item["name"]] = features.feature_values(item, namespace, len(x))
        except features.FORMULA_ERRORS as exc:
            old_values[item["name"]] = str(exc)

    trials = []
    best = None
    for candidate, column in options:
        # A trigger changes at most one feature per constraint response.
        slots = range(len(old)) if len(old) == evolution.FEATURE_CAP else (None,)
        for slot in slots:
            retained = [item for i, item in enumerate(old) if i != slot]
            dropped = old[slot]["name"] if slot is not None else None
            proposed = old.copy()
            if slot is None:
                proposed.append(candidate)
            else:
                proposed[slot] = candidate
            trial = {"added": candidate["name"], "removed": dropped, "feature_set": proposed}
            trials.append(trial)
            invalid = [item["name"] for item in retained if isinstance(old_values[item["name"]], str)]
            if invalid:
                trial.update(status="invalid_retained_features", reasons={name: old_values[name] for name in invalid})
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
                score = additive.loo_response(x, y, matrix, response, scale, restarts=restarts, seed=seed)
            except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
                trial.update(status="fit_failed", reason=str(exc))
                continue
            gain = features.relative_gain(baseline["nrmse"], score["nrmse"])
            trial.update(status="scored", loo=score, relative_gain=gain)
            if best is None or score["nrmse"] < best["loo"]["nrmse"]:
                best = trial

    adopted = (best is not None and baseline["nrmse"] > 0
               and best["loo"]["nrmse"] <= baseline["nrmse"] * (1 - min_gain))
    decision = "adopted" if adopted else "below_threshold" if best is not None else "no_eligible_candidate"
    return {
        "decision": decision, "shadow_features": old,
        "active_features": best["feature_set"] if adopted else old,
        "shadow_loo": baseline, "best_trial": best, "redundant_proposals": redundant, "trials": trials,
    }


def decide(x, y, namespace, active, options, baselines, scales, *, min_gain, restarts, seed):
    return {
        response: compare_response(
            x, y[:, j], namespace, active[response], options[response], baselines[response], response,
            float(scales[j]), min_gain=min_gain, restarts=restarts, seed=seed + j,
        )
        for j, response in enumerate(gp.RESPONSES)
    }
