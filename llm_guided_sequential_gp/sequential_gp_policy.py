"""Select a validated subset of up to three LLM length-scale prior proposals.

The sealed M_t prediction is verified against its snapshot before any new fit.
Only the requested truth row joins the revealed prefix. The unchanged rules and
all nonempty proposal subsets use identical fit settings and the same LOO folds.
Independent test data never participate in selection.
"""

from __future__ import annotations

import argparse
import contextlib
import copy
import io
import itertools
import json
import math
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path

try:
    from . import sequential_gp_diagnostics as diagnostics
    from . import sequential_gp_model as gp_model
    from . import sequential_gp_prompt as prompt_io
except ImportError:
    import sequential_gp_diagnostics as diagnostics
    import sequential_gp_model as gp_model
    import sequential_gp_prompt as prompt_io


HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "outputs" / "semiparametric_guarded_20_20261109"


@dataclass(frozen=True)
class GateSettings:
    """Fixed selection thresholds shared by every subset in an experiment."""

    min_absolute_loo_gain: float = 0.001
    min_relative_loo_gain: float = 0.01
    max_response_relative_worsening: float = 0.02

    def validate(self) -> None:
        for name, value in asdict(self).items():
            if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
        if self.min_relative_loo_gain >= 1 or self.max_response_relative_worsening >= 1:
            raise ValueError("relative LOO thresholds must be below 1")


def read_validated_answer(request_dir: Path) -> tuple[dict, dict, dict]:
    manifest, context, _ = prompt_io.load_request(request_dir)
    if manifest["stage"] != "sequential":
        raise ValueError("policy requires a sequential-stage prompt request")
    validation = gp_model.read_json(request_dir / "validation_result.json")
    raw = (request_dir / "llm_response_raw.txt").read_bytes()
    expected = {
        "request_id": manifest["request_id"],
        "stage": "sequential",
        "context_sha256": manifest["context_sha256"],
        "prompt_sha256": manifest["prompt_sha256"],
        "response_sha256": prompt_io.digest(raw),
    }
    if any(validation.get(key) != value for key, value in expected.items()):
        raise ValueError("LLM validation record does not match the saved request and raw answer")
    if validation.get("valid") is not True:
        raise ValueError(f"LLM answer was not validated: {validation.get('error')}")
    try:
        answer = json.loads(raw.decode("utf-8-sig"),
                            object_pairs_hook=prompt_io.unique_pairs,
                            parse_constant=prompt_io.reject_constant)
    except (UnicodeError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot decode the saved LLM answer: {error}") from error
    prompt_io.validate_sequential(answer, context)
    if validation.get("answer") != answer:
        raise ValueError("validated answer differs from the saved raw LLM answer")
    if answer["status"] != "OK":
        raise ValueError(f"LLM reported INVALID_INPUT: {answer['error']}")
    return manifest, context, answer


def read_ledger(path: Path | None, target: int) -> list[dict]:
    if path is None:
        if target > 21:
            raise ValueError("later sequential steps require the previous intervention ledger")
        return []
    ledger = gp_model.read_json(path)
    entries = ledger.get("entries")
    if not isinstance(entries, list):
        raise ValueError("ledger must contain an entries array")
    if ledger.get("schema_version") == 2 and target == 21 and entries == []:
        return []
    if ledger.get("schema_version") != 2:
        raise ValueError("completed steps require a version-2 intervention ledger")
    if set(ledger) != {"schema_version", "entries", "gate_settings",
                       "fit_settings", "experiment_settings"}:
        raise ValueError("version-2 ledger has unexpected fields")
    indices = []
    candidate_keys = {"response", "variable", "from_category", "to_category", "reason"}
    for entry in entries:
        if not isinstance(entry, dict) or set(entry) != {
            "sample_index", "design_id", "adopted", "proposed_candidates",
            "adopted_candidates", "outcome",
        }:
            raise ValueError("ledger entry has unexpected fields")
        index = entry["sample_index"]
        if type(index) is not int or not 21 <= index < target:
            raise ValueError("ledger entry must identify an earlier sequential sample")
        if type(entry["design_id"]) is not str or not entry["design_id"].strip():
            raise ValueError("ledger design_id must be nonempty")
        proposals, adopted = entry["proposed_candidates"], entry["adopted_candidates"]
        if (not isinstance(proposals, list) or len(proposals) > 3
                or not isinstance(adopted, list) or len(adopted) > len(proposals)
                or any(not isinstance(item, dict) or set(item) != candidate_keys
                       for item in proposals + adopted)):
            raise ValueError("ledger candidate arrays are malformed")
        if (any(any(type(item[key]) is not str or not item[key].strip()
                    for key in candidate_keys) for item in proposals + adopted)
                or any(item["response"] not in gp_model.RESPONSES
                       or item["variable"] not in gp_model.VARIABLES
                       or item["from_category"] not in gp_model.CATEGORY_LOG_SHIFT
                       or item["to_category"] not in gp_model.CATEGORY_LOG_SHIFT
                       or item["from_category"] == item["to_category"]
                       for item in proposals + adopted)):
            raise ValueError("ledger contains an invalid prior adjustment")
        if any(item not in proposals for item in adopted) or len({
            (item["response"], item["variable"]) for item in proposals
        }) != len(proposals):
            raise ValueError("ledger adopted candidates must be a unique proposal subset")
        outcome = entry["outcome"]
        expected = ("kept" if not proposals else "adopted" if adopted else "rolled_back")
        if outcome != expected or type(entry["adopted"]) is not bool or entry["adopted"] != bool(adopted):
            raise ValueError("ledger outcome and adopted candidates disagree")
        indices.append(index)
    if indices != list(range(21, target)):
        raise ValueError("ledger must contain every completed sequential step from 21 through t")
    return entries


def candidate_subsets(candidates: list[dict]) -> list[tuple[int, ...]]:
    if len(candidates) > 3:
        raise ValueError("at most three prior adjustments may be proposed")
    return [indices for size in range(1, len(candidates) + 1)
            for indices in itertools.combinations(range(len(candidates)), size)]


def propose_rules(current: dict, changes: list[dict]) -> dict:
    proposed = copy.deepcopy(gp_model.validate_rules(current))
    seen = set()
    for change in changes:
        response, variable = change["response"], change["variable"]
        pair = response, variable
        if pair in seen:
            raise ValueError(f"duplicate proposed relationship: {response}/{variable}")
        seen.add(pair)
        old = proposed["lengthscale_prior_categories"][response][variable]
        if change["from_category"] != old or change["to_category"] == old:
            raise ValueError(f"proposal disagrees with current prior: {response}/{variable}")
        proposed["lengthscale_prior_categories"][response][variable] = change["to_category"]
    return gp_model.validate_rules(proposed)


def loo_scores(report: dict) -> dict[str, float]:
    metrics = report.get("loo_metrics")
    if not isinstance(metrics, dict) or set(metrics) != set(gp_model.RESPONSES):
        raise ValueError("every rule set needs LOO metrics for all three responses")
    scores = {}
    for response in gp_model.RESPONSES:
        value = metrics[response].get("loo_nrmse")
        if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
            raise ValueError(f"invalid LOO-NRMSE for {response}")
        scores[response] = float(value)
    return scores


def compare_loo(baseline: dict, candidate: dict,
                settings: GateSettings) -> tuple[bool, dict]:
    old, new = loo_scores(baseline), loo_scores(candidate)
    old_mean = sum(old.values()) / len(old)
    new_mean = sum(new.values()) / len(new)
    gain = old_mean - new_mean
    required = max(settings.min_absolute_loo_gain,
                   settings.min_relative_loo_gain * old_mean)
    changes = {response: new[response] - old[response] for response in gp_model.RESPONSES}
    tolerance = {
        response: settings.max_response_relative_worsening * max(old[response], 1e-12)
        for response in gp_model.RESPONSES
    }
    passed = gain > 0 and gain >= required and all(
        changes[response] <= tolerance[response] for response in gp_model.RESPONSES
    )
    return passed, {
        "baseline_mean_loo_nrmse": old_mean,
        "candidate_mean_loo_nrmse": new_mean,
        "absolute_mean_gain": gain,
        "required_mean_gain": required,
        "per_response_nrmse_change": changes,
        "per_response_worsening_tolerance": tolerance,
        "passed": passed,
    }


def write_new_prefix(path: Path, train: Path, designs: Path, truth: Path,
                     target: int) -> None:
    rows = gp_model.read_csv(train, ("sample_index", "design_id", *gp_model.UNIT_FIELDS,
                                     *gp_model.RESPONSES))
    design = diagnostics.csv_rows(designs, ("sample_index", "design_id", *gp_model.UNIT_FIELDS),
                                  target=target)[0]
    revealed = diagnostics.csv_rows(truth, ("sample_index", "design_id", *gp_model.RESPONSES),
                                    target=target)[0]
    fields = list(rows[0])
    if set(gp_model.RESPONSES) & set(design):
        raise ValueError("design-only CSV must not expose response columns")
    new_row = {}
    for field in fields:
        if field in gp_model.RESPONSES:
            new_row[field] = revealed[field]
        elif field in design:
            new_row[field] = design[field]
        else:
            raise ValueError(f"new design does not provide training column {field}")
    if (new_row["sample_index"] != revealed["sample_index"]
            or new_row["design_id"] != revealed["design_id"]):
        raise ValueError("new design and revealed truth identify different samples")
    gp_model.write_csv(path, fields, [*rows, new_row])


def selected_final_index(designs: Path, sequential_count: int | None) -> int:
    rows = gp_model.read_csv(designs, ("sample_index", "design_id", *gp_model.UNIT_FIELDS))
    gp_model.parse_rows(rows, training=False)
    indices = [int(row["sample_index"]) for row in rows]
    if indices != list(range(gp_model.INITIAL_COUNT + 1,
                             gp_model.INITIAL_COUNT + len(rows) + 1)):
        raise ValueError("sequential design indices must be contiguous from sample 21")
    count = len(rows) if sequential_count is None else sequential_count
    if type(count) is not int or not 1 <= count <= len(rows):
        raise ValueError(f"sequential_count must be between 1 and {len(rows)} available designs")
    return gp_model.INITIAL_COUNT + count


def prior_config_from_report(report: dict) -> gp_model.PriorConfig:
    fit = report.get("fit_settings", {})
    if report.get("schema_version") != 2 or fit.get("lengthscale_prior_mode") != "soft_log_normal_map":
        raise ValueError("sealed model must use the version-2 soft prior")
    length_bounds = fit["length_scale_optimizer_bounds"]
    noise_bounds = fit["nugget_optimizer_bounds_standardized_variance"]
    config = gp_model.PriorConfig(
        base_factor=fit["lengthscale_prior_base_factor"],
        log_std=fit["lengthscale_prior_log_std"],
        length_scale_lower=length_bounds[0],
        length_scale_upper=length_bounds[1],
        noise_lower=noise_bounds[0],
        noise_upper=noise_bounds[1],
        noise_initial=fit["nugget_optimizer_initial_standardized_variance"],
    )
    config.validate(len(gp_model.VARIABLES))
    return config


def run(*, train: Path, reference: Path, designs: Path, truth: Path,
        prior_context: Path, model_dir: Path, request_dir: Path,
        ledger_path: Path | None, output: Path, settings: GateSettings,
        restarts: int = 3, seed: int = 2026,
        sequential_count: int | None = None) -> dict:
    settings.validate()
    if type(restarts) is not int or restarts < 0:
        raise ValueError("restarts must be a nonnegative integer")
    manifest, context, answer = read_validated_answer(request_dir)
    target = context["sample_index"]
    final = selected_final_index(designs, sequential_count)
    if target > final:
        raise ValueError(f"sample {target} exceeds the selected final sample {final}")
    experiment_settings = {
        "initial_count": gp_model.INITIAL_COUNT,
        "sequential_count": final - gp_model.INITIAL_COUNT,
        "final_sample_index": final,
    }
    entries = read_ledger(ledger_path, target)
    with tempfile.TemporaryDirectory(prefix="sequential_gp_policy_verify_") as temp:
        with contextlib.redirect_stdout(io.StringIO()):
            verified = diagnostics.run(train, designs, truth, model_dir, Path(temp),
                                       target, ledger_path, prior_context)
    if verified != context:
        raise ValueError("saved LLM context differs from the verified current-step diagnostic")
    current_rules = gp_model.validate_rules(gp_model.read_json(model_dir / "rules_used.json"))
    previous_report = gp_model.read_json(model_dir / "fit_report.json")
    prior_config = prior_config_from_report(previous_report)
    if (previous_report["fit_settings"]["optimizer_restarts"] != restarts
            or previous_report["fit_settings"]["random_seed"] != seed):
        raise ValueError("fit settings differ from the sealed prior model")
    fit_settings = {"restarts": restarts, "seed": seed, "prior_config": asdict(prior_config)}
    if ledger_path is not None and target > 21:
        ledger = gp_model.read_json(ledger_path)
        if (ledger.get("gate_settings") != asdict(settings)
                or ledger.get("fit_settings") != fit_settings
                or ledger.get("experiment_settings") != experiment_settings):
            raise ValueError("experiment, gate, and fit settings differ from the previous ledger")
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"policy output is not empty; choose a fresh --output: {output}")

    signature = gp_model.canonical_sha256({
        "request_manifest_sha256": gp_model.file_sha256(request_dir / "request_manifest.json"),
        "validation_sha256": gp_model.file_sha256(request_dir / "validation_result.json"),
        "response_sha256": gp_model.file_sha256(request_dir / "llm_response_raw.txt"),
        "train_sha256": gp_model.file_sha256(train),
        "reference_sha256": gp_model.file_sha256(reference),
        "designs_sha256": gp_model.file_sha256(designs),
        "prior_context_sha256": gp_model.file_sha256(prior_context),
        "previous_model_sha256": gp_model.file_sha256(model_dir / "gp_models.pkl"),
        "previous_rules_sha256": gp_model.file_sha256(model_dir / "rules_used.json"),
        "revealed_row": diagnostics.csv_rows(
            truth, ("sample_index", "design_id", *gp_model.RESPONSES), target=target
        )[0],
        "diagnostic_context_sha256": gp_model.canonical_sha256(context),
        "ledger_sha256": gp_model.file_sha256(ledger_path) if ledger_path else None,
        "gate_settings": asdict(settings), "fit_settings": fit_settings,
        "experiment_settings": experiment_settings,
    })
    output.mkdir(parents=True, exist_ok=True)
    new_train = output / f"revealed_{target:04d}.csv"
    write_new_prefix(new_train, train, designs, truth, target)
    next_designs = designs if target < final else None
    baseline_dir = output / "baseline_model"
    baseline = gp_model.run(
        new_train, reference, baseline_dir, rules_path=model_dir / "rules_used.json",
        initial_decision=model_dir / "rules_used.json", predict_path=next_designs,
        predict_index=target + 1 if next_designs else None,
        restarts=restarts, seed=seed, loo=True, prior_config=prior_config,
    )
    selected = None
    comparisons = []
    proposals = answer["candidates"]
    for subset in candidate_subsets(proposals):
        changes = [proposals[i] for i in subset]
        rules = propose_rules(current_rules, changes)
        label = "_".join(str(i + 1) for i in subset)
        rules_path = output / f"candidate_{label}_rules.json"
        gp_model.write_json(rules_path, rules)
        candidate_dir = output / f"candidate_{label}_model"
        candidate = gp_model.run(
            new_train, reference, candidate_dir, rules_path=rules_path,
            initial_decision=rules_path, predict_path=next_designs,
            predict_index=target + 1 if next_designs else None,
            restarts=restarts, seed=seed, loo=True, prior_config=prior_config,
        )
        for key in ("train_sha256", "reference_initial_sha256",
                    "blind_prediction_sample_index", "fit_settings"):
            if baseline[key] != candidate[key]:
                raise ValueError(f"candidate {label} differs in shared fit setting {key}")
        passed, gate = compare_loo(baseline, candidate, settings)
        comparison = {
            "candidate_indices": [i + 1 for i in subset],
            "candidate_model_dir": str(candidate_dir.resolve()),
            "rules_sha256": gp_model.canonical_sha256(rules),
            "loo_metrics": candidate["loo_metrics"],
            "gate": gate,
        }
        comparisons.append(comparison)
        if passed:
            rank = (gate["candidate_mean_loo_nrmse"], len(subset), subset)
            if selected is None or rank < selected[0]:
                selected = rank, comparison, rules, changes

    chosen_dir = (Path(selected[1]["candidate_model_dir"]) if selected else baseline_dir)
    final_rules = selected[2] if selected else current_rules
    adopted = selected[3] if selected else []
    outcome = "adopted" if adopted else "rolled_back" if proposals else "kept"
    entry = {
        "sample_index": target, "design_id": context["design_id"],
        "adopted": bool(adopted), "proposed_candidates": proposals,
        "adopted_candidates": adopted, "outcome": outcome,
    }
    ledger_file = output / "intervention_ledger.json"
    gp_model.write_json(ledger_file, {
        "schema_version": 2, "experiment_settings": experiment_settings,
        "gate_settings": asdict(settings), "fit_settings": fit_settings,
        "entries": [*entries, entry],
    })
    record = {
        "schema_version": 2,
        "sample_index": target, "design_id": context["design_id"],
        "input_signature_sha256": signature, "request_id": manifest["request_id"],
        "llm_decision": answer["decision"], "llm_reason": answer["reason"],
        "proposed_candidates": proposals, "adopted_candidates": adopted,
        "outcome": outcome, "gate_settings": asdict(settings),
        "fit_settings": fit_settings, "experiment_settings": experiment_settings,
        "new_training_file": str(new_train.resolve()),
        "baseline_model_dir": str(baseline_dir.resolve()),
        "candidate_comparisons": comparisons,
        "adopted_model_dir": str(chosen_dir.resolve()),
        "adopted_rules_sha256": gp_model.canonical_sha256(final_rules),
        "baseline_loo_metrics": baseline["loo_metrics"],
        "next_blind_sample_index": target + 1 if next_designs else None,
        "ledger_file": str(ledger_file.resolve()),
        "scope": "Revealed-prefix LOO selects prior rules; future truth and independent tests are excluded.",
    }
    gp_model.write_json(output / "decision_record.json", record)
    print(f"Sample {target}: LLM={answer['decision']}, outcome={outcome}, "
          f"adopted {len(adopted)} of {len(proposals)} proposed changes")
    print(f"Adopted model: {chosen_dir}")
    print(f"Decision record: {output / 'decision_record.json'}")
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train", type=Path, default=DATA_DIR / "initial20.csv")
    parser.add_argument("--reference-initial", type=Path, default=DATA_DIR / "initial20.csv")
    parser.add_argument("--designs", type=Path, default=DATA_DIR / "sequential_designs.csv")
    parser.add_argument("--truth", type=Path, default=DATA_DIR / "sequential_truth.csv")
    parser.add_argument("--prior-context", type=Path, default=DATA_DIR / "initial20_context.json")
    parser.add_argument("--model-dir", type=Path)
    parser.add_argument("--request", type=Path)
    parser.add_argument("--ledger", type=Path)
    parser.add_argument("--sequential-count", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--restarts", type=int, default=3)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--min-absolute-loo-gain", type=float, default=0.001)
    parser.add_argument("--min-relative-loo-gain", type=float, default=0.01)
    parser.add_argument("--max-response-relative-worsening", type=float, default=0.02)
    args = parser.parse_args()
    try:
        train = args.train.resolve()
        count = len(gp_model.read_csv(
            train, ("sample_index", "design_id", *gp_model.UNIT_FIELDS, *gp_model.RESPONSES)
        ))
        target = count + 1
        settings = GateSettings(
            min_absolute_loo_gain=args.min_absolute_loo_gain,
            min_relative_loo_gain=args.min_relative_loo_gain,
            max_response_relative_worsening=args.max_response_relative_worsening,
        )
        run(
            train=train, reference=args.reference_initial.resolve(),
            designs=args.designs.resolve(), truth=args.truth.resolve(),
            prior_context=args.prior_context.resolve(),
            model_dir=(args.model_dir or HERE / "outputs" / "models" /
                       f"samples_{count:04d}").resolve(),
            request_dir=(args.request or HERE / "outputs" / "prompt_requests" /
                         f"sequential_{target:04d}").resolve(),
            ledger_path=args.ledger.resolve() if args.ledger else None,
            output=(args.output or args.designs.resolve().parent / "policy" /
                    f"sample_{target:04d}").resolve(),
            settings=settings, restarts=args.restarts, seed=args.seed,
            sequential_count=args.sequential_count,
        )
    except (ValueError, OSError, KeyError, IndexError, TypeError,
            gp_model.np.linalg.LinAlgError) as error:
        parser.exit(2, f"Policy step failed: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
