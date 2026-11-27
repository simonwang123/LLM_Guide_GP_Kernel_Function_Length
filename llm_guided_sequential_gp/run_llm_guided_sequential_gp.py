"""Run a 20+N LLM-guided sequential GP experiment, one revealed point at a time.

Example from the project root:
  .venv/Scripts/python.exe llm_guided_sequential_gp/run_llm_guided_sequential_gp.py \
      --initial-count 20 --sequential-steps 10 \
      --provider openrouter --model PROVIDER/MODEL

The run directory records each LLM request, response, diagnosis, fit, policy
decision and ledger. A rerun with identical settings resumes completed steps
and preserves partial local fits as separate attempts. It never silently sends
an API request again after a failed or invalid call. An optional independent
40-point test set is prepared and scored only after the final online step.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import math
import pickle
import tempfile
from dataclasses import asdict
from pathlib import Path

try:
    from . import sequential_gp_data as data_io
    from . import sequential_gp_diagnostics as diagnostics
    from . import sequential_gp_llm as llm_io
    from . import sequential_gp_model as gp_model
    from . import sequential_gp_policy as policy
    from . import sequential_gp_prompt as prompt_io
except ImportError:  # Direct execution from this directory.
    import sequential_gp_data as data_io
    import sequential_gp_diagnostics as diagnostics
    import sequential_gp_llm as llm_io
    import sequential_gp_model as gp_model
    import sequential_gp_policy as policy
    import sequential_gp_prompt as prompt_io


HERE = Path(__file__).resolve().parent


def valid_json_marker(path: Path) -> bool:
    if not path.is_file():
        return False
    try:
        gp_model.read_json(path)
    except ValueError:
        return False
    return True


def artifact_directory(parent: Path, name: str, marker: str) -> Path:
    """Reuse one completed artifact, or preserve partial attempts and choose a fresh directory."""
    candidates = [parent / name, *sorted(parent.glob(f"{name}_attempt_*"))]
    completed = [path for path in candidates if valid_json_marker(path / marker)]
    if len(completed) > 1:
        raise ValueError(f"multiple completed {name} artifacts exist in {parent}")
    if completed:
        return completed[0]
    for path in candidates:
        if not path.exists() or (path.is_dir() and not any(path.iterdir())):
            return path
    index = 1
    while (parent / f"{name}_attempt_{index:04d}").exists():
        index += 1
    return parent / f"{name}_attempt_{index:04d}"


def request_directory(parent: Path) -> Path:
    """Retry incomplete prompt rendering, but never conceal evidence of an API attempt."""
    candidates = [parent / "request", *sorted(parent.glob("request_attempt_*"))]
    for path in candidates:
        if not valid_json_marker(path / "request_manifest.json") and any(
            (path / name).exists() for name in (
                "api_call.json", "provider_response.json", "llm_response_raw.txt",
                "validation_result.json",
            )
        ):
            raise ValueError(f"incomplete request has API or answer files; inspect {path}")
    return artifact_directory(parent, "request", "request_manifest.json")


def write_json_atomic(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def ensure_run_config(root: Path, config: dict) -> None:
    path = root / "run_config.json"
    temporary = path.with_name(path.name + ".tmp")
    if path.exists():
        if gp_model.read_json(path) != config:
            raise ValueError(f"existing run has different settings; choose another --output-root: {root}")
    elif root.exists() and any(root.iterdir()):
        if set(root.iterdir()) == {temporary}:
            if valid_json_marker(temporary) and gp_model.read_json(temporary) != config:
                raise ValueError(f"interrupted run configuration differs: {root}")
            write_json_atomic(path, config)
        else:
            raise ValueError(f"run directory is nonempty without run_config.json: {root}")
    else:
        write_json_atomic(path, config)


def ensure_prepared_data(run: Path, output: Path, steps: int,
                         test_dataset: Path | None = None) -> dict:
    summary_path = output / "data_summary.json"
    if not summary_path.exists():
        if output.exists() and any(output.iterdir()):
            raise ValueError(f"prepared data directory is incomplete: {output}")
        return data_io.prepare(run, output, test_dataset, steps)
    summary = gp_model.read_json(summary_path)
    counts = summary.get("sample_counts", {})
    if (counts.get("initial") != data_io.INITIAL_COUNT
            or counts.get("sequential") != steps
            or Path(summary.get("source_run", "")).resolve() != run
            or summary.get("dataset_sha256") != data_io.sha256(run / "dataset.csv")
            or summary.get("history_sha256") != data_io.sha256(run / "optimization_history.csv")
            or counts.get("independent_test") != (data_io.TEST_COUNT if test_dataset else 0)
            or summary.get("test_dataset_sha256") != (
                data_io.sha256(test_dataset) if test_dataset else None
            )):
        raise ValueError("prepared data do not match the requested source and sample count")
    names = ["initial20.csv", "sequential_designs.csv", "sequential_truth.csv",
             "initial20_context.json"]
    if test_dataset:
        names.append("holdout/test_holdout.csv")
    for name in names:
        if not (output / name).is_file():
            raise ValueError(f"prepared data are missing {name}")
    if policy.selected_final_index(output / "sequential_designs.csv", None) != 20 + steps:
        raise ValueError("prepared design count differs from run configuration")
    prompt_io.validate_context("initial", prompt_io.read_json_file(output / "initial20_context.json"))
    with tempfile.TemporaryDirectory(prefix="sequential_gp_data_verify_") as temporary:
        expected = Path(temporary)
        with contextlib.redirect_stdout(io.StringIO()):
            data_io.prepare(run, expected, test_dataset, steps)
        if any(gp_model.file_sha256(output / name) != gp_model.file_sha256(expected / name)
               for name in names):
            raise ValueError("prepared CSV or initial context differs from the source data")
    return summary


def ensure_request(stage: str, context_path: Path, request_dir: Path) -> dict:
    context = prompt_io.validate_context(stage, prompt_io.read_json_file(context_path))
    manifest_path = request_dir / "request_manifest.json"
    if manifest_path.exists():
        manifest, saved_context, _ = prompt_io.load_request(request_dir)
        if manifest["stage"] != stage or saved_context != context:
            raise ValueError(f"saved request belongs to different {stage} context: {request_dir}")
        return manifest
    if request_dir.exists() and any(request_dir.iterdir()):
        raise ValueError(f"request directory is incomplete: {request_dir}")
    return prompt_io.render(stage, context_path, request_dir, None)


def verify_saved_answer(request_dir: Path) -> dict:
    manifest, context, _ = prompt_io.load_request(request_dir)
    validation = gp_model.read_json(request_dir / "validation_result.json")
    raw = (request_dir / "llm_response_raw.txt").read_bytes()
    expected = {
        "request_id": manifest["request_id"],
        "stage": manifest["stage"],
        "context_sha256": manifest["context_sha256"],
        "prompt_sha256": manifest["prompt_sha256"],
        "response_sha256": prompt_io.digest(raw),
    }
    if any(validation.get(name) != value for name, value in expected.items()):
        raise ValueError(f"LLM validation does not match saved request: {request_dir}")
    if validation.get("valid") is not True:
        raise ValueError(f"LLM output is invalid: {validation.get('error')}")
    try:
        answer = prompt_io.parse_llm_json(raw)
    except (UnicodeError, ValueError) as error:
        raise ValueError(f"cannot decode saved LLM answer: {error}") from error
    if manifest["stage"] == "initial":
        prompt_io.validate_initial(answer)
    else:
        prompt_io.validate_sequential(answer, context)
    if validation.get("answer") != answer:
        raise ValueError("validated LLM answer differs from raw response")
    if answer["status"] != "OK":
        raise ValueError(f"LLM reported INVALID_INPUT: {answer['error']}")
    return answer


def ensure_llm_answer(request_dir: Path, llm_settings: dict) -> dict:
    validation_path = request_dir / "validation_result.json"
    raw_path = request_dir / "llm_response_raw.txt"
    if validation_path.is_file():
        if not raw_path.is_file():
            raise ValueError(f"validation exists without raw LLM response: {request_dir}")
        return verify_saved_answer(request_dir)
    if raw_path.is_file():
        result = prompt_io.validate_response(request_dir, raw_path)
        if result["valid"] is not True:
            raise ValueError(f"saved LLM output is invalid: {result['error']}")
        return verify_saved_answer(request_dir)
    if any((request_dir / name).exists() for name in ("api_call.json", "provider_response.json")):
        raise ValueError(
            f"prior API attempt did not yield a valid raw response: {request_dir}; "
            "inspect it and use a new run directory for a fresh call"
        )
    result = llm_io.run(request_dir, **llm_settings, dry_run=False)
    if result.get("status") != "validated":
        raise ValueError(f"LLM API call did not yield a valid answer: {result.get('status')}; {request_dir}")
    return verify_saved_answer(request_dir)


def ensure_initial_model(data_dir: Path, request_dir: Path, model_dir: Path,
                         restarts: int, seed: int,
                         prior_config: gp_model.PriorConfig = gp_model.DEFAULT_PRIOR_CONFIG) -> dict:
    train = data_dir / "initial20.csv"
    designs = data_dir / "sequential_designs.csv"
    decision = request_dir / "validation_result.json"
    report_path = model_dir / "fit_report.json"
    if not report_path.exists():
        if model_dir.exists() and any(model_dir.iterdir()):
            raise ValueError(f"cold-start model directory is incomplete: {model_dir}")
        return gp_model.run(
            train, train, model_dir, rules_path=None, initial_decision=decision,
            predict_path=designs, predict_index=21, restarts=restarts, seed=seed, validate=False,
            prior_config=prior_config,
        )
    report = gp_model.read_json(report_path)
    rules = gp_model.rules_from_initial_decision(
        decision, gp_model.read_csv(train, ("sample_index", "design_id", *gp_model.UNIT_FIELDS,
                                            *gp_model.RESPONSES))
    )
    if (report.get("schema_version") not in (2, 3)
            or report.get("training_sample_count") != 20
            or report.get("train_sha256") != gp_model.file_sha256(train)
            or report.get("reference_initial_sha256") != gp_model.file_sha256(train)
            or report.get("rules_sha256") != gp_model.canonical_sha256(rules)
            or report.get("blind_prediction_source") != str(designs.resolve())
            or report.get("blind_prediction_sample_index") != 21
            or report.get("fit_settings", {}).get("optimizer_restarts") != restarts
            or report.get("fit_settings", {}).get("random_seed") != seed
            or asdict(policy.prior_config_from_report(report)) != asdict(prior_config)
            or gp_model.read_json(model_dir / "rules_used.json") != rules
            or not (model_dir / "gp_models.pkl").is_file()
            or not (model_dir / "rules_used.json").is_file()
            or not (model_dir / "blind_prediction.csv").is_file()):
        raise ValueError("saved cold-start model does not match this experiment")
    return report


def ensure_initial_ledger(initial_dir: Path) -> Path:
    path = initial_dir / "intervention_ledger.json"
    expected = {"schema_version": 2, "entries": []}
    if path.exists():
        if gp_model.read_json(path) != expected:
            raise ValueError("initial intervention ledger is not empty")
    else:
        write_json_atomic(path, expected)
    policy.read_ledger(path, 21)
    return path


def verify_sealed_prediction(train: Path, designs: Path, model_dir: Path,
                             target: int) -> None:
    """Check the M_t prediction and identities before opening the truth CSV."""
    report = gp_model.read_json(model_dir / "fit_report.json")
    train_rows = gp_model.read_csv(
        train, ("sample_index", "design_id", *gp_model.UNIT_FIELDS, *gp_model.RESPONSES)
    )
    design = diagnostics.csv_rows(
        designs, ("sample_index", "design_id", *gp_model.UNIT_FIELDS), target=target
    )[0]
    blind_rows = gp_model.read_csv(
        model_dir / "blind_prediction.csv", ("sample_index", "design_id")
    )
    if (len(train_rows) != target - 1 or len(blind_rows) != 1
            or report.get("training_sample_count") != target - 1
            or report.get("train_sha256") != gp_model.file_sha256(train)
            or report.get("blind_prediction_source") != str(designs.resolve())
            or report.get("blind_prediction_sample_index") != target
            or int(blind_rows[0]["sample_index"]) != target
            or blind_rows[0]["design_id"] != design["design_id"]):
        raise ValueError(f"sealed M_{target - 1} prediction does not match sample {target}")
    for response in gp_model.RESPONSES:
        for suffix in ("mean_m", "observed_std_m", "latent_std_m"):
            value = float(blind_rows[0][f"{response}_{suffix}"])
            if not math.isfinite(value) or (suffix != "mean_m" and value < 0):
                raise ValueError(f"sealed prediction has invalid {response}/{suffix}")
        if float(blind_rows[0][f"{response}_observed_std_m"]) <= 0:
            raise ValueError(f"sealed prediction has nonpositive observed std for {response}")
    if not (model_dir / "gp_models.pkl").is_file():
        raise ValueError("sealed model snapshot is missing")
    rules = gp_model.read_json(model_dir / "rules_used.json")
    if report.get("rules_sha256") != gp_model.canonical_sha256(rules):
        raise ValueError("sealed model rules disagree with the fit report")
    # The snapshot must be a locally trusted output of sequential_gp_model.py.
    try:
        with (model_dir / "gp_models.pkl").open("rb") as stream:
            snapshot = pickle.load(stream)
    except (ImportError, ModuleNotFoundError) as error:
        raise ValueError("cannot load sealed GP snapshot in this Python environment") from error
    if (snapshot.get("training_sample_count") != target - 1
            or snapshot.get("train_sha256") != gp_model.file_sha256(train)
            or snapshot.get("rules") != rules
            or tuple(snapshot.get("responses", ())) != gp_model.RESPONSES
            or tuple(snapshot.get("variables", ())) != gp_model.VARIABLES):
        raise ValueError("sealed GP snapshot disagrees with training data or rules")
    fitted = snapshot.get("fitted")
    if not isinstance(fitted, list) or len(fitted) != len(gp_model.RESPONSES):
        raise ValueError("sealed GP snapshot does not contain three fitted responses")
    x = diagnostics.unit_design(design, "next design")
    for response, fitted_response in zip(gp_model.RESPONSES, fitted, strict=True):
        mean, observed_std, latent_std = diagnostics.model_prediction(fitted_response, x)
        for suffix, calculated in (("mean_m", mean), ("observed_std_m", observed_std),
                                   ("latent_std_m", latent_std)):
            sealed = float(blind_rows[0][f"{response}_{suffix}"])
            if not math.isclose(sealed, calculated, rel_tol=1e-8, abs_tol=1e-10):
                raise ValueError(f"sealed prediction differs from M_{target - 1}: {response}/{suffix}")


def ensure_policy_record(policy_dir: Path, request_dir: Path, train: Path,
                         model_dir: Path, designs: Path, truth: Path,
                         target: int, final: int, settings: policy.GateSettings,
                         restarts: int, seed: int) -> dict:
    path = policy_dir / "decision_record.json"
    if not path.exists():
        if policy_dir.exists() and any(policy_dir.iterdir()):
            raise ValueError(f"policy output is incomplete: {policy_dir}")
        return {}
    record = gp_model.read_json(path)
    manifest, _, answer = policy.read_validated_answer(request_dir)
    new_train = policy_dir / f"revealed_{target:04d}.csv"
    chosen = Path(record.get("adopted_model_dir", "")).resolve()
    baseline = (policy_dir / "baseline_model").resolve()
    proposals = answer["candidates"]
    comparisons = record.get("candidate_comparisons")
    subsets = policy.candidate_subsets(proposals)
    expected_candidate_dirs = {
        tuple(index + 1 for index in subset): (
            policy_dir / f"candidate_{'_'.join(str(index + 1) for index in subset)}_model"
        ).resolve()
        for subset in subsets
    }
    adopted = record.get("adopted_candidates")
    adopted_indices = tuple(index + 1 for index, item in enumerate(proposals)
                            if isinstance(adopted, list) and item in adopted)
    expected_chosen = expected_candidate_dirs.get(adopted_indices, baseline)
    previous_report = gp_model.read_json(model_dir / "fit_report.json")
    expected_fit_settings = {
        "restarts": restarts, "seed": seed,
        "prior_config": policy.asdict(policy.prior_config_from_report(previous_report)),
    }
    if (record.get("sample_index") != target
            or record.get("schema_version") != 2
            or record.get("request_id") != manifest["request_id"]
            or record.get("llm_decision") != answer["decision"]
            or record.get("proposed_candidates") != proposals
            or not isinstance(adopted, list)
            or any(item not in proposals for item in adopted)
            or len(adopted_indices) != len(adopted)
            or record.get("outcome") != (
                "adopted" if adopted else "rolled_back" if proposals else "kept"
            )
            or record.get("gate_settings") != policy.asdict(settings)
            or record.get("fit_settings") != expected_fit_settings
            or record.get("experiment_settings") != {
                "initial_count": 20, "sequential_count": final - 20,
                "final_sample_index": final,
            }
            or Path(record.get("new_training_file", "")).resolve() != new_train.resolve()
            or Path(record.get("baseline_model_dir", "")).resolve() != baseline
            or chosen != expected_chosen
            or not isinstance(comparisons, list)
            or len(comparisons) != len(subsets)
            or record.get("next_blind_sample_index") != (target + 1 if target < final else None)):
        raise ValueError(f"saved policy decision does not match this step: {path}")
    baseline_report = gp_model.read_json(baseline / "fit_report.json")
    if (baseline_report.get("train_sha256") != gp_model.file_sha256(new_train)
            or baseline_report.get("rules_sha256") != previous_report.get("rules_sha256")
            or baseline_report.get("validation_metrics") != record.get("baseline_validation_metrics")):
        raise ValueError("saved baseline fit does not match this policy step")
    passing = []
    for comparison, subset in zip(comparisons, subsets, strict=True):
        indices = tuple(index + 1 for index in subset)
        candidate_dir = expected_candidate_dirs[indices]
        if (comparison.get("candidate_indices") != list(indices)
                or Path(comparison.get("candidate_model_dir", "")).resolve()
                != candidate_dir):
            raise ValueError("saved candidate comparison does not match the proposed subset")
        candidate_report = gp_model.read_json(candidate_dir / "fit_report.json")
        passed, gate = policy.compare_validation(baseline_report, candidate_report, settings)
        if (candidate_report.get("train_sha256") != baseline_report["train_sha256"]
                or candidate_report.get("fit_settings") != baseline_report.get("fit_settings")
                or candidate_report.get("rules_sha256") != comparison.get("rules_sha256")
                or candidate_report.get("validation_metrics") != comparison.get("validation_metrics")
                or gate != comparison.get("gate")):
            raise ValueError("saved candidate comparison differs from its fitted model")
        if passed:
            passing.append((gate["candidate_mean_composite_nrmse"], len(indices),
                            tuple(index - 1 for index in indices), candidate_dir))
    winner = min(passing) if passing else None
    if (chosen != (winner[3] if winner else baseline)
            or adopted != ([proposals[index] for index in winner[2]] if winner else [])):
        raise ValueError("saved adopted candidate subset is not the best passing fit")
    previous = gp_model.read_csv(
        train, ("sample_index", "design_id", *gp_model.UNIT_FIELDS, *gp_model.RESPONSES)
    )
    updated = gp_model.read_csv(
        new_train, ("sample_index", "design_id", *gp_model.UNIT_FIELDS, *gp_model.RESPONSES)
    )
    with tempfile.TemporaryDirectory(prefix="sequential_gp_prefix_verify_") as temporary:
        expected_train = Path(temporary) / "revealed.csv"
        policy.write_new_prefix(expected_train, train, designs, truth, target)
        if gp_model.file_sha256(new_train) != gp_model.file_sha256(expected_train):
            raise ValueError("saved training prefix differs from the revealed design and truth")
    fit = gp_model.read_json(chosen / "fit_report.json")
    if (len(updated) != target or updated[:-1] != previous
            or fit.get("train_sha256") != gp_model.file_sha256(new_train)
            or fit.get("rules_sha256") != record.get("adopted_rules_sha256")
            or fit.get("blind_prediction_sample_index") != record["next_blind_sample_index"]):
        raise ValueError(f"saved adopted model does not match the revealed prefix: {chosen}")
    ledger_file = (policy_dir / "intervention_ledger.json").resolve()
    if Path(record.get("ledger_file", "")).resolve() != ledger_file:
        raise ValueError("saved policy ledger path is inconsistent")
    entries = policy.read_ledger(ledger_file, target + 1)
    if (len(entries) != target - 20 or entries[-1]["sample_index"] != target
            or entries[-1]["proposed_candidates"] != proposals
            or entries[-1]["adopted_candidates"] != adopted
            or entries[-1]["outcome"] != record["outcome"]):
        raise ValueError("saved policy ledger does not include this step")
    if target < final and not (chosen / "blind_prediction.csv").is_file():
        raise ValueError("adopted model has no sealed prediction for the next step")
    return record


def prediction_metrics(predictions: list[dict], reference_initial: Path) -> dict:
    """Score sealed predictions with the initial-20 response scale."""
    initial = gp_model.read_csv(
        reference_initial, ("sample_index", "design_id", *gp_model.UNIT_FIELDS,
                            *gp_model.RESPONSES),
    )
    scores = {}
    for response in gp_model.RESPONSES:
        scale = float(gp_model.np.std(
            [float(row[response]) for row in initial], ddof=1,
        ))
        if not math.isfinite(scale) or scale <= 0:
            raise ValueError(f"initial response scale is invalid for {response}")
        errors = [item[response]["truth_m"] - item[response]["mean_m"]
                  for item in predictions]
        observed_stds = [item[response]["observed_std_m"] for item in predictions]
        if (not predictions or any(not math.isfinite(value) for value in errors)
                or any(not math.isfinite(value) or value <= 0
                       for value in observed_stds)):
            raise ValueError(f"predictions or standard deviations are invalid for {response}")
        rmse = math.sqrt(sum(error**2 for error in errors) / len(errors))
        scores[response] = {
            "count": len(errors),
            "rmse_m": rmse,
            "mae_m": sum(abs(error) for error in errors) / len(errors),
            "nrmse_initial20_std": rmse / scale,
            "initial20_std_m": scale,
            "nominal_95_coverage": sum(
                abs(error) <= 1.96 * std
                for error, std in zip(errors, observed_stds, strict=True)
            ) / len(errors),
        }
    return scores


def evaluate_independent_test(model_dir: Path, holdout: Path,
                              reference_initial: Path, output: Path) -> dict:
    """Read the independent test truths only after all online steps are complete."""
    report = gp_model.read_json(model_dir / "fit_report.json")
    rules = gp_model.validate_rules(gp_model.read_json(model_dir / "rules_used.json"))
    if (report.get("schema_version") not in (2, 3)
            or report.get("rules_sha256") != gp_model.canonical_sha256(rules)):
        raise ValueError("final model report and prior rules disagree")
    with (model_dir / "gp_models.pkl").open("rb") as stream:
        snapshot = pickle.load(stream)
    if (snapshot.get("schema_version") != 2
            or snapshot.get("rules") != rules
            or snapshot.get("train_sha256") != report.get("train_sha256")
            or snapshot.get("training_sample_count") != report.get("training_sample_count")
            or tuple(snapshot.get("responses", ())) != gp_model.RESPONSES
            or tuple(snapshot.get("variables", ())) != gp_model.VARIABLES):
        raise ValueError("final model snapshot does not match the saved fit")
    rows = gp_model.read_csv(
        holdout, ("test_index", "design_id", *gp_model.UNIT_FIELDS,
                  *gp_model.RESPONSES),
    )
    if len(rows) != data_io.TEST_COUNT:
        raise ValueError(f"independent test set must contain {data_io.TEST_COUNT} rows")
    x = gp_model.np.asarray([
        [float(row[field]) for field in gp_model.UNIT_FIELDS] for row in rows
    ], dtype=float)
    if not gp_model.np.isfinite(x).all() or gp_model.np.any((x < 0) | (x > 1)):
        raise ValueError("independent test designs must be finite normalized values")
    data = [{"test_index": int(row["test_index"]), "design_id": row["design_id"]}
            for row in rows]
    predictions = [{} for _ in rows]
    for response, fitted in zip(gp_model.RESPONSES, snapshot["fitted"], strict=True):
        mean, observed_std, latent_std = gp_model.predict_one(fitted, x)
        if (not gp_model.np.isfinite(mean).all()
                or not gp_model.np.isfinite(observed_std).all()
                or not gp_model.np.isfinite(latent_std).all()
                or gp_model.np.any(observed_std <= 0)):
            raise ValueError(f"independent test predictions are invalid for {response}")
        for index, row in enumerate(rows):
            truth = float(row[response])
            if not math.isfinite(truth):
                raise ValueError(f"independent test truth is invalid for {response}")
            predictions[index][response] = {
                "truth_m": truth, "mean_m": float(mean[index]),
                "observed_std_m": float(observed_std[index]),
            }
            data[index].update({
                f"{response}_truth": truth,
                f"{response}_mean": float(mean[index]),
                f"{response}_observed_std": float(observed_std[index]),
                f"{response}_latent_std": float(latent_std[index]),
            })
    metrics = prediction_metrics(predictions, reference_initial)
    with tempfile.TemporaryDirectory(prefix="sequential_gp_test_verify_") as temporary:
        expected_csv = Path(temporary) / "test_predictions.csv"
        gp_model.write_csv(expected_csv, list(data[0]), data)
        prediction_hash = gp_model.file_sha256(expected_csv)
        result = {
            "schema_version": 2,
            "scope": "Independent test set, read only after the final sequential model was selected.",
            "final_model_dir": str(model_dir.resolve()),
            "final_model_sha256": gp_model.file_sha256(model_dir / "gp_models.pkl"),
            "holdout_sha256": gp_model.file_sha256(holdout),
            "reference_initial_sha256": gp_model.file_sha256(reference_initial),
            "test_predictions_file": str((output / "test_predictions.csv").resolve()),
            "test_predictions_sha256": prediction_hash,
            "metrics": metrics,
        }
        report_path = output / "test_metrics.json"
        if report_path.is_file():
            if (gp_model.read_json(report_path) != result
                    or gp_model.file_sha256(output / "test_predictions.csv") != prediction_hash):
                raise ValueError("saved independent test evaluation differs from this final model")
        elif output.exists() and any(output.iterdir()):
            raise ValueError(f"independent test output is incomplete: {output}")
        else:
            output.mkdir(parents=True, exist_ok=True)
            (output / "test_predictions.csv").write_bytes(expected_csv.read_bytes())
            write_json_atomic(report_path, result)
    return result


def run_experiment(*, source_run: Path, output_root: Path,
                   initial_count: int, sequential_steps: int,
                   llm_settings: dict, restarts: int, seed: int,
                   gate_settings: policy.GateSettings,
                   prior_config: gp_model.PriorConfig = gp_model.DEFAULT_PRIOR_CONFIG,
                   test_dataset: Path | None = None) -> dict:
    if initial_count != data_io.INITIAL_COUNT:
        raise ValueError("current prompts and GP model require --initial-count 20")
    if type(sequential_steps) is not int or sequential_steps < 1:
        raise ValueError("--sequential-steps must be a positive integer")
    if restarts < 0:
        raise ValueError("--restarts must be nonnegative")
    gate_settings.validate()
    prior_config.validate(len(gp_model.VARIABLES))
    source_run = source_run.resolve()
    output_root = output_root.resolve()
    test_dataset = test_dataset.resolve() if test_dataset else None
    config = {
        "schema_version": 4,
        "source_run": str(source_run),
        "source_dataset_sha256": gp_model.file_sha256(source_run / "dataset.csv"),
        "source_history_sha256": gp_model.file_sha256(source_run / "optimization_history.csv"),
        "independent_test_source": str(test_dataset) if test_dataset else None,
        "independent_test_source_sha256": (
            gp_model.file_sha256(test_dataset) if test_dataset else None
        ),
        "initial_count": initial_count,
        "sequential_steps": sequential_steps,
        "llm": {
            key: value for key, value in llm_settings.items()
            if key != "extra_body"
        } | {"extra_body_sha256": gp_model.canonical_sha256(llm_settings["extra_body"])},
        "fit": {"restarts": restarts, "seed": seed,
                "prior_config": asdict(prior_config)},
        "gate_settings": policy.asdict(gate_settings),
        "implementation_sha256": {
            name: gp_model.file_sha256(HERE / name) for name in (
                "sequential_gp_data.py", "sequential_gp_prompt.py",
                "sequential_gp_llm.py", "sequential_gp_model.py",
                "sequential_gp_prior.py",
                "sequential_gp_diagnostics.py", "sequential_gp_policy.py",
                "sequential_gp_validation.py", "run_llm_guided_sequential_gp.py",
                "prompts/gp_common_rules.md", "prompts/gp_initial_tiers.md",
                "prompts/gp_sequential_diagnosis.md",
            )
        },
    }
    ensure_run_config(output_root, config)
    data_dir = artifact_directory(output_root, "data", "data_summary.json")
    ensure_prepared_data(source_run, data_dir, sequential_steps)
    train = data_dir / "initial20.csv"
    designs = data_dir / "sequential_designs.csv"
    truth = data_dir / "sequential_truth.csv"
    initial_request = request_directory(output_root / "initial")
    ensure_request("initial", data_dir / "initial20_context.json", initial_request)
    ensure_llm_answer(initial_request, llm_settings)
    initial_model = artifact_directory(output_root / "initial", "model", "fit_report.json")
    ensure_initial_model(data_dir, initial_request, initial_model, restarts, seed,
                         prior_config)
    initial_ledger = ensure_initial_ledger(output_root / "initial")
    summary_path = output_root / "run_summary.json"
    summary = {
        "schema_version": 2,
        "status": "running",
        "initial_count": 20,
        "sequential_steps": sequential_steps,
        "final_sample_index": 20 + sequential_steps,
        "data_dir": str(data_dir),
        "initial_model_dir": str(initial_model),
        "initial_ledger_file": str(initial_ledger),
        "prior_config": asdict(prior_config),
        "independent_test_source": str(test_dataset) if test_dataset else None,
        "steps": [],
    }
    write_json_atomic(summary_path, summary)
    model_dir = initial_model
    ledger_path = initial_ledger
    final = 20 + sequential_steps
    for target in range(21, final + 1):
        print(f"Sequential sample {target}/{final}: checking sealed prediction")
        verify_sealed_prediction(train, designs, model_dir, target)
        step_dir = output_root / "steps" / f"sample_{target:04d}"
        diagnostic_dir = artifact_directory(step_dir, "diagnostics",
                                            "diagnostic_manifest.json")
        request_dir = request_directory(step_dir)
        policy_dir = artifact_directory(step_dir, "policy", "decision_record.json")
        # This is the first read of the current truth row in the online loop.
        diagnostic = diagnostics.run(
            train, designs, truth, model_dir, diagnostic_dir, target,
            ledger_path, data_dir / "initial20_context.json",
        )
        ensure_request("sequential", diagnostic_dir / "diagnostic_context.json", request_dir)
        ensure_llm_answer(request_dir, llm_settings)
        record = ensure_policy_record(
            policy_dir, request_dir, train, model_dir, designs, truth, target, final,
            gate_settings, restarts, seed,
        )
        if not record:
            record = policy.run(
                train=train, reference=data_dir / "initial20.csv",
                designs=designs, truth=truth,
                prior_context=data_dir / "initial20_context.json",
                model_dir=model_dir, request_dir=request_dir,
                ledger_path=ledger_path, output=policy_dir,
                settings=gate_settings, restarts=restarts, seed=seed,
                sequential_count=sequential_steps,
            )
        adopted_report = gp_model.read_json(Path(record["adopted_model_dir"]) /
                                            "fit_report.json")
        adopted_rules = gp_model.validate_rules(gp_model.read_json(
            Path(record["adopted_model_dir"]) / "rules_used.json"
        ))
        sealed_prediction = {
            response: {
                "truth_m": diagnostic["responses"][response]["truth_m"],
                "mean_m": diagnostic["responses"][response]["blind_mean_m"],
                "observed_std_m": diagnostic["responses"][response]["blind_observed_std_m"],
            }
            for response in gp_model.RESPONSES
        }
        summary["steps"].append({
            "sample_index": target,
            "design_id": diagnostic["design_id"],
            "llm_decision": record["llm_decision"],
            "outcome": record["outcome"],
            "proposed_candidates": record["proposed_candidates"],
            "adopted_candidates": record["adopted_candidates"],
            "baseline_validation_metrics": record["baseline_validation_metrics"],
            "adopted_validation_metrics": adopted_report["validation_metrics"],
            "adopted_prior_categories": adopted_rules["lengthscale_prior_categories"],
            "sealed_prediction": sealed_prediction,
            "decision_record": str(policy_dir / "decision_record.json"),
            "adopted_model_dir": record["adopted_model_dir"],
            "ledger_file": record["ledger_file"],
        })
        write_json_atomic(summary_path, summary)
        train = Path(record["new_training_file"])
        model_dir = Path(record["adopted_model_dir"])
        ledger_path = Path(record["ledger_file"])
    summary["online_metrics"] = prediction_metrics(
        [step["sealed_prediction"] for step in summary["steps"]],
        data_dir / "initial20.csv",
    )
    summary["outcome_counts"] = {
        outcome: sum(step["outcome"] == outcome for step in summary["steps"])
        for outcome in ("kept", "adopted", "rolled_back")
    }
    summary["proposed_adjustment_count"] = sum(
        len(step["proposed_candidates"]) for step in summary["steps"]
    )
    summary["adopted_adjustment_count"] = sum(
        len(step["adopted_candidates"]) for step in summary["steps"]
    )
    summary["final_model_dir"] = str(model_dir)
    summary["final_ledger_file"] = str(ledger_path)
    if test_dataset:
        holdout_data_dir = artifact_directory(
            output_root, "independent_test_data", "data_summary.json"
        )
        ensure_prepared_data(source_run, holdout_data_dir, sequential_steps,
                             test_dataset)
        for name in ("initial20.csv", "sequential_designs.csv", "sequential_truth.csv",
                     "initial20_context.json"):
            if (gp_model.file_sha256(holdout_data_dir / name)
                    != gp_model.file_sha256(data_dir / name)):
                raise ValueError("independent test preparation changed the online data")
        evaluation_dir = artifact_directory(
            output_root, "independent_test", "test_metrics.json"
        )
        evaluation = evaluate_independent_test(
            model_dir, holdout_data_dir / "holdout" / "test_holdout.csv",
            data_dir / "initial20.csv", evaluation_dir,
        )
        summary["independent_test_metrics"] = evaluation["metrics"]
        summary["independent_test_report"] = str(evaluation_dir / "test_metrics.json")
        summary["independent_test_status"] = "complete"
    else:
        summary["independent_test_status"] = "not_requested"
    summary["status"] = "complete"
    write_json_atomic(summary_path, summary)
    print(f"Completed {initial_count}+{sequential_steps}; final model: {model_dir}")
    print(f"Summary: {summary_path}")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=data_io.DEFAULT_RUN,
                        help="source run with dataset.csv and optimization_history.csv")
    parser.add_argument("--test-dataset", type=Path,
                        help="optional independent 40-point dataset; evaluated only after the final step")
    parser.add_argument("--output-root", type=Path,
                        help="isolated experiment directory; rerun with same settings to resume")
    parser.add_argument("--initial-count", type=int, default=20)
    parser.add_argument("--sequential-steps", type=int, default=10)
    parser.add_argument("--model", required=True, help="LLM provider model ID")
    parser.add_argument("--provider", choices=("openai", "openrouter", "custom"),
                        default="openai")
    parser.add_argument("--base-url", help="API base URL; required for custom provider")
    parser.add_argument("--api-style", choices=("responses", "chat-completions"))
    parser.add_argument("--api-key-env", help="environment variable containing the API key")
    parser.add_argument("--extra-body", type=Path,
                        help="JSON object of provider-specific request parameters")
    parser.add_argument("--http-referer", help="OpenRouter HTTP-Referer header")
    parser.add_argument("--app-title", help="OpenRouter X-OpenRouter-Title header")
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument("--retries", type=int, default=2)
    parser.add_argument("--restarts", type=int, default=3)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--prior-base-factor", type=float,
                        default=gp_model.DEFAULT_PRIOR_CONFIG.base_factor)
    parser.add_argument("--prior-log-std", type=float,
                        default=gp_model.DEFAULT_PRIOR_CONFIG.log_std)
    parser.add_argument("--min-absolute-nrmse-gain", type=float, default=0.0005)
    parser.add_argument("--min-relative-nrmse-gain", type=float, default=0.003)
    parser.add_argument("--max-response-relative-worsening", type=float, default=0.02)
    parser.add_argument("--max-nlpd-increase", type=float, default=0.02)
    parser.add_argument("--min-nlpd-gain-for-tie", type=float, default=0.005)
    parser.add_argument("--min-delta-map", type=float, default=-0.25)
    args = parser.parse_args()
    try:
        if args.provider == "custom":
            if not args.base_url or not args.api_key_env:
                raise ValueError("custom provider requires --base-url and --api-key-env")
            defaults = (args.base_url, "chat-completions", args.api_key_env)
        else:
            defaults = llm_io.PROVIDER_DEFAULTS[args.provider]
        if args.provider != "openrouter" and (args.http_referer or args.app_title):
            raise ValueError("OpenRouter headers require --provider openrouter")
        base_url = args.base_url or defaults[0]
        api_style = args.api_style or defaults[1]
        key_env = args.api_key_env or defaults[2]
        llm_io.endpoint(base_url, api_style)
        extra_body = prompt_io.read_json_file(args.extra_body) if args.extra_body else {}
        if not isinstance(extra_body, dict):
            raise ValueError("--extra-body must contain a JSON object")
        extra_headers = {}
        if args.http_referer:
            extra_headers["HTTP-Referer"] = args.http_referer
        if args.app_title:
            extra_headers["X-OpenRouter-Title"] = args.app_title
        llm_settings = {
            "model": args.model, "provider": args.provider,
            "base_url": base_url, "api_style": api_style,
            "key_env": key_env, "timeout": args.timeout,
            "retries": args.retries, "extra_body": extra_body,
            "extra_headers": extra_headers,
        }
        gate_settings = policy.GateSettings(
            min_absolute_nrmse_gain=args.min_absolute_nrmse_gain,
            min_relative_nrmse_gain=args.min_relative_nrmse_gain,
            max_response_relative_worsening=args.max_response_relative_worsening,
            max_nlpd_increase=args.max_nlpd_increase,
            min_nlpd_gain_for_tie=args.min_nlpd_gain_for_tie,
            min_delta_map=args.min_delta_map,
        )
        prior_config = gp_model.PriorConfig(
            base_factor=args.prior_base_factor,
            log_std=args.prior_log_std,
        )
        source_run = args.run.resolve()
        output_root = (args.output_root or HERE / "outputs" / "experiments" /
                       f"{source_run.name}_20_plus_{args.sequential_steps}_summary_taylor_v4").resolve()
        run_experiment(
            source_run=source_run, output_root=output_root,
            initial_count=args.initial_count, sequential_steps=args.sequential_steps,
            llm_settings=llm_settings, restarts=args.restarts,
            seed=args.seed, gate_settings=gate_settings,
            prior_config=prior_config,
            test_dataset=args.test_dataset,
        )
    except (ValueError, OSError, KeyError, IndexError, TypeError,
            gp_model.np.linalg.LinAlgError) as error:
        parser.exit(2, f"Sequential experiment failed: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
