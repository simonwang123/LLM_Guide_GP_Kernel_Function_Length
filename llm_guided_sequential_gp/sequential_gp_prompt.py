"""Render GP prompts and strictly validate cold-start or sequential LLM JSON.

This standalone experimental script does not call an LLM or fit a GP. Running it
without arguments renders the cold-start prompt using the current 20-point data
output. For a sequential prompt, supply the current diagnostic context.

The LLM never has to echo a hash. A request directory stores the exact prompt,
context snapshot, and local hashes alongside the raw answer and validation.

The answer contracts use schema version 2. This script checks structure and
current prior categories; evidence quality and model acceptance belong to the
later decision pipeline.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
PROMPTS = HERE / "prompts"
DEFAULT_INITIAL_CONTEXT = (
    HERE / "outputs" / "semiparametric_guarded_20_20261109_20_plus_40_summary_taylor_v4"
    / "initial20_context.json"
)
RESPONSES = (
    "crest_max_settlement_m",
    "dam_max_settlement_m",
    "dam_max_horizontal_displacement_m",
)
VARIABLES = (
    "upstream_slope",
    "downstream_slope",
    "crest_width",
    "modulus_zone_height",
    "secondary_top_width",
    "secondary_height",
    "secondary_bottom_width",
    "secondary_right_slope_distance",
)
PRIOR_CATEGORIES = frozenset({"short", "base", "long"})
INITIAL_KEYS = {"schema_version", "status", "assignments", "error"}
SEQUENTIAL_KEYS = {
    "schema_version", "status", "sample_index", "design_id", "decision",
    "candidates", "reason", "error",
}
CANDIDATE_KEYS = {"response", "variable", "from_category", "to_category", "reason"}
INITIAL_CONTEXT_KEYS = {
    "schema_version", "initial_count", "initial_observations_sha256",
    "design_variables", "design_bounds_physical", "input_scale", "responses",
    "prior_sign", "prior_tier", "prior_tier_meaning", "statistics",
}
SEQUENTIAL_CONTEXT_KEYS = {
    "schema_version", "stage", "sample_index", "design_id",
    "training_sample_count", "statistics", "coverage", "responses", "current_model",
    "engineering_prior", "recent_interventions", "recent_interventions_status",
    "diagnostic_limit",
}


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_json(value: object) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def read_json_file(path: Path) -> object:
    try:
        return json.loads(path.read_text(encoding="utf-8-sig"),
                          object_pairs_hook=unique_pairs, parse_constant=reject_constant)
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(f"cannot read valid JSON from {path}: {error}") from error


def unique_pairs(pairs: list[tuple[str, object]]) -> dict:
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON number: {value}")


_JSON_CODE_FENCE_RE = re.compile(
    r"```[ \t]*(?:json\b)?[ \t]*(?:\r?\n)?([\s\S]*?)\s*```",
    flags=re.IGNORECASE,
)


def strip_json_code_fence(text: str) -> str:
    """Unwrap one outer Markdown ``` or ```json fence when the reply is wrapped as a code block."""
    stripped = text.strip()
    if stripped.startswith("```"):
        match = _JSON_CODE_FENCE_RE.fullmatch(stripped)
        if match is not None:
            return match.group(1).strip()
    return stripped


def parse_llm_json(raw: bytes | str) -> object:
    """Decode UTF-8-SIG LLM output, unwrap an outer JSON code fence, and parse strictly."""
    if isinstance(raw, (bytes, bytearray)):
        text = bytes(raw).decode("utf-8-sig")
    elif isinstance(raw, str):
        text = raw.lstrip("\ufeff")
    else:
        raise TypeError("LLM response must be bytes or str")
    return json.loads(
        strip_json_code_fence(text),
        object_pairs_hook=unique_pairs,
        parse_constant=reject_constant,
    )


def exact_keys(value: object, expected: set[str], label: str) -> dict:
    if not isinstance(value, dict) or set(value) != expected:
        actual = set(value) if isinstance(value, dict) else type(value).__name__
        raise ValueError(f"{label} keys must be {sorted(expected)}; got {actual}")
    return value


def nonempty_string(value: object, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a nonempty string")
    return value


def finite_number(value: object, label: str) -> float:
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f"{label} must be a finite number")
    return float(value)


def named_numbers(value: object, names: tuple[str, ...], label: str,
                  *, unit: bool = False, positive: bool = False) -> dict:
    values = exact_keys(value, set(names), label)
    for name in names:
        number = finite_number(values[name], f"{label}/{name}")
        if unit and not -1e-10 <= number <= 1 + 1e-10:
            raise ValueError(f"{label}/{name} must be within [0, 1]")
        if positive and number <= 0:
            raise ValueError(f"{label}/{name} must be positive")
    return values


def validate_prior_matrix(value: object, label: str) -> None:
    if (not isinstance(value, list) or len(value) != len(VARIABLES)
            or any(not isinstance(row, list) or len(row) != len(RESPONSES)
                   for row in value)):
        raise ValueError(f"{label} must have eight variable rows and three response columns")
    for i, row in enumerate(value):
        for j, number in enumerate(row):
            finite_number(number, f"{label}/{i}/{j}")


def validate_statistics(value: object) -> None:
    stats = exact_keys(value, {
        "partial_spearman_design_response", "delta_loo_r2_quad_univariate", "interpretation",
    }, "statistics")
    nonempty_string(stats["interpretation"], "statistics interpretation")
    for key in ("partial_spearman_design_response", "delta_loo_r2_quad_univariate"):
        matrix = stats[key]
        if (not isinstance(matrix, list) or len(matrix) != len(VARIABLES)
                or any(not isinstance(row, list) or len(row) != len(RESPONSES)
                       for row in matrix)):
            raise ValueError(f"{key} must be an 8 by 3 matrix")
        for i, row in enumerate(matrix):
            for j, value in enumerate(row):
                if value is not None:
                    number = finite_number(value, f"{key}/{i}/{j}")
                    if key == "partial_spearman_design_response" and not -1 <= number <= 1:
                        raise ValueError(f"{key}/{i}/{j} must be within [-1, 1]")


def validate_taylor(value: object, label: str) -> None:
    result = exact_keys(value, {
        "method", "neighbor_distances_unit_space", "mean_absolute_contribution_total_m",
        "no_dominant_contribution", "top3", "approximation_checks",
    }, label)
    nonempty_string(result["method"], f"{label} method")
    distances = result["neighbor_distances_unit_space"]
    if not isinstance(distances, list) or len(distances) != 3:
        raise ValueError(f"{label} requires three neighbor distances")
    for distance in distances:
        if finite_number(distance, f"{label} distance") < 0:
            raise ValueError(f"{label} distance must be nonnegative")
    total = finite_number(result["mean_absolute_contribution_total_m"], f"{label} total")
    if total < 0 or type(result["no_dominant_contribution"]) is not bool:
        raise ValueError(f"{label} has invalid total or dominance flag")
    if result["no_dominant_contribution"] != (total <= 1e-12):
        raise ValueError(f"{label} dominance flag disagrees with total")
    top3 = result["top3"]
    if not isinstance(top3, list) or len(top3) != (0 if result["no_dominant_contribution"] else 3):
        raise ValueError(f"{label} must have zero or three ranked variables")
    seen, previous = set(), math.inf
    for row in top3:
        item = exact_keys(row, {
            "variable", "mean_absolute_contribution_m", "share_of_all_absolute_contributions",
            "mean_signed_contribution_m", "direction_counts",
        }, f"{label} top3 item")
        variable = item["variable"]
        if variable not in VARIABLES or variable in seen:
            raise ValueError(f"{label} repeats or misnames a variable")
        seen.add(variable)
        magnitude = finite_number(item["mean_absolute_contribution_m"], f"{label} magnitude")
        share = finite_number(item["share_of_all_absolute_contributions"], f"{label} share")
        finite_number(item["mean_signed_contribution_m"], f"{label} signed mean")
        if magnitude < 0 or magnitude > previous or not 0 <= share <= 1:
            raise ValueError(f"{label} has invalid ranking or share")
        if not math.isclose(share * total, magnitude, rel_tol=1e-8, abs_tol=1e-12):
            raise ValueError(f"{label} share disagrees with magnitude")
        previous = magnitude
        counts = exact_keys(item["direction_counts"], {"positive", "negative", "zero"},
                            f"{label} direction counts")
        if any(type(count) is not int or count < 0 for count in counts.values()) or sum(counts.values()) != 3:
            raise ValueError(f"{label} direction counts must sum to three")
    checks = result["approximation_checks"]
    if not isinstance(checks, list) or len(checks) != 3:
        raise ValueError(f"{label} requires three approximation checks")
    for rank, raw in enumerate(checks, 1):
        check = exact_keys(raw, {
            "neighbor_rank", "predicted_increment_m", "first_order_increment_m",
            "remainder_m", "sum_absolute_contributions_m", "poor_approximation",
        }, f"{label} approximation")
        if check["neighbor_rank"] != rank or type(check["poor_approximation"]) is not bool:
            raise ValueError(f"{label} approximation rank or flag is invalid")
        for key in ("predicted_increment_m", "first_order_increment_m", "remainder_m",
                    "sum_absolute_contributions_m"):
            finite_number(check[key], f"{label} {key}")


def validate_context(stage: str, context: object) -> dict:
    if not isinstance(context, dict):
        raise ValueError("context must be a JSON object")
    if stage == "initial":
        exact_keys(context, INITIAL_CONTEXT_KEYS, "initial context")
        if type(context.get("schema_version")) is not int or context["schema_version"] != 2:
            raise ValueError("initial context must have schema_version = 2")
        if type(context.get("initial_count")) is not int or context["initial_count"] != 20:
            raise ValueError("initial context must contain exactly 20 initial samples")
        if context.get("responses") != list(RESPONSES):
            raise ValueError("initial context response order does not match the template")
        if context.get("design_variables") != list(VARIABLES):
            raise ValueError("initial context variable order does not match the template")
        value = context.get("initial_observations_sha256")
        if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
            raise ValueError("initial context has no valid observations hash")
        for name in ("prior_sign", "prior_tier", "statistics"):
            if name not in context:
                raise ValueError(f"initial context is missing {name}")
        validate_prior_matrix(context["prior_sign"], "prior_sign")
        validate_prior_matrix(context["prior_tier"], "prior_tier")
        bounds = exact_keys(context["design_bounds_physical"], set(VARIABLES),
                            "design_bounds_physical")
        for variable, pair in bounds.items():
            if (not isinstance(pair, list) or len(pair) != 2
                    or finite_number(pair[0], f"{variable} lower bound") >=
                    finite_number(pair[1], f"{variable} upper bound")):
                raise ValueError(f"{variable} must have increasing physical bounds")
        validate_statistics(context["statistics"])
    elif stage == "sequential":
        exact_keys(context, SEQUENTIAL_CONTEXT_KEYS, "sequential context")
        if context.get("schema_version") != 3 or context.get("stage") != "sequential":
            raise ValueError("sequential context must have schema_version = 3 and stage = sequential")
        sample = context.get("sample_index")
        if type(sample) is not int or sample <= 20:
            raise ValueError("sequential context needs an integer sample_index > 20")
        design_id = nonempty_string(context.get("design_id"), "sequential context design_id")
        if type(context.get("training_sample_count")) is not int or context["training_sample_count"] != sample - 1:
            raise ValueError("sequential training_sample_count must be sample_index - 1")
        validate_statistics(context.get("statistics"))
        model = context.get("current_model")
        if not isinstance(model, dict) or model.get("lengthscale_rule_mode") != "soft_log_normal":
            raise ValueError("sequential context requires a fitted soft_log_normal prior model")
        response_states = exact_keys(context.get("responses"), set(RESPONSES), "sequential responses")
        for response in RESPONSES:
            state = response_states[response]
            if not isinstance(state, dict):
                raise ValueError(f"{response} diagnostic must be an object")
            categories = exact_keys(state.get("lengthscale_prior_categories"),
                                    set(VARIABLES), f"{response} prior categories")
            if any(type(category) is not str or category not in PRIOR_CATEGORIES
                   for category in categories.values()):
                raise ValueError(f"{response} has an invalid length-scale prior category")
            named_numbers(state.get("fitted_length_scales_unit_space"), VARIABLES,
                          f"{response} fitted length scales", positive=True)
            finite_number(state.get("truth_m"), f"{response} truth_m")
            finite_number(state.get("blind_mean_m"), f"{response} blind_mean_m")
            if finite_number(state.get("blind_observed_std_m"),
                             f"{response} blind_observed_std_m") <= 0:
                raise ValueError(f"{response} blind_observed_std_m must be positive")
            validate_taylor(state.get("taylor_attribution"), f"{response} Taylor attribution")
        if not isinstance(context.get("coverage"), dict):
            raise ValueError("sequential context is missing coverage")
        prior = context.get("engineering_prior")
        if not isinstance(prior, dict) or "prior_sign" not in prior or "prior_tier" not in prior:
            raise ValueError("sequential context is missing engineering_prior")
        validate_prior_matrix(prior["prior_sign"], "engineering_prior/prior_sign")
        validate_prior_matrix(prior["prior_tier"], "engineering_prior/prior_tier")
        status = context.get("recent_interventions_status")
        interventions = context.get("recent_interventions")
        if type(status) is not str or status not in {"provided", "not_provided"}:
            raise ValueError("sequential context has invalid recent_interventions_status")
        if status == "provided" and not isinstance(interventions, list):
            raise ValueError("provided recent_interventions must be an array")
        if status == "not_provided" and interventions is not None:
            raise ValueError("recent_interventions must be null when not provided")
    else:
        raise ValueError(f"unknown stage: {stage}")
    return context


def render(stage: str, context_file: Path, output: Path, template_file: Path | None) -> dict:
    context = validate_context(stage, read_json_file(context_file))
    common_path = PROMPTS / "gp_common_rules.md"
    stage_path = template_file or PROMPTS / (
        "gp_initial_tiers.md" if stage == "initial" else "gp_sequential_diagnosis.md"
    )
    try:
        common = common_path.read_text(encoding="utf-8")
        stage_template = stage_path.read_text(encoding="utf-8")
    except OSError as error:
        raise ValueError(f"prompt template is unavailable: {error}") from error
    placeholder = "{{INITIAL20_CONTEXT_JSON}}" if stage == "initial" else "{{SEQUENTIAL_CONTEXT_JSON}}"
    markers = re.findall(r"\{\{[A-Z][A-Z0-9_]*\}\}", stage_template)
    if markers != [placeholder]:
        raise ValueError(f"stage template must contain only one {placeholder}; found {markers}")
    context_text = json.dumps(context, ensure_ascii=False, indent=2, allow_nan=False)
    prompt = common.rstrip() + "\n\n---\n\n" + stage_template.replace(
        placeholder, context_text
    ).rstrip() + "\n"
    context_hash = digest(canonical_json(context))
    prompt_hash = digest(prompt.encode("utf-8"))
    manifest = {
        "schema_version": 1,
        "stage": stage,
        "request_id": f"{stage}-{prompt_hash[:16]}",
        "context_sha256": context_hash,
        "prompt_sha256": prompt_hash,
        "common_template_sha256": digest(common.encode("utf-8")),
        "stage_template_sha256": digest(stage_template.encode("utf-8")),
        "context_source": str(context_file.resolve()),
        "stage_template_source": str(stage_path.resolve()),
    }
    manifest_path = output / "request_manifest.json"
    if output.exists():
        if manifest_path.exists():
            old_manifest = read_json_file(manifest_path)
            if not isinstance(old_manifest, dict) or old_manifest.get("prompt_sha256") != prompt_hash:
                raise ValueError(f"output contains a different request; choose another --output: {output}")
        elif any(output.iterdir()):
            raise ValueError(f"output directory is not an empty request directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    (output / "context.json").write_text(context_text + "\n", encoding="utf-8")
    # Write exact UTF-8 bytes: text-mode newline conversion on Windows would
    # otherwise change the prompt hash after it was recorded in the manifest.
    (output / "prompt.md").write_bytes(prompt.encode("utf-8"))
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def validate_initial(answer: dict) -> None:
    exact_keys(answer, INITIAL_KEYS, "initial answer")
    if type(answer["schema_version"]) is not int or answer["schema_version"] != 2:
        raise ValueError("initial answer schema_version must be integer 2")
    if answer["status"] == "INVALID_INPUT":
        if answer["assignments"] != []:
            raise ValueError("INVALID_INPUT must have an empty assignments array")
        nonempty_string(answer["error"], "initial error")
        return
    if answer["status"] != "OK" or answer["error"] is not None:
        raise ValueError("initial status must be OK with error = null, or INVALID_INPUT")
    rows = answer["assignments"]
    expected = [(response, variable) for response in RESPONSES for variable in VARIABLES]
    if not isinstance(rows, list) or len(rows) != len(expected):
        raise ValueError("initial assignments must contain exactly 24 records")
    for index, ((response, variable), row) in enumerate(zip(expected, rows, strict=True), 1):
        exact_keys(row, {"response", "variable", "prior_category", "reason"},
                   f"assignment {index}")
        if (row["response"], row["variable"]) != (response, variable):
            raise ValueError(f"assignment {index} must be {response}/{variable}")
        if type(row["prior_category"]) is not str or row["prior_category"] not in PRIOR_CATEGORIES:
            raise ValueError(f"assignment {index} has an illegal prior category")
        nonempty_string(row["reason"], f"assignment {index} reason")


def validate_sequential(answer: dict, context: dict) -> None:
    exact_keys(answer, SEQUENTIAL_KEYS, "sequential answer")
    if type(answer["schema_version"]) is not int or answer["schema_version"] != 2:
        raise ValueError("sequential answer schema_version must be integer 2")
    if type(answer["sample_index"]) is not int or answer["sample_index"] != context["sample_index"]:
        raise ValueError("sequential answer sample_index does not match this request")
    if answer["design_id"] != context["design_id"]:
        raise ValueError("sequential answer design_id does not match this request")
    if answer["status"] == "INVALID_INPUT":
        if answer["decision"] is not None or answer["candidates"] != [] or answer["reason"] is not None:
            raise ValueError("INVALID_INPUT must have null decision/reason and empty candidates")
        nonempty_string(answer["error"], "sequential error")
        return
    if answer["status"] != "OK" or answer["error"] is not None:
        raise ValueError("sequential status must be OK with error = null, or INVALID_INPUT")
    nonempty_string(answer["reason"], "sequential reason")
    if answer["decision"] == "KEEP":
        if answer["candidates"] != []:
            raise ValueError("KEEP must have no candidates")
        return
    if answer["decision"] != "MODIFY":
        raise ValueError("sequential decision must be KEEP or MODIFY")
    candidates = answer["candidates"]
    if not isinstance(candidates, list) or not 1 <= len(candidates) <= 3:
        raise ValueError("MODIFY must contain one to three candidates")
    seen = set()
    for index, raw in enumerate(candidates, 1):
        item = exact_keys(raw, CANDIDATE_KEYS, f"candidate {index}")
        response, variable = item["response"], item["variable"]
        if (type(response) is not str or type(variable) is not str
                or response not in RESPONSES or variable not in VARIABLES):
            raise ValueError(f"candidate {index} has an unknown response or variable")
        pair = response, variable
        if pair in seen:
            raise ValueError(f"candidate {index} repeats {response}/{variable}")
        seen.add(pair)
        current = context["responses"][response]["lengthscale_prior_categories"][variable]
        if item["from_category"] != current:
            raise ValueError(f"candidate {index} from_category differs from current prior")
        if type(item["to_category"]) is not str or item["to_category"] not in PRIOR_CATEGORIES:
            raise ValueError(f"candidate {index} has an illegal to_category")
        if item["to_category"] == current:
            raise ValueError(f"candidate {index} does not change the prior category")
        nonempty_string(item["reason"], f"candidate {index} reason")


def load_request(request_dir: Path) -> tuple[dict, dict, bytes]:
    """Read a rendered request and verify its saved context and prompt bytes."""
    manifest = read_json_file(request_dir / "request_manifest.json")
    exact_keys(manifest, {
        "schema_version", "stage", "request_id", "context_sha256", "prompt_sha256",
        "common_template_sha256", "stage_template_sha256", "context_source",
        "stage_template_source",
    }, "request manifest")
    if manifest["schema_version"] != 1 or manifest["stage"] not in {"initial", "sequential"}:
        raise ValueError("request manifest has no valid stage")
    stage = manifest["stage"]
    context = validate_context(stage, read_json_file(request_dir / "context.json"))
    if digest(canonical_json(context)) != manifest.get("context_sha256"):
        raise ValueError("saved context does not match the request manifest")
    try:
        prompt = (request_dir / "prompt.md").read_bytes()
    except OSError as error:
        raise ValueError(f"cannot read saved prompt: {error}") from error
    if digest(prompt) != manifest.get("prompt_sha256"):
        raise ValueError("saved prompt does not match the request manifest")
    return manifest, context, prompt


def validate_response(request_dir: Path, response_file: Path) -> dict:
    manifest, context, _ = load_request(request_dir)
    stage = manifest["stage"]
    try:
        raw = response_file.read_bytes()
        text = raw.decode("utf-8-sig")
    except (OSError, UnicodeError) as error:
        raise ValueError(f"cannot read LLM response: {error}") from error
    (request_dir / "llm_response_raw.txt").write_bytes(raw)
    result = {
        "schema_version": 1,
        "request_id": manifest["request_id"],
        "stage": stage,
        "context_sha256": manifest["context_sha256"],
        "prompt_sha256": manifest["prompt_sha256"],
        "response_sha256": digest(raw),
        "valid": False,
        "answer": None,
        "error": None,
    }
    try:
        answer = parse_llm_json(text)
        if not isinstance(answer, dict):
            raise ValueError("LLM response must be one JSON object")
        if stage == "initial":
            validate_initial(answer)
        else:
            validate_sequential(answer, context)
        result["valid"] = True
        result["answer"] = answer
    except (json.JSONDecodeError, ValueError, TypeError) as error:
        result["error"] = str(error)
    (request_dir / "validation_result.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command")
    render_parser = sub.add_parser("render", help="assemble and save a prompt request")
    render_parser.add_argument("--stage", choices=("initial", "sequential"), default="initial")
    render_parser.add_argument("--context", type=Path)
    render_parser.add_argument("--template", type=Path,
                               help="stage template override; useful before the sequential template is written")
    render_parser.add_argument("--output", type=Path)
    check_parser = sub.add_parser("validate", help="validate a raw LLM response for a saved request")
    check_parser.add_argument("--request", required=True, type=Path)
    check_parser.add_argument("--response", required=True, type=Path)
    args = parser.parse_args(sys.argv[1:] or ["render"])
    try:
        if args.command == "render":
            context = args.context or DEFAULT_INITIAL_CONTEXT
            if args.stage == "sequential" and args.context is None:
                raise ValueError("--context is required for the sequential stage")
            parsed = validate_context(args.stage, read_json_file(context))
            name = "initial" if args.stage == "initial" else f"sequential_{parsed['sample_index']:04d}"
            output = (args.output or HERE / "outputs" / "prompt_requests" / name).resolve()
            manifest = render(args.stage, context.resolve(), output,
                              args.template.resolve() if args.template else None)
            print(f"Request: {manifest['request_id']}")
            print(f"Prompt:  {output / 'prompt.md'}")
            print(f"Context: {output / 'context.json'}")
            print(f"Manifest: {output / 'request_manifest.json'}")
            return 0
        if args.command == "validate":
            result = validate_response(args.request.resolve(), args.response.resolve())
            print(f"Stage: {result['stage']}; valid: {result['valid']}")
            print(f"Result: {args.request.resolve() / 'validation_result.json'}")
            if result["error"]:
                print(f"Error: {result['error']}")
            return 0 if result["valid"] else 2
        parser.print_help()
        return 2
    except ValueError as error:
        parser.exit(2, f"Prompt processing failed: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
