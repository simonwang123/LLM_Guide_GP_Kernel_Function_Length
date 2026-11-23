"""Build residual-based physical feature requests and validate their replies."""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np

if __package__:
    from . import sequential_gp_additive as additive
    from . import sequential_gp_feature_candidates as candidates
    from . import sequential_gp_feature_selection as features
    from . import sequential_gp_llm as llm
    from . import sequential_gp_model as gp
    from . import sequential_gp_primitives as primitives
else:
    import sequential_gp_additive as additive
    import sequential_gp_feature_candidates as candidates
    import sequential_gp_feature_selection as features
    import sequential_gp_llm as llm
    import sequential_gp_model as gp
    import sequential_gp_primitives as primitives


PROMPT_TEMPLATE = Path(__file__).resolve().parent / "prompts" / "gp_feature_evolution.md"
CONTEXT_TOKEN = "{{CONTEXT_JSON}}"
TOP_COUNT = 5
PROPOSAL_COUNT = 2
FEATURE_CAP = 3


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


def response_context(rows, x, y, namespace, formulas, baseline):
    residual = np.asarray(baseline["residual_m"])
    order = np.argsort(-np.abs(residual), kind="stable")
    hard = order[:TOP_COUNT]
    # Draw distinct controls from the lower-error half, then match in the current GP's feature space.
    low_error = np.argsort(np.abs(residual), kind="stable")
    controls = set(low_error[: max(TOP_COUNT, len(rows) // 2)]) - set(hard)
    geometry = x
    feature_values = {}
    for item in formulas:
        try:
            feature_values[item["name"]] = features.feature_values(item, namespace, len(x))
        except features.FORMULA_ERRORS:
            feature_values[item["name"]] = None
    available = [values for values in feature_values.values() if values is not None]
    if available:
        values = np.column_stack(available)
        keep, lo, span = additive.scale_physical(values)
        if keep.any():
            geometry = (values[:, keep] - lo) / span

    def item(index):
        return {
            "sample_index": int(rows[index]["sample_index"]), "design_id": rows[index]["design_id"],
            "original_variables": {name: float(namespace[name][index]) for name in gp.VARIABLES},
            "geometric_primitives": {name: float(namespace[name][index]) for name in primitives.PRIMITIVE_METADATA},
            "current_feature_values": {
                name: float(values[index]) if values is not None else None for name, values in feature_values.items()
            },
            "loo_prediction_m": baseline["predicted_m"][index], "observed_m": float(y[index]),
            "observed_minus_loo_m": float(residual[index]),
        }

    pairs = []
    if len(controls) < TOP_COUNT:
        raise ValueError("insufficient distinct low-residual controls")
    for index in hard:
        control = min(controls, key=lambda i: (np.sum((geometry[i] - geometry[index]) ** 2), i))
        controls.remove(control)
        pairs.append({"high_residual": item(index), "matched_low_residual": item(control)})
    signs = residual[hard]
    under = int(np.count_nonzero(signs > 0))
    over = int(np.count_nonzero(signs < 0))
    direction = "mostly_underestimated" if under >= 4 else "mostly_overestimated" if over >= 4 else "mixed_or_tied"
    return {
        "current_formulas": formulas, "loo_nrmse": baseline["nrmse"],
        "loo_formula_fallback": baseline.get("formula_fallback"),
        "top5_signs": {"underestimated": under, "overestimated": over, "zero": TOP_COUNT - under - over},
        "top5_direction": direction, "matched_pairs": pairs,
    }


def make_context(rows, x, y, namespace, active, baselines):
    return {
        "sample_count": len(rows), "response_unit": "m", "residual_definition": "observed minus LOO prediction",
        "original_variable_units": {name: "m" if candidates.SYMBOL_UNITS[name] == 1 else "1"
                                    for name in gp.VARIABLES},
        "aliases": candidates.ALIAS_TO_CANONICAL,
        "primitive_catalog": primitives.PRIMITIVE_METADATA,
        "responses": {
            response: response_context(rows, x, y[:, j], namespace, active[response], baselines[response])
            for j, response in enumerate(gp.RESPONSES)
        },
    }


def render_request(path: Path, context: dict) -> None:
    template = PROMPT_TEMPLATE.read_text(encoding="utf-8")
    if template.count(CONTEXT_TOKEN) != 1:
        raise ValueError("evolution prompt must contain exactly one context placeholder")
    payload = json.dumps(context, ensure_ascii=False, indent=2, allow_nan=False)
    prompt = template.replace(CONTEXT_TOKEN, payload).encode("utf-8")
    write_json(path / "context.json", context)
    write_bytes(path / "prompt.md", prompt)
    write_json(path / "request_manifest.json", {
        "schema_version": 1, "stage": "feature_evolution", "sample_count": context["sample_count"],
        "context_sha256": gp.canonical_sha256(context), "prompt_sha256": gp.file_sha256(path / "prompt.md"),
    })


def validate_reply(raw: bytes, x, namespace, active):
    options = {response: [] for response in gp.RESPONSES}
    audit = {response: [] for response in gp.RESPONSES}
    try:
        reply = json.loads(raw.decode("utf-8-sig"), object_pairs_hook=gp.unique_json_pairs,
                           parse_constant=gp.reject_json_constant)
        if not isinstance(reply, dict) or set(reply) != set(gp.RESPONSES):
            raise ValueError("reply must contain exactly the three response keys in the prompt")
    except (UnicodeError, ValueError) as exc:
        return options, audit, str(exc)

    basis = np.column_stack([np.ones(len(x)), x - x.mean(axis=0)])
    for response in gp.RESPONSES:
        rows = reply[response]
        if not isinstance(rows, list) or len(rows) != PROPOSAL_COUNT:
            audit[response].append({"status": "rejected", "reason": "expected exactly two candidate formulas"})
            continue
        names = {item["name"] for item in active[response]}
        for position, proposed in enumerate(rows, 1):
            check = {"position": position, "name": proposed.get("name") if isinstance(proposed, dict) else None}
            try:
                feature, _ = candidates.validate_candidate(proposed, namespace, basis)
                column = np.asarray(feature["values"])
                if not np.isfinite(column.std()) or column.std() <= features.std_floor(column):
                    raise ValueError("candidate is near constant on the revealed prefix")
                if feature["name"] in names:
                    raise ValueError("candidate name duplicates an existing or proposed feature")
                names.add(feature["name"])
                formula = {field: feature[field] for field in features.FORMULA_FIELDS}
                options[response].append((formula, column))
                check.update(status="valid", formula=formula, unit=feature["unit"])
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
