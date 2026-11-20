"""Send one saved GP prompt to an LLM API, then validate the returned JSON.

Examples (run from the project root):
  uv run python llm_guided_sequential_gp/sequential_gp_llm.py --request llm_guided_sequential_gp/outputs/prompt_requests/initial --model YOUR_MODEL --dry-run

  uv run python llm_guided_sequential_gp/sequential_gp_llm.py --request llm_guided_sequential_gp/outputs/prompt_requests/initial --model YOUR_MODEL

Use --api-style chat-completions and --base-url for a compatible provider.
For OpenRouter, use --provider openrouter --model PROVIDER/MODEL and set
OPENROUTER_API_KEY in the current process environment.
The API key is read only from the named environment variable. This script does
not create prompts, fit models, or silently regenerate invalid LLM answers.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

try:
    from . import sequential_gp_prompt as prompt_io
except ImportError:  # Direct execution from this directory.
    import sequential_gp_prompt as prompt_io


RETRYABLE_HTTP = {408, 409, 429, 500, 502, 503, 504}
RESERVED_BODY_KEYS = {
    "model", "input", "messages", "stream", "n", "tools", "tool_choice",
    "background", "conversation", "previous_response_id",
}
PROVIDER_DEFAULTS = {
    "openai": ("https://api.openai.com/v1", "responses", "OPENAI_API_KEY"),
    "openrouter": (
        "https://openrouter.ai/api/v1", "chat-completions", "OPENROUTER_API_KEY"
    ),
}


class APIError(Exception):
    def __init__(self, message: str, *, status: int | None = None, retryable: bool = False):
        super().__init__(message)
        self.status = status
        self.retryable = retryable


def write_json(path: Path, value: dict) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def endpoint(base_url: str, api_style: str) -> str:
    base = base_url.strip().rstrip("/")
    parsed = urlsplit(base)
    if parsed.username or parsed.password or parsed.query or parsed.fragment:
        raise ValueError("base URL must not contain credentials, query, or fragment")
    is_local = parsed.hostname in {"localhost", "127.0.0.1", "::1"}
    if parsed.scheme != "https" and not (parsed.scheme == "http" and is_local):
        raise ValueError("base URL must use HTTPS, except for a local test server")
    if not parsed.hostname:
        raise ValueError("base URL has no host")
    path = "/responses" if api_style == "responses" else "/chat/completions"
    if parsed.path.rstrip("/").endswith(("/responses", "/chat/completions")):
        raise ValueError("pass the API base URL, not a full endpoint URL")
    return base + path


def request_body(model: str, prompt: str, style: str, extra: dict) -> dict:
    body = {"model": model}
    if style == "responses":
        body["input"] = [{"role": "user", "content": prompt}]
    else:
        body["messages"] = [{"role": "user", "content": prompt}]
    body.update(extra)
    return body


def post_json(
    url: str, body: dict, api_key: str, timeout: float,
    extra_headers: dict[str, str] | None = None,
) -> tuple[dict, int]:
    payload = json.dumps(body, ensure_ascii=False, allow_nan=False).encode("utf-8")
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "Accept": "application/json",
    }
    headers.update(extra_headers or {})
    request = Request(
        url,
        data=payload,
        method="POST",
        headers=headers,
    )
    try:
        with urlopen(request, timeout=timeout) as response:
            raw = response.read()
            status = response.status
    except HTTPError as error:
        raw = error.read(2048)
        detail = raw.decode("utf-8", errors="replace").replace(api_key, "[REDACTED]")
        raise APIError(
            f"HTTP {error.code}: {detail[:500]}",
            status=error.code,
            retryable=error.code in RETRYABLE_HTTP,
        ) from error
    except (URLError, TimeoutError, OSError) as error:
        raise APIError(f"network error: {error}", retryable=True) from error
    try:
        result = json.loads(raw.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as error:
        raise APIError(f"HTTP {status} returned invalid JSON: {error}", status=status) from error
    if not isinstance(result, dict):
        raise APIError(f"HTTP {status} returned a non-object JSON value", status=status)
    return result, status


def extract_answer(payload: dict, style: str) -> tuple[str, str | None, object, str | None]:
    """Return assistant text, provider ID, usage object, and provider model."""
    if style == "responses":
        if payload.get("status") != "completed":
            raise APIError(f"Responses API status is {payload.get('status')!r}, not completed")
        outputs = payload.get("output")
        if not isinstance(outputs, list):
            raise APIError("Responses API output is missing")
        parts = []
        for item in outputs:
            if not isinstance(item, dict) or item.get("type") != "message":
                continue
            for content in item.get("content") or []:
                if isinstance(content, dict) and content.get("type") == "output_text":
                    if isinstance(content.get("text"), str):
                        parts.append(content["text"])
        answer = "".join(parts)
    else:
        choices = payload.get("choices")
        if not isinstance(choices, list) or len(choices) != 1:
            raise APIError("Chat Completions must return exactly one choice")
        choice = choices[0]
        if not isinstance(choice, dict):
            raise APIError("Chat Completions choice must be an object")
        if choice.get("finish_reason") != "stop":
            raise APIError(f"Chat Completions finish_reason is {choice.get('finish_reason')!r}")
        message = choice.get("message")
        answer = message.get("content") if isinstance(message, dict) else None
    if not isinstance(answer, str) or not answer.strip():
        raise APIError("the API returned no assistant text to validate")
    return answer, payload.get("id"), payload.get("usage"), payload.get("model")


def run(
    request_dir: Path, *, model: str, base_url: str, api_style: str,
    key_env: str, timeout: float, retries: int, extra_body: dict,
    dry_run: bool, provider: str = "openai",
    extra_headers: dict[str, str] | None = None,
) -> dict:
    manifest, _, prompt_bytes = prompt_io.load_request(request_dir)
    if not model.strip():
        raise ValueError("--model must be a nonempty model ID")
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", key_env):
        raise ValueError("--api-key-env must name one environment variable")
    if timeout <= 0 or retries < 0:
        raise ValueError("--timeout must be positive and --retries must be nonnegative")
    if provider not in {"openai", "openrouter", "custom"}:
        raise ValueError(f"unknown provider: {provider}")
    for name, value in (extra_headers or {}).items():
        if not isinstance(value, str) or not value.isascii() or "\r" in value or "\n" in value:
            raise ValueError(f"{name} header must be ASCII text without line breaks")
    url = endpoint(base_url, api_style)
    if RESERVED_BODY_KEYS & set(extra_body):
        raise ValueError(f"extra body cannot override: {sorted(RESERVED_BODY_KEYS & set(extra_body))}")
    prompt = prompt_bytes.decode("utf-8")
    body = request_body(model, prompt, api_style, extra_body)
    body_hash = prompt_io.digest(prompt_io.canonical_json(body))
    if dry_run:
        return {
            "request_id": manifest["request_id"], "stage": manifest["stage"],
            "endpoint": url, "provider": provider, "api_style": api_style, "model": model,
            "request_body_sha256": body_hash, "dry_run": True,
        }

    key = os.environ.get(key_env, "")
    if not key:
        raise ValueError(f"set the {key_env} environment variable before calling the API")
    for name in ("api_call.json", "provider_response.json", "llm_response_raw.txt"):
        if (request_dir / name).exists():
            raise ValueError(
                f"{name} already exists; use a fresh request directory to avoid a duplicate API call"
            )

    attempts = []
    started = time.monotonic()
    call = {
        "schema_version": 1,
        "request_id": manifest["request_id"],
        "stage": manifest["stage"],
        "prompt_sha256": manifest["prompt_sha256"],
        "context_sha256": manifest["context_sha256"],
        "request_body_sha256": body_hash,
        "api_style": api_style,
        "provider": provider,
        "endpoint": url,
        "requested_model": model,
        "api_key_environment_variable": key_env,
        "started_at_utc": datetime.now(UTC).isoformat(),
        "attempts": attempts,
        "status": "pending",
    }
    for number in range(1, retries + 2):
        attempt_started = time.monotonic()
        try:
            provider_response, http_status = post_json(
                url, body, key, timeout, extra_headers
            )
            attempts.append({"number": number, "http_status": http_status,
                             "elapsed_seconds": round(time.monotonic() - attempt_started, 3)})
            # Save the full provider object before extraction or schema checks.
            write_json(request_dir / "provider_response.json", provider_response)
            text, provider_id, usage, provider_model = extract_answer(
                provider_response, api_style
            )
            response_path = request_dir / "llm_response_raw.txt"
            response_path.write_bytes(text.encode("utf-8"))
            validation = prompt_io.validate_response(request_dir, response_path)
            call.update(
                status="validated" if validation["valid"] else "invalid_llm_output",
                http_status=http_status,
                provider_id=provider_id,
                provider_model=provider_model,
                usage=usage,
                validation_valid=validation["valid"],
                validation_error=validation["error"],
            )
            break
        except APIError as error:
            attempts.append({"number": number, "http_status": error.status,
                             "elapsed_seconds": round(time.monotonic() - attempt_started, 3),
                             "error": str(error)})
            if not error.retryable or number > retries:
                call.update(status="api_error", error=str(error), http_status=error.status)
                break
            time.sleep(min(2 ** (number - 1), 8))
        except ValueError as error:
            call.update(status="local_validation_error", error=str(error))
            break
    call["elapsed_seconds"] = round(time.monotonic() - started, 3)
    call["completed_at_utc"] = datetime.now(UTC).isoformat()
    write_json(request_dir / "api_call.json", call)
    return call


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Send one saved GP prompt to an LLM API and validate its JSON answer.",
        epilog=(
            "OpenRouter: set OPENROUTER_API_KEY, then pass "
            "--provider openrouter --model PROVIDER/MODEL. "
            "Use --dry-run to inspect the endpoint before sending."
        ),
    )
    parser.add_argument("--request", required=True, type=Path,
                        help="directory created by sequential_gp_prompt.py render")
    parser.add_argument("--model", required=True, help="provider model ID")
    parser.add_argument("--provider", choices=("openai", "openrouter", "custom"),
                        default="openai", help="sets default URL, API style, and key variable")
    parser.add_argument("--base-url", help="API root; overrides the provider default")
    parser.add_argument("--api-style", choices=("responses", "chat-completions"),
                        help="overrides the provider default")
    parser.add_argument("--api-key-env",
                        help="overrides the provider key variable; value is never logged")
    parser.add_argument("--http-referer", help="optional OpenRouter site URL header")
    parser.add_argument("--app-title", help="optional OpenRouter app title header")
    parser.add_argument("--timeout", type=float, default=120.0, help="request timeout in seconds")
    parser.add_argument("--retries", type=int, default=2,
                        help="retries for temporary HTTP/network errors only")
    parser.add_argument("--extra-body", type=Path,
                        help="optional JSON object of provider-specific request parameters")
    parser.add_argument("--dry-run", action="store_true",
                        help="verify request and show endpoint without network access or API key")
    args = parser.parse_args()
    try:
        if args.provider == "custom":
            if not args.base_url or not args.api_key_env:
                raise ValueError("custom provider requires --base-url and --api-key-env")
            defaults = (args.base_url, "chat-completions", args.api_key_env)
        else:
            defaults = PROVIDER_DEFAULTS[args.provider]
        base_url = args.base_url or defaults[0]
        api_style = args.api_style or defaults[1]
        key_env = args.api_key_env or defaults[2]
        if args.provider != "openrouter" and (args.http_referer or args.app_title):
            raise ValueError("--http-referer and --app-title require --provider openrouter")
        extra_headers = {}
        if args.http_referer:
            extra_headers["HTTP-Referer"] = args.http_referer
        if args.app_title:
            extra_headers["X-OpenRouter-Title"] = args.app_title
        extra = prompt_io.read_json_file(args.extra_body) if args.extra_body else {}
        if not isinstance(extra, dict):
            raise ValueError("--extra-body must contain a JSON object")
        result = run(
            args.request.resolve(), model=args.model, base_url=base_url,
            api_style=api_style, key_env=key_env,
            timeout=args.timeout, retries=args.retries, extra_body=extra,
            dry_run=args.dry_run, provider=args.provider, extra_headers=extra_headers,
        )
    except (ValueError, UnicodeError) as error:
        parser.exit(2, f"LLM call failed before sending: {error}\n")
    if result.get("dry_run"):
        print(f"Dry run: {result['stage']} request {result['request_id']}")
        print(f"POST {result['endpoint']} | provider={result['provider']} | model={result['model']}")
        print(f"Request body SHA-256: {result['request_body_sha256']}")
        return 0
    print(f"API status: {result['status']}; attempts: {len(result['attempts'])}")
    print(f"Call record: {args.request.resolve() / 'api_call.json'}")
    if result["status"] in {"validated", "invalid_llm_output"}:
        print(f"LLM text: {args.request.resolve() / 'llm_response_raw.txt'}")
        print(f"Validation: {args.request.resolve() / 'validation_result.json'}")
    if result.get("error") or result.get("validation_error"):
        print(f"Error: {result.get('error') or result.get('validation_error')}")
    return 0 if result["status"] == "validated" else 2


if __name__ == "__main__":
    raise SystemExit(main())
