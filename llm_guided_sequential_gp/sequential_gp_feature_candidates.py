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
    config = primitives.DamGeometryConfig(geometry_template) if geometry_template else None
    namespace = primitives.build_formula_namespace(x, config=config)
    namespace = {
        name: values for name, values in namespace.items()
        if name in SYMBOL_UNITS or name in ALIAS_TO_CANONICAL
    }

    summary = [
        "统计仅来自样本 1–20。下列 Spearman 相关用于提出候选机制。",
        "",
        "| 响应 | 均值 (mm) | 样本标准差 (mm) | 观测范围 (mm) |",
        "|---|---:|---:|---:|",
    ]
    for key, values in zip(PROMPT_RESPONSES, 1000 * y.T, strict=True):
        summary.append(f"| `{key}` | {values.mean():.3f} | {values.std(ddof=1):.3f} | [{values.min():.3f}, {values.max():.3f}] |")
    summary.extend([
        "",
        "| 原语 | 单位 | 观测范围 | 定义 | ρ(C1) | ρ(C2) | ρ(C3) |",
        "|---|---|---:|---|---:|---:|---:|",
    ])
    response_ranks = pd.DataFrame(y).rank()
    for name, meta in primitives.PRIMITIVE_METADATA.items():
        values = namespace[name]
        ranks = pd.Series(values).rank()
        correlations = response_ranks.corrwith(ranks) if ranks.nunique() > 1 else [math.nan] * 3
        rho = " | ".join(f"{v:+.3f}" if math.isfinite(v) else "N/A" for v in correlations)
        summary.append(f"| `{name}` | {meta['unit']} | [{values.min():.3f}, {values.max():.3f}] | {meta['desc']} | {rho} |")

    template = TEMPLATE.read_text(encoding="utf-8")
    if template.count(PLACEHOLDER) != 1:
        raise ValueError("物理特征模板应包含一个原语统计占位符")
    prompt = template.replace(PLACEHOLDER, "\n".join(summary)).encode("utf-8")
    # 保存本地求值所需的几何快照。API 请求只使用渲染后的统计摘要。
    context = {"design_ids": design_ids, "namespace": {name: values.tolist() for name, values in namespace.items()}}
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
    context = prompt_io.exact_keys(prompt_io.read_json_file(request / "context.json"), {
        "design_ids", "namespace",
    }, "物理特征上下文")
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

    evaluator = primitives.SafeFormulaEvaluator(tuple(namespace), min_std=0)
    for node in ast.walk(tree):
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
            with np.errstate(all="ignore"):
                denominator = evaluator.evaluate("denominator", ast.unparse(node.right), namespace, n)
            if (not denominator.is_valid or denominator.values is None
                    or np.min(denominator.values) <= 0 <= np.max(denominator.values)
                    or (isinstance(node.right, ast.Name) and node.right.id in {"X_c_sec", "d_us_min"})):
                raise ValueError("分母在已观测样本上变号、触零或含未经保护的敏感净距/偏心量")
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


def validate_candidate(candidate: object, namespace: dict, basis: np.ndarray) -> tuple[dict, np.ndarray]:
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
    return feature, standardized


def equivalent_feature(a: np.ndarray, b: np.ndarray, scaled_a: np.ndarray, scaled_b: np.ndarray) -> bool:
    if (np.allclose(scaled_a, scaled_b, atol=1e-9, rtol=0)
            or np.allclose(scaled_a, -scaled_b, atol=1e-9, rtol=0)):
        return True
    # 倒数或经常数缩放的倒数在两个特征均有定义时，逐点乘积为非零常数。
    product = (a / np.max(np.abs(a))) * (b / np.max(np.abs(b)))
    return bool(product[0] != 0 and np.allclose(product, product[0], rtol=1e-9, atol=0))


def validate(request: Path, raw_reply: bytes) -> tuple[dict, Path]:
    manifest, context, _ = load_request(request)
    response_hash = prompt_io.digest(raw_reply)
    output = request / "responses" / response_hash
    output.mkdir(parents=True, exist_ok=True)
    (output / "llm_response_raw.txt").write_bytes(raw_reply)
    report = {"schema_version": 1, "prompt_sha256": manifest["prompt_sha256"],
              "response_sha256": response_hash, "valid": False, "errors": [], "checks": {}}
    candidates = {}
    try:
        parsed = json.loads(raw_reply.decode("utf-8-sig"), object_pairs_hook=prompt_io.unique_pairs,
                            parse_constant=prompt_io.reject_constant)
        answer = prompt_io.exact_keys(parsed, set(PROMPT_RESPONSES), "物理特征回复")
        namespace = {name: np.array(values, dtype=float) for name, values in context["namespace"].items()}
        x = np.column_stack([namespace[name] for name in DESIGN_VARIABLES])
        basis = np.column_stack([np.ones(len(x)), x - x.mean(axis=0)])
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
                    feature, standardized = validate_candidate(item, namespace, basis)
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
        report, output = validate(args.request, reply.encode("utf-8"))
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
    send = commands.add_parser("call", help="向提供商发送已渲染的提示词")
    send.add_argument("--request", type=Path, required=True)
    send.add_argument("--provider", default="openai", choices=("openai", "openrouter", "custom"))
    send.add_argument("--model", required=True)
    send.add_argument("--base-url")
    send.add_argument("--api-style", choices=("responses", "chat-completions"))
    send.add_argument("--api-key-env")
    send.add_argument("--timeout", type=float, default=120.0)
    send.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    try:
        if args.command == "render":
            result = render(args.initial_csv, args.output, args.geometry_template)
            print(f"{args.output / 'prompt.md'} | SHA-256 {result['prompt_sha256']}")
            return 0
        if args.command == "validate":
            result, output = validate(args.request, args.response.read_bytes())
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
