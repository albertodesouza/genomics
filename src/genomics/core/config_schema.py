from __future__ import annotations

import ast
import importlib
import inspect
import json
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Type, get_args, get_origin


@dataclass(frozen=True)
class ConfigSchemaSpec:
    kind: str
    title: str
    model_module: str
    model_name: str
    loader_module: str
    loader_name: str
    path_markers: tuple[str, ...]


CONFIG_SCHEMA_SPECS: Dict[str, ConfigSchemaSpec] = {
    "genotype": ConfigSchemaSpec(
        kind="genotype",
        title="Genotype-Based Predictor",
        model_module="genomics.predictors.genotype_based.config",
        model_name="PipelineConfig",
        loader_module="genomics.predictors.genotype_based.config",
        loader_name="load_config",
        path_markers=("configs/predictors/genotype_based/", "/genotype_based/"),
    ),
    "variant": ConfigSchemaSpec(
        kind="variant",
        title="Variant Transformer Predictor",
        model_module="genomics.predictors.variant_transformer.config",
        model_name="PipelineConfig",
        loader_module="genomics.predictors.variant_transformer.config",
        loader_name="load_config",
        path_markers=("configs/predictors/variant_transformer/", "/variant_transformer/"),
    ),
}


def schema_kinds() -> List[str]:
    return sorted(CONFIG_SCHEMA_SPECS)


def get_schema_spec(kind: str) -> ConfigSchemaSpec:
    try:
        return CONFIG_SCHEMA_SPECS[kind]
    except KeyError as exc:
        raise ValueError(f"Unknown config kind: {kind}") from exc


def infer_schema_kind(path: Path) -> Optional[str]:
    normalized = str(Path(path)).replace("\\", "/")
    for spec in CONFIG_SCHEMA_SPECS.values():
        if any(marker in normalized for marker in spec.path_markers):
            return spec.kind
    return None


def load_schema_model(kind: str) -> Type[Any]:
    spec = get_schema_spec(kind)
    module = importlib.import_module(spec.model_module)
    return getattr(module, spec.model_name)


def load_typed_config(kind: str, path: Path) -> Any:
    spec = get_schema_spec(kind)
    module = importlib.import_module(spec.loader_module)
    loader = getattr(module, spec.loader_name)
    return loader(Path(path))


def json_schema_for_kind(kind: str) -> Dict[str, Any]:
    model = load_schema_model(kind)
    return model.model_json_schema()


def describe_config_schema(kind: str) -> List[Dict[str, Any]]:
    model = load_schema_model(kind)
    return list(_describe_model(model))


def _describe_model(model: Type[Any], prefix: str = "") -> Iterable[Dict[str, Any]]:
    field_docs = _field_docstrings(model)
    for name, field in model.model_fields.items():
        path = f"{prefix}.{name}" if prefix else name
        annotation = field.annotation
        nested = _nested_model_type(annotation)
        description = field.description or field_docs.get(name) or ""
        yield {
            "path": path,
            "type": _type_label(annotation),
            "required": bool(field.is_required()),
            "default": _default_label(field),
            "description": description,
        }
        if nested is not None:
            yield from _describe_model(nested, path)


def _field_docstrings(model: Type[Any]) -> Dict[str, str]:
    try:
        source = textwrap.dedent(inspect.getsource(model))
    except (OSError, TypeError):
        return {}
    try:
        module_ast = ast.parse(source)
    except SyntaxError:
        return {}
    class_node = next((node for node in module_ast.body if isinstance(node, ast.ClassDef)), None)
    if class_node is None:
        return {}

    docs: Dict[str, str] = {}
    previous_field: Optional[str] = None
    for node in class_node.body:
        field_name = _assigned_field_name(node)
        if field_name:
            previous_field = field_name
            continue
        if previous_field and isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
            docs[previous_field] = " ".join(node.value.value.split())
            previous_field = None
            continue
        previous_field = None
    return docs


def _assigned_field_name(node: ast.AST) -> Optional[str]:
    if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
        return node.target.id
    if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
        return node.targets[0].id
    return None


def _nested_model_type(annotation: Any) -> Optional[Type[Any]]:
    try:
        from pydantic import BaseModel
    except Exception:
        return None
    if inspect.isclass(annotation) and issubclass(annotation, BaseModel):
        return annotation
    return None


def _type_label(annotation: Any) -> str:
    origin = get_origin(annotation)
    args = get_args(annotation)
    if origin is None:
        return getattr(annotation, "__name__", str(annotation).replace("typing.", ""))
    if str(origin) == "typing.Literal":
        return " | ".join(json.dumps(arg, ensure_ascii=False) for arg in args)
    origin_name = getattr(origin, "__name__", str(origin).replace("typing.", ""))
    if origin_name == "Union":
        return " | ".join(_type_label(arg) for arg in args)
    if args:
        return f"{origin_name}[{', '.join(_type_label(arg) for arg in args)}]"
    return origin_name


def _default_label(field: Any) -> str:
    if field.is_required():
        return "<required>"
    if field.default_factory is not None:
        return "<generated>"
    if field.default is None:
        return "null"
    if isinstance(field.default, (str, int, float, bool, list, dict)):
        return json.dumps(field.default, ensure_ascii=False, default=str)
    return str(field.default)
