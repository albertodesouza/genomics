from __future__ import annotations

from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence


def _get_attr_or_key(value: Any, key: str, default: Any = None) -> Any:
    if isinstance(value, Mapping):
        return value.get(key, default)
    return getattr(value, key, default)


def derived_target_value(metadata: Mapping[str, Any], config: Any) -> Optional[str]:
    source_field = _get_attr_or_key(config, "source_field")
    class_map = _get_attr_or_key(config, "class_map", {}) or {}
    exclude_unmapped = bool(_get_attr_or_key(config, "exclude_unmapped", False))
    if not source_field:
        return None
    source_value = metadata.get(source_field)
    if source_value is None:
        return None
    for class_name, source_values in class_map.items():
        if source_value in source_values:
            return str(class_name)
    return None if exclude_unmapped else str(source_value)


def target_value(metadata: Mapping[str, Any], target: str, derived_targets: Optional[Mapping[str, Any]] = None) -> Optional[str]:
    if target == "frog_likelihood":
        return None
    if derived_targets and target in derived_targets:
        return derived_target_value(metadata, derived_targets[target])
    value = metadata.get(target)
    if value is None:
        return None
    return str(value)


def classes_from_known(known_classes: Optional[Sequence[str]]) -> Optional[List[str]]:
    if not known_classes:
        return None
    return [str(item) for item in known_classes]


def classes_from_distribution(metadata: Mapping[str, Any], target: str) -> Optional[List[str]]:
    key = f"{target}_distribution"
    distribution = metadata.get(key, {})
    if isinstance(distribution, Mapping) and distribution:
        return sorted(str(item) for item in distribution.keys())
    return None


def classes_from_records(records: Iterable[Mapping[str, Any]], target: str, derived_targets: Optional[Mapping[str, Any]] = None) -> List[str]:
    values = {
        value
        for value in (target_value(record, target, derived_targets) for record in records)
        if value is not None
    }
    return sorted(values)


def build_class_maps(classes: Sequence[str]) -> tuple[Dict[str, int], Dict[int, str]]:
    target_to_idx = {str(name): idx for idx, name in enumerate(classes)}
    idx_to_target = {idx: name for name, idx in target_to_idx.items()}
    return target_to_idx, idx_to_target
