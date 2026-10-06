from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple


YELLOW = "\033[33;1m"
RED = "\033[31;1m"
RESET = "\033[0m"


def _load_aligned_tsv(path: Path) -> Tuple[List[str], List[Tuple[str, str, str]]]:
    comments: List[str] = []
    rows: List[Tuple[str, str, str]] = []

    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                comments.append(line)
                continue
            parts = line.split("\t")
            if parts[0] == "sample_id":
                continue
            if len(parts) != 3:
                raise ValueError(f"Linha invalida em {path}: esperado 3 colunas, recebeu {len(parts)}")
            rows.append((parts[0], parts[1], parts[2]))

    if not rows:
        raise ValueError(f"Nenhuma sequencia encontrada em {path}")
    return comments, rows


def _format_base(base: str, ref_base: str, is_reference: bool, color: bool) -> str:
    if not color:
        return base
    if base == "X":
        return f"{YELLOW}{base}{RESET}"
    if not is_reference and base != ref_base:
        return f"{RED}{base}{RESET}"
    return base


def convert_aligned_tsv_to_columns(
    input_path: Path,
    output_path: Path,
    include_position: bool = True,
    color: bool = False,
) -> Path:
    comments, rows = _load_aligned_tsv(input_path)
    expected_length = len(rows[0][1])
    for sample_id, h1, h2 in rows:
        if len(h1) != expected_length or len(h2) != expected_length:
            raise ValueError(
                f"Comprimento inconsistente para {sample_id}: "
                f"H1={len(h1)}, H2={len(h2)}, esperado={expected_length}"
            )

    columns: Dict[str, str] = {}
    for sample_id, h1, h2 in rows:
        if sample_id == "REF":
            columns["REF"] = h1
        else:
            columns[f"{sample_id}_H1"] = h1
            columns[f"{sample_id}_H2"] = h2

    column_names = list(columns.keys())
    if include_position:
        header = ["pos"] + column_names
    else:
        header = column_names

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as out:
        for comment in comments:
            out.write(f"{comment}\n")
        out.write("\t".join(header) + "\n")

        ref_sequence = columns.get("REF")
        if ref_sequence is None:
            raise ValueError("Coluna REF ausente no arquivo de entrada")

        for idx in range(expected_length):
            ref_base = ref_sequence[idx]
            values = [
                _format_base(seq[idx], ref_base, column_name == "REF", color)
                for column_name, seq in columns.items()
            ]
            if include_position:
                out.write("\t".join([str(idx + 1)] + values) + "\n")
            else:
                out.write("\t".join(values) + "\n")

    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert aligned DNA TSV to position-by-column TXT")
    parser.add_argument("input_path", type=Path)
    parser.add_argument("output_path", type=Path)
    parser.add_argument("--no-position", action="store_true", help="omit the 1-based aligned position column")
    parser.add_argument("--color", action="store_true", help="highlight X in yellow and bases different from REF in red")
    args = parser.parse_args()

    output_path = convert_aligned_tsv_to_columns(
        args.input_path.resolve(),
        args.output_path.resolve(),
        include_position=not args.no_position,
        color=args.color,
    )
    print(output_path)


if __name__ == "__main__":
    main()
