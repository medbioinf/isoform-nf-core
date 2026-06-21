#!/usr/bin/env python3
"""Convert multi-FASTA IUPred2A output to the block format imported by ISAR."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert IUPred2A output from a multi-FASTA run into the header-block "
            "format expected by IsoformSwitchAnalyzeR::analyzeIUPred2A()."
        )
    )
    parser.add_argument("input_iupred", type=Path)
    parser.add_argument("input_fasta", type=Path)
    parser.add_argument("output_iupred", type=Path)
    return parser.parse_args()


def read_fasta_lengths(path: Path) -> list[tuple[str, int]]:
    records: list[tuple[str, int]] = []
    current_id: str | None = None
    current_length = 0

    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                records.append((current_id, current_length))
            current_id = line[1:].split()[0]
            current_length = 0
        else:
            current_length += len(line)

    if current_id is not None:
        records.append((current_id, current_length))

    return records


def read_prediction_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 4:
            continue
        rows.append(parts[:4])
    return rows


def main() -> None:
    args = parse_args()
    records = read_fasta_lengths(args.input_fasta)
    rows = read_prediction_rows(args.input_iupred)
    expected_rows = sum(length for _, length in records)

    if expected_rows != len(rows):
        raise SystemExit(
            f"Cannot split IUPred2A rows by FASTA lengths: expected {expected_rows}, found {len(rows)}"
        )

    output_lines: list[str] = []
    offset = 0
    for isoform_id, length in records:
        output_lines.append(f">{isoform_id}")
        for index, row in enumerate(rows[offset : offset + length], start=1):
            _, residue, iupred_score, anchor_score = row
            output_lines.append(f"{index}\t{residue}\t{iupred_score}\t{anchor_score}")
        offset += length

    args.output_iupred.parent.mkdir(parents=True, exist_ok=True)
    args.output_iupred.write_text("\n".join(output_lines) + "\n")


if __name__ == "__main__":
    main()
