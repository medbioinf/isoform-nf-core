#!/usr/bin/env python3
"""Convert DeepLoc2 output into the table shape imported by ISAR."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ISAR_COLUMNS = [
    "Protein_ID",
    "Localizations",
    "Signals",
    "Cytoplasm",
    "Nucleus",
    "Extracellular",
    "Cell membrane",
    "Mitochondrion",
    "Plastid",
    "Endoplasmic reticulum",
    "Lysosome/Vacuole",
    "Golgi apparatus",
    "Peroxisome",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert DeepLoc2 CSV output to the 13-column format expected by "
            "IsoformSwitchAnalyzeR::analyzeDeepLoc2(). DeepLoc2.1 adds membrane "
            "type columns that ISAR 2.10 does not import directly."
        )
    )
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    with args.input_csv.open(newline="") as input_handle:
        reader = csv.DictReader(input_handle)
        if reader.fieldnames is None:
            raise SystemExit(f"Input CSV is empty: {args.input_csv}")

        missing = [column for column in ISAR_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise SystemExit(
                "Input does not look like DeepLoc2 output. Missing columns: "
                + ", ".join(missing)
            )

        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.output_csv.open("w", newline="") as output_handle:
            writer = csv.DictWriter(output_handle, fieldnames=ISAR_COLUMNS)
            writer.writeheader()
            for row in reader:
                writer.writerow({column: row[column] for column in ISAR_COLUMNS})


if __name__ == "__main__":
    main()
