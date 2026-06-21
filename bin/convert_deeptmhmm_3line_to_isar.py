#!/usr/bin/env python3
"""Convert DeepTMHMM 3-line output into the tabular format ISAR imports."""

from __future__ import annotations

import argparse
from pathlib import Path


REGION_NAMES = {
    "I": "inside",
    "O": "outside",
    "M": "TMhelix",
    "S": "signal",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert DeepTMHMM predicted_topologies.3line output to a four-column "
            "tab-separated file for IsoformSwitchAnalyzeR::analyzeDeepTMHMM()."
        )
    )
    parser.add_argument("input_3line", type=Path)
    parser.add_argument("output_tsv", type=Path)
    return parser.parse_args()


def topology_regions(topology: str):
    if not topology:
        return

    start = 1
    current = topology[0]

    for position, state in enumerate(topology[1:], start=2):
        if state != current:
            yield REGION_NAMES.get(current, current), start, position - 1
            start = position
            current = state

    yield REGION_NAMES.get(current, current), start, len(topology)


def main() -> None:
    args = parse_args()
    lines = args.input_3line.read_text().splitlines()

    rows: list[str] = []
    for index in range(0, len(lines), 3):
        if index + 2 >= len(lines) or not lines[index].startswith(">"):
            continue

        isoform_id = lines[index][1:].split()[0]
        topology = lines[index + 2].strip()
        for region_type, start, end in topology_regions(topology):
            rows.append(f"{isoform_id}\t{region_type}\t{start}\t{end}")

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.output_tsv.write_text("\n".join(rows) + ("\n" if rows else ""))


if __name__ == "__main__":
    main()
