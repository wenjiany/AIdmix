#!/usr/bin/env python3
"""Combine complete one-sample AIdmix TSV files with header validation."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    output = Path(args.output).resolve()
    paths = sorted(
        path for path in Path(args.input_dir).glob("*.tsv")
        if path.resolve() != output
    )
    if not paths:
        raise SystemExit(f"No TSV files found in {args.input_dir}")

    fields: list[str] | None = None
    rows: list[dict[str, str]] = []
    for path in paths:
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            current_fields = list(reader.fieldnames or [])
            current_rows = list(reader)
        if len(current_rows) != 1:
            raise SystemExit(f"Expected one row in {path}; found {len(current_rows)}")
        if fields is None:
            fields = current_fields
        elif fields != current_fields:
            raise SystemExit(f"Header differs in {path}")
        rows.append(current_rows[0])

    rows.sort(key=lambda row: row.get("sample", ""))
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + f".tmp.{os.getpid()}")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, output)
    print(f"Combined {len(rows)} samples into {output}")


if __name__ == "__main__":
    main()
