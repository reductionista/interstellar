#!/usr/bin/env python3
"""
Combine multiple heliolinx-format detection CSVs into a single deduplicated file.

Usage:
    python python_scripts/combine_detections.py file1.csv file2.csv ... -o output.csv

Deduplicates by obsid (unique diaSource identifier). Output is sorted by MJD.
"""

import argparse
import sys
import pandas as pd
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Combine detection CSV files")
    parser.add_argument("inputs", nargs="+", help="Input CSV files")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    args = parser.parse_args()

    dfs = []
    reference_cols = None
    for path in args.inputs:
        df = pd.read_csv(path)
        if reference_cols is None:
            reference_cols = set(df.columns)
        elif set(df.columns) != reference_cols:
            missing = reference_cols - set(df.columns)
            extra = set(df.columns) - reference_cols
            print(f"ERROR: {path} has mismatched columns", file=sys.stderr)
            if missing:
                print(f"  Missing: {sorted(missing)}", file=sys.stderr)
            if extra:
                print(f"  Extra:   {sorted(extra)}", file=sys.stderr)
            sys.exit(1)
        print(f"  {path}: {len(df)} rows, MJD {df['mjd'].min():.1f}–{df['mjd'].max():.1f}")
        dfs.append(df)

    combined = pd.concat(dfs, ignore_index=True)
    before = len(combined)
    combined = combined.drop_duplicates(subset=["obsid"]).sort_values("mjd").reset_index(drop=True)
    dupes = before - len(combined)

    print(f"\nCombined: {before} rows → {len(combined)} after removing {dupes} duplicates")
    print(f"MJD range: {combined['mjd'].min():.1f}–{combined['mjd'].max():.1f}")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out, index=False)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
