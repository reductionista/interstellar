#!/usr/bin/env python3
"""
Extract unlinked moving source detections from Fink LSST parquet files.

Reads parquet files from a fink_datatransfer output directory, filters for
unlinked sources (ssObjectId == 0) with positive psfFlux, and writes a CSV
suitable for heliolinx input.

Usage:
    python data_prep/fink_extract_detections.py \\
        --indir fink/ftransfer_lsst_2026-05-30_895376 \\
        --output data/fink_detections.csv

    # Include negative-flux detections (darker than template):
    python data_prep/fink_extract_detections.py \\
        --indir fink/ftransfer_lsst_2026-05-30_895376 \\
        --output data/fink_detections.csv \\
        --include-negative

Output CSV columns:
    mjd, ra_deg, dec_deg, mag, psfFlux, band, obscode, obsid

psfFlux is in nJy. mag = 31.4 - 2.5*log10(psfFlux) for positive flux; 99.9 otherwise.
ssObjectId == 0 means unlinked (Fink uses 0, not NULL, as the sentinel value).
"""

import argparse
import csv
import glob
import math
import sys
from pathlib import Path

import pandas as pd

OBSCODE = 'X05'   # Vera Rubin Observatory, Cerro Pachón
FLUX_ZP = 31.4    # AB mag zero point for Rubin nJy fluxes

FIELDNAMES = ['mjd', 'ra_deg', 'dec_deg', 'mag', 'psfFlux', 'band', 'obscode', 'obsid']


def flux_to_mag(flux_njy):
    if flux_njy is None or flux_njy <= 0:
        return None
    return FLUX_ZP - 2.5 * math.log10(flux_njy)


def extract_file(path, seen_ids, writer, include_negative):
    """Process one parquet file. Returns (n_written, n_skipped_linked, n_skipped_negative)."""
    df = pd.read_parquet(path)
    n_written = n_linked = n_negative = 0

    for _, row in df.iterrows():
        ds = row['diaSource']
        if ds is None:
            continue

        src_id = ds.get('diaSourceId')
        if src_id in seen_ids:
            continue

        # Skip sources linked to known solar system objects
        sso_id = ds.get('ssObjectId')
        if sso_id is not None and sso_id != 0:
            n_linked += 1
            continue

        flux = ds.get('psfFlux')
        is_negative = (flux is not None and flux <= 0)

        if is_negative and not include_negative:
            n_negative += 1
            continue

        mjd = ds.get('midpointMjdTai')
        ra  = ds.get('ra')
        dec = ds.get('dec')   # Fink uses 'dec', not 'decl'
        if None in (mjd, ra, dec):
            continue

        mag = flux_to_mag(flux)
        writer.writerow({
            'mjd':     mjd,
            'ra_deg':  ra,
            'dec_deg': dec,
            'mag':     round(mag, 4) if mag is not None else 99.9,
            'psfFlux': round(flux, 4) if flux is not None else '',
            'band':    ds.get('band', ''),
            'obscode': OBSCODE,
            'obsid':   src_id,
        })
        seen_ids.add(src_id)
        n_written += 1

    return n_written, n_linked, n_negative


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--indir',            required=True,
                        help='Directory containing Fink parquet files')
    parser.add_argument('--output',           default='data/fink_detections.csv')
    parser.add_argument('--include-negative', action='store_true',
                        help='Include detections with negative psfFlux (darker than template)')
    args = parser.parse_args()

    indir   = Path(args.indir)
    outpath = Path(args.output)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    files = sorted(indir.glob('**/*.parquet'))
    if not files:
        print(f"No parquet files found in {indir}", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(files)} parquet files in {indir}")

    # Pre-load seen IDs if output already exists
    seen_ids = set()
    if outpath.exists() and outpath.stat().st_size > 0:
        with open(outpath) as f:
            for row in __import__('csv').DictReader(f):
                try:
                    seen_ids.add(int(row['obsid']))
                except (KeyError, ValueError):
                    pass
        print(f"Loaded {len(seen_ids)} existing source IDs")

    write_header = not outpath.exists() or outpath.stat().st_size == 0
    outfile = open(outpath, 'a', newline='')
    writer  = csv.DictWriter(outfile, fieldnames=FIELDNAMES)
    if write_header:
        writer.writeheader()

    total_written = total_linked = total_negative = 0
    for i, f in enumerate(files):
        n_w, n_l, n_n = extract_file(f, seen_ids, writer, args.include_negative)
        total_written  += n_w
        total_linked   += n_l
        total_negative += n_n

        if (i + 1) % 500 == 0:
            outfile.flush()
            print(f"  {i+1}/{len(files)} files — "
                  f"{total_written} written, {total_linked} linked (skipped), "
                  f"{total_negative} negative (skipped)")

    outfile.flush()
    outfile.close()

    print(f"\nDone:")
    print(f"  {total_written} unlinked detections → {outpath}")
    print(f"  {total_linked} skipped (linked to known SSO)")
    if not args.include_negative:
        print(f"  {total_negative} skipped (negative psfFlux; rerun with --include-negative to keep)")


if __name__ == '__main__':
    main()
