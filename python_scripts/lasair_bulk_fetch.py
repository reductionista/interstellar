#!/usr/bin/env python3
"""
Bulk historical detection fetcher from Lasair LSST.

Queries the Lasair objects table for moving source candidates and fetches
per-epoch detection data for each, writing to a CSV suitable for heliolinx.
Checkpoints progress so it can be stopped and resumed without re-fetching.

Usage:
    # Test: last 2 weeks, small sky patch, stay within 100/hr public limit
    python data_prep/lasair_bulk_fetch.py \\
        --mjd-min 61145 --mjd-max 61160 \\
        --ra-min 60 --ra-max 120 \\
        --rate-limit 90

    # Full sky, all post-Rubin data (after heavy-user approval):
    python data_prep/lasair_bulk_fetch.py --rate-limit 9000

Output CSV columns: mjd, ra_deg, dec_deg, mag, band, obscode, obsid
"""

import argparse
import collections
import csv
import json
import math
import sys
import time
from pathlib import Path

import settings
from lasair import LasairError, lasair_client as lasair

ENDPOINT  = 'https://api.lasair.lsst.ac.uk/api'
OBSCODE   = 'X05'   # Vera Rubin Observatory, Cerro Pachón
FLUX_ZP   = 31.4    # AB mag zero point for Rubin nJy fluxes
PAGE_SIZE = 1000    # objects per query() call

DEFAULT_OUTPUT     = Path(__file__).parent.parent / 'data' / 'lasair_detections.csv'
DEFAULT_CHECKPOINT = Path(__file__).parent.parent / 'data' / 'lasair_bulk_fetch_checkpoint.json'

FIELDNAMES = ['mjd', 'ra_deg', 'dec_deg', 'mag', 'psfFlux', 'band', 'obscode', 'obsid']


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def flux_to_mag(flux_njy):
    if flux_njy is None or flux_njy <= 0:
        return None
    return FLUX_ZP - 2.5 * math.log10(flux_njy)


def build_conditions(args):
    parts = [f'nDiaSources <= {args.ndia_max}',
             f'firstDiaSourceMjdTai > {args.mjd_min}']
    if args.mjd_max is not None:
        parts.append(f'firstDiaSourceMjdTai < {args.mjd_max}')
    if args.ra_min is not None:
        parts.append(f'ra >= {args.ra_min}')
    if args.ra_max is not None:
        parts.append(f'ra < {args.ra_max}')
    if args.dec_min is not None:
        parts.append(f'decl >= {args.dec_min}')
    if args.dec_max is not None:
        parts.append(f'decl < {args.dec_max}')
    return ' AND '.join(parts)


class RateLimiter:
    """Rolling-window rate limiter: allows at most max_calls calls per 3600 s."""
    def __init__(self, max_calls_per_hour):
        self.max_calls = max_calls_per_hour
        self.timestamps = collections.deque()

    def wait(self):
        now = time.monotonic()
        # Drop calls older than 1 hour
        while self.timestamps and now - self.timestamps[0] > 3600:
            self.timestamps.popleft()
        if len(self.timestamps) >= self.max_calls:
            sleep_for = 3600 - (now - self.timestamps[0]) + 0.1
            print(f"  [rate limit] sleeping {sleep_for:.1f}s ...")
            time.sleep(sleep_for)
        self.timestamps.append(time.monotonic())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output',      default=str(DEFAULT_OUTPUT))
    parser.add_argument('--checkpoint',  default=str(DEFAULT_CHECKPOINT))
    parser.add_argument('--rate-limit',  type=int,   default=90,
                        help='Max API calls per hour (default 90, safe for public tier)')
    parser.add_argument('--mjd-min',     type=float, default=61097.0,
                        help='Minimum firstDiaSourceMjdTai (default: Rubin start)')
    parser.add_argument('--mjd-max',     type=float, default=None,
                        help='Maximum firstDiaSourceMjdTai')
    parser.add_argument('--ra-min',      type=float, default=None)
    parser.add_argument('--ra-max',      type=float, default=None)
    parser.add_argument('--dec-min',     type=float, default=None)
    parser.add_argument('--dec-max',     type=float, default=None)
    parser.add_argument('--ndia-max',    type=int,   default=5,
                        help='Max nDiaSources per object (default 5)')
    parser.add_argument('--reset',       action='store_true',
                        help='Ignore existing checkpoint and start from scratch')
    args = parser.parse_args()

    outpath   = Path(args.output)
    ckpt_path = Path(args.checkpoint)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    conditions = build_conditions(args)
    print(f"Filter conditions: {conditions}")

    # Load checkpoint
    offset = 0
    n_processed = 0
    n_detections = 0
    if not args.reset and ckpt_path.exists():
        ckpt = json.loads(ckpt_path.read_text())
        if ckpt.get('conditions') == conditions:
            offset       = ckpt.get('offset', 0)
            n_processed  = ckpt.get('n_processed', 0)
            n_detections = ckpt.get('n_detections', 0)
            print(f"Resuming from offset {offset} ({n_processed} objects already done)")
        else:
            print("Checkpoint conditions differ — starting fresh")

    # Load already-seen source IDs from output CSV (prevents row duplication on reset/overlap)
    seen_source_ids = set()
    if outpath.exists() and outpath.stat().st_size > 0:
        with open(outpath) as f:
            for row in csv.DictReader(f):
                try:
                    seen_source_ids.add(int(row['obsid']))
                except (KeyError, ValueError):
                    pass
        print(f"Loaded {len(seen_source_ids)} existing source IDs")

    write_header = not outpath.exists() or outpath.stat().st_size == 0
    outfile = open(outpath, 'a', newline='')
    writer  = csv.DictWriter(outfile, fieldnames=FIELDNAMES)
    if write_header:
        writer.writeheader()

    L       = lasair(settings.lasair_token, endpoint=ENDPOINT)
    rl      = RateLimiter(args.rate_limit)
    t_start = time.time()

    print(f"Rate limit: {args.rate_limit} calls/hour  |  Output: {outpath}")
    print("Ctrl-C to stop (progress is checkpointed).\n")

    try:
        while True:
            # --- fetch one page of object IDs ---
            rl.wait()
            try:
                page = L.query('diaObjectId', 'objects', conditions,
                               limit=PAGE_SIZE, offset=offset)
            except LasairError as e:
                print(f"query() error at offset {offset}: {e}", file=sys.stderr)
                break

            if not page:
                print("No more objects — fetch complete.")
                break

            object_ids = [row['diaObjectId'] for row in page]

            # --- fetch per-epoch detections for each object ---
            for oid in object_ids:
                rl.wait()
                try:
                    result = L.object(oid, lite=True)
                except LasairError as e:
                    print(f"  object({oid}) error: {e}", file=sys.stderr)
                    continue

                rows_written = 0
                for src in result.get('diaSourcesList', []):
                    src_id = src.get('diaSourceId')
                    if src_id in seen_source_ids:
                        continue
                    if src.get('ssObjectId') is not None:
                        continue
                    mjd = src.get('midpointMjdTai')
                    ra  = src.get('ra')
                    dec = src.get('decl')
                    if None in (mjd, ra, dec):
                        continue
                    flux = src.get('psfFlux')
                    mag  = flux_to_mag(flux)
                    writer.writerow({
                        'mjd':     mjd,
                        'ra_deg':  ra,
                        'dec_deg': dec,
                        'mag':     round(mag, 4) if mag is not None else 99.9,
                        'psfFlux': round(flux, 4) if flux is not None else '',
                        'band':    src.get('band', ''),
                        'obscode': OBSCODE,
                        'obsid':   src_id,
                    })
                    seen_source_ids.add(src_id)
                    rows_written += 1

                n_detections += rows_written
                n_processed  += 1

            offset += len(page)

            # Checkpoint and progress report after every page
            ckpt_path.write_text(json.dumps({
                'conditions':   conditions,
                'offset':       offset,
                'n_processed':  n_processed,
                'n_detections': n_detections,
            }, indent=2))
            outfile.flush()

            elapsed    = time.time() - t_start
            rate_actual = n_processed / (elapsed / 3600) if elapsed > 0 else 0
            print(f"  offset={offset}  processed={n_processed}  "
                  f"detections={n_detections}  rate={rate_actual:.0f}/hr  "
                  f"elapsed={elapsed/3600:.2f}hr")

            if len(page) < PAGE_SIZE:
                print("Last page reached — fetch complete.")
                break

    except KeyboardInterrupt:
        print("\nStopped. Progress saved to checkpoint.")
    finally:
        outfile.flush()
        outfile.close()

    print(f"\nDone: {n_processed} objects, {n_detections} detections → {outpath}")


if __name__ == '__main__':
    main()
