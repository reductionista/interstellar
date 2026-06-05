#!/usr/bin/env python3
"""
Lasair Kafka subscriber for LSST unlinked moving source detections.

Subscribes to the Lasair filter topic and accumulates per-epoch detections
into a CSV file suitable for heliolinx input. Designed to run continuously
(or in batches via --max-messages) and append to the same output file across
multiple runs without duplicating sources.

Usage:
    python python_scripts/lasair_kafka_subscriber.py
    python python_scripts/lasair_kafka_subscriber.py --max-messages 500
    python python_scripts/lasair_kafka_subscriber.py --output data/my_detections.csv

Output CSV columns:
    mjd, ra_deg, dec_deg, mag, psfFlux, band, obscode, obsid

psfFlux is the raw value in nJy (empty=no template/null, negative=darker than template).
mag is derived from psfFlux where positive; 99.9 otherwise.
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

from lasair import lasair_consumer

KAFKA_SERVER = 'lasair-lsst-kafka.lsst.ac.uk:9092'
TOPIC        = 'lasair_965nDiaSources5'
GROUP_ID     = 'interstellar_subscriber_01'
OBSCODE      = 'X05'   # Vera Rubin Observatory, Cerro Pachón

DEFAULT_OUTPUT = Path(__file__).parent.parent / 'data' / 'lasair_detections.csv'

FIELDNAMES = ['mjd', 'ra_deg', 'dec_deg', 'mag', 'psfFlux', 'band', 'obscode', 'obsid']

# AB magnitude zero point: Rubin psfFlux is in nJy, so m = 31.4 - 2.5*log10(flux_nJy)
FLUX_ZP = 31.4


def flux_to_mag(flux_njy):
    if flux_njy is None or flux_njy <= 0:
        return None
    return FLUX_ZP - 2.5 * math.log10(flux_njy)


def extract_detections(alert, seen_ids):
    """Return list of detection rows from a full-alert packet, skipping duplicates and linked sources."""
    # Kafka full-alert packets nest the source list under alert['alert']
    inner = alert.get('alert', alert)
    rows = []
    for src in inner.get('diaSourcesList', []):
        src_id = src.get('diaSourceId')
        if src_id in seen_ids:
            continue
        if src.get('ssObjectId') is not None:
            continue
        # Skip reliability filter — the ROB classifier is trained on transients,
        # not moving objects, so low scores don't mean spurious for our use case.
        mjd  = src.get('midpointMjdTai')
        ra   = src.get('ra')
        dec  = src.get('decl')
        if None in (mjd, ra, dec):
            continue
        flux = src.get('psfFlux')
        mag  = flux_to_mag(flux)
        rows.append({
            'mjd':     mjd,
            'ra_deg':  ra,
            'dec_deg': dec,
            'mag':     round(mag, 4) if mag is not None else 99.9,
            'psfFlux': round(flux, 4) if flux is not None else '',
            'band':    src.get('band', ''),
            'obscode': OBSCODE,
            'obsid':   src_id,
        })
        seen_ids.add(src_id)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output',       default=str(DEFAULT_OUTPUT), help='Output CSV path (appended)')
    parser.add_argument('--max-messages', type=int, default=None,      help='Stop after this many Kafka messages')
    parser.add_argument('--group-id',     default=GROUP_ID,            help='Kafka consumer group ID')
    args = parser.parse_args()

    outpath = Path(args.output)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    # Pre-load already-seen source IDs so re-runs don't duplicate rows
    seen_ids = set()
    if outpath.exists() and outpath.stat().st_size > 0:
        with open(outpath) as f:
            for row in csv.DictReader(f):
                try:
                    seen_ids.add(int(row['obsid']))
                except (KeyError, ValueError):
                    pass
        print(f"Loaded {len(seen_ids)} existing source IDs from {outpath}")

    write_header = not outpath.exists() or outpath.stat().st_size == 0
    outfile = open(outpath, 'a', newline='')
    writer  = csv.DictWriter(outfile, fieldnames=FIELDNAMES)
    if write_header:
        writer.writeheader()

    consumer = lasair_consumer(KAFKA_SERVER, args.group_id, TOPIC)
    print(f"Subscribed to {TOPIC} (group={args.group_id})")
    print("Ctrl-C to stop.\n")

    n_msg = n_det = 0
    try:
        while True:
            msg = consumer.poll(timeout=10)
            if msg is None:
                print(f"  [{n_msg} msgs, {n_det} detections] waiting for new alerts...")
                if args.max_messages and n_msg >= args.max_messages:
                    break
                continue
            if msg.error():
                print(f"Kafka error: {msg.error()}", file=sys.stderr)
                break

            try:
                alert = json.loads(msg.value())
            except json.JSONDecodeError as e:
                print(f"JSON error: {e}", file=sys.stderr)
                continue

            rows = extract_detections(alert, seen_ids)
            writer.writerows(rows)
            n_det += len(rows)
            n_msg += 1

            if n_msg % 100 == 0:
                outfile.flush()
                print(f"  {n_msg} messages, {n_det} detections written")

            if args.max_messages and n_msg >= args.max_messages:
                break

    except KeyboardInterrupt:
        print("\nStopped by user.")
    finally:
        outfile.flush()
        outfile.close()

    print(f"\nDone: {n_msg} messages, {n_det} new detections → {outpath}")


if __name__ == '__main__':
    main()
