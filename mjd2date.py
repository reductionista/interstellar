#!/usr/bin/env python3
"""Convert MJD values to human-readable UTC dates."""

import sys
from datetime import datetime, timezone, timedelta

MJD_EPOCH = datetime(1858, 11, 17, tzinfo=timezone.utc)


def mjd_to_datetime(mjd: float) -> datetime:
    return MJD_EPOCH + timedelta(days=mjd)


def main():
    args = sys.argv[1:]
    if not args:
        print("Usage: mjd2date.py <mjd> [<mjd> ...]")
        sys.exit(1)
    for arg in args:
        mjd = float(arg)
        dt = mjd_to_datetime(mjd)
        print(f"MJD {mjd:.1f}  →  {dt.strftime('%Y-%m-%d %H:%M UTC')}  ({dt.strftime('%b %-d, %Y')})")


if __name__ == "__main__":
    main()
