#!/usr/bin/env python3
"""
Fetch and parse MPC 80-column astrometry for 3I/ATLAS (C/2025 N1) and convert
to a simple CSV format for heliolinc/make_tracklets.

Output CSV columns (1-indexed):
  1: obsid  2: MJD  3: RA_deg  4: Dec_deg  5: mag  6: band  7: obscode

Also writes the matching column format file for make_tracklets.

MPC 80-column format:
  Col 1-5:   packed minor planet number (or blank)
  Col 6:     column 5 designation suffix / discovery flag
  Col 8-12:  packed provisional designation
  Col 10-11: note1, note2
  Col 12:    note3
  Col 13:    type (C = CCD, S/s = satellite pair member)
  Col 15-32: date (YYYY MM DD.dddddd)
  Col 33-44: RA  HH MM SS.ddd
  Col 45:    sign of Dec
  Col 46-56: Dec DD MM SS.dd
  Col 66-70: magnitude
  Col 71-72: band
  Col 78-80: observatory code
"""

import urllib.request
import urllib.error
import math
import sys
import os

MPC_URL   = "https://www.minorplanetcenter.net/tmp2/3I.txt"
_TESTS = os.path.join(os.path.dirname(__file__), "../heliolinc2_src/tests")
OUTCSV    = os.path.join(_TESTS, "3I_detections.csv")
COLFMT    = os.path.join(_TESTS, "colformat_3I_01.txt")


def parse_mpc80_date(date_str):
    """Parse 'YYYY MM DD.ddddd' → MJD (UTC approximation)."""
    parts = date_str.split()
    year  = int(parts[0])
    month = int(parts[1])
    day_f = float(parts[2])
    day   = int(day_f)
    frac  = day_f - day

    # Julian Day Number for integer Y M D
    a = (14 - month) // 12
    y = year + 4800 - a
    m = month + 12 * a - 3
    jdn = day + (153*m + 2)//5 + 365*y + y//4 - y//100 + y//400 - 32045
    jd  = jdn - 0.5 + frac   # JD at 0h UT = JDN - 0.5; add fraction of day
    return jd - 2400000.5    # MJD


def parse_mpc80_ra(ra_str):
    """Parse 'HH MM SS.ddd' → decimal degrees."""
    parts = ra_str.split()
    h, m, s = int(parts[0]), int(parts[1]), float(parts[2])
    return (h + m/60.0 + s/3600.0) * 15.0


def parse_mpc80_dec(sign, dec_str):
    """Parse sign char + 'DD MM SS.dd' → decimal degrees."""
    parts = dec_str.split()
    d, m, s = int(parts[0]), int(parts[1]), float(parts[2])
    val = d + m/60.0 + s/3600.0
    return -val if sign == '-' else val


def parse_mpc80_line(line):
    """Parse one 80-column MPC observation line. Returns dict or None."""
    if len(line) < 80:
        return None
    # Skip if not a valid obs type (skip blank lines, header, etc.)
    obs_type = line[14]  # column 15 (0-indexed 14): C, S, s, etc.
    if obs_type not in ('C', 'S', 's', 'A', 'B', 'T', 'M', 'P', 'E', 'H'):
        return None

    try:
        date_str = line[15:32].strip()   # cols 16-32
        ra_str   = line[32:44].strip()   # cols 33-44
        sign     = line[44]              # col 45
        dec_str  = line[45:56].strip()   # cols 46-56
        mag_str  = line[65:70].strip()   # cols 66-70
        band     = line[70:72].strip()   # cols 71-72
        obscode  = line[77:80].strip()   # cols 78-80

        if not date_str or not ra_str or not dec_str:
            return None

        mjd   = parse_mpc80_date(date_str)
        ra    = parse_mpc80_ra(ra_str)
        dec   = parse_mpc80_dec(sign, dec_str)
        mag   = float(mag_str) if mag_str else 99.9
        band  = band if band else "V"

        return {
            "mjd": mjd, "ra": ra, "dec": dec,
            "mag": mag, "band": band, "obscode": obscode
        }
    except (ValueError, IndexError):
        return None


print(f"Fetching 3I/ATLAS MPC astrometry from {MPC_URL} ...")
try:
    with urllib.request.urlopen(MPC_URL, timeout=60) as resp:
        raw = resp.read().decode("utf-8", errors="replace")
except urllib.error.URLError as e:
    print(f"ERROR fetching MPC data: {e}", file=sys.stderr)
    sys.exit(1)

lines = raw.splitlines()
print(f"  Downloaded {len(lines)} lines")

obs_list = []
for i, line in enumerate(lines):
    obs = parse_mpc80_line(line)
    if obs:
        obs_list.append(obs)

print(f"  Parsed {len(obs_list)} valid observations")

if not obs_list:
    print("ERROR: no observations parsed — check MPC format", file=sys.stderr)
    sys.exit(1)

# Sort by MJD
obs_list.sort(key=lambda o: o["mjd"])

# Write CSV
with open(OUTCSV, "w") as f:
    f.write("obsid,mjd,ra_deg,dec_deg,mag,band,obscode\n")
    for i, obs in enumerate(obs_list):
        f.write(f"{i+1},{obs['mjd']:.6f},{obs['ra']:.6f},{obs['dec']:.6f},"
                f"{obs['mag']:.2f},{obs['band']},{obs['obscode']}\n")

print(f"Wrote {OUTCSV}")

# Write column format file (1-indexed column numbers matching the CSV)
with open(COLFMT, "w") as f:
    f.write("IDCOL 1\n")
    f.write("MJDCOL 2\n")
    f.write("RACOL 3\n")
    f.write("DECCOL 4\n")
    f.write("MAGCOL 5\n")
    f.write("BANDCOL 6\n")
    f.write("OBSCODECOL 7\n")

print(f"Wrote {COLFMT}")

# Print date range
mjds = [o["mjd"] for o in obs_list]
print(f"MJD range: {min(mjds):.2f} – {max(mjds):.2f}")
print(f"RA range:  {min(o['ra'] for o in obs_list):.3f}° – {max(o['ra'] for o in obs_list):.3f}°")
print(f"Dec range: {min(o['dec'] for o in obs_list):.3f}° – {max(o['dec'] for o in obs_list):.3f}°")
print("Done.")
