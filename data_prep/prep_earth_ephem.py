#!/usr/bin/env python3
"""
Fetch Earth state vectors from JPL Horizons for the 3I/ATLAS observation period
and save in the 4-line block format expected by heliolinc's read_horizons_fileLD.

Format expected:
  $$SOE
  <JD> = A.D. <date> TDB
   X = <x km> Y = <y km> Z = <z km>
   VX= <vx km/s> VY= <vy km/s> VZ= <vz km/s>
   LT= <lt> RG= <rg> RR= <rr>
  $$EOE
"""

import os
import urllib.request
import urllib.parse
import sys

START = "2025-May-01"
STOP  = "2026-Jul-01"
OUTFILE = os.path.join(os.path.dirname(__file__), "../heliolinc2_src/tests/Earth1day2025.csv")

params = {
    "format":      "text",
    "COMMAND":     "'399'",
    "OBJ_DATA":    "'NO'",
    "MAKE_EPHEM":  "'YES'",
    "EPHEM_TYPE":  "'VECTORS'",
    "CENTER":      "'500@10'",
    "START_TIME":  f"'{START}'",
    "STOP_TIME":   f"'{STOP}'",
    "STEP_SIZE":   "'1 d'",
    "OUT_UNITS":   "'KM-S'",
    "VEC_TABLE":   "'3'",
    "VEC_CORR":    "'NONE'",
    "CSV_FORMAT":  "'YES'",
    "REF_PLANE":   "'ECLIPTIC'",
    "REF_SYSTEM":  "'J2000'",
}

url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urllib.parse.urlencode(params)

print(f"Fetching Earth ephemeris from JPL Horizons ({START} to {STOP})...")
print(f"URL: {url[:80]}...")

try:
    with urllib.request.urlopen(url, timeout=60) as resp:
        data = resp.read().decode("utf-8")
except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)

# The API returns JSON-wrapped text; strip the outer wrapper if present
if data.startswith("{"):
    import json
    obj = json.loads(data)
    text = obj.get("result", "")
else:
    text = data

if "$$SOE" not in text and "$$EOE" not in text:
    print("ERROR: Response does not contain $$SOE/$$EOE markers. Printing first 500 chars:")
    print(text[:500])
    sys.exit(1)

with open(OUTFILE, "w") as f:
    f.write(text)

# Count lines in SOE block
lines = text.splitlines()
in_block = False
n_epochs = 0
for line in lines:
    if "$$SOE" in line:
        in_block = True
        continue
    if "$$EOE" in line:
        break
    if in_block and line.strip() and "=" in line and not line.startswith(" X") and not line.startswith(" VX") and not line.startswith(" LT"):
        n_epochs += 1

print(f"Saved {OUTFILE} — {n_epochs} epochs")
print("Done.")
