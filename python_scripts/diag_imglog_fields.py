#!/usr/bin/env python3
"""Check hlimage field names and whether observer position is actually populated."""
import numpy as np
from pathlib import Path
import heliolinx
from heliolinx import solarsyst_dyn_geo as sdg

ROOT      = Path(__file__).parent.parent
TEST_DATA = ROOT / "test_data"

earth_raw = sdg.load_earth_ephemerides(TEST_DATA / "Earth1day2025.csv")
earthpos  = heliolinx.create_EarthState(len(earth_raw))
earthpos['MJD'] = earth_raw['MJD']
for src, dst in [('X','x'),('Y','y'),('Z','z'),('VX','vx'),('VY','vy'),('VZ','vz')]:
    earthpos[dst] = earth_raw[src]

obs_lookup = {}
with open(TEST_DATA / "ObsCodes_clean.txt") as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 4:
            try:
                obs_lookup[parts[0]] = (float(parts[1]), float(parts[2]), float(parts[3]))
            except ValueError:
                pass

# Create a single test image entry
img = heliolinx.create_hlimage(1)
print("hlimage dtype fields:", img.dtype.names)
print()

# Populate it for a known date with Rubin observatory (X05)
mjd_test = 61150.0
code = 'X05'
Long, pcos, psin = obs_lookup[code]
obsvec = heliolinx.observer_vel(mjd_test, Long, pcos, psin, earthpos)
print(f"observer_vel() returns (len={len(obsvec)}):")
print(f"  obsvec = {obsvec}")
print()

img[0]['MJD'] = mjd_test
img[0]['RA']  = 150.0
img[0]['Dec'] = 2.0
img[0]['obscode'] = b'X05'
img[0]['X']  = obsvec[0]; img[0]['Y']  = obsvec[1]; img[0]['Z']  = obsvec[2]
img[0]['VX'] = obsvec[3]; img[0]['VY'] = obsvec[4]; img[0]['VZ'] = obsvec[5]
img[0]['exptime'] = -1.0

print("After assignment, reading back:")
print(f"  img['X']  = {img['X']}")
print(f"  img['Y']  = {img['Y']}")
print(f"  img['Z']  = {img['Z']}")
print(f"  img['VX'] = {img['VX']}")
print(f"  img['VY'] = {img['VY']}")
print(f"  img['VZ'] = {img['VZ']}")
print()

# Try lowercase too
try:
    print(f"  img['x']  = {img['x']}")
    print(f"  img['vx'] = {img['vx']}")
except Exception as e:
    print(f"  lowercase fields not present: {e}")
