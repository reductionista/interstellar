#!/usr/bin/env python3
"""
Generate a hypothesis grid for heliolinc_interstellar covering hyperbolic inbound orbits.

Grid axes:
  r    = heliocentric distance (AU): range typical for near-Earth detections
  rdot = radial velocity (AU/day): negative = inbound, magnitude > v_escape → hyperbolic

For a bound orbit:  |rdot| < v_escape(r) = sqrt(2*GM_sun / r)
For hyperbolic:     rdot < -v_escape(r)  (inbound)

v_escape(r) in AU/day:
  GM_sun = 132712440041.279 km^3/s^2 = 132712440041.279 / (149597870.7)^3 AU^3/s^2
         = 2.9591220828 × 10^-4 AU^3/day^2  (using 1 day = 86400 s)
  v_esc(r) = sqrt(2 * 2.9591220828e-4 / r) AU/day

The third column is the mean radial acceleration in units of GM_sun/r^2.
For a purely hyperbolic (unperturbed) trajectory ≈ 0 (the code uses it as a perturbation
term around the Keplerian; set to 1.0 as per existing test files).

3I/ATLAS parameters for reference:
  q ≈ 1.36 AU (perihelion distance)
  v_inf ≈ 58 km/s = 58/1731.457 AU/day ≈ 0.0335 AU/day
  e ≈ 6.14 (very hyperbolic)
  At r = 5 AU inbound: v_r ≈ -sqrt(v_inf^2 + 2*GM/r) ≈ -sqrt(0.0335^2 + 2*2.96e-4/5)
    = -sqrt(0.001122 + 0.0001184) ≈ -sqrt(0.001241) ≈ -0.0352 AU/day

For Oumuamua (v_inf ≈ 26 km/s ≈ 0.015 AU/day) at r=5 AU:
  v_r ≈ -sqrt(0.015^2 + 0.0001184) ≈ -sqrt(0.000225 + 0.000118) ≈ -0.0185 AU/day

Strategy: grid covers r = 0.5 to 15 AU, rdot from -(v_esc + Δv) to -(v_esc + v_inf_max)
where v_inf_max ~ 200 km/s (galactic velocity dispersion upper bound) = 0.1155 AU/day.
"""

import math
from pathlib import Path

GM_SUN  = 2.9591220828e-4   # AU^3 / day^2
TEST_DATA = Path(__file__).parent.parent / "tests" / "data"

def v_esc(r):
    """Escape velocity at r (AU) in AU/day."""
    return math.sqrt(2.0 * GM_SUN / r)

def rdot_hyperbolic(r, v_inf_au_day, inbound=True):
    """Total radial speed at r for given v_inf (AU/day). Negative if inbound."""
    speed = math.sqrt(v_inf_au_day**2 + 2.0 * GM_SUN / r)
    return -speed if inbound else speed

# r grid: fine near the Sun, coarser at large distances.
# 2–6 AU uses 0.25 AU steps (was 0.5) — this is the prime detection zone for
# interstellars observed pre/post-perihelion at typical survey depths.
r_values = []
r = 0.5
while r <= 1.5:
    r_values.append(round(r, 3))
    r += 0.25
r = 2.0
while r <= 6.0:
    r_values.append(round(r, 3))
    r += 0.25
r = 7.0
while r <= 15.0:
    r_values.append(round(r, 3))
    r += 1.0

# v_inf grid (AU/day): from very slow (~5 km/s) to fast (~200 km/s).
# The 20–75 km/s range covers all three known interstellars
# ('Oumuamua 26, Borisov 32, 3I/ATLAS 58) and uses 5 km/s steps (was 10–15).
# Outside that sweet spot the original spacing is kept.
# 1 km/s = 1/1731.457 AU/day
KM_PER_S_TO_AU_PER_DAY = 1.0 / 1731.457

v_inf_kms_values = [
     5, 10, 15,                              # slow end: 5 km/s steps
    20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75,  # sweet spot: 5 km/s steps
    87, 100, 115, 130, 145, 160, 180, 200,   # fast end: ~15 km/s steps
]
v_inf_au_values  = [v * KM_PER_S_TO_AU_PER_DAY for v in v_inf_kms_values]

rows = []
for r in r_values:
    # Also include a barely-hyperbolic entry (v_inf = 1 km/s) to probe the boundary
    for v_inf in [0.001 * KM_PER_S_TO_AU_PER_DAY] + v_inf_au_values:
        for inbound in (True, False):
            rdot = rdot_hyperbolic(r, v_inf, inbound=inbound)
            # mean_accel = 1.0 means purely Keplerian centripetal (no extra perturbation)
            rows.append((r, rdot, 1.0))

# Deduplicate and sort
rows = sorted(set((round(r,4), round(rd,6), a) for r,rd,a in rows))

# Write inbound-only grid (backward-compatible)
outfile_01 = TEST_DATA / "heliohyp_interstellar01.txt"
inbound_rows = [(r, rd, a) for r, rd, a in rows if rd < 0]
with open(outfile_01, "w") as f:
    f.write("#r(AU) rdot(AU/day) mean_accel\n")
    for r, rdot, a in inbound_rows:
        f.write(f"{r:.4f} {rdot:.6f} {a:.1f}\n")
print(f"Wrote {len(inbound_rows)} inbound hypothesis rows to {outfile_01}")

# Write full inbound+outbound grid
outfile_02 = TEST_DATA / "heliohyp_interstellar02.txt"
with open(outfile_02, "w") as f:
    f.write("#r(AU) rdot(AU/day) mean_accel\n")
    for r, rdot, a in rows:
        f.write(f"{r:.4f} {rdot:.6f} {a:.1f}\n")
print(f"Wrote {len(rows)} inbound+outbound hypothesis rows to {outfile_02}")

# Summary
r_ref = 3.0  # AU at which 3I was well-detected
ve_ref = v_esc(r_ref)
rdot_3I = rdot_hyperbolic(r_ref, 58 * KM_PER_S_TO_AU_PER_DAY)
print(f"\nReference: r={r_ref} AU, v_esc={ve_ref:.4f} AU/day ({ve_ref*1731.457:.1f} km/s)")
print(f"  3I/ATLAS (v_inf≈58 km/s): rdot ≈ {rdot_3I:.4f} AU/day ({rdot_3I*1731.457:.1f} km/s)")
print(f"  Oumuamua (v_inf≈26 km/s): rdot ≈ {rdot_hyperbolic(r_ref, 26*KM_PER_S_TO_AU_PER_DAY):.4f} AU/day")
