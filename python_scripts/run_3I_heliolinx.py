#!/usr/bin/env python3
"""
heliolinx Python-API port of run_3I_test.sh.

Runs the full makeTracklets → heliolinc → linkPurify pipeline on the
3I/ATLAS inbound-window detections using the heliolinx Python bindings,
with no CLI subprocess calls and no intermediate files.

Inputs  (from test_data/):
  3I_detections_inbound.csv   MPC astrometry filtered to MJD 60810–60870
  Earth1day2025.csv           JPL Horizons Earth ephemeris
  ObsCodes_clean.txt          MPC observatory codes
  heliohyp_interstellar01.txt Hypothesis grid (r, rdot, mean_accel)

Outputs (to test_data/heliolinx/):
  3I_clusters.csv             Raw heliolinc clusters
  3I_refined.csv              linkPurify-refined candidates
"""

import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

import heliolinx
from heliolinx import solarsyst_dyn_geo as sdg

# ---------------------------------------------------------------------------
# CLI args
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="3I/ATLAS heliolinx test pipeline")
parser.add_argument("--clustrad", type=float, default=15_000_000.0,
                    help="Clustering radius in km (default: 15000000)")
parser.add_argument("--maxrms", type=float, default=50_000_000.0,
                    help="linkPurify maxrms in km (default: 50000000, effectively disabled)")
parser.add_argument("--output-dir", type=str, default=None,
                    help="Output directory (default: test_data/heliolinx/)")
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------
TEST_DATA = Path(__file__).parent.parent / "test_data"
OUT_DIR   = Path(args.output_dir) if args.output_dir else TEST_DATA / "heliolinx"

MJDREF         = 60840.0
MJD_MIN        = 60810.0
MJD_MAX        = 60870.0
IMAGETIMETOL   = 10.0 / 86400.0   # 10 seconds in days

print(f"=== 3I/ATLAS clustrad sweep: clustrad={args.clustrad/1e3:.0f}K km  maxrms={args.maxrms/1e3:.0f}K km  out={OUT_DIR} ===")

# ---------------------------------------------------------------------------
# Step 1: Earth ephemeris + ObsCodes
# ---------------------------------------------------------------------------
print("=== Step 1: Load Earth ephemeris and ObsCodes ===")

earth_raw = sdg.load_earth_ephemerides(TEST_DATA / "Earth1day2025.csv")

# observer_vel and heliolinc both require EarthState dtype (lowercase x,y,z…),
# but load_earth_ephemerides returns uppercase X,Y,Z,VX,VY,VZ.
earthpos = heliolinx.create_EarthState(len(earth_raw))
earthpos['MJD'] = earth_raw['MJD']
for src, dst in [('X','x'), ('Y','y'), ('Z','z'), ('VX','vx'), ('VY','vy'), ('VZ','vz')]:
    earthpos[dst] = earth_raw[src]

print(f"  Loaded {len(earthpos)} Earth ephemeris epochs")
print(f"  MJD range: {earthpos['MJD'].min():.1f} – {earthpos['MJD'].max():.1f}")

# Build obscode lookup: code → (Longitude, pcos, psin).
# sdg.parse_ObsCodes uses fixed column offsets that only work for the compact
# heliolinx bundled format; our ObsCodes_clean.txt is space-separated with
# longer longitude values, so we parse it directly with split().
obs_lookup = {}
with open(TEST_DATA / "ObsCodes_clean.txt") as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 4:
            try:
                obs_lookup[parts[0]] = (float(parts[1]), float(parts[2]), float(parts[3]))
            except ValueError:
                pass
print(f"  Loaded {len(obs_lookup)} observatory codes")

# ---------------------------------------------------------------------------
# Step 2: Load detections → hldet array
# ---------------------------------------------------------------------------
print("\n=== Step 2: Load detections ===")

df = pd.read_csv(TEST_DATA / "3I_detections_inbound.csv")
df = df[(df['mjd'] >= MJD_MIN) & (df['mjd'] <= MJD_MAX)].copy()

# Deduplicate by (mjd, obscode) — duplicate-MJD pairs cause istimedup rejection
# in linkPurify.  Keep the first occurrence (arbitrary but deterministic).
n_before = len(df)
df = df.drop_duplicates(subset=['mjd', 'obscode']).reset_index(drop=True)
print(f"  Dropped {n_before - len(df)} duplicate (mjd, obscode) rows")

df.sort_values(['mjd', 'obscode', 'ra_deg'], inplace=True, ignore_index=True)

# Filter out any obscodes not in our lookup (space telescopes, etc.)
known_mask = df['obscode'].isin(obs_lookup)
n_unknown = (~known_mask).sum()
if n_unknown > 0:
    unknown = df.loc[~known_mask, 'obscode'].unique()
    print(f"  Dropping {n_unknown} detections with unknown obscodes: {unknown}")
    df = df[known_mask].reset_index(drop=True)

print(f"  {len(df)} detections from {df['obscode'].nunique()} observatories")

detvec = heliolinx.create_hldet(len(df))
detvec['MJD']      = df['mjd'].values
detvec['RA']       = df['ra_deg'].values
detvec['Dec']      = df['dec_deg'].values
detvec['mag']      = df['mag'].values
detvec['band']     = [s.encode('ascii')[:5] for s in df['band']]
detvec['obscode']  = [s.encode('ascii')[:5] for s in df['obscode']]
detvec['idstring'] = [str(x).encode('ascii')[:20] for x in df['obsid']]

# ---------------------------------------------------------------------------
# Step 3: Build image log with observer positions
# ---------------------------------------------------------------------------
print("\n=== Step 3: Build image log ===")

# Group consecutive detections at the same (MJD ± IMAGETIMETOL, obscode) into one image.
# sorted_df is already sorted by (mjd, obscode).
image_rows = []
sorted_df = df.sort_values(['mjd', 'obscode']).reset_index(drop=True)
i = 0
while i < len(sorted_df):
    row_i = sorted_df.iloc[i]
    j = i + 1
    while j < len(sorted_df) and \
          abs(sorted_df.iloc[j]['mjd'] - row_i['mjd']) <= IMAGETIMETOL and \
          sorted_df.iloc[j]['obscode'] == row_i['obscode']:
        j += 1
    group = sorted_df.iloc[i:j]
    image_rows.append({
        'MJD':     group['mjd'].mean(),
        'RA':      group['ra_deg'].mean(),
        'Dec':     group['dec_deg'].mean(),
        'obscode': row_i['obscode'],
    })
    i = j

print(f"  {len(image_rows)} images from {len(df)} detections")

imglog = heliolinx.create_hlimage(len(image_rows))
for k, img in enumerate(image_rows):
    code = img['obscode']
    Long, pcos, psin = obs_lookup[code]
    # observer_vel returns [x, y, z, vx, vy, vz] in km, km/s (heliocentric)
    obsvec = heliolinx.observer_vel(img['MJD'], Long, pcos, psin, earthpos)
    imglog[k]['MJD']     = img['MJD']
    imglog[k]['RA']      = img['RA']
    imglog[k]['Dec']     = img['Dec']
    imglog[k]['obscode'] = code.encode('ascii')[:5]
    imglog[k]['X']       = obsvec[0]
    imglog[k]['Y']       = obsvec[1]
    imglog[k]['Z']       = obsvec[2]
    imglog[k]['VX']      = obsvec[3]
    imglog[k]['VY']      = obsvec[4]
    imglog[k]['VZ']      = obsvec[5]
    imglog[k]['exptime'] = -1.0

# ---------------------------------------------------------------------------
# Step 4: Load hypothesis grid → hlradhyp
# ---------------------------------------------------------------------------
print("\n=== Step 4: Load hypothesis grid ===")

hyp_raw = np.loadtxt(TEST_DATA / "heliohyp_interstellar01.txt", comments='#')
radhyp  = heliolinx.create_hlradhyp(len(hyp_raw))
radhyp['HelioRad'] = hyp_raw[:, 0]   # AU
radhyp['R_dot']    = hyp_raw[:, 1]   # AU/day (negative = inbound)
radhyp['R_dubdot'] = hyp_raw[:, 2]   # dimensionless mean_accel (multiplier of −GM/r²)

print(f"  {len(radhyp)} hypothesis grid points")
print(f"  r range:    {radhyp['HelioRad'].min():.2f} – {radhyp['HelioRad'].max():.2f} AU")
print(f"  rdot range: {radhyp['R_dot'].min():.4f} – {radhyp['R_dot'].max():.4f} AU/day")

# ---------------------------------------------------------------------------
# Step 5: makeTracklets
# ---------------------------------------------------------------------------
print("\n=== Step 5: makeTracklets ===")

mt_config = heliolinx.MakeTrackletsConfig()
mt_config.imagerad     = 5.0              # deg — image grouping radius
mt_config.maxvel       = 15.0             # deg/day
mt_config.maxtime      = 2.0 / 24.0       # days (2 hr intra-night window)
mt_config.minarc       = 0.5              # arcsec minimum arc
mt_config.maxgcr       = 0.5             # arcsec max great-circle residual
mt_config.imagetimetol = IMAGETIMETOL
mt_config.forcerun     = 1

pairdets, tracklets, trk2det = heliolinx.makeTracklets(mt_config, detvec, imglog)
print(f"  {len(tracklets)} tracklets, {len(pairdets)} paired detections")

if len(tracklets) == 0:
    print("ERROR: no tracklets formed — cannot continue", file=sys.stderr)
    sys.exit(1)

# ---------------------------------------------------------------------------
# Step 6: heliolinc
# ---------------------------------------------------------------------------
print("\n=== Step 6: heliolinc ===")

hl_config = heliolinx.HeliolincConfig()
hl_config.MJDref       = MJDREF
hl_config.clustrad     = args.clustrad
hl_config.dbscan_npt   = 2
hl_config.minobsnights = 2
hl_config.mintimespan  = 0.5            # days
hl_config.mingeodist   = 0.05           # AU
hl_config.maxgeodist   = 20.0           # AU
hl_config.geologstep   = 1.5
hl_config.max_v_inf    = 200.0          # km/s — allow interstellar speeds
# use_univar 9 = NotKepler mode (8 + 1): correct for high-v∞ unbound orbits
hl_config.use_univar   = 9

clusters, clust2det = heliolinx.heliolinc(
    hl_config, imglog, pairdets, tracklets, trk2det, radhyp, earthpos)
print(f"  {len(clusters)} raw clusters")

if len(clusters) == 0:
    print("WARNING: heliolinc found no clusters")

# ---------------------------------------------------------------------------
# Step 7: linkPurify
# ---------------------------------------------------------------------------
print("\n=== Step 7: linkPurify ===")

lp_config = heliolinx.LinkPurifyConfig()
lp_config.maxrms       = args.maxrms
lp_config.minobsnights = 2
lp_config.minpointnum  = 4

refined, refined2det = heliolinx.linkPurify(
    lp_config, imglog, pairdets, clusters, clust2det)
print(f"  {len(refined)} refined candidates")

# ---------------------------------------------------------------------------
# Step 8: Save outputs
# ---------------------------------------------------------------------------
print(f"\n=== Step 8: Save outputs to {OUT_DIR} ===")

OUT_DIR.mkdir(parents=True, exist_ok=True)

tag = f"clustrad{int(args.clustrad/1000)}K"
clusters_file = OUT_DIR / f"3I_clusters_{tag}.csv"
refined_file  = OUT_DIR / f"3I_refined_{tag}.csv"

clusters_df = pd.DataFrame(clusters)
clusters_df.to_csv(clusters_file, index=False)
print(f"  Wrote {clusters_file.name} ({len(clusters_df)} rows)")

refined_df = pd.DataFrame(refined)
refined_df.to_csv(refined_file, index=False)
print(f"  Wrote {refined_file.name} ({len(refined_df)} rows)")

# Quick sanity check: look for hyperbolic candidates
if len(refined_df) > 0 and 'orbit_e' in refined_df.columns:
    hyp = refined_df[refined_df['orbit_e'] > 1.0]
    print(f"\n  Hyperbolic candidates (e > 1): {len(hyp)}")
    if len(hyp) > 0:
        show_cols = [c for c in ['orbit_a', 'orbit_e', 'orbit_i', 'orbit_incl', 'orbit_inc',
                                  'obsnights', 'timespan'] if c in refined_df.columns]
        print(hyp[show_cols].head(5).to_string())
    print(f"\n  All columns: {list(refined_df.columns)}")

print("\n=== Pipeline complete ===")
