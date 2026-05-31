#!/usr/bin/env python3
"""
Run the full heliolinx pipeline on real Rubin/LSST detections from Lasair.

Inputs (from data/ and test_data/):
  data/lasair_detections.csv          Unlinked Rubin detections (from Lasair)
  test_data/Earth1day2025.csv         JPL Horizons Earth ephemeris
  test_data/ObsCodes_clean.txt        MPC observatory codes
  test_data/heliohyp_interstellar02.txt  Inbound+outbound hypothesis grid

Outputs (to data/heliolinx/):
  lasair_clusters.csv                 Raw heliolinc clusters
  lasair_refined.csv                  linkPurify-refined candidates
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

import heliolinx
from heliolinx import solarsyst_dyn_geo as sdg


def _flux_consistent(flux_i, flux_j, band_i, band_j):
    """
    Return True if two detections have flux properties consistent with being
    the same stationary source rather than two unrelated objects.

    Checks:
      - Same flux sign (both positive or both negative)
      - If same band: flux ratio within 10x
      - If different bands: flux ratio within 100x (stellar colours vary)
    If either flux is None/zero, only the sign check is skipped.
    """
    # Sign check: stationary sources stay in the same positive/negative regime
    if flux_i is not None and flux_j is not None:
        if (flux_i > 0) != (flux_j > 0):
            return False   # one positive, one negative → likely different objects
        # Flux ratio check (only meaningful for positive fluxes)
        if flux_i > 0 and flux_j > 0:
            ratio = max(flux_i, flux_j) / min(flux_i, flux_j)
            limit = 10.0 if band_i == band_j else 100.0
            if ratio > limit:
                return False  # flux too different to be same source
    return True


def filter_stationary(df, threshold_arcsec=5.0, min_other_nights=2):
    """
    Remove detections that appear at essentially the same sky position on
    multiple different nights AND whose flux/band properties are consistent
    with being the same stationary source.

    A detection is flagged as stationary only when BOTH conditions hold:
      1. It has cross-night positional neighbors within threshold_arcsec on
         at least min_other_nights distinct nights.
      2. Each such neighbor passes _flux_consistent() — same flux sign,
         and flux ratio within 10x (same band) or 100x (cross-band).

    Conservative defaults:
      threshold_arcsec=5   — well below the ~22 arcsec/day minimum motion
                             of the slowest object in our hypothesis grid
                             (r=15 AU, v_inf=5 km/s → ~22 arcsec/day transverse)
      min_other_nights=2   — requires flux-consistent repeats on 2+ other nights
    """
    if threshold_arcsec <= 0:
        return df, 0

    ra  = np.radians(df['ra_deg'].values)
    dec = np.radians(df['dec_deg'].values)
    xyz = np.column_stack([
        np.cos(dec) * np.cos(ra),
        np.cos(dec) * np.sin(ra),
        np.sin(dec),
    ])
    nights  = ((df['mjd'].values - 0.5)).astype(int)
    fluxes  = df['psfFlux'].values if 'psfFlux' in df.columns else np.full(len(df), None)
    bands   = df['band'].values   if 'band'   in df.columns else np.full(len(df), '')

    chord = 2.0 * np.sin(np.radians(threshold_arcsec / 3600.0) / 2.0)

    try:
        from scipy.spatial import KDTree
        tree = KDTree(xyz)
        neighbor_lists = tree.query_ball_tree(tree, chord)
    except ImportError:
        dists = np.sqrt(((xyz[:, None] - xyz[None, :]) ** 2).sum(axis=2))
        neighbor_lists = [list(np.where(dists[i] <= chord)[0]) for i in range(len(df))]

    stationary = np.zeros(len(df), dtype=bool)
    for i, neighbors in enumerate(neighbor_lists):
        confirming_nights = set()
        for j in neighbors:
            if j == i or nights[j] == nights[i]:
                continue
            fi = fluxes[i] if not (isinstance(fluxes[i], float) and np.isnan(fluxes[i])) else None
            fj = fluxes[j] if not (isinstance(fluxes[j], float) and np.isnan(fluxes[j])) else None
            if _flux_consistent(fi, fj, bands[i], bands[j]):
                confirming_nights.add(nights[j])
        if len(confirming_nights) >= min_other_nights:
            stationary[i] = True

    n_removed = stationary.sum()
    return df[~stationary].reset_index(drop=True), n_removed

ROOT      = Path(__file__).parent.parent
TEST_DATA = ROOT / "test_data"
DATA      = ROOT / "data"
OUT_DIR   = DATA / "heliolinx"

MJDREF       = 61138.5
MJD_MIN      = 61090.0
MJD_MAX      = 61188.0
IMAGETIMETOL = 10.0 / 86400.0   # 10 seconds in days

# ---------------------------------------------------------------------------
print("=== Step 1: Load Earth ephemeris and ObsCodes ===")

earth_raw = sdg.load_earth_ephemerides(TEST_DATA / "Earth1day2025.csv")

earthpos = heliolinx.create_EarthState(len(earth_raw))
earthpos['MJD'] = earth_raw['MJD']
for src, dst in [('X','x'), ('Y','y'), ('Z','z'), ('VX','vx'), ('VY','vy'), ('VZ','vz')]:
    earthpos[dst] = earth_raw[src]

print(f"  Loaded {len(earthpos)} Earth ephemeris epochs")

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
print("\n=== Step 2: Load detections ===")

df = pd.read_csv(DATA / "lasair_detections.csv")
df = df[(df['mjd'] >= MJD_MIN) & (df['mjd'] <= MJD_MAX)].copy()

df.sort_values(['mjd', 'obscode', 'ra_deg'], inplace=True, ignore_index=True)

known_mask = df['obscode'].isin(obs_lookup)
n_unknown = (~known_mask).sum()
if n_unknown > 0:
    print(f"  Dropping {n_unknown} detections with unknown obscodes: {df.loc[~known_mask, 'obscode'].unique()}")
    df = df[known_mask].reset_index(drop=True)

print(f"  {len(df)} detections from {df['obscode'].nunique()} observatories")

print("\n=== Step 2b: Filter stationary sources ===")
df, n_stationary = filter_stationary(df, threshold_arcsec=5.0, min_other_nights=2)
print(f"  Removed {n_stationary} stationary detections → {len(df)} remaining")

detvec = heliolinx.create_hldet(len(df))
detvec['MJD']      = df['mjd'].values
detvec['RA']       = df['ra_deg'].values
detvec['Dec']      = df['dec_deg'].values
detvec['mag']      = df['mag'].values
detvec['band']     = [s.encode('ascii')[:5] for s in df['band']]
detvec['obscode']  = [s.encode('ascii')[:5] for s in df['obscode']]
detvec['idstring'] = [str(x).encode('ascii')[:20] for x in df['obsid']]

# ---------------------------------------------------------------------------
print("\n=== Step 3: Build image log ===")

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
print("\n=== Step 4: Load hypothesis grid ===")

hyp_raw = np.loadtxt(TEST_DATA / "heliohyp_interstellar02.txt", comments='#')
radhyp  = heliolinx.create_hlradhyp(len(hyp_raw))
radhyp['HelioRad'] = hyp_raw[:, 0]
radhyp['R_dot']    = hyp_raw[:, 1]
radhyp['R_dubdot'] = hyp_raw[:, 2]

print(f"  {len(radhyp)} hypothesis grid points")
print(f"  r range:    {radhyp['HelioRad'].min():.2f} – {radhyp['HelioRad'].max():.2f} AU")
print(f"  rdot range: {radhyp['R_dot'].min():.4f} – {radhyp['R_dot'].max():.4f} AU/day")

# ---------------------------------------------------------------------------
print("\n=== Step 5: makeTracklets ===")

mt_config = heliolinx.MakeTrackletsConfig()
mt_config.imagerad     = 5.0
mt_config.maxvel       = 15.0
mt_config.maxtime      = 2.0 / 24.0
mt_config.minarc       = 0.5
mt_config.maxgcr       = 0.5
mt_config.imagetimetol = IMAGETIMETOL
mt_config.forcerun     = 1

pairdets, tracklets, trk2det = heliolinx.makeTracklets(mt_config, detvec, imglog)
print(f"  {len(tracklets)} tracklets, {len(pairdets)} paired detections")

if len(tracklets) == 0:
    print("No tracklets formed — data too sparse for heliolinc. Exiting.")
    sys.exit(0)

# ---------------------------------------------------------------------------
print("\n=== Step 6: heliolinc ===")

hl_config = heliolinx.HeliolincConfig()
hl_config.MJDref       = MJDREF
hl_config.clustrad     = 15_000_000.0
hl_config.dbscan_npt   = 2
hl_config.minobsnights = 2
hl_config.mintimespan  = 0.5
hl_config.mingeodist   = 0.05
hl_config.maxgeodist   = 20.0
hl_config.geologstep   = 1.5
hl_config.max_v_inf    = 200.0
hl_config.use_univar   = 9

clusters, clust2det = heliolinx.heliolinc(
    hl_config, imglog, pairdets, tracklets, trk2det, radhyp, earthpos)
print(f"  {len(clusters)} raw clusters")

# ---------------------------------------------------------------------------
print("\n=== Step 7: linkPurify ===")

lp_config = heliolinx.LinkPurifyConfig()
lp_config.maxrms       = 50_000_000.0
lp_config.minobsnights = 2
lp_config.minpointnum  = 4

refined, refined2det = heliolinx.linkPurify(
    lp_config, imglog, pairdets, clusters, clust2det)
print(f"  {len(refined)} refined candidates")

# ---------------------------------------------------------------------------
print(f"\n=== Step 8: Save outputs to {OUT_DIR} ===")

OUT_DIR.mkdir(exist_ok=True)

pd.DataFrame(clusters).to_csv(OUT_DIR / "lasair_clusters.csv", index=False)
print(f"  Wrote lasair_clusters.csv ({len(clusters)} rows)")

refined_df = pd.DataFrame(refined)
refined_df.to_csv(OUT_DIR / "lasair_refined.csv", index=False)
print(f"  Wrote lasair_refined.csv ({len(refined_df)} rows)")

if len(refined_df) > 0 and 'orbit_e' in refined_df.columns:
    hyp = refined_df[refined_df['orbit_e'] > 1.0]
    print(f"\n  Hyperbolic candidates (e > 1): {len(hyp)}")
    if len(hyp) > 0:
        print(refined_df[['orbit_a', 'orbit_e', 'orbit_incl', 'obsnights', 'timespan']].to_string())

print("\n=== Pipeline complete ===")
