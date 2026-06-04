#!/usr/bin/env python3
"""
Run the full heliolinx pipeline on unlinked Rubin/LSST detections.

Reads a detection CSV (from lasair_kafka_subscriber.py, lasair_bulk_fetch.py,
or fink_extract_detections.py) and runs makeTracklets → heliolinc → linkPurify,
writing cluster and refined-candidate CSVs to data/heliolinx/.

Usage:
    python data_prep/run_heliolinx.py \\
        --input data/fink_detections_March1-April31.csv \\
        --output-prefix fink_apr_may \\
        --mjd-min 61122 --mjd-max 61213

    python data_prep/run_heliolinx.py \\
        --input data/lasair_detections.csv \\
        --output-prefix lasair

The output prefix is used to name the output files:
    data/heliolinx/<prefix>_clusters.csv
    data/heliolinx/<prefix>_refined.csv
"""

import argparse
import sys
import numpy as np
import pandas as pd
from pathlib import Path

import heliolinx
from heliolinx import solarsyst_dyn_geo as sdg

ROOT      = Path(__file__).parent.parent
TEST_DATA = ROOT / "test_data"
DATA      = ROOT / "data"

IMAGETIMETOL = 10.0 / 86400.0   # 10 seconds in days


# ---------------------------------------------------------------------------
# Stationary source pre-filter
# ---------------------------------------------------------------------------

def _flux_consistent(flux_i, flux_j, band_i, band_j):
    """
    Return True if two detections have flux properties consistent with being
    the same stationary source rather than two unrelated objects.

    Checks:
      - Same flux sign (both positive or both negative)
      - If same band: flux ratio within 10x
      - If different bands: flux ratio within 100x (stellar colours vary)
    """
    if flux_i is not None and flux_j is not None:
        if (flux_i > 0) != (flux_j > 0):
            return False
        if flux_i > 0 and flux_j > 0:
            ratio = max(flux_i, flux_j) / min(flux_i, flux_j)
            if ratio > (10.0 if band_i == band_j else 100.0):
                return False
    return True


def filter_stationary(df, threshold_arcsec=5.0, min_other_nights=2):
    """
    Remove detections that appear at the same sky position on multiple different
    nights with flux properties consistent with a stationary source.

    Conservative defaults keep all objects moving faster than ~22 arcsec/day
    (the slowest object in our hypothesis grid: r=15 AU, v_inf=5 km/s).
    Requiring flux-consistent matches on min_other_nights=2 separate nights
    prevents two unrelated moving objects from being confused with one stationary.
    """
    if threshold_arcsec <= 0:
        return df, 0

    ra  = np.radians(df['ra_deg'].values)
    dec = np.radians(df['dec_deg'].values)
    xyz = np.column_stack([np.cos(dec)*np.cos(ra),
                           np.cos(dec)*np.sin(ra),
                           np.sin(dec)])
    nights = (df['mjd'].values - 0.5).astype(int)
    fluxes = df['psfFlux'].values if 'psfFlux' in df.columns else np.full(len(df), np.nan)
    bands  = df['band'].values    if 'band'   in df.columns else np.full(len(df), '')

    chord = 2.0 * np.sin(np.radians(threshold_arcsec / 3600.0) / 2.0)

    try:
        from scipy.spatial import KDTree
        neighbor_lists = KDTree(xyz).query_ball_tree(KDTree(xyz), chord)
    except ImportError:
        dists = np.sqrt(((xyz[:, None] - xyz[None, :]) ** 2).sum(axis=2))
        neighbor_lists = [list(np.where(dists[i] <= chord)[0]) for i in range(len(df))]

    stationary = np.zeros(len(df), dtype=bool)
    for i, neighbors in enumerate(neighbor_lists):
        confirming_nights = set()
        for j in neighbors:
            if j == i or nights[j] == nights[i]:
                continue
            fi = None if np.isnan(fluxes[i]) else float(fluxes[i])
            fj = None if np.isnan(fluxes[j]) else float(fluxes[j])
            if _flux_consistent(fi, fj, bands[i], bands[j]):
                confirming_nights.add(nights[j])
        if len(confirming_nights) >= min_other_nights:
            stationary[i] = True

    return df[~stationary].reset_index(drop=True), int(stationary.sum())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input',         required=True,
                        help='Detection CSV (output of lasair/fink extraction scripts)')
    parser.add_argument('--output-prefix', required=True,
                        help='Prefix for output filenames (e.g. "fink_apr_may")')
    parser.add_argument('--mjd-min',  type=float, default=None,
                        help='Filter detections to MJD >= this value')
    parser.add_argument('--mjd-max',  type=float, default=None,
                        help='Filter detections to MJD <= this value')
    parser.add_argument('--mjd-ref',  type=float, default=None,
                        help='Reference MJD for heliolinc (default: midpoint of data)')
    parser.add_argument('--hyp-grid', default=str(TEST_DATA / 'heliohyp_interstellar02.txt'),
                        help='Hypothesis grid file (default: heliohyp_interstellar02.txt)')
    parser.add_argument('--stationary-threshold', type=float, default=5.0,
                        help='Arcsec radius for stationary-source filter (0 to disable)')
    parser.add_argument('--stationary-nights', type=int, default=2,
                        help='Min cross-night matches required to flag as stationary')
    args = parser.parse_args()

    out_dir = DATA / 'heliolinx'
    out_dir.mkdir(exist_ok=True)

    # --- Step 1: Earth ephemeris + ObsCodes ---
    print("=== Step 1: Load Earth ephemeris and ObsCodes ===")
    earth_raw = sdg.load_earth_ephemerides(TEST_DATA / "Earth1day2025.csv")
    earthpos  = heliolinx.create_EarthState(len(earth_raw))
    earthpos['MJD'] = earth_raw['MJD']
    for src, dst in [('X','x'),('Y','y'),('Z','z'),('VX','vx'),('VY','vy'),('VZ','vz')]:
        earthpos[dst] = earth_raw[src]
    print(f"  {len(earthpos)} ephemeris epochs, "
          f"MJD {earthpos['MJD'].min():.1f}–{earthpos['MJD'].max():.1f}")

    obs_lookup = {}
    with open(TEST_DATA / "ObsCodes_clean.txt") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    obs_lookup[parts[0]] = (float(parts[1]), float(parts[2]), float(parts[3]))
                except ValueError:
                    pass
    print(f"  {len(obs_lookup)} observatory codes")

    # --- Step 2: Load detections ---
    print("\n=== Step 2: Load detections ===")
    df = pd.read_csv(args.input)
    if args.mjd_min is not None:
        df = df[df['mjd'] >= args.mjd_min]
    if args.mjd_max is not None:
        df = df[df['mjd'] <= args.mjd_max]
    df = df.copy()

    df.sort_values(['mjd', 'obscode', 'ra_deg'], inplace=True, ignore_index=True)

    known_mask = df['obscode'].isin(obs_lookup)
    if (~known_mask).sum() > 0:
        print(f"  Dropping {(~known_mask).sum()} detections with unknown obscodes: "
              f"{df.loc[~known_mask,'obscode'].unique()}")
        df = df[known_mask].reset_index(drop=True)

    print(f"  {len(df)} detections from {df['obscode'].nunique()} observatories, "
          f"MJD {df['mjd'].min():.1f}–{df['mjd'].max():.1f}")

    # --- Step 2b: Filter stationary sources ---
    print("\n=== Step 2b: Filter stationary sources ===")
    df, n_stat = filter_stationary(df, args.stationary_threshold, args.stationary_nights)
    print(f"  Removed {n_stat} stationary detections → {len(df)} remaining")

    mjd_ref = args.mjd_ref if args.mjd_ref else (df['mjd'].min() + df['mjd'].max()) / 2
    print(f"  Using MJDREF = {mjd_ref:.2f}")

    detvec = heliolinx.create_hldet(len(df))
    detvec['MJD']      = df['mjd'].values
    detvec['RA']       = df['ra_deg'].values
    detvec['Dec']      = df['dec_deg'].values
    detvec['mag']      = df['mag'].values
    detvec['band']     = [s.encode('ascii')[:5] for s in df['band']]
    detvec['obscode']  = [s.encode('ascii')[:5] for s in df['obscode']]
    detvec['idstring'] = [str(x).encode('ascii')[:20] for x in df['obsid']]

    # --- Step 3: Build image log ---
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
        image_rows.append({'MJD': group['mjd'].mean(), 'RA': group['ra_deg'].mean(),
                           'Dec': group['dec_deg'].mean(), 'obscode': row_i['obscode']})
        i = j
    print(f"  {len(image_rows)} images")

    imglog = heliolinx.create_hlimage(len(image_rows))
    for k, img in enumerate(image_rows):
        code = img['obscode']
        Long, pcos, psin = obs_lookup[code]
        obsvec = heliolinx.observer_vel(img['MJD'], Long, pcos, psin, earthpos)
        imglog[k]['MJD']     = img['MJD']
        imglog[k]['RA']      = img['RA']
        imglog[k]['Dec']     = img['Dec']
        imglog[k]['obscode'] = code.encode('ascii')[:5]
        imglog[k]['X']  = obsvec[0]; imglog[k]['Y']  = obsvec[1]; imglog[k]['Z']  = obsvec[2]
        imglog[k]['VX'] = obsvec[3]; imglog[k]['VY'] = obsvec[4]; imglog[k]['VZ'] = obsvec[5]
        imglog[k]['exptime'] = -1.0

    # --- Step 4: Hypothesis grid ---
    print("\n=== Step 4: Load hypothesis grid ===")
    hyp_raw = np.loadtxt(args.hyp_grid, comments='#')
    radhyp  = heliolinx.create_hlradhyp(len(hyp_raw))
    radhyp['HelioRad'] = hyp_raw[:, 0]
    radhyp['R_dot']    = hyp_raw[:, 1]
    radhyp['R_dubdot'] = hyp_raw[:, 2]
    print(f"  {len(radhyp)} hypotheses, "
          f"r {radhyp['HelioRad'].min():.1f}–{radhyp['HelioRad'].max():.1f} AU, "
          f"rdot {radhyp['R_dot'].min():.4f}–{radhyp['R_dot'].max():.4f} AU/day")

    # --- Step 5: makeTracklets ---
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
        print("No tracklets formed — data too sparse. Exiting.")
        sys.exit(0)

    # --- Step 6: heliolinc ---
    print("\n=== Step 6: heliolinc ===")
    hl_config = heliolinx.HeliolincConfig()
    hl_config.MJDref       = mjd_ref
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

    # --- Step 7: linkPurify ---
    print("\n=== Step 7: linkPurify ===")
    lp_config = heliolinx.LinkPurifyConfig()
    lp_config.maxrms       = 50_000_000.0
    lp_config.minobsnights = 2
    lp_config.minpointnum  = 4

    refined, refined2det = heliolinx.linkPurify(
        lp_config, imglog, pairdets, clusters, clust2det)
    print(f"  {len(refined)} refined candidates")

    # --- Step 8: Save ---
    print(f"\n=== Step 8: Save outputs to {out_dir} ===")
    prefix = args.output_prefix

    pd.DataFrame(clusters).to_csv(out_dir / f"{prefix}_clusters.csv", index=False)
    print(f"  Wrote {prefix}_clusters.csv ({len(clusters)} rows)")

    refined_df = pd.DataFrame(refined)
    refined_df.to_csv(out_dir / f"{prefix}_refined.csv", index=False)
    print(f"  Wrote {prefix}_refined.csv ({len(refined_df)} rows)")

    if len(refined_df) > 0 and 'orbit_e' in refined_df.columns:
        hyp = refined_df[refined_df['orbit_e'] > 1.0]
        print(f"\n  Hyperbolic candidates (e > 1): {len(hyp)}")
        if len(hyp) > 0:
            print(refined_df[['orbit_a','orbit_e','obsnights','timespan','astromRMS']].to_string())

    print("\n=== Pipeline complete ===")


if __name__ == '__main__':
    main()
