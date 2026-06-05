#!/usr/bin/env python3
"""
Run the full heliolinx pipeline on unlinked Rubin/LSST detections.

Reads a detection CSV (from lasair_kafka_subscriber.py, lasair_bulk_fetch.py,
or fink_extract_detections.py) and runs makeTracklets → heliolinc → linkPurify,
writing cluster and refined-candidate CSVs to data/heliolinx/.

Usage:
    python python_scripts/run_heliolinx.py \\
        --input data/fink_detections_March1-April31.csv \\
        --output-prefix fink_apr_may \\
        --mjd-min 61122 --mjd-max 61213

    python python_scripts/run_heliolinx.py \\
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

DEG2RAD = np.pi / 180.0
# Ecliptic north pole in equatorial ICRS Cartesian (J2000 obliquity 23.4392°).
# Using h · ECLIPTIC_NORTH / |h| gives cos(inclination relative to ecliptic).
_OBL = 23.4392 * DEG2RAD
ECLIPTIC_NORTH = np.array([0.0, -np.sin(_OBL), np.cos(_OBL)])

# observer_vel() returns position in km and velocity in km/s.
# The inclination formula uses hyp_r in AU and hyp_rdot in AU/day, so we
# convert obs_pos/obs_vel to the same units.
AU_KM      = 1.495978707e8   # km per AU
S_PER_DAY  = 86400.0         # seconds per day


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
# Per-hypothesis inclination filter
# ---------------------------------------------------------------------------

def compute_tracklet_inclinations(tracklets, imglog, hyp_r, hyp_rdot):
    """
    For each tracklet, compute the orbital inclination (degrees, ecliptic J2000)
    implied by the given hypothesis (r in AU, rdot in AU/day).

    Strategy: project each tracklet to a heliocentric 3D state vector using the
    hypothesis distance + the observed angular velocity, then derive inclination
    from the angular momentum vector h = r × v.

    The coordinate system is equatorial ICRS (same as heliolinx internals).
    Inclination is measured against the ecliptic north pole, not the celestial
    north pole — the two differ by the obliquity (~23.4°).
    """
    i1 = tracklets['Img1']    # image indices into imglog
    i2 = tracklets['Img2']

    # Midpoint sky position (radians)
    ra_mid  = (tracklets['RA1']  + tracklets['RA2'])  * 0.5 * DEG2RAD
    dec_mid = (tracklets['Dec1'] + tracklets['Dec2']) * 0.5 * DEG2RAD

    # Angular velocity (rad/day) from the two endpoints
    dt = imglog['MJD'][i2] - imglog['MJD'][i1]          # (N,) days
    dt = np.where(np.abs(dt) < 1e-10, 1e-10, dt)        # guard against zero
    dra_dt  = (tracklets['RA2']  - tracklets['RA1'])  * DEG2RAD / dt
    ddec_dt = (tracklets['Dec2'] - tracklets['Dec1']) * DEG2RAD / dt

    # Observer heliocentric state at tracklet midpoint.
    # imglog X/Y/Z are in km and VX/VY/VZ in km/s (from observer_vel()).
    # Convert to AU and AU/day to match hyp_r (AU) and hyp_rdot (AU/day).
    obs_pos = np.column_stack([
        (imglog['X'][i1]  + imglog['X'][i2])  * 0.5,
        (imglog['Y'][i1]  + imglog['Y'][i2])  * 0.5,
        (imglog['Z'][i1]  + imglog['Z'][i2])  * 0.5,
    ]) / AU_KM
    obs_vel = np.column_stack([
        (imglog['VX'][i1] + imglog['VX'][i2]) * 0.5,
        (imglog['VY'][i1] + imglog['VY'][i2]) * 0.5,
        (imglog['VZ'][i1] + imglog['VZ'][i2]) * 0.5,
    ]) * (S_PER_DAY / AU_KM)

    # Unit vector toward object (N, 3)
    cd, sd = np.cos(dec_mid), np.sin(dec_mid)
    cr, sr = np.cos(ra_mid),  np.sin(ra_mid)
    u = np.column_stack([cd * cr, cd * sr, sd])

    # Time-derivative of unit vector (N, 3), units 1/day (radians dimensionless)
    du_dt = np.column_stack([
        -cd * sr * dra_dt - sd * cr * ddec_dt,
         cd * cr * dra_dt - sd * sr * ddec_dt,
         cd * ddec_dt,
    ])

    # Heliocentric position and velocity (AU, AU/day)
    r_obj = obs_pos + hyp_r    * u
    v_obj = obs_vel + hyp_rdot * u + hyp_r * du_dt

    # Angular momentum h = r × v (N, 3), AU²/day
    h = np.cross(r_obj, v_obj)
    h_mag = np.linalg.norm(h, axis=1)
    h_mag = np.where(h_mag < 1e-30, 1e-30, h_mag)   # guard against degenerate radial motion

    # Inclination relative to ecliptic J2000
    cos_i = np.clip(h @ ECLIPTIC_NORTH / h_mag, -1.0, 1.0)
    return np.degrees(np.arccos(cos_i))


def apply_incl_filter(tracklets, trk2det, incl, min_incl):
    """
    Keep only tracklets with implied inclination >= min_incl degrees.
    Remaps tracklet indices in trk2det to match the filtered array.
    Returns (filtered_tracklets, filtered_trk2det).
    """
    keep = incl >= min_incl
    if keep.all():
        return tracklets, trk2det

    new_idx = np.full(len(tracklets), -1, dtype=np.int64)
    new_idx[keep] = np.arange(int(keep.sum()))

    valid = keep[trk2det['i1']]
    ft2d = trk2det[valid].copy()
    ft2d['i1'] = new_idx[trk2det['i1'][valid]]

    return tracklets[keep], ft2d


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
    parser.add_argument('--clustrad', type=float, default=15_000_000.0,
                        help='heliolinc DBSCAN clustering radius in KM (not AU). Controls how '
                             'close two state vectors must be in phase space to join the same '
                             'cluster. Effective radius scales with geocentric distance, so the '
                             'same value is more permissive at larger r. Too large → degenerate '
                             'megaclusters that absorb unrelated detections (default 15e6).')
    parser.add_argument('--maxrms', type=float, default=50_000_000.0,
                        help='linkPurify max RMS residual in KM for accepting a purified '
                             'linkage (default 50e6).')
    parser.add_argument('--hyp-batch-size', type=int, default=50,
                        help='Process the hypothesis grid in batches of this many rows, '
                             'running heliolinc+linkPurify per batch and merging the refined '
                             'candidates. Bounds peak RAM (heliolinc accumulates all linkages '
                             'across hypotheses before purifying, which OOMs on large grids). '
                             'Use 0 or a value >= grid size to process all at once.')
    parser.add_argument('--min-incl', type=float, default=0.0,
                        help='Per-hypothesis inclination pre-filter (degrees). For each '
                             'hypothesis, compute the orbital inclination implied by each '
                             'tracklet and drop tracklets with inclination < this threshold '
                             'before passing to heliolinc. Interstellar objects have '
                             'isotropically distributed inclinations (P(i<15°) ≈ 3.4%%), '
                             'while prograde solar-system objects cluster near 0°. '
                             'Suggested value: 15. Forces one hypothesis per heliolinc call. '
                             'Default 0 = disabled (existing behaviour).')
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

    # --- Steps 6+7: heliolinc + linkPurify, batched over the hypothesis grid ---
    #
    # heliolinc accumulates the candidate linkages from *every* hypothesis in RAM
    # before linkPurify runs. On large grids (644 hyps × ~3300 linkages each, with
    # thousands of detections per linkage) this reaches billions of detection
    # references and OOMs. Hypotheses are clustered independently, so we process the
    # grid in batches: heliolinc+linkPurify per batch, keep only the small purified
    # output, and free the giant cluster set before the next batch. Peak RAM is then
    # bounded by one batch instead of the whole grid.
    print("\n=== Steps 6+7: heliolinc + linkPurify (batched) ===")
    prefix = args.output_prefix

    hl_config = heliolinx.HeliolincConfig()
    hl_config.MJDref       = mjd_ref
    hl_config.clustrad     = args.clustrad
    hl_config.dbscan_npt   = 4
    hl_config.minobsnights = 2
    hl_config.mintimespan  = 0.5
    hl_config.mingeodist   = 0.05
    hl_config.maxgeodist   = 20.0
    hl_config.geologstep   = 1.5
    hl_config.max_v_inf    = 200.0
    hl_config.use_univar   = 9

    lp_config = heliolinx.LinkPurifyConfig()
    lp_config.maxrms       = args.maxrms
    lp_config.minobsnights = 2
    lp_config.minpointnum  = 4

    n_hyp = len(radhyp)
    if args.min_incl > 0:
        # Inclination filter requires per-hypothesis tracklet sets — force batch=1.
        batch_size = 1
        print(f"  Inclination pre-filter: min_incl={args.min_incl}°  "
              f"(forcing batch_size=1, one hypothesis per heliolinc call)")
    else:
        batch_size = args.hyp_batch_size if args.hyp_batch_size > 0 else n_hyp
        batch_size = min(batch_size, n_hyp)
    n_batches = (n_hyp + batch_size - 1) // batch_size
    print(f"  {n_hyp} hypotheses in {n_batches} batch(es) of up to {batch_size}")

    clusters_path = out_dir / f"{prefix}_clusters.csv"
    refined_path  = out_dir / f"{prefix}_refined.csv"
    if clusters_path.exists():
        clusters_path.unlink()          # don't append onto a stale run
    if refined_path.exists():
        refined_path.unlink()
    refined_written = 0
    total_clusters = 0

    for b in range(n_batches):
        lo, hi = b * batch_size, min((b + 1) * batch_size, n_hyp)

        # --- Optional per-hypothesis inclination filter ---
        if args.min_incl > 0:
            # Single hypothesis per iteration when filtering is active.
            hyp_r    = float(radhyp['HelioRad'][lo])
            hyp_rdot = float(radhyp['R_dot'][lo])
            incl = compute_tracklet_inclinations(tracklets, imglog, hyp_r, hyp_rdot)
            ft, ft2d = apply_incl_filter(tracklets, trk2det, incl, args.min_incl)
            n_kept = len(ft)
            n_drop = len(tracklets) - n_kept
            if n_kept < hl_config.dbscan_npt:
                # Fewer surviving tracklets than the minimum cluster size — skip.
                continue
            run_tracklets, run_trk2det = ft, ft2d
            incl_note = f" ({n_kept}/{len(tracklets)} tracklets after incl filter, dropped {n_drop})"
        else:
            run_tracklets, run_trk2det = tracklets, trk2det
            incl_note = ""

        # Copy the hypothesis slice into a fresh contiguous array.
        sub = heliolinx.create_hlradhyp(hi - lo)
        sub['HelioRad'] = radhyp['HelioRad'][lo:hi]
        sub['R_dot']    = radhyp['R_dot'][lo:hi]
        sub['R_dubdot'] = radhyp['R_dubdot'][lo:hi]

        print(f"\n  --- Batch {b + 1}/{n_batches}: hypotheses {lo}–{hi - 1}{incl_note} ---")
        clusters, clust2det = heliolinx.heliolinc(
            hl_config, imglog, pairdets, run_tracklets, run_trk2det, sub, earthpos)
        print(f"    {len(clusters)} raw clusters")
        total_clusters += len(clusters)

        # Append raw clusters to disk so we never hold the full set in RAM.
        pd.DataFrame(clusters).to_csv(
            clusters_path, mode='a', header=(b == 0), index=False)

        refined, refined2det = heliolinx.linkPurify(
            lp_config, imglog, pairdets, clusters, clust2det)
        print(f"    {len(refined)} refined candidates")
        if len(refined) > 0:
            rdf = pd.DataFrame(refined)
            rdf.to_csv(refined_path, mode='a', header=(refined_written == 0), index=False)
            refined_written += len(rdf)

        # Free the big per-batch structures before the next iteration.
        del clusters, clust2det, refined, refined2det

    print(f"\n  Total raw clusters across all batches: {total_clusters}")
    print(f"  Wrote {prefix}_clusters.csv")

    # --- Step 8: Merge refined candidates and save ---
    print(f"\n=== Step 8: Save outputs to {out_dir} ===")
    if refined_written > 0:
        refined_df = pd.read_csv(refined_path)
        # The same object can be purified under adjacent hypotheses in different
        # batches; cross-batch dedup isn't done by linkPurify, so drop exact
        # duplicate rows. (Near-duplicates of a real candidate are harmless — we
        # inspect candidates individually — but exact dupes just add noise.)
        before = len(refined_df)
        refined_df = refined_df.drop_duplicates().reset_index(drop=True)
        if len(refined_df) < before:
            print(f"  Dropped {before - len(refined_df)} exact-duplicate refined rows")
    else:
        refined_df = pd.DataFrame()

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
