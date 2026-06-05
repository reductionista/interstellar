#!/usr/bin/env python3
"""
Diagnostic: print inclination distributions for sample hypotheses.
Helps understand why --min-incl 15 is dropping 0 tracklets.

Usage:
    python data_prep/diag_inclinations.py
"""

import numpy as np
import pandas as pd
from pathlib import Path

import heliolinx
from heliolinx import solarsyst_dyn_geo as sdg

ROOT      = Path(__file__).parent.parent
TEST_DATA = ROOT / "test_data"
DATA      = ROOT / "data"
IMAGETIMETOL = 10.0 / 86400.0
DEG2RAD = np.pi / 180.0
_OBL = 23.4392 * DEG2RAD
ECLIPTIC_NORTH = np.array([0.0, -np.sin(_OBL), np.cos(_OBL)])
AU_KM     = 1.495978707e8
S_PER_DAY = 86400.0

# Sample hypotheses: (r AU, rdot AU/day, label)
# rdot conversions: 1 AU/day ≈ 1731 km/s; 0.01 AU/day ≈ 17 km/s
SAMPLE_HYPS = [
    (1.0,  -0.010, "r=1 AU, v=-17 km/s (near inbound)"),
    (1.0,   0.010, "r=1 AU, v=+17 km/s (near outbound)"),
    (3.0,  -0.030, "r=3 AU, v=-52 km/s"),
    (5.0,  -0.050, "r=5 AU, v=-87 km/s"),
    (5.0,  -0.116, "r=5 AU, v=-200 km/s (fast inbound)"),
    (1.0,   0.000, "r=1 AU, rdot=0 (transverse motion only)"),
]

def compute_inclinations(tracklets, imglog, hyp_r, hyp_rdot):
    i1 = tracklets['Img1']
    i2 = tracklets['Img2']
    ra_mid  = (tracklets['RA1']  + tracklets['RA2'])  * 0.5 * DEG2RAD
    dec_mid = (tracklets['Dec1'] + tracklets['Dec2']) * 0.5 * DEG2RAD
    dt = imglog['MJD'][i2] - imglog['MJD'][i1]
    dt = np.where(np.abs(dt) < 1e-10, 1e-10, dt)
    dra_dt  = (tracklets['RA2']  - tracklets['RA1'])  * DEG2RAD / dt
    ddec_dt = (tracklets['Dec2'] - tracklets['Dec1']) * DEG2RAD / dt
    obs_pos = np.column_stack([
        (imglog['X'][i1] + imglog['X'][i2]) * 0.5,
        (imglog['Y'][i1] + imglog['Y'][i2]) * 0.5,
        (imglog['Z'][i1] + imglog['Z'][i2]) * 0.5,
    ]) / AU_KM
    obs_vel = np.column_stack([
        (imglog['VX'][i1] + imglog['VX'][i2]) * 0.5,
        (imglog['VY'][i1] + imglog['VY'][i2]) * 0.5,
        (imglog['VZ'][i1] + imglog['VZ'][i2]) * 0.5,
    ]) * (S_PER_DAY / AU_KM)
    cd, sd = np.cos(dec_mid), np.sin(dec_mid)
    cr, sr = np.cos(ra_mid),  np.sin(ra_mid)
    u = np.column_stack([cd*cr, cd*sr, sd])
    du_dt = np.column_stack([
        -cd*sr*dra_dt - sd*cr*ddec_dt,
         cd*cr*dra_dt - sd*sr*ddec_dt,
         cd*ddec_dt,
    ])
    r_obj = obs_pos + hyp_r    * u
    v_obj = obs_vel + hyp_rdot * u + hyp_r * du_dt
    h = np.cross(r_obj, v_obj)
    h_mag = np.linalg.norm(h, axis=1)
    h_mag = np.where(h_mag < 1e-30, 1e-30, h_mag)
    cos_i = np.clip(h @ ECLIPTIC_NORTH / h_mag, -1.0, 1.0)
    return np.degrees(np.arccos(cos_i))

# --- Load earth ephemeris and obscodes ---
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

# --- Load a sample of the combined detections ---
print("Loading detections...")
df = pd.read_csv(DATA / "fink_detections_combined.csv")
known_mask = df['obscode'].isin(obs_lookup)
df = df[known_mask].reset_index(drop=True)
print(f"  {len(df)} detections, MJD {df['mjd'].min():.1f}–{df['mjd'].max():.1f}")
print(f"  RA range: {df['ra_deg'].min():.1f}–{df['ra_deg'].max():.1f}")
print(f"  Dec range: {df['dec_deg'].min():.1f}–{df['dec_deg'].max():.1f}")

# --- Build detvec and imglog ---
df.sort_values(['mjd', 'obscode', 'ra_deg'], inplace=True, ignore_index=True)
detvec = heliolinx.create_hldet(len(df))
detvec['MJD']      = df['mjd'].values
detvec['RA']       = df['ra_deg'].values
detvec['Dec']      = df['dec_deg'].values
detvec['mag']      = df['mag'].values
detvec['band']     = [s.encode('ascii')[:5] for s in df['band']]
detvec['obscode']  = [s.encode('ascii')[:5] for s in df['obscode']]
detvec['idstring'] = [str(x).encode('ascii')[:20] for x in df['obsid']]

sorted_df = df.sort_values(['mjd', 'obscode']).reset_index(drop=True)
image_rows = []
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

# --- makeTracklets ---
print("Building tracklets...")
mt_config = heliolinx.MakeTrackletsConfig()
mt_config.imagerad     = 5.0
mt_config.maxvel       = 15.0
mt_config.maxtime      = 2.0 / 24.0
mt_config.minarc       = 0.5
mt_config.maxgcr       = 0.5
mt_config.imagetimetol = IMAGETIMETOL
mt_config.forcerun     = 1
pairdets, tracklets, trk2det = heliolinx.makeTracklets(mt_config, detvec, imglog)
print(f"  {len(tracklets)} tracklets formed")

# --- Also show angular velocity distribution of the tracklets ---
i1 = tracklets['Img1']
i2 = tracklets['Img2']
dt = imglog['MJD'][i2] - imglog['MJD'][i1]
dra  = (tracklets['RA2']  - tracklets['RA1'])
ddec = (tracklets['Dec2'] - tracklets['Dec1'])
angvel = np.sqrt(dra**2 + ddec**2) / np.where(np.abs(dt) < 1e-10, 1e-10, dt)  # deg/day
print(f"\nTracklet angular velocity (deg/day):")
for pct in [0, 1, 5, 25, 50, 75, 95, 99, 100]:
    print(f"  {pct:3d}th pct: {np.percentile(angvel, pct):.4f}")

# --- Inclination distributions for sample hypotheses ---
print("\n" + "="*70)
print("Inclination distributions (degrees, ecliptic J2000):")
print("="*70)
for hyp_r, hyp_rdot, label in SAMPLE_HYPS:
    incl = compute_inclinations(tracklets, imglog, hyp_r, hyp_rdot)
    n_below_15 = (incl < 15).sum()
    n_below_5  = (incl < 5).sum()
    print(f"\n{label}")
    print(f"  below 15°: {n_below_15}/{len(incl)} ({100*n_below_15/len(incl):.1f}%)")
    print(f"  below  5°: {n_below_5}/{len(incl)} ({100*n_below_5/len(incl):.1f}%)")
    for pct in [0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100]:
        print(f"  {pct:3d}th pct: {np.percentile(incl, pct):.1f}°")
