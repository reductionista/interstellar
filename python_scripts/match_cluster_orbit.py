#!/usr/bin/env python3
"""
Recover the member detections of a heliolinx cluster by propagating its
Herget-fit orbit and matching against the input detection catalog.

heliolinx (run_heliolinx.py) does not persist the clust2det / refined2det
mapping, so the only way to recover which detections built a given refined
cluster — after the fact — is to take the fitted orbit (orbitX..orbitVZ at
orbit_MJD, ecliptic heliocentric km / km/s), propagate it across the data
timespan, and find detections whose RA/Dec lie on the predicted track.

Frame notes:
  * The heliolinx state vectors share the frame of the Earth ephemeris, which
    is ECLIPTIC J2000 (Earth Z-component ~1e-5 of |r|, i.e. Earth lies in the
    X-Y plane). RA/Dec in the detection catalog are EQUATORIAL (ICRS), so the
    geocentric vector must be rotated ecliptic -> equatorial before computing
    RA/Dec.
  * Positions are heliocentric; geocentric = obj_helio - earth_helio.
  * Light-time is corrected (one iteration): a 118 km/s mover shifts ~80" over
    the ~1000 s light travel time at ~2 AU, which matters at the arcsec level.

Usage:
  python python_scripts/match_cluster_orbit.py \
      --refined data/heliolinx/sweep_c1000000_refined.csv \
      --row 149 \
      --detections data/fink_detections_March1-April31.csv \
      --radius 60
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path

AU_KM   = 1.495978707e8
GM_SUN  = 1.32712440018e11          # km^3/s^2
C_KM_S  = 299792.458                # speed of light, km/s
EPS     = np.radians(23.43929111)   # obliquity of the ecliptic (J2000)
COS_E, SIN_E = np.cos(EPS), np.sin(EPS)


def ecl_to_eq(v):
    """Rotate an ecliptic-frame cartesian vector to the equatorial frame."""
    x, y, z = v
    return np.array([x, y * COS_E - z * SIN_E, y * SIN_E + z * COS_E])


def load_earth(path):
    """Parse a JPL Horizons vector table (ecliptic, heliocentric, km)."""
    lines = open(path).readlines()
    start = next(i for i, l in enumerate(lines) if '$$SOE' in l) + 1
    end   = next(i for i, l in enumerate(lines) if '$$EOE' in l)
    rows = []
    for l in lines[start:end]:
        p = [s.strip() for s in l.strip().rstrip(',').split(',')]
        if len(p) >= 8:
            rows.append([float(p[0]) - 2400000.5] + [float(x) for x in p[2:8]])
    return pd.DataFrame(rows, columns=['MJD', 'X', 'Y', 'Z', 'VX', 'VY', 'VZ'])


def earth_pos(earth, mjd):
    """Linear-interpolated heliocentric Earth position (km) at mjd."""
    i = int(np.clip(np.searchsorted(earth['MJD'].values, mjd), 1, len(earth) - 1))
    t0, t1 = earth['MJD'].iloc[i - 1], earth['MJD'].iloc[i]
    f = (mjd - t0) / (t1 - t0)
    return np.array([earth[c].iloc[i - 1] + f * (earth[c].iloc[i] - earth[c].iloc[i - 1])
                     for c in ['X', 'Y', 'Z']])


def kepler_accel(r):
    return -GM_SUN * r / np.linalg.norm(r) ** 3


def rk4(r, v, dt):
    k1r, k1v = v, kepler_accel(r)
    k2r, k2v = v + 0.5 * dt * k1v, kepler_accel(r + 0.5 * dt * k1r)
    k3r, k3v = v + 0.5 * dt * k2v, kepler_accel(r + 0.5 * dt * k2r)
    k4r, k4v = v + dt * k3v, kepler_accel(r + dt * k3r)
    return (r + dt * (k1r + 2 * k2r + 2 * k3r + k4r) / 6,
            v + dt * (k1v + 2 * k2v + 2 * k3v + k4v) / 6)


def make_propagator(r0, v0, t0, t_lo, t_hi, dt_grid=0.005):
    """
    Integrate the 2-body orbit ONCE onto a uniform time grid spanning
    [t_lo, t_hi] (days), then return f(mjd) -> heliocentric position (km) by
    linear interpolation. Avoids re-propagating from t0 for every detection.
    """
    pad = 1.0
    grid = np.arange(t_lo - pad, t_hi + pad + dt_grid, dt_grid)

    # State at the lowest grid time, reached by one continuous RK4 from t0.
    def integrate_to(target):
        dt_s = (target - t0) * 86400.0
        n = max(1, int(abs(dt_s) / 600))
        r, v = r0.copy(), v0.copy()
        step = dt_s / n
        for _ in range(n):
            r, v = rk4(r, v, step)
        return r, v

    r, v = integrate_to(grid[0])
    pos = np.empty((len(grid), 3))
    pos[0] = r
    step = dt_grid * 86400.0
    for k in range(1, len(grid)):
        r, v = rk4(r, v, step)
        pos[k] = r

    def prop(mjd):
        return np.array([np.interp(mjd, grid, pos[:, j]) for j in range(3)])
    return prop


def predict_radec(prop, earth, mjd, light_time=True):
    """Apparent (RA, Dec) in degrees, with optional one-step light-time correction."""
    obs = earth_pos(earth, mjd)
    t_emit = mjd
    for _ in range(2 if light_time else 1):
        geo_ecl = prop(t_emit) - obs
        if not light_time:
            break
        lt_days = np.linalg.norm(geo_ecl) / C_KM_S / 86400.0
        t_emit = mjd - lt_days
    geo_eq = ecl_to_eq(geo_ecl)
    ra  = np.degrees(np.arctan2(geo_eq[1], geo_eq[0])) % 360.0
    dec = np.degrees(np.arcsin(geo_eq[2] / np.linalg.norm(geo_eq)))
    return ra, dec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--refined', required=True)
    ap.add_argument('--row', type=int, required=True, help='clusternum / row index')
    ap.add_argument('--detections', required=True)
    ap.add_argument('--earth', default=str(Path(__file__).parent.parent /
                                           'tests' / 'data' / 'Earth1day2025.csv'))
    ap.add_argument('--radius', type=float, default=60.0, help='match radius, arcsec')
    args = ap.parse_args()

    ref = pd.read_csv(args.refined)
    c = ref[ref['clusternum'] == args.row].iloc[0]
    r0 = np.array([c['orbitX'], c['orbitY'], c['orbitZ']])
    v0 = np.array([c['orbitVX'], c['orbitVY'], c['orbitVZ']])
    t0 = c['orbit_MJD']

    print(f"Cluster row {args.row}: e={c['orbit_e']:.3f}  a={c['orbit_a']:.4f} AU  "
          f"astromRMS={c['astromRMS']:.3f}\"  obsnights={int(c['obsnights'])}  "
          f"timespan={c['timespan']:.2f} d")
    print(f"  Herget epoch MJD {t0:.4f}, r={np.linalg.norm(r0)/AU_KM:.2f} AU, "
          f"|v|={np.linalg.norm(v0):.1f} km/s\n")

    earth = load_earth(args.earth)
    dets = pd.read_csv(args.detections)
    print(f"Detection catalog: {len(dets)} rows, "
          f"MJD {dets['mjd'].min():.1f}-{dets['mjd'].max():.1f}")

    prop = make_propagator(r0, v0, t0, dets['mjd'].min(), dets['mjd'].max())

    # Predicted track across the data timespan
    print("\nPredicted apparent track (light-time corrected):")
    for mjd in np.linspace(dets['mjd'].min(), dets['mjd'].max(), 9):
        ra, dec = predict_radec(prop, earth, mjd)
        rh = np.linalg.norm(prop(mjd)) / AU_KM
        print(f"  MJD {mjd:8.2f}: RA={ra:7.3f}  Dec={dec:7.3f}  r_helio={rh:.2f} AU")

    # Match every detection
    rows = []
    for _, d in dets.iterrows():
        ra_p, dec_p = predict_radec(prop, earth, d['mjd'])
        dra  = (d['ra_deg'] - ra_p) * np.cos(np.radians(dec_p)) * 3600.0
        ddec = (d['dec_deg'] - dec_p) * 3600.0
        sep  = np.hypot(dra, ddec)
        if sep < args.radius:
            rows.append({'obsid': d['obsid'], 'mjd': round(d['mjd'], 4),
                         'ra': round(d['ra_deg'], 5), 'dec': round(d['dec_deg'], 5),
                         'sep_arcsec': round(sep, 2), 'band': d['band']})

    if not rows:
        print(f"\nNo detections within {args.radius}\".")
        return
    m = pd.DataFrame(rows).sort_values('mjd').reset_index(drop=True)
    print(f"\n{len(m)} detections within {args.radius}\" of the predicted track:\n")
    print(m.to_string(index=False))
    print(f"\nDistinct nights: {sorted(set(m['mjd'].astype(int)))}")


if __name__ == '__main__':
    main()
