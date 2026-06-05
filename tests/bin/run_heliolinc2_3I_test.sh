#!/usr/bin/env bash
# Runs heliolinc_interstellar pipeline from legacy C++ heliolinc2 repo (not heliolinx!) on 3I/ATLAS MPC astrometry.
#
# Run from data_prep/ (or anywhere) — all paths are relative to heliolinc2_src/tests/.
# Prerequisites:
#   - prep_earth_ephem.py and prep_3I_detections.py have already been run
#   - heliolinc_interstellar and make_tracklets are in PATH
#   - make_hyp_grid.py has already been run

set -e

TESTS_DIR="$(cd "$(dirname "$0")/../heliolinc2_src/tests" && pwd)"
cd "$TESTS_DIR"

echo "=== Running from $TESTS_DIR ==="

# Reference MJD: center of inbound window (MJD 60810–60870).
# At MJD 60840, 3I is at r≈4.067 AU, near the r=4.0 grid point (gap 0.067 AU).
# rdot≈-0.0331 AU/day, closest grid hypothesis rdot≈-0.031 AU/day (gap 0.002 AU/day).
# Expected position error ~0.08 AU, well within clustrad=0.3 AU.
MJDREF=60840

# === Step 1: make_tracklets ===
# - maxvel 15 deg/day: 3I was moving ~1-10 deg/day depending on phase
# - maxtime 2.0 hr: allow 2-hour window for intra-night pairs
# - imrad 5.0 deg: image radius (group nearby simultaneous detections into one image)
# - maxGCR 0.5 arcsec: max great circle residual for multi-point tracklets
# - minarc 0.5 arcsec: min arc length for valid tracklet

echo ""
echo "=== Step 1: make_tracklets ==="
time make_tracklets \
  -dets       3I_detections_inbound.csv \
  -outimgs    3I_images.txt \
  -pairs      3I_pairs.txt \
  -pairdets   3I_pairdets.csv \
  -colformat  colformat_3I_01.txt \
  -maxvel     15.0 \
  -imrad      5.0 \
  -maxtime    2.0 \
  -maxGCR     0.5 \
  -minarc     0.5 \
  -earth      Earth1day2025.txt \
  -obscode    ObsCodes_clean.txt \
  -forcerun

echo "make_tracklets done."
echo "Pairs file lines: $(wc -l < 3I_pairs.txt)"
echo "Pairdets lines:   $(wc -l < 3I_pairdets.csv)"

# === Step 2: heliolinc_interstellar ===
# - clustrad 0.2 AU: clustering radius in heliocentric space
# - npt 2: minimum 2 points for a DBSCAN cluster (pairs)
# - minobsnights 2: require at least 2 nights
# - mintimespan 0.5: require at least 0.5-day baseline
# - mingeodist 0.05: ignore objects closer than 0.05 AU geocentric
# - maxgeodist 20.0: ignore objects farther than 20 AU
# - geologstep 1.5: log spacing for geocentric distance bins

echo ""
echo "=== Step 2: heliolinc_interstellar ==="
time heliolinc_interstellar \
  -dets         3I_pairdets.csv \
  -pairs        3I_pairs.txt \
  -mjd          $MJDREF \
  -obspos       Earth1day2025.txt \
  -heliodist    heliohyp_interstellar01.txt \
  -clustrad     15000000 \
  -npt          2 \
  -minobsnights 2 \
  -mintimespan  0.5 \
  -mingeodist   0.05 \
  -maxgeodist   20.0 \
  -geologstep   1.5 \
  -out          3I_clusters.csv \
  -outrms       3I_clusters_rms.csv

echo "heliolinc_interstellar done."
echo "Clusters found (RMS file): $(( $(wc -l < 3I_clusters_rms.csv) - 1 ))"
echo "Total detection records in cluster file: $(( $(wc -l < 3I_clusters.csv) - 1 ))"

# === Step 3: link_refine_Herget ===
echo ""
echo "=== Step 3: link_refine_Herget ==="
printf "3I_clusters.csv 3I_clusters_rms.csv\n" > 3I_clusterlist.txt

time link_refine_Herget \
  -pairdet  3I_pairdets.csv \
  -lflist   3I_clusterlist.txt \
  -mjd      $MJDREF \
  -maxrms   50000000 \
  -outfile  3I_refined.csv \
  -outrms   3I_refined_rms.csv

echo "link_refine_Herget done."
echo "Refined candidates: $(wc -l < 3I_refined.csv)"
echo ""
echo "=== Pipeline complete ==="
echo "Key outputs in $TESTS_DIR:"
echo "  3I_clusters.csv      — raw heliolinc clusters (one per row)"
echo "  3I_clusters_rms.csv  — cluster RMS/quality metrics"
echo "  3I_refined.csv       — Herget-refined orbit fits"
echo "  3I_refined_rms.csv   — refined orbit RMS"
