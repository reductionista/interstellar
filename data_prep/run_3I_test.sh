#!/usr/bin/env bash
# Run heliolinc_interstellar pipeline on 3I/ATLAS MPC astrometry.
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

# Reference MJD: near center of obs window (Jun–Oct 2025; perihelion Oct 29 = MJD ~60947)
MJDREF=60920

# === Step 1: make_tracklets ===
# - maxvel 15 deg/day: 3I was moving ~1-10 deg/day depending on phase
# - maxtime 2.0 hr: allow 2-hour window for intra-night pairs
# - imrad 5.0 deg: image radius (group nearby simultaneous detections into one image)
# - maxGCR 0.5 arcsec: max great circle residual for multi-point tracklets
# - minarc 0.5 arcsec: min arc length for valid tracklet

echo ""
echo "=== Step 1: make_tracklets ==="
time make_tracklets \
  -dets       3I_detections.csv \
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
  -clustrad     0.2 \
  -npt          2 \
  -minobsnights 2 \
  -mintimespan  0.5 \
  -mingeodist   0.05 \
  -maxgeodist   20.0 \
  -geologstep   1.5 \
  -out          3I_clusters.csv \
  -outrms       3I_clusters_rms.csv

echo "heliolinc_interstellar done."
echo "Cluster candidates: $(wc -l < 3I_clusters.csv)"

# === Step 3: link_refine_Herget ===
echo ""
echo "=== Step 3: link_refine_Herget ==="
printf "3I_clusters.csv 3I_clusters_rms.csv\n" > 3I_clusterlist.txt

time link_refine_Herget \
  -pairdet  3I_pairdets.csv \
  -lflist   3I_clusterlist.txt \
  -mjd      $MJDREF \
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
