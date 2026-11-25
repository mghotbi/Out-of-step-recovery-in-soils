#!/usr/bin/env bash
# Regenerate every figure and check every number quoted in the manuscript.
# Usage (from the repository root):  bash run_all.sh
set -euo pipefail
cd "$(dirname "$0")/python"
for s in fig1_concept.py fig5_evidence.py fig6_test.py figS2_predictions.py figS3_capacity_inventory.py figS4_reanalysis_detail.py; do
  echo "== $s"; python "$s"
done
if python -c "import cv2" >/dev/null 2>&1; then
  python figS1a_relabel_artwork.py
else
  test -s ../figures/Supplementary_Figure_1a_artwork.png
  echo "OpenCV unavailable: using the bundled final Supplementary Figure 1a artwork (relabeling was not rerun)."
fi
python figS1_worked_example.py b
python figS1_compose.py
python la_integral_audit.py
python verify_reported_values.py   # must print ALL PASS
python verify_yu_source.py         # current Figure 5d and Supplementary Figure 3
if command -v Rscript >/dev/null 2>&1; then
  cd ../R
  for s in verify_reported_values.R verify_theory.R independent_math_check.R verify_weak_escape.R; do
    echo "== $s"; Rscript "$s"
  done
else
  echo "Rscript not found: skipping the independent R checks."
fi
echo "Done."
