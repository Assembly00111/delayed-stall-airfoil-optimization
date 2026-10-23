#!/bin/bash
# =============================================================================
# run_domain_sweep.sh
#
# Domain-independence sweep for the airfoil2d_transient OpenFOAM case.
# Keeps snappyHexMeshDict (surface/region refinement) fixed and only changes
# the far-field box size in blockMeshDict, scaling the background block
# divisions so the far-field cell density stays roughly constant across runs.
#
# USAGE
#   1. Edit BASE_CASE below to point at a CLEAN copy of the case (0/, constant/
#      minus polyMesh, system/ -- i.e. before any blockMesh/snappyHexMesh has
#      been run on it).
#   2. Edit the DOMAINS array if you want different sizes.
#   3. Put this script, blockMeshDict.template and snappyHexMeshDict.template
#      in the same directory.
#   4. Run:  bash run_domain_sweep.sh
#
# OUTPUT
#   ./sweep_domain/<tag>/   -- one full case per domain size
#   ./sweep_domain/<tag>/log.*         -- logs of each meshing/solving step
#   ./sweep_domain/<tag>/postProcessing/forces/0/forceCoeffs.dat
#
# Then post-process with:
#   python3 post_process.py sweep_domain/U1_x-3 sweep_domain/U2_x-5 \
#           sweep_domain/U3_x-7 --out domain_sweep
# =============================================================================
set -euo pipefail

# ---- USER SETTINGS ---------------------------------------------------------

BASE_CASE="./airfoil2d_transient"      # <-- path to a clean template case
WORKDIR="./sweep_domain"
BLOCK_TEMPLATE="./blockMeshDict.template"
SNAPPY_TEMPLATE="./snappyHexMeshDict.template"

# Refinement is held FIXED at these (baseline / "medium") levels for the
# whole domain sweep -- only the far-field box size changes.
WLMIN=4; WLMAX=5; B1=3; B2=3; B3=4

# Background cell density of the ORIGINAL case (cells per unit length),
# used to scale NX, NY with domain size so far-field resolution is held
# constant. Original: XMIN=-5 XMAX=27 NX=25 ; YMIN=-8 YMAX=8 NY=12
DENSITY_X="0.78125"   # 25 / (27 - (-5))
DENSITY_Y="0.75"      # 12 / (8  - (-8))

# Upstream-distance sweep: XMAX/YMIN/YMAX held fixed, only XMIN (i.e. how far
# upstream the inlet sits) is varied.
#
# tag          XMIN  XMAX  YMIN  YMAX
DOMAINS=(
  "U1_x-3       -3   27    -8    8"
  "U2_x-5       -5   27    -8    8"
  "U3_x-7       -7   27    -8    8"
)

# ---- END USER SETTINGS -----------------------------------------------------

mkdir -p "$WORKDIR"

reset_case () {
  # Strip a case back to 0/, constant/{geometry,physicalProperties,momentumTransport},
  # system/ only -- removes any mesh/results left over from a previous run.
  local dir="$1"
  rm -rf "$dir/constant/polyMesh"
  rm -rf "$dir/postProcessing"
  rm -f  "$dir"/log.*
  find "$dir" -maxdepth 1 -regextype posix-extended -regex '.*/[0-9]+(\.[0-9]+)?' ! -name '0' -exec rm -rf {} + 2>/dev/null || true
}

for entry in "${DOMAINS[@]}"; do
  read -r TAG XMIN XMAX YMIN YMAX <<< "$entry"
  CASE_DIR="$WORKDIR/$TAG"

  echo "=================================================================="
  echo " Domain case: $TAG   x=[$XMIN,$XMAX]  y=[$YMIN,$YMAX]"
  echo "=================================================================="

  rm -rf "$CASE_DIR"
  cp -r "$BASE_CASE" "$CASE_DIR"
  reset_case "$CASE_DIR"

  NX=$(awk -v a="$XMIN" -v b="$XMAX" -v d="$DENSITY_X" 'BEGIN{printf "%d", (b-a)*d + 0.5}')
  NY=$(awk -v a="$YMIN" -v b="$YMAX" -v d="$DENSITY_Y" 'BEGIN{printf "%d", (b-a)*d + 0.5}')
  echo "  -> background divisions NX=$NX NY=$NY"

  sed -e "s/@XMIN@/${XMIN}/g" -e "s/@XMAX@/${XMAX}/g" -e "s/@YMIN@/${YMIN}/g" -e "s/@YMAX@/${YMAX}/g" -e "s/@NX@/${NX}/g" -e "s/@NY@/${NY}/g" "$BLOCK_TEMPLATE" > "$CASE_DIR/system/blockMeshDict"

  # locationInMesh must sit inside THIS domain, safely upstream of the wing --
  # 1 unit in from the inlet face is enough margin since the wing sits near x=0..1.
  LOCX=$(awk -v a="$XMIN" 'BEGIN{printf "%.3f", a + 1}')
  LOCY=0
  # refinementBox1 starts exactly at this case's own inlet (no buffer) so
  # every domain in the sweep gets the SAME mesh treatment at the inlet --
  # no variable-length stretch of coarse background cells that changes
  # with domain size.
  BOX1_XMIN=$(awk -v a="$XMIN" 'BEGIN{printf "%.3f", a}')
  sed -e "s/@WING_LMIN@/${WLMIN}/g" -e "s/@WING_LMAX@/${WLMAX}/g" -e "s/@BOX1_LEVEL@/${B1}/g" -e "s/@BOX2_LEVEL@/${B2}/g" -e "s/@BOX3_LEVEL@/${B3}/g" -e "s/@LOCX@/${LOCX}/g" -e "s/@LOCY@/${LOCY}/g" -e "s/@BOX1_XMIN@/${BOX1_XMIN}/g" "$SNAPPY_TEMPLATE" > "$CASE_DIR/system/snappyHexMeshDict"

  (
    cd "$CASE_DIR"
    blockMesh                    | tee log.blockMesh
    snappyHexMesh -overwrite     | tee log.snappyHexMesh
    extrudeMesh                  | tee log.extrudeMesh
    createPatch -overwrite       | tee log.createPatch
    checkMesh -latestTime        | tee log.checkMesh
    pimpleFoam                   | tee log.pimpleFoam
  )

  echo "  -> done: $CASE_DIR"
done

echo "All domain-sweep cases finished. Run post_process.py on $WORKDIR/* next."
