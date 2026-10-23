#!/bin/bash
# =============================================================================
# run_mesh_sweep.sh
#
# Mesh-independence sweep for the airfoil2d_transient OpenFOAM case.
# Keeps the domain size FIXED (set it to whatever run_domain_sweep.sh showed
# to be domain-independent) and only changes the snappyHexMesh refinement
# levels near the wing and in the three nested refinement boxes.
#
# USAGE
#   1. Edit BASE_CASE below (clean template case, same as for the domain sweep).
#   2. Edit XMIN/XMAX/YMIN/YMAX/NX/NY below to the domain-independent size
#      you found with run_domain_sweep.sh.
#   3. Edit the MESHES array if you want different refinement combinations.
#   4. Put this script, blockMeshDict.template and snappyHexMeshDict.template
#      in the same directory.
#   5. Run:  bash run_mesh_sweep.sh
#
# OUTPUT
#   ./sweep_mesh/<tag>/   -- one full case per mesh density
#   ./sweep_mesh/<tag>/log.checkMesh    -- contains total cell count
#   ./sweep_mesh/<tag>/postProcessing/forces/0/forceCoeffs.dat
#
# Then post-process with:
#   python3 post_process.py sweep_mesh/M1_coarse sweep_mesh/M2_medium \
#           sweep_mesh/M3_fine sweep_mesh/M4_xfine --out mesh_sweep
# =============================================================================
set -euo pipefail

# ---- USER SETTINGS ---------------------------------------------------------

BASE_CASE="./airfoil2d_transient"
WORKDIR="./sweep_mesh"
BLOCK_TEMPLATE="./blockMeshDict.template"
SNAPPY_TEMPLATE="./snappyHexMeshDict.template"

# FIXED domain -- set these to the domain-independent size you found first.
# Defaults below match the original case (D2_baseline).
XMIN=-5; XMAX=27; YMIN=-8; YMAX=8; NX=25; NY=12

# tag           wing_Lmin  wing_Lmax  box1  box2  box3
MESHES=(
  "M1_coarse     3          4          2     2     3"
  "M2_medium     4          5          3     3     4"
  "M3_fine       5          6          4     4     5"
  "M4_xfine      6          7          5     5     6"
)

# ---- END USER SETTINGS -----------------------------------------------------

mkdir -p "$WORKDIR"

reset_case () {
  local dir="$1"
  rm -rf "$dir/constant/polyMesh"
  rm -rf "$dir/postProcessing"
  rm -f  "$dir"/log.*
  find "$dir" -maxdepth 1 -regextype posix-extended -regex '.*/[0-9]+(\.[0-9]+)?' ! -name '0' -exec rm -rf {} + 2>/dev/null || true
}

for entry in "${MESHES[@]}"; do
  read -r TAG WLMIN WLMAX B1 B2 B3 <<< "$entry"
  CASE_DIR="$WORKDIR/$TAG"

  echo "=================================================================="
  echo " Mesh case: $TAG   wing level=($WLMIN $WLMAX)  boxes=($B1,$B2,$B3)"
  echo "=================================================================="

  rm -rf "$CASE_DIR"
  cp -r "$BASE_CASE" "$CASE_DIR"
  reset_case "$CASE_DIR"

  # fixed domain, same for every mesh in this sweep
  sed -e "s/@XMIN@/${XMIN}/g" -e "s/@XMAX@/${XMAX}/g" -e "s/@YMIN@/${YMIN}/g" -e "s/@YMAX@/${YMAX}/g" -e "s/@NX@/${NX}/g" -e "s/@NY@/${NY}/g" "$BLOCK_TEMPLATE" > "$CASE_DIR/system/blockMeshDict"

  # swept refinement levels; locationInMesh and box1 fixed relative to
  # the (fixed) domain used throughout the mesh sweep
  LOCX=$(awk -v a="$XMIN" 'BEGIN{printf "%.3f", a + 1}')
  LOCY=0
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

echo "All mesh-sweep cases finished. Run post_process.py on $WORKDIR/* next."
