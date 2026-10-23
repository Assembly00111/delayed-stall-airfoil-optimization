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
#   2. Edit the DOMAINS array if you want different sizes, and AOA_LIST if
#      you want more than just AoA=0.
#   3. Put this script, blockMeshDict.template, snappyHexMeshDict.template
#      AND rotate_wing.py in the same directory.
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
ROTATE_SCRIPT="./rotate_wing.py"       # from compute_centroid / rotate_wing.py

# Refinement is held FIXED at these (baseline / "medium") levels for the
# whole domain sweep -- only the far-field box size changes.
WLMIN=4; WLMAX=5; B1=3; B2=3; B3=4

# Background cell density of the ORIGINAL case (cells per unit length),
# used to scale NX, NY with domain size so far-field resolution is held
# constant. Original: XMIN=-5 XMAX=27 NX=25 ; YMIN=-8 YMAX=8 NY=12
DENSITY_X="0.78125"   # 25 / (27 - (-5))
DENSITY_Y="0.75"      # 12 / (8  - (-8))

# Angle-of-attack values to sweep. Wing is rotated about its own area
# centroid (see rotate_wing.py --find-centroid) so its position doesn't
# shift, only its orientation. AoA=0 needs no rotation -- the original
# geometry is copied unchanged. Start with just (0) until the AoA=0 domain
# sweep itself is confirmed clean; add 10/20/30 once that's settled.
AOA_LIST=(5 10 15 20 25)
# Pivot = area centroid (center of mass) of THIS airfoil. "auto" computes it
# from the wing file via rotate_wing.py --find-centroid, so a different airfoil
# gets its own COM automatically. Or set both to numbers to override.
PIVOT_X="auto"
PIVOT_Y="auto"

# Upstream-distance sweep: XMAX/YMIN/YMAX held fixed, only XMIN (i.e. how far
# upstream the inlet sits) is varied. NOTE: at higher AoA, blockage grows, so
# you may need to widen this per-AoA once you get to nonzero angles (see
# CASE_AND_SCRIPTS_DOCUMENTATION.md Part 4) -- for now the same distances are
# used at every angle so you can directly compare "does AoA change how much
# upstream distance is needed".
#
# tag          XMIN  XMAX  YMIN  YMAX
DOMAINS=(
  "U1_x-3       -3   27    -8    8"
  "U2_x-5       -5   27    -8    8"
  "U3_x-7       -7   27    -8    8"
)

# ---- END USER SETTINGS -----------------------------------------------------

mkdir -p "$WORKDIR"

if [ "$PIVOT_X" = "auto" ] || [ "$PIVOT_Y" = "auto" ]; then
  CEN=$(python3 "$ROTATE_SCRIPT" --find-centroid "$BASE_CASE/constant/geometry/cylinder_2d.obj" | grep "Area centroid" | tr -d "()," | awk "{print \$3, \$4}")
  PIVOT_X=$(echo "$CEN" | awk "{print \$1}")
  PIVOT_Y=$(echo "$CEN" | awk "{print \$2}")
fi
echo "Rotation pivot (COM): ($PIVOT_X, $PIVOT_Y)"

# Pre-generate one rotated wing file per AoA (once, not once per domain case)
ROT_DIR="$WORKDIR/_rotated_wings"
mkdir -p "$ROT_DIR"
for AOA in "${AOA_LIST[@]}"; do
  ROT_OBJ="$ROT_DIR/cylinder_2d_aoa${AOA}.obj"
  if [ "$AOA" -eq 0 ]; then
    cp "$BASE_CASE/constant/geometry/cylinder_2d.obj" "$ROT_OBJ"
  else
    python3 "$ROTATE_SCRIPT" "$BASE_CASE/constant/geometry/cylinder_2d.obj" "$ROT_OBJ" \
        --pivot "$PIVOT_X" "$PIVOT_Y" --aoa "$AOA"
  fi
done

reset_case () {
  # Strip a case back to 0/, constant/{geometry,physicalProperties,momentumTransport},
  # system/ only -- removes any mesh/results left over from a previous run.
  local dir="$1"
  rm -rf "$dir/constant/polyMesh"
  rm -rf "$dir/postProcessing"
  rm -f  "$dir"/log.*
  find "$dir" -maxdepth 1 -regextype posix-extended -regex '.*/[0-9]+(\.[0-9]+)?' ! -name '0' -exec rm -rf {} + 2>/dev/null || true
}

# step <cmd> [args]: run an OpenFOAM command, log stdout+stderr to log.<cmd>,
# and return the COMMAND's exit status (pipefail), not tee's.
step () { local name="$1"; "$@" 2>&1 | tee "log.$name"; }

FAILED_CASES=()

for AOA in "${AOA_LIST[@]}"; do
for entry in "${DOMAINS[@]}"; do
  read -r TAG XMIN XMAX YMIN YMAX <<< "$entry"
  CASE_DIR="$WORKDIR/AoA${AOA}_${TAG}"

  echo "=================================================================="
  echo " AoA=${AOA}  Domain case: $TAG   x=[$XMIN,$XMAX]  y=[$YMIN,$YMAX]"
  echo "=================================================================="

  rm -rf "$CASE_DIR"
  cp -r "$BASE_CASE" "$CASE_DIR"
  reset_case "$CASE_DIR"

  # drop in the wing geometry for THIS angle -- filename stays cylinder_2d.obj
  # so snappyHexMeshDict.geometry.wing.file never needs to change per case
  cp "$ROT_DIR/cylinder_2d_aoa${AOA}.obj" "$CASE_DIR/constant/geometry/cylinder_2d.obj"

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

  # Run the pipeline. step() sends BOTH stdout and stderr to log.<command>
  # (OpenFOAM fatal errors go to stderr), and the && chain stops this case at
  # the first failing step instead of silently carrying on.
  if ( cd "$CASE_DIR" \
       && step blockMesh \
       && step snappyHexMesh -overwrite \
       && step extrudeMesh \
       && step createPatch -overwrite \
       && step checkMesh -latestTime \
       && step pimpleFoam ); then
    if [ -f "$CASE_DIR/postProcessing/forces/0/forceCoeffs.dat" ]; then
      echo "  -> done: $CASE_DIR"
    else
      echo "  !! pimpleFoam finished but forceCoeffs.dat is missing in $CASE_DIR"
      FAILED_CASES+=("$CASE_DIR (no forceCoeffs.dat)")
    fi
  else
    echo "  !! FAILED: $CASE_DIR -- last lines of its newest log:"
    LAST_LOG=$(ls -t "$CASE_DIR"/log.* 2>/dev/null | head -1)
    if [ -n "$LAST_LOG" ]; then echo "     ($LAST_LOG)"; tail -n 15 "$LAST_LOG"; else echo "     (no log files were written)"; fi
    FAILED_CASES+=("$CASE_DIR")
  fi
done
done

if [ ${#FAILED_CASES[@]} -gt 0 ]; then
  echo "=================================================================="
  echo "SWEEP FINISHED WITH FAILURES:"
  for c in "${FAILED_CASES[@]}"; do echo "  - $c"; done
  exit 1
fi
echo "All domain-sweep cases finished. Run post_process.py on $WORKDIR/AoA0_* first."
