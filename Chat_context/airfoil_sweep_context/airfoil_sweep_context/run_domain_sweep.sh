#!/bin/bash
# =============================================================================
# run_domain_sweep.sh   (unified sweep driver)
#
# One script for every independence study on the airfoil2d_transient case.
# Pick the study with SWEEP, the wing angles with AOA, both as environment
# variables (no file editing needed):
#
#   SWEEP=upstream   (default) vary the inlet distance      XMIN = -3 ... -15
#   SWEEP=downstream vary the outlet position               XMAX = 15 ... 45
#   SWEEP=lateral    vary the half-height (top/bottom)      YMAX = -YMIN = 4 ... 16
#   SWEEP=mesh       vary snappyHexMesh refinement levels   coarse ... xfine
#
# Each study holds everything else fixed. The fixed values come from the
# studies before it, so you must state them explicitly (the script refuses to
# guess):
#
#   upstream   : nothing needed
#   downstream : UP_FIXED   (converged XMIN from the upstream study, e.g. -10)
#   lateral    : UP_FIXED  DOWN_FIXED
#   mesh       : UP_FIXED  DOWN_FIXED  HALF_H
#
# EXAMPLES
#   AOA="0" bash run_domain_sweep.sh
#   AOA="0 15 25" SWEEP=downstream UP_FIXED=-10 bash run_domain_sweep.sh
#   AOA="0 15 25" SWEEP=lateral UP_FIXED=-10 DOWN_FIXED=27 bash run_domain_sweep.sh
#   AOA="0 15 25" SWEEP=mesh UP_FIXED=-10 DOWN_FIXED=27 HALF_H=8 bash run_domain_sweep.sh
#
# OTHER SWITCHES (environment variables)
#   DRY=1         print the planned cases and stop (check before a long run)
#   MESH_ONLY=1   build + checkMesh only, no solve (minutes; inspect in ParaView)
#   FORCE=1       redo cases that are already complete (default: skip them)
#   WORKDIR=...   output folder (default depends on SWEEP, see below)
#   AOA="0 5 10"  integer angles in degrees (default "0"); positive = nose up,
#                 wing rotated about its own area centroid (center of mass)
#
# FILES NEEDED IN THE SAME FOLDER
#   blockMeshDict.template  snappyHexMeshDict.template  rotate_wing.py
#   and a CLEAN copy of the case in ./airfoil2d_transient (see BASE_CASE)
#
# OUTPUT (per case)   $WORKDIR/AoA<angle>_<tag>/
#   log.<command>  one per step, stdout AND stderr        mesh_backup/  mesh copy
#   postProcessing/forces/0/forceCoeffs.dat               kept time folders
#
# Then:   python3 post_process.py $WORKDIR/AoA* --avg-fraction 0.33 --out <name>
# =============================================================================
set -euo pipefail

# ---- SETTINGS --------------------------------------------------------------

BASE_CASE="./airfoil2d_transient"      # clean case: no constant/polyMesh, no results
BLOCK_TEMPLATE="./blockMeshDict.template"
SNAPPY_TEMPLATE="./snappyHexMeshDict.template"
ROTATE_SCRIPT="./rotate_wing.py"

# Time-step saving (for ParaView). forceCoeffs.dat is written every step
# regardless. Fields are written every WRITE_INTERVAL time units (purgeWrite 0),
# then, after a case finishes OK, time folders OUTSIDE [KEEP_T0, KEEP_T1] are
# deleted. 0/ and the final time are always kept. PRUNE=0 keeps everything.
WRITE_INTERVAL=0.5
KEEP_T0=60
KEEP_T1=120
PRUNE=1

# Background mesh density of the original case (cells per unit length), so the
# far-field cell SIZE stays the same however big the domain gets.
# Original: x from -5 to 27 with NX=25 ; y from -8 to 8 with NY=12
DENSITY_X="0.78125"
DENSITY_Y="0.75"

# Baseline refinement (wing level min/max, then box1/box2/box3 levels).
BASE_LEVELS="4 5 3 3 4"

# Pivot for the angle-of-attack rotation: "auto" = area centroid of the wing.
PIVOT_X="auto"
PIVOT_Y="auto"

# ---- STUDY DEFINITIONS -----------------------------------------------------
# Every case is one line:  tag XMIN XMAX YMIN YMAX WLMIN WLMAX B1 B2 B3

SWEEP="${SWEEP:-upstream}"
MESH_ONLY="${MESH_ONLY:-0}"
FORCE="${FORCE:-0}"
DRY="${DRY:-0}"
read -ra AOA_LIST <<< "${AOA:-0}"

case "$SWEEP" in
  upstream)
    WORKDIR_DEFAULT="./sweep_domain"
    SPECS=(
      "U1_x-3    -3  27 -8 8 $BASE_LEVELS"
      "U2_x-5    -5  27 -8 8 $BASE_LEVELS"
      "U3_x-7    -7  27 -8 8 $BASE_LEVELS"
      "U4_x-10  -10  27 -8 8 $BASE_LEVELS"
      "U5_x-15  -15  27 -8 8 $BASE_LEVELS"
    ) ;;
  downstream)
    : "${UP_FIXED:?SWEEP=downstream needs UP_FIXED, the converged XMIN from the upstream study, e.g. UP_FIXED=-10}"
    WORKDIR_DEFAULT="./sweep_downstream"
    SPECS=(
      "D1_xd15  $UP_FIXED 15 -8 8 $BASE_LEVELS"
      "D2_xd20  $UP_FIXED 20 -8 8 $BASE_LEVELS"
      "D3_xd27  $UP_FIXED 27 -8 8 $BASE_LEVELS"
      "D4_xd35  $UP_FIXED 35 -8 8 $BASE_LEVELS"
      "D5_xd45  $UP_FIXED 45 -8 8 $BASE_LEVELS"
    ) ;;
  lateral)
    : "${UP_FIXED:?SWEEP=lateral needs UP_FIXED (converged upstream XMIN)}"
    : "${DOWN_FIXED:?SWEEP=lateral needs DOWN_FIXED (converged downstream XMAX)}"
    WORKDIR_DEFAULT="./sweep_lateral"
    SPECS=(
      "L1_h4   $UP_FIXED $DOWN_FIXED  -4  4 $BASE_LEVELS"
      "L2_h6   $UP_FIXED $DOWN_FIXED  -6  6 $BASE_LEVELS"
      "L3_h8   $UP_FIXED $DOWN_FIXED  -8  8 $BASE_LEVELS"
      "L4_h12  $UP_FIXED $DOWN_FIXED -12 12 $BASE_LEVELS"
      "L5_h16  $UP_FIXED $DOWN_FIXED -16 16 $BASE_LEVELS"
    ) ;;
  mesh)
    : "${UP_FIXED:?SWEEP=mesh needs UP_FIXED (converged upstream XMIN)}"
    : "${DOWN_FIXED:?SWEEP=mesh needs DOWN_FIXED (converged downstream XMAX)}"
    : "${HALF_H:?SWEEP=mesh needs HALF_H (converged half-height, e.g. 8)}"
    WORKDIR_DEFAULT="./sweep_mesh"
    SPECS=(
      "M1_coarse  $UP_FIXED $DOWN_FIXED -$HALF_H $HALF_H 3 4 2 2 3"
      "M2_medium  $UP_FIXED $DOWN_FIXED -$HALF_H $HALF_H 4 5 3 3 4"
      "M3_fine    $UP_FIXED $DOWN_FIXED -$HALF_H $HALF_H 5 6 4 4 5"
      "M4_xfine   $UP_FIXED $DOWN_FIXED -$HALF_H $HALF_H 6 7 5 5 6"
    ) ;;
  *)
    echo "Unknown SWEEP='$SWEEP' (use upstream, downstream, lateral or mesh)" >&2; exit 2 ;;
esac
WORKDIR="${WORKDIR:-$WORKDIR_DEFAULT}"

for A in "${AOA_LIST[@]}"; do
  [[ "$A" =~ ^-?[0-9]+$ ]] || { echo "AOA values must be integers, got '$A'" >&2; exit 2; }
done

echo "=================================================================="
echo " SWEEP=$SWEEP   AoA list: ${AOA_LIST[*]}   -> ${#SPECS[@]} cases per angle"
echo " WORKDIR=$WORKDIR   MESH_ONLY=$MESH_ONLY   FORCE=$FORCE"
echo "=================================================================="
if [ "$DRY" = "1" ]; then
  for A in "${AOA_LIST[@]}"; do for s in "${SPECS[@]}"; do echo "  AoA${A}_${s}"; done; done
  echo "(DRY=1: nothing was run)"; exit 0
fi

mkdir -p "$WORKDIR"

# ---- ROTATED WINGS ---------------------------------------------------------

if [ "$PIVOT_X" = "auto" ] || [ "$PIVOT_Y" = "auto" ]; then
  CEN=$(python3 "$ROTATE_SCRIPT" --find-centroid "$BASE_CASE/constant/geometry/cylinder_2d.obj" | grep "Area centroid" | tr -d "(),")
  PIVOT_X=$(echo "$CEN" | awk '{print $3}')
  PIVOT_Y=$(echo "$CEN" | awk '{print $4}')
fi
echo "Rotation pivot (COM): ($PIVOT_X, $PIVOT_Y)"

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

# ---- HELPERS ---------------------------------------------------------------

reset_case () {
  # strip a case back to 0/, constant/{geometry,...}, system/
  local dir="$1"
  rm -rf "$dir/constant/polyMesh" "$dir/postProcessing"
  rm -f  "$dir"/log.*
  find "$dir" -maxdepth 1 -regextype posix-extended -regex '.*/[0-9]+(\.[0-9]+)?' ! -name '0' -exec rm -rf {} + 2>/dev/null || true
}

# step <cmd> [args]: run an OpenFOAM command, log stdout+stderr to log.<cmd>,
# and return the COMMAND's exit status (pipefail), not tee's.
step () { local name="$1"; "$@" 2>&1 | tee "log.$name"; }

# complete = solver reached its end AND wrote the forces
is_complete () {
  [ -f "$1/postProcessing/forces/0/forceCoeffs.dat" ] && grep -q "^End" "$1/log.pimpleFoam" 2>/dev/null
}

# delete numeric time folders outside [t0, t1]; always keep 0 and the last time
prune_times () {
  local dir="$1" t0="$2" t1="$3" last d
  last=$(ls -1 "$dir" | grep -E "^[0-9]+(\.[0-9]+)?$" | sort -g | tail -1)
  for d in $(ls -1 "$dir" | grep -E "^[0-9]+(\.[0-9]+)?$"); do
    [ "$d" = "0" ] && continue
    [ "$d" = "$last" ] && continue
    if awk -v t="$d" -v a="$t0" -v b="$t1" 'BEGIN{exit !(t>=a && t<=b)}'; then continue; fi
    rm -rf "$dir/$d"
  done
}

# ---- MAIN LOOP -------------------------------------------------------------

FAILED_CASES=()

for AOA in "${AOA_LIST[@]}"; do
for entry in "${SPECS[@]}"; do
  read -r TAG XMIN XMAX YMIN YMAX WLMIN WLMAX B1 B2 B3 <<< "$entry"
  CASE_DIR="$WORKDIR/AoA${AOA}_${TAG}"

  echo "=================================================================="
  echo " AoA=${AOA}  $TAG   x=[$XMIN,$XMAX] y=[$YMIN,$YMAX]  levels: wing($WLMIN $WLMAX) boxes($B1 $B2 $B3)"
  echo "=================================================================="

  if [ "$FORCE" != "1" ] && is_complete "$CASE_DIR"; then
    echo "  -> already complete, skipping (FORCE=1 to redo): $CASE_DIR"
    continue
  fi

  rm -rf "$CASE_DIR"
  cp -r "$BASE_CASE" "$CASE_DIR"
  reset_case "$CASE_DIR"

  # the case's own run.sh starts with foamCleanCase, which deletes the mesh,
  # logs, results and time folders -- remove it from sweep copies
  rm -f "$CASE_DIR/run.sh" "$CASE_DIR/run_preprocessing.sh"

  # write controls: keep every written time (pruned to the viz window later)
  sed -i -e "s/^writeInterval[[:space:]].*/writeInterval   ${WRITE_INTERVAL};/" -e "s/^purgeWrite[[:space:]].*/purgeWrite      0;/" "$CASE_DIR/system/controlDict"
  grep -q "^purgeWrite      0;" "$CASE_DIR/system/controlDict" || echo "  !! WARNING: could not set purgeWrite in $CASE_DIR/system/controlDict"

  # wing for THIS angle; the filename never changes so snappyHexMeshDict needs no edit
  cp "$ROT_DIR/cylinder_2d_aoa${AOA}.obj" "$CASE_DIR/constant/geometry/cylinder_2d.obj"

  NX=$(awk -v a="$XMIN" -v b="$XMAX" -v d="$DENSITY_X" 'BEGIN{printf "%d", (b-a)*d + 0.5}')
  NY=$(awk -v a="$YMIN" -v b="$YMAX" -v d="$DENSITY_Y" 'BEGIN{printf "%d", (b-a)*d + 0.5}')
  echo "  -> background divisions NX=$NX NY=$NY"

  sed -e "s/@XMIN@/${XMIN}/g" -e "s/@XMAX@/${XMAX}/g" -e "s/@YMIN@/${YMIN}/g" -e "s/@YMAX@/${YMAX}/g" -e "s/@NX@/${NX}/g" -e "s/@NY@/${NY}/g" "$BLOCK_TEMPLATE" > "$CASE_DIR/system/blockMeshDict"

  # locationInMesh: 1 unit inside this case's inlet (must be in the fluid, in the domain)
  # refinementBox1 starts exactly at this case's inlet (no buffer)
  LOCX=$(awk -v a="$XMIN" 'BEGIN{printf "%.3f", a + 1}')
  BOX1_XMIN=$(awk -v a="$XMIN" 'BEGIN{printf "%.3f", a}')
  sed -e "s/@WING_LMIN@/${WLMIN}/g" -e "s/@WING_LMAX@/${WLMAX}/g" -e "s/@BOX1_LEVEL@/${B1}/g" -e "s/@BOX2_LEVEL@/${B2}/g" -e "s/@BOX3_LEVEL@/${B3}/g" -e "s/@LOCX@/${LOCX}/g" -e "s/@LOCY@/0/g" -e "s/@BOX1_XMIN@/${BOX1_XMIN}/g" "$SNAPPY_TEMPLATE" > "$CASE_DIR/system/snappyHexMeshDict"

  # the && chain stops this case at the first failing step
  if ( cd "$CASE_DIR" \
       && step blockMesh \
       && step snappyHexMesh -overwrite \
       && step extrudeMesh \
       && step createPatch -overwrite \
       && step checkMesh -latestTime \
       && { if [ -f constant/polyMesh/points ]; then rm -rf mesh_backup && cp -r constant/polyMesh mesh_backup; else echo "  !! WARNING: no constant/polyMesh/points after meshing. Mesh files found at:"; find . -name points; fi; } \
       && { [ "$MESH_ONLY" = "1" ] || step pimpleFoam; } ); then
    if [ "$MESH_ONLY" = "1" ]; then
      echo "  -> mesh built (MESH_ONLY): $CASE_DIR  (backup in $CASE_DIR/mesh_backup)"
    elif [ -f "$CASE_DIR/postProcessing/forces/0/forceCoeffs.dat" ]; then
      if [ "$PRUNE" = "1" ]; then prune_times "$CASE_DIR" "$KEEP_T0" "$KEEP_T1"; fi
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
echo "All $SWEEP cases finished. Next: python3 post_process.py $WORKDIR/AoA* --avg-fraction 0.33 --out ${SWEEP}_sweep"
