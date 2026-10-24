# airfoil2d_transient — Case & Script Documentation

Reference for the OpenFOAM domain- and mesh-independence study on the 2D airfoil case, kept in step with the scripts in this package. The shorter hand-off summary is `PROJECT_CONTEXT.md`; this file is the detailed reference.

## Contents

| Part | What it covers |
|---|---|
| Quick start | The whole workflow in a few commands |
| 1 | What every file in the base case does |
| 2 | The two templates (`blockMeshDict.template`, `snappyHexMeshDict.template`) |
| 3 | `run_domain_sweep.sh`, the one script that runs all four studies |
| 4 | `post_process.py` |
| 5 | Angle of attack: method and `rotate_wing.py` |
| 6 | Saving time steps and looking at the results in ParaView |
| 7 | Running the studies: order, and how to judge convergence |
| 8 | Troubleshooting: symptoms, causes, fixes |
| Appendix A | Reference tables (environment variables, settings, outputs) |
| Appendix B | Change history of the scripts |

---

## Quick start

1. **Folder layout.** One folder containing `run_domain_sweep.sh`, `run_mesh_sweep.sh`, `blockMeshDict.template`, `snappyHexMeshDict.template`, `rotate_wing.py`, `post_process.py`, and a **clean** copy of the case named `airfoil2d_transient/` (no `constant/polyMesh`, no time folders except `0`, no `postProcessing`, no `log.*`). Set `endTime` in that copy's `system/controlDict` (original 500; the study uses 150).
2. **OpenFOAM 12 sourced** in the terminal (`blockMesh -help` should work). Python 3 is needed for `rotate_wing.py` and `post_process.py`.
3. **See the plan first:** `AOA="0 15 25" DRY=1 bash run_domain_sweep.sh` lists every case and starts nothing.
4. **Build the meshes only** (minutes, no solver): `AOA="0 15 25" MESH_ONLY=1 bash run_domain_sweep.sh`. Open one case in ParaView and check the wing, the refinement and, for AoA above 0, the tilt direction (Part 6).
5. **Run for real**, detached so closing the terminal does not kill it:
   `AOA="0 15 25" nohup bash run_domain_sweep.sh > sweep.log 2>&1 &`
   then `tail -f sweep.log`. One case takes about 35 minutes at `endTime 150`.
6. **Post-process:** `python3 post_process.py sweep_domain/AoA* --avg-fraction 0.33 --out upstream_sweep` (Part 4).
7. **Move to the next study** with the converged values (Part 7), for example `SWEEP=downstream UP_FIXED=-10`.

---

## Part 1 — What each file in the base case does

```
airfoil2d_transient/
├── 0/                          initial & boundary conditions (t=0)
├── constant/
│   ├── geometry/cylinder_2d.obj   the wing surface (triangulated/OBJ)
│   ├── physicalProperties         viscosity -> sets Reynolds number
│   ├── momentumTransport          laminar / turbulence model choice
│   └── polyMesh/                  the MESH ITSELF (generated, not edited)
├── system/
│   ├── blockMeshDict              outer domain box + background mesh
│   ├── snappyHexMeshDict          local refinement around the wing
│   ├── extrudeMeshDict            2D "pseudo-3D" thin-slab fix-up
│   ├── createPatchDict            turns Front/Back into empty (true 2D)
│   ├── controlDict                solver, time step, endTime, forces output
│   ├── fvSchemes                  discretization schemes
│   ├── fvSolution                 linear solvers, PIMPLE settings
│   ├── decomposeParDict           parallel run settings (unused unless -parallel)
│   └── residuals                  extra residual-monitoring function object
├── run.sh, run_preprocessing.sh, README.FIRST   the case's own run scripts
│                                  (run.sh starts with foamCleanCase: see below)
├── plot_cl, plot_probesX, plot_probesY          gnuplot helpers at the top level
└── scripts0/plot_cl, plot_cd, plot_coeffs        gnuplot helper scripts
```
A run adds: time folders, `postProcessing/` (force output) and `log.*` files. The case was prepared with OpenFOAM **12** (Foundation).

### `0/U` — initial & boundary velocity
```
internalField   uniform (1 0 0);
boundaryField
{
    Top     { type slip; }
    Bottom  { type slip; }
    Inlet   { type fixedValue; value uniform (1 0 0); }
    Outlet  { type inletOutlet; inletValue uniform (0 0 0); value uniform (1 0 0); }
    wing    { type fixedValue; value uniform (0 0 0); }   // no-slip
    front/back { type empty; }                             // 2D
}
```
- Free-stream speed is `(1 0 0)` — magnitude 1, **along +x**. This is where the flow direction is set.
- `Top`/`Bottom` are `slip` (zero normal velocity, zero-gradient tangential) — they behave like frictionless walls, **not** open far-field boundaries. This matters a lot for AoA (see Part 5).
- `wing` is no-slip (viscous body surface).

### `0/p` — initial & boundary pressure
Standard incompressible setup: `zeroGradient` everywhere except a `fixedValue 0` at the `Outlet` (pressure reference) — nothing here needs to change for AoA or domain sweeps.

### `constant/geometry/cylinder_2d.obj` — the wing surface
This is the actual body geometry `snappyHexMesh` snaps to. The file stores 41 unique (x,y) points, each duplicated at `z=±0.5` (a thin extrusion for the pseudo-2D workflow — see `extrudeMeshDict` below). Reading the points **in file order** (rather than just their bounding box) shows this is a genuine closed, cambered thin-foil cross-section — chord from `x=0` to `x=1`, a blunt trailing edge (~0.006 gap), camber line rising to a maximum of about `y≈0.05` around mid-chord. Its true area centroid (computed by integrating the closed polygon, not just averaging the digitized points — see `rotate_wing.py --find-centroid`) is:
```
Area centroid ≈ (0.41735, 0.04068)
```
This is the pivot used for the angle-of-attack rotation in Part 5.

### `constant/physicalProperties` — sets the Reynolds number
```
viscosityModel  constant;
nu              [0 2 -1 0 0 0 0] 0.01;
```
With chord `C = 1` and `U∞ = 1` (from `0/U`), `Re = U·C/ν = 1×1/0.01 = 100`. To change `Re`, edit `nu` here.

### `constant/momentumTransport`
```
simulationType laminar;
```
No turbulence model — fully laminar, consistent with Re = 100.

### `constant/polyMesh/`
This is **output**, not input — it's the actual generated mesh (points/faces/owner/neighbour/boundary + `cellLevel`/`pointLevel` refinement history from `snappyHexMesh`). The sweep script deletes this directory and regenerates it fresh for every case (`reset_case()`), which is why a "clean" `BASE_CASE` copy shouldn't need it at all. After meshing, the script also saves a copy of the finished mesh in `<case>/mesh_backup/`; if a case ever loses its mesh files (ParaView then says "contains no meshes"), restore it with `cp -r mesh_backup/. constant/polyMesh/` (see Part 8).

### `system/blockMeshDict` — the outer domain box
Defines 8 corner vertices of a rectangular box and a background hex block with divisions `(NX NY 1)`. `blockMeshDict.template` parameterizes this as `@XMIN@ @XMAX@ @YMIN@ @YMAX@ @NX@ @NY@` (see Part 2). The z-extent (-0.05 to 0.05) is fixed. The 6 boundary patches (`Top`, `Bottom`, `Inlet`, `Outlet`, `Front`, `Back`) are named here and referenced by `0/U`, `0/p`.

### `system/snappyHexMeshDict` — local mesh refinement
Key entries:
- `geometry.wing` — points to `constant/geometry/cylinder_2d.obj`.
- `geometry.refinementBox1/2/3` — three nested boxes (decreasing size, increasing refinement level) that locally refine the mesh near the wing and its wake. **Their coordinates are hard-coded absolute numbers in the original file** — in a domain sweep that is a trap for `refinementBox1`, whose start now follows the inlet (see Part 2).
- `castellatedMeshControls.refinementSurfaces.wing.level (4 5)` — min/max refinement level right at the wing surface.
- `castellatedMeshControls.locationInMesh (-5 0 0)` — a single point snappyHexMesh uses to identify "this side is fluid, the wing is solid." **Must be inside the domain and outside the wing** for every case, which is why the template sets it per case (see Part 2).
- `addLayersControls` — prism boundary-layer settings at the wing surface (`nSurfaceLayers 3`, `finalLayerThickness`, etc.).

### `system/extrudeMeshDict` — 2D thin-slab fix-up
`blockMesh`+`snappyHexMesh` build a thin 3D slab (z from −0.05 to −0.15 after this step — yes, it's asymmetric, that's intentional in this workflow). `extrudeMesh` re-extrudes the `Front` patch by a fixed thickness to produce a single clean layer of hex cells in z, replacing whatever the 3D snap produced there. This is a standard trick for building genuinely 2D OpenFOAM cases through the inherently-3D `snappyHexMesh`.

### `system/createPatchDict`
Converts the `Front`/`Back` patches (created by `extrudeMesh`) into `type empty` patches — this is what makes the case run as 2D rather than 3D.

### `system/controlDict` — solver, time control, and force output
```
application  pimpleFoam;
endTime      500;
deltaT       0.005;
```
Also contains the function objects that write your data:
```
forceCoeffsIncompressible(patches=(wing), magUInf=1, lRef=1, Aref=0.1,
    CofR=(0 0 0), liftDir=(0 1 0), dragDir=(1 0 0), pitchAxis=(0 0 1))
```
Also in this file: `writeControl runTime;`, `writeInterval 0.1;` and `purgeWrite 10;` at the top level (only the last 10 written time folders are kept). The original case has `endTime 500`; **the study has used `endTime 150`** because the results had converged by about t = 100. The sweep script does not touch `endTime` (set it in the clean base case), but it overrides `writeInterval` and `purgeWrite` in each sweep copy (see Part 6).

The `forceCoeffsIncompressible` object writes `postProcessing/forces/0/forceCoeffs.dat` (columns `Time Cm Cd Cl Cl(f) Cl(r)`), which `post_process.py` parses. Note **`liftDir`/`dragDir` are fixed lab-frame directions** (y / x), independent of the wing's orientation — this is exactly why rotating the *geometry* for AoA (rather than the inflow) is convenient: drag stays "force along free-stream" and lift stays "force perpendicular to free-stream" with zero changes needed here.

### `system/fvSchemes` — discretization
```
divSchemes { div(phi,U) Gauss linear; }   // <- central differencing, unbounded
```
Worth knowing: `Gauss linear` for convection has no numerical diffusion. It's more accurate than upwind schemes but has zero built-in damping of grid-scale noise — on a locally poor-quality mesh cell, this scheme is more likely to let spurious oscillations grow than a bounded/upwind-biased scheme would. It is one possible contributor to grid-scale noise on a poor-quality cell. It is not something you need to change, but keep it in mind if a case shows irregular oscillations (the cause of the one such case in the upstream study was never isolated; see Part 8).

### `system/fvSolution` — linear solvers & PIMPLE
`nCorrectors 2`, `nNonOrthogonalCorrectors 4`, `momentumPredictor yes` — standard PIMPLE transient settings. `relaxationFactors` (U: 0.9, p: 0.8) apply under-relaxation even though this is a transient run (common for extra stability with `nOuterCorrectors 1`).

### `system/decomposeParDict`
Only relevant if you ever run in parallel (`mpirun ... pimpleFoam -parallel`); unused by the sweep scripts as written (they run each case serially).

### `system/residuals`
An extra function object that logs the initial residual of `p`, `U`, `C` per time step — useful for debugging solver convergence, separate from the forces output.

### `run.sh`, `run_preprocessing.sh`, `README.FIRST`
These are the case's own original run instructions — `run_domain_sweep.sh` replicates exactly this pipeline (`blockMesh → snappyHexMesh → extrudeMesh → createPatch → checkMesh → pimpleFoam`) for each swept case. **Do not run the case's own `run.sh` inside a finished case:** it begins with `foamCleanTutorials; foamCleanCase`, which deletes the mesh, logs, `postProcessing` and time folders. The sweep script removes `run.sh` and `run_preprocessing.sh` from every sweep copy for this reason.

### `scripts0/plot_cl`, `plot_cd`, `plot_coeffs`
Gnuplot scripts for quick-look plotting of the force coefficients directly from `postProcessing/forces/0/forceCoeffs.dat` — `post_process.py` supersedes these for sweep comparisons, but they're handy for a single case.

---

## Part 2 — The two templates

`blockMeshDict.template` and `snappyHexMeshDict.template` are copies of the case's own `blockMeshDict` and `snappyHexMeshDict` in which every value that changes between cases is replaced by a placeholder: a capitalised word wrapped in at-signs. For each case the sweep script runs `sed` on the template and writes the result into that case's `system/` folder. You never edit a sweep case's dictionaries by hand, and you never edit the placeholders in the templates.

| Placeholder | File | Meaning | Value the script puts in |
|---|---|---|---|
| `@XMIN@` `@XMAX@` `@YMIN@` `@YMAX@` | blockMesh | Corners of the domain box | From the case line of the study |
| `@NX@` `@NY@` | blockMesh | Background cells in x and y | `round((XMAX-XMIN) * 0.78125)` and `round((YMAX-YMIN) * 0.75)` |
| `@WING_LMIN@` `@WING_LMAX@` | snappy | Refinement level range at the wing surface | 4 and 5 (baseline; the mesh study changes them) |
| `@BOX1_LEVEL@` `@BOX2_LEVEL@` `@BOX3_LEVEL@` | snappy | Levels of the three nested refinement boxes | 3, 3, 4 (baseline) |
| `@LOCX@` `@LOCY@` | snappy | The `locationInMesh` seed point | `XMIN + 1` and `0` |
| `@BOX1_XMIN@` | snappy | Where `refinementBox1` starts in x | `XMIN` (the inlet) |

**Why the density matters.** `NX` and `NY` are computed from cells per unit length (25 cells over the original 32 units, 12 over 16), so a longer domain gets more cells of the same size, not the same number of bigger cells. Otherwise a domain study would also change the mesh coarseness.

**Why two of the placeholders exist (both were real failures):**
- **`locationInMesh`.** The original file has the fixed point `(-5 0 0)`, correct only for a domain starting at x = -5. For a smaller domain it falls outside the mesh and `snappyHexMesh` stops with "Point (-5 0 0) is not inside the mesh". It is now one unit inside each case's inlet, upstream of the wing and in the fluid.
- **`refinementBox1` start.** It was fixed at x = -3 while the inlet moved. For the smallest domain (inlet at x = -3) it sat exactly on the inlet; for the widest of the three cases then run (inlet at x = -7) it sat 4 units away. That x = -7 case was the one that showed irregular lift oscillations and continuity errors. Starting the box at the inlet in every case removed the inconsistency and the oscillations did not return. This is the most likely explanation, but the cause was not isolated by an experiment (Part 8).

**What is still fixed** (copied from the original case): box1 ends at x = 25 with y limits ±4; box2 covers x from -2 to 10, y ±3; box3 covers x from -1 to 5, y ±2; the boundary-layer settings; the z-extent ±0.05. These matter for the downstream and lateral studies (Part 7).

**Adding another templated value:** put a new `@NAME@` in the template, add `-e "s/@NAME@/${VALUE}/g"` to the matching `sed` line in `run_domain_sweep.sh`, then confirm `grep -c "@" sweep_*/AoA*/system/*Dict` prints 0 for every file (no placeholder left unfilled).

---

## Part 3 — `run_domain_sweep.sh` (the unified sweep driver), step by step

One script runs all four studies. The study and angles are chosen with environment variables, so you rarely need to edit the file:

```bash
AOA="0 15 25" SWEEP=upstream   bash run_domain_sweep.sh
AOA="0 15 25" SWEEP=downstream UP_FIXED=-10 bash run_domain_sweep.sh
AOA="0 15 25" SWEEP=lateral    UP_FIXED=-10 DOWN_FIXED=27 bash run_domain_sweep.sh
AOA="0 15 25" SWEEP=mesh       UP_FIXED=-10 DOWN_FIXED=27 HALF_H=8 bash run_domain_sweep.sh
```
`run_mesh_sweep.sh` is now a 3-line wrapper that sets `SWEEP=mesh`. Extra switches: `DRY=1` (list the cases and stop), `MESH_ONLY=1` (build and check meshes, no solve), `FORCE=1` (redo finished cases), `WORKDIR=...`.

**1. Header and `set -euo pipefail`.** Usage notes, then the safety flags: stop on any error, error on undefined variables, and make a pipeline fail if the command in it fails (not just `tee`).

**2. SETTINGS block.** `BASE_CASE` (clean copy of the case), the two templates, `rotate_wing.py`; `WRITE_INTERVAL=0.5`, `KEEP_T0=60`, `KEEP_T1=120`, `PRUNE=1` (time-folder saving); `DENSITY_X/Y` (background cells per unit length of the original case: 25 cells over 32 units, 12 over 16); `BASE_LEVELS="4 5 3 3 4"` (wing level min/max, then the three box levels); `PIVOT_X/Y="auto"`.

**3. STUDY DEFINITIONS.** A `case "$SWEEP"` block builds `SPECS`, one line per case: `tag XMIN XMAX YMIN YMAX WLMIN WLMAX B1 B2 B3`.

| SWEEP | varies | tags | needs |
|---|---|---|---|
| upstream | XMIN = -3, -5, -7, -10, -15 | `U1_x-3` … | nothing |
| downstream | XMAX = 15, 20, 27, 35, 45 | `D1_xd15` … | `UP_FIXED` |
| lateral | half-height 4, 6, 8, 12, 16 | `L1_h4` … | `UP_FIXED DOWN_FIXED` |
| mesh | refinement levels coarse … xfine | `M1_coarse` … | `UP_FIXED DOWN_FIXED HALF_H` |

The `: "${UP_FIXED:?...}"` lines make the script refuse to start if a fixed value is missing, so it can never silently reuse a wrong default. Angles must be integers.

**4. Summary and `DRY`.** Prints the plan; with `DRY=1` it lists every `AoA<angle>_<tag>` case and exits.

**5. Rotated wings.** With `PIVOT_X/Y="auto"` the pivot is the area centroid from `rotate_wing.py --find-centroid`. One rotated copy of `cylinder_2d.obj` is made per angle in `$WORKDIR/_rotated_wings/` (angle 0 is a plain copy).

**6. Helpers.**
- `reset_case`: deletes `constant/polyMesh`, `postProcessing`, `log.*` and every numeric time folder except `0/` in a freshly copied case.
- `step <cmd>`: runs an OpenFOAM command and writes stdout AND stderr to `log.<cmd>`. OpenFOAM prints fatal errors on stderr, so the old script left empty logs on failure.
- `is_complete`: a case counts as finished only if `log.pimpleFoam` contains `End` and `forceCoeffs.dat` exists.
- `prune_times`: deletes time folders outside [`KEEP_T0`, `KEEP_T1`]; always keeps `0` and the last time.

**7. Main loop**, for each angle and each case:
1. Skip it if `is_complete` (unless `FORCE=1`). This is what stops a re-run from wiping finished results.
2. Delete any old copy, copy `BASE_CASE`, run `reset_case`.
3. Delete `run.sh` and `run_preprocessing.sh` from the copy (the case's `run.sh` begins with `foamCleanCase`, which deletes the mesh and results).
4. `sed` the top-level `writeInterval` and `purgeWrite` in `controlDict` (the anchor `^` skips the indented `forceCoeffs` one) and warn if `purgeWrite` could not be set.
5. Copy that angle's rotated wing over `constant/geometry/cylinder_2d.obj` (the file name never changes, so `snappyHexMeshDict` needs no edit).
6. Compute `NX = (XMAX-XMIN) * DENSITY_X` and `NY = (YMAX-YMIN) * DENSITY_Y` (rounded), so the far-field cell size stays constant, and fill `blockMeshDict.template`.
7. Fill `snappyHexMeshDict.template`: the refinement levels, `locationInMesh` at `(XMIN+1, 0, 0)` (must lie in the fluid and inside the domain), and `refinementBox1` starting exactly at `XMIN`.
8. Run `blockMesh → snappyHexMesh -overwrite → extrudeMesh → createPatch -overwrite → checkMesh -latestTime` chained with `&&`, so a case stops at its first failing step. Then copy `constant/polyMesh` to `mesh_backup/` (with a warning that lists where the mesh files are if `points` is missing), then run `pimpleFoam` unless `MESH_ONLY=1`.
9. On success: prune time folders and print `done`; if `forceCoeffs.dat` is missing, record a failure. On failure: print the last 15 lines of the newest log and record it.

**8. End.** Prints every failed case and exits with an error, or points you to `post_process.py`.

---

## Part 4 — `post_process.py`

```bash
python3 post_process.py sweep_domain/AoA* --avg-fraction 0.33 --out upstream_sweep
```
Cases are grouped by angle and ordered by the swept value (upstream distance, outlet position, half-height; or cell count for mesh studies). Per case it reports `Cd_mean`, `Cl_mean`, `Cl_rms`, `St`, `state`, and two change columns: `d%prev` (versus the previous case) and `d%max` (versus the last case, your best estimate).

- `state` is `steady` (Cl barely fluctuates, means do not drift), `DRIFTING` (steady-looking but still moving, so run longer), `osc` (periodic shedding, settled) or `osc-GROWING`. `St` is shown only for oscillating cases (`Cl_rms` above `--osc-threshold`, default 1e-3); for steady flows a frequency analysis only finds noise.
- The "Independence check" line per angle compares only the last two cases against `--tol` (default 1%). Also read the `d%max` column: one small last step can be a coincidence.
- Outputs: `<out>.csv`, `<out>.png` (Cd, Cl and Cd-difference versus the swept value, one line per angle) and `<out>_timeseries.png` (Cd(t) and Cl(t), one row per angle).
- `--avg-fraction 0.33` averages over the last third of the run (t of about 100 to 150 for `endTime 150`), after the start-up transient.

---

## Part 5 — Angle of attack: the method and `rotate_wing.py`

Your reasoning is correct and matches standard CFD practice: **blockage effects grow with angle of attack.** At AoA, the body's *frontal projected height* (roughly `C·sin(AoA) + t·cos(AoA)`, where `t` is thickness) increases, so for a domain of fixed height the *blockage ratio* (body height ÷ domain height) grows — the slip walls "feel" the body more, and the flow needs more room, both upstream and possibly laterally, to relax back to the true free-stream condition before it reaches the inlet's influence. So yes: a domain that's independent at AoA = 0° is not guaranteed to still be independent at AoA = 30°, and re-running the domain study at a few AoA values (not just AoA = 0°) is the right instinct if you want a single domain size that's genuinely safe across your whole AoA range.

### The key design decision: rotate the geometry, not the inflow

Because `Top`/`Bottom` are `slip` walls **aligned with the x-axis**, and `liftDir=(0 1 0)`/`dragDir=(1 0 0)` in `controlDict` are fixed lab-frame directions, the clean way to change AoA here is to **rotate the wing surface** (`cylinder_2d.obj`) about its center of mass and keep the free-stream (`0/U`) and slip walls exactly as they are now. If you instead tilted the inlet velocity, the flow would have a nonzero y-component right at the `Top`/`Bottom` slip walls, which enforce zero normal velocity — that would create artificial, non-physical shear right at your outer boundary, exactly what a domain-independence study is trying to avoid. Rotating the geometry keeps every boundary condition and every `forceCoeffs` setting in `controlDict` (lift/drag directions, `CofR`, `pitchAxis`) unchanged — only the triangulated surface file itself changes.

### `rotate_wing.py` — rotates the wing about its center of mass

Rather than composing OpenFOAM's `surfaceTransformPoints` (which only rotates about the origin, and this wing isn't centered on the origin), use the included `rotate_wing.py`, which does a proper translate → rotate → translate-back about an arbitrary pivot, in one step, with a built-in visual sanity check.

**1. Find the pivot** (already done for you, but this is how, and how to redo it if you ever change geometry):
```bash
python3 rotate_wing.py --find-centroid constant/geometry/cylinder_2d.obj
#   Area centroid:  (0.41735, 0.04068)
```
This is the true area centroid of the closed 2D cross-section (computed by polygon integration, not a naive average of the digitized points) — the correct pivot for a "center of mass" rotation of a uniform-density thin body.

**2. Generate the rotated geometry**, with a sanity-check plot:
```bash
python3 rotate_wing.py constant/geometry/cylinder_2d.obj \
    constant/geometry/cylinder_2d_aoa10.obj \
    --pivot 0.41735 0.04068 --aoa 10 --plot aoa10_check.png
```
**Sign convention (verified, not assumed):** I tested this at AoA=+15° and plotted the before/after outline — the leading edge (near x=0) rotates *up* and the trailing edge (near x=1) rotates *down*, which is the standard "nose up" positive angle of attack into a free stream moving in +x. The pivot point itself stays exactly fixed in space (true rotation, not an accidental translation) — confirmed visually, not just by the math. **Always regenerate `--plot` and glance at it whenever you change the pivot or geometry file**, since the sign of a rotation depends on the point-ordering convention baked into the original `.obj`, which isn't guaranteed to be the same for a different geometry file.

**3. Generate the whole family in one loop**, e.g. for 10° steps from 0° to 30°:
```bash
for a in 0 10 20 30; do
  python3 rotate_wing.py constant/geometry/cylinder_2d.obj \
      constant/geometry/cylinder_2d_aoa${a}.obj \
      --pivot 0.41735 0.04068 --aoa $a
done
```

### What has to change per AoA case
- Only the wing file itself: the sweep script copies that angle's rotated wing over `constant/geometry/cylinder_2d.obj`, so `snappyHexMeshDict` keeps pointing at the same file name.
- **Nothing else** needs to change in `0/U`, `0/p`, or `controlDict` — that's the payoff of rotating the geometry instead of the flow.
- Worth double-checking, not necessarily changing: the three `refinementBox` extents in `snappyHexMeshDict` currently have a fixed, fairly narrow lateral (`y`) half-width (box1: ±4, box2: ±3, box3: ±2). At higher AoA the wing and its wake will sit further off the `y=0` centerline — check that the rotated wing (and its expected wake) still comfortably fits inside these boxes; if not, they'll need widening in `y` for the higher-AoA cases.
- Per your own point: at higher AoA you may want the **domain height** (`YMIN`/`YMAX`), not just the upstream distance, to grow — worth sweeping both, or at minimum checking blockage ratio (projected frontal height ÷ domain height) stays comparably small across your AoA range.

### How it's wired into `run_domain_sweep.sh`

- The angles come from the environment: `AOA="0 5 10 15 20 25" bash run_domain_sweep.sh` (integers; default is just `0`).
- The pivot is computed automatically: `PIVOT_X/PIVOT_Y="auto"` runs `rotate_wing.py --find-centroid` on the wing file, so a different airfoil gets its own centre of mass. Set both to numbers to override.
- Before the main loop, one rotated wing per angle is written to `$WORKDIR/_rotated_wings/`; the main loop is nested angle-outside, case-inside, and case folders are named `AoA<N>_<tag>`, for example `AoA10_U2_x-5`.
- Each case gets its angle's rotated wing copied over `constant/geometry/cylinder_2d.obj`. The file name never changes, so `snappyHexMeshDict` needs no edit.

Tested with stubbed OpenFOAM commands and the real `rotate_wing.py`: the angle-0 wing is byte-identical to the original, and re-measuring the centroid of a rotated wing returns the same point, so the body turns about its centre of mass without drifting.

**Practical order:** run `MESH_ONLY=1` first for the angles you plan to use, open one rotated case in ParaView to check the wing tilts the right way (positive = leading edge up), then run the solver.

---

## Part 6 — Saving time steps, and looking at the results in ParaView

### What gets saved
Forces are written every time step into `postProcessing/forces/0/forceCoeffs.dat` whatever the settings below, so `post_process.py` never depends on the saved fields.

The original `controlDict` writes fields every 0.1 time units and keeps only the last 10 (`purgeWrite 10`), so only t = 149.1 to 150 survive, which is no use for animation. The sweep script therefore changes each sweep copy as follows:
1. `writeInterval` becomes `WRITE_INTERVAL` (default **0.5**) and `purgeWrite` becomes 0 (keep everything). Only the top-level lines are changed; the indented `writeInterval` of the force function object is left alone.
2. While a case runs, all its time folders exist (300 for `endTime 150`).
3. After the case finishes successfully, time folders **outside `KEEP_T0` to `KEEP_T1`** (default **60 to 120**) are deleted. Folder `0` and the final time are always kept, so you end up with 0, 60, 60.5 … 120 and 150.

Change these at the top of `run_domain_sweep.sh`: `WRITE_INTERVAL` (0.25 gives smoother animations, more disk), `KEEP_T0`, `KEEP_T1`, and `PRUNE=0` to keep every written folder. Folder sizes were not measured; after the first case check `du -sh` on it before running many.

**Cases finished before this feature existed** (the first AoA 0 to 25 upstream runs) only have the last 10 time folders. The skip-finished rule will not redo them; use `FORCE=1` if you want animations of them.

At the angles run so far the flow settled to a steady state, so one snapshot shows everything. If shedding appears at a higher angle or a different Reynolds number, the 60 to 120 window is where to look.

### Opening a case in ParaView
1. In the case folder run `touch case.foam` (any name ending in `.foam`).
2. Open it in ParaView; in Properties tick `internalMesh` and the patches `wing`, `Inlet`, `Outlet`, `Top`, `Bottom` (not `Front`/`Back`, which are `empty`). Under Cell Arrays tick `U` and `p`. Press **Apply**.
3. Press the **+Z** view button to look straight at the x–y plane.
4. **Mesh:** representation "Surface With Edges", zoom to the wing. Check the wing outline and tilt, that refinement is finest at the wing and wake, that the refined region starts at the inlet, and that the 3 boundary layers are continuous.
5. **Flow:** go to the last time (or any time from 60 to 120), colour by `U` (Magnitude) or `p`. For **vorticity** apply Filters → Gradient on `U` with *Compute Vorticity* ticked and colour by `Vorticity` (Z component, a symmetric range such as ±10 works well). For **streamlines** use Filters → Stream Tracer with vectors `U` and a line seed across the flow just upstream of the wing, for example from (-1, -1, 0) to (-1, 1, 0). Show the `wing` patch in a solid colour.
6. **Mesh quality:** run `checkMesh -latestTime -writeAllFields` in the case folder, reload in ParaView and colour by fields such as `nonOrthoAngle`, `cellAspectRatio` or `skewness` to see where the worst cells are.

### If ParaView says "contains no meshes"
`constant/polyMesh` is missing `points`, `faces`, `owner`, `neighbour`, `boundary`. First try `cp -r mesh_backup/. constant/polyMesh/`. If there is no backup, rebuild the mesh in a **copy** (this does not touch your results):
```bash
cp -r AoA0_U1_x-3 AoA0_U1_x-3_viz && cd AoA0_U1_x-3_viz
rm -rf constant/polyMesh
blockMesh && snappyHexMesh -overwrite && extrudeMesh && createPatch -overwrite
checkMesh -latestTime | grep -E "points:|cells:"      # must match the original log.checkMesh
```
If the cell and point counts match the case's original `log.checkMesh`, the saved fields line up with the rebuilt mesh; if not, do not trust the pictures. Do not use the case's `run.sh` for this.

---

## Part 7 — Running the studies: order, and how to judge convergence

**Order.** Upstream, then downstream, then lateral (recommended, see below), then mesh. Each study holds everything else fixed, and takes the fixed values from the studies before it, which is why the script demands them (`UP_FIXED`, `DOWN_FIXED`, `HALF_H`).

| Study | Command | Case values |
|---|---|---|
| Upstream | `AOA="…" SWEEP=upstream bash run_domain_sweep.sh` | XMIN = -3, -5, -7, -10, -15 (XMAX 27, y ±8) |
| Downstream | `… SWEEP=downstream UP_FIXED=<xmin> …` | XMAX = 15, 20, 27, 35, 45 |
| Lateral | `… SWEEP=lateral UP_FIXED=<xmin> DOWN_FIXED=<xmax> …` | half-height 4, 6, 8, 12, 16 |
| Mesh | `… SWEEP=mesh UP_FIXED=… DOWN_FIXED=… HALF_H=<h> …` | levels 3 4 2 2 3 / 4 5 3 3 4 / 5 6 4 4 5 / 6 7 5 5 6 |

**Why angle matters.** Blockage grows with angle, so a domain that is large enough at 0° may not be at 25°. Run each study at several angles (for example 0, 15, 25) and take the **largest** requirement across them as the value you fix for all angles. To save time, the later studies can use fewer angles than the upstream one.

**Judging a result.** For each angle `post_process.py` prints the change between consecutive cases (`d%prev`), the difference from the last case (`d%max`), and an "Independence check" comparing the last two cases to a tolerance (default 1%; the tolerance is your decision, change it with `--tol`).
- Read both columns, not only the verdict: one small last step can be a coincidence, and `d%max` shows whether the earlier cases were close to the best estimate.
- A study has plateaued only when the values stop changing as the swept parameter grows. In the first upstream results (x = -3, -5, -7) Cd was still falling by 3 to 10% between the last two cases at five of the six angles (0, 5, 10, 15, 25), so those distances were not yet enough; at 20° it was not even monotonic, which is a reason to look at its meshes in ParaView before trusting it.
- Before comparing numbers, check the `state` column: cases should be `steady` (or `osc`, settled). `DRIFTING` and `osc-GROWING` mean the run needs to be longer. Compare **Cl_mean** as well as Cd, especially above 0°.
- Look at `*_timeseries.png` too; a flat line after the start-up spike is what a converged steady case looks like.

**Things these studies do not cover:**
- Only the **upstream** distance has been varied so far. The half-height (slip walls at ±8, which act like a channel) and the downstream length are fixed, so upstream results say nothing about them.
- **Downstream:** `refinementBox1` ends at x = 25, so for XMAX = 35 or 45 the extra wake lies on the coarse background mesh (cells about 1.3 chords). That may affect the result; consider extending the box with XMAX. For XMAX below 25 the box simply extends past the domain, which is fine.
- **Lateral:** the refinement boxes have fixed y limits (±4, ±3, ±2); at high angles check the wing and wake still sit well inside them.
- **Mesh:** the refinement levels were identical (4 5 3 3 4) in all domain studies. The background cell size is also held constant; only the `snappyHexMesh` levels change in the mesh study.
- **Time step and averaging window:** `deltaT 0.005` and the `--avg-fraction` window were not studied.

---

## Part 8 — Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| bash: `syntax error near unexpected token '>'` on a `sed`/`find` line | Windows line endings (CRLF) breaking backslash line continuations | `sed -i 's/\r$//' file`. The script now keeps those commands on single lines; edit it only with an editor that keeps Unix line endings. |
| `FOAM FATAL ERROR: Point (-5 0 0) is not inside the mesh` | Hard-coded `locationInMesh` outside a smaller domain | Already fixed: the template sets it to `XMIN + 1`. If you edit the case's own dict, keep the point inside the domain and outside the wing. |
| One case shows huge, irregular `Cl` swings and continuity errors; its checkMesh has higher non-orthogonality and skewness than the others | Seen only in the x = -7 case (the widest then) while `refinementBox1` was fixed at x = -3 | `refinementBox1` now starts at the inlet and the problem did not return. Cause not proven by an isolated experiment. Compare `checkMesh` numbers across cases if it reappears. |
| A log file is empty or missing, or `post_process.py` cannot find `postProcessing` | The run died early or the case folder was wiped. The old script logged only standard output, so error messages never reached the logs | The script now logs standard output and error to `log.<command>`, stops that case at the first failing step, prints the end of the newest log, and lists failed cases at the end. If a command says "command not found", OpenFOAM is not sourced in that terminal. |
| ParaView: "case.foam contains no meshes" | `constant/polyMesh` lost its mesh files. This happened once; the cause was not found | Restore from `mesh_backup/` or rebuild in a copy (Part 6). The script skips finished cases, deletes `run.sh` from sweep copies and keeps a backup to make a repeat unlikely. If it recurs, note what was run in that folder just before. |
| Re-running the script did not redo a case | The case is already complete (`log.pimpleFoam` ends with `End` and the forces file exists) | `FORCE=1 bash run_domain_sweep.sh` |
| Only time folders 149.1 to 150 exist | `purgeWrite 10` in a run made before the saving feature was added | Re-run those cases with `FORCE=1` if you need the fields (Part 6). |
| `St` around 40, "still transient" on a flat signal, or large % changes between different angles | The first version of `post_process.py` analysed noise and compared rows in alphabetical order | Fixed in the current `post_process.py` (Part 4). |
| `SWEEP=downstream needs UP_FIXED …` and the script exits | Intended: a study will not guess the fixed values from earlier studies | Pass the value, for example `UP_FIXED=-10`. |
| `AOA values must be integers` | Angles like 7.5 are not supported | Use whole degrees. |
| First thought for a flat Cl was "shedding takes long to start" | That reading was wrong for this case: low-Re flow around this foil is steady at the angles run | `state` in `post_process.py` tells steady from oscillating, so judge by that. |

---

## Appendix A — Reference tables

**Environment variables of `run_domain_sweep.sh`**

| Variable | Meaning | Default |
|---|---|---|
| `SWEEP` | `upstream`, `downstream`, `lateral` or `mesh` | `upstream` |
| `AOA` | Space-separated integer angles in degrees; positive = leading edge up | `0` |
| `UP_FIXED` | XMIN held fixed (downstream, lateral, mesh) | required for those studies |
| `DOWN_FIXED` | XMAX held fixed (lateral, mesh) | required for those studies |
| `HALF_H` | Half-height held fixed (mesh) | required for the mesh study |
| `DRY=1` | List the cases and stop | off |
| `MESH_ONLY=1` | Build and check meshes, no solve | off |
| `FORCE=1` | Redo cases that are already complete | off |
| `WORKDIR` | Output folder | `sweep_domain`, `sweep_downstream`, `sweep_lateral` or `sweep_mesh` |

**Settings at the top of the script**

| Setting | Value | Meaning |
|---|---|---|
| `BASE_CASE` | `./airfoil2d_transient` | Clean case that is copied for every sweep case |
| `WRITE_INTERVAL` | 0.5 | Field write interval for sweep copies |
| `KEEP_T0`, `KEEP_T1` | 60, 120 | Time window kept after a case finishes |
| `PRUNE` | 1 | 0 keeps every written time folder |
| `DENSITY_X`, `DENSITY_Y` | 0.78125, 0.75 | Background cells per unit length |
| `BASE_LEVELS` | `4 5 3 3 4` | Wing level min/max, then box1/box2/box3 levels |
| `PIVOT_X`, `PIVOT_Y` | `auto` | Rotation pivot; `auto` = area centroid of the wing |

**What a case folder contains** (`<WORKDIR>/AoA<angle>_<tag>/`): `log.blockMesh`, `log.snappyHexMesh`, `log.extrudeMesh`, `log.createPatch`, `log.checkMesh`, `log.pimpleFoam` (standard output and error); `mesh_backup/`; `constant/polyMesh/`; `postProcessing/forces/0/forceCoeffs.dat`; time folders `0`, 60 to 120 and the final time. Rotated wings are in `<WORKDIR>/_rotated_wings/`.

**`post_process.py` options:** `--avg-fraction` (default 0.5; use 0.33 for the last third), `--osc-threshold` (default 1e-3, `Cl_rms` above which a case counts as oscillating), `--drift-tol` (default 0.2%), `--tol` (default 1%), `--out`.

**`post_process.py` CSV columns:** `case, aoa, x_value, N_cells, Cd_mean, Cd_pct_prev, Cd_pct_vs_largest, Cl_mean, Cl_pct_prev, Cl_pct_vs_largest, Cl_rms, St, state, Cd_drift_pct, Cl_drift_pct`. ("vs_largest" means versus the last case of that angle: the largest domain or finest mesh.)

**Case tags understood by `post_process.py`:** `U1_x-3` (upstream), `D2_xd20` (downstream), `L3_h8` (half-height), `M2_medium` (mesh, ordered by cell count).

---

## Appendix B — Change history of the scripts

| # | Change | Why |
|---|---|---|
| 1 | First scripts: a domain sweep, a mesh sweep and `post_process.py`, driven by the two templates; background cells scaled with domain size | Make the studies repeatable and keep far-field cell size constant |
| 2 | Multi-line `sed`/`find` commands put on single lines | Bash syntax error on the user's WSL setup (line endings) |
| 3 | `locationInMesh` templated as `XMIN + 1`; `snappyHexMeshDict` generated for the domain sweep as well | The fixed point fell outside smaller domains |
| 4 | Domain sweep changed to upstream distance only (XMIN = -3, -5, -7) | The user wanted to isolate the upstream effect |
| 5 | `post_process.py`: saturation check and raw Cd/Cl time-history plot | To see whether a case had settled before comparing numbers |
| 6 | `refinementBox1` start templated: first a 2-unit buffer from the inlet, then no buffer (start = XMIN) | Box start was fixed while the inlet moved; the 2-unit buffer was a stop-gap with no rigorous basis and was dropped on the user's suggestion |
| 7 | Case documentation written | Explain every case file and script |
| 8 | Angle of attack: `rotate_wing.py` (centroid pivot, sign check plot); angle loop in the sweep script; later `PIVOT_X/Y=auto` | Study blockage at several angles; each airfoil gets its own centre of mass |
| 9 | Standard error logged; `&&` chain per case; check that `forceCoeffs.dat` exists; failed cases listed at the end | Failures left empty logs and silent gaps |
| 10 | `purgeWrite 0` with `WRITE_INTERVAL`, pruning to t = 60 to 120; skip finished cases; `run.sh` removed from sweep copies; `mesh_backup/`; `MESH_ONLY` | Animations in ParaView; protect results and meshes from accidental wipes |
| 11 | `post_process.py` rewritten: grouped by angle, `St` only for oscillating flow, steady/drifting/oscillating states, `Cl_mean`, independence check, per-angle plots | The first version produced false alarms on steady flows and compared rows across angles |
| 12 | One unified script with `SWEEP=upstream/downstream/lateral/mesh`, environment variables, `DRY`, `FORCE`; `run_mesh_sweep.sh` reduced to a wrapper; `post_process.py` reads all case tags; angle list now `AOA` instead of an `AOA_LIST` array in the file | The mesh sweep had not received the safety features; avoid two diverging scripts |
| 13 | Template header comments rewritten (they still pointed to arrays that no longer exist); banner shows OpenFOAM 12 | Keep the comments true; comments only, real content unchanged and re-tested |
| 14 | `PROJECT_CONTEXT.md` and the zip of all files | A compact hand-off for a new conversation |

All script changes were tested with stubbed OpenFOAM commands (dry runs, every study mode, skipping, failures, pruning) and with the real Python scripts; real solver runs were done by the user.

---
