# Project context: domain & mesh independence study (airfoil2d_transient, OpenFOAM)

Read this first. It replaces the full earlier chat. Everything marked **[FILL IN]** is a value the user decided or measured that was never reported back in the chat; ask the user for it.

---

## 1. Goal

Run **domain-independence studies (upstream, downstream, and probably lateral) and then a mesh-independence study** on a 2D OpenFOAM airfoil case, **at several angles of attack (AoA)**. The user's reasoning: at higher AoA the body blocks more flow, so a domain that is large enough at 0° may not be large enough at 25°. The order the user wants: upstream, then downstream, then mesh.

Reference paper (provided by the user): Shukla, Sharma, Agrawal, Bhardwaj, *Hydrodynamics and propulsion of a hydrofoil undergoing leading-edge pitching and traveling wave-based surface undulation*, Phys. Fluids 36, 043605 (2024). It uses an in-house sharp-interface immersed-boundary code, a NACA0012, Re = 5000, and a 30C × 20C domain (10C upstream of the leading edge, 10C above and below, ~20C downstream). **The paper does not describe its mesh/domain tests** (validation was "not reproduced for brevity"), so it served only as a setup reference. The user's case is a different, simpler one (below).

## 2. Environment

- OpenFOAM **12** (Foundation, build 12-86e126a7bc4d) inside **WSL**; case folders live on a Windows drive (`/mnt/d/...`). Python 3, gnuplot and ParaView (Linux build, opened from WSL) are available to the user.
- The assistant's sandbox has **no OpenFOAM**. All scripts were tested only with **stubbed OpenFOAM commands** (fake `blockMesh` etc.) plus the real Python scripts. The user ran the real cases and reported results. Treat anything about solver behaviour as unverified unless the user showed it.
- Windows line endings (CRLF) once broke a script; if bash reports odd syntax errors, run `sed -i 's/\r$//' file`.

## 3. The base case (what is in the user's `airfoil2d_transient/`)

- Solver `pimpleFoam`, **laminar**, 2D (single cell in z, `empty` front/back). `nu = 0.01`, inflow `U = (1 0 0)`, chord = 1, so **Re = 100**. `deltaT = 0.005`. The user changed `endTime` from 500 to **150** because results had converged by t ≈ 100. One case takes about **35 min** at endTime 150.
- Boundaries: Inlet fixedValue (1 0 0); Outlet inletOutlet (p fixed 0); **Top/Bottom are `slip` walls** (not open far-field); wing no-slip.
- Force output: `forceCoeffsIncompressible` (patch `wing`, magUInf 1, lRef 1, Aref 0.1, CofR (0 0 0), liftDir (0 1 0), dragDir (1 0 0)) writes `postProcessing/forces/0/forceCoeffs.dat`, columns `Time Cm Cd Cl Cl(f) Cl(r)`. Lift and drag directions are fixed lab-frame directions, which is why AoA is made by rotating the **geometry**, not the inflow.
- Original domain: x from -5 to 27, y from -8 to 8, background mesh `NX=25, NY=12`. Refinement (`snappyHexMeshDict`): wing level (4 5); three nested boxes: box1 x[-3,25] y±4 level 3, box2 x[-2,10] y±3 level 3, box3 x[-1,5] y±2 level 4; 3 boundary layers.
- Mesh pipeline: `blockMesh → snappyHexMesh -overwrite → extrudeMesh → createPatch -overwrite → checkMesh → pimpleFoam`. About 11.9k–13.4k cells per case in the upstream study.
- `fvSchemes`: `div(phi,U) Gauss linear` (central, unbounded). Writing: originally `writeInterval 0.1`, `purgeWrite 10` (only the last 10 time folders survive).
- Wing geometry `constant/geometry/cylinder_2d.obj`: 41 points of a 2D profile extruded in z (82 vertices). Chord 1, y in [0, 0.0847]. **It is a thin cambered profile lying entirely above y = 0, not a symmetric NACA0012.** Area 0.049837, **area centroid (center of mass) (0.41735, 0.04068)**. The vertex average (0.408, 0.034) is not the centroid.

## 4. Files delivered (all in this zip, same folder layout needed to run)

| File | Purpose |
|---|---|
| `run_domain_sweep.sh` | **Unified sweep driver** for all four studies (details in section 5). |
| `run_mesh_sweep.sh` | 3-line wrapper that runs `SWEEP=mesh`. |
| `blockMeshDict.template` | Domain box with placeholders `@XMIN@ @XMAX@ @YMIN@ @YMAX@ @NX@ @NY@`. |
| `snappyHexMeshDict.template` | Placeholders for wing/box refinement levels, `@LOCX@ @LOCY@` (`locationInMesh`) and `@BOX1_XMIN@`. |
| `rotate_wing.py` | Rotates the wing about a pivot; `--find-centroid` computes the area centroid; `--plot` draws a check. Positive AoA lifts the leading edge. |
| `post_process.py` | Reads `forceCoeffs.dat`, groups by AoA, reports Cd/Cl/St/state and convergence, makes plots and a CSV. |
| `CASE_AND_SCRIPTS_DOCUMENTATION.md` | Long reference, kept in step with the scripts: quick start; every case file (Part 1); the two templates (2); the sweep script step by step (3); `post_process.py` (4); the AoA method and `rotate_wing.py` (5); saving time steps and ParaView (6); order of studies and judging convergence (7); troubleshooting (8); reference tables and the change history of the scripts (appendices). |
| `aoa15_sign_check.png` | Picture confirming the rotation direction and pivot. |

To use: put these files next to a **clean copy** of the case named `airfoil2d_transient/` (no `constant/polyMesh`, no results). The case itself is not in the zip (it is ~7 MB; the user has it).

## 5. How the sweep script works

```bash
AOA="0 15 25" SWEEP=upstream   bash run_domain_sweep.sh                          # XMIN = -3,-5,-7,-10,-15
AOA="0 15 25" SWEEP=downstream UP_FIXED=<xmin> bash run_domain_sweep.sh           # XMAX = 15,20,27,35,45
AOA="0 15 25" SWEEP=lateral    UP_FIXED=<xmin> DOWN_FIXED=<xmax> bash run_domain_sweep.sh   # half-height 4,6,8,12,16
AOA="0 15 25" SWEEP=mesh       UP_FIXED=<xmin> DOWN_FIXED=<xmax> HALF_H=<h> bash run_domain_sweep.sh
```
- Cases are folders `AoA<angle>_<tag>` (tags `U1_x-3`, `D2_xd20`, `L3_h8`, `M2_medium`) under `sweep_domain/`, `sweep_downstream/`, `sweep_lateral/` or `sweep_mesh/`. Each study refuses to start without the fixed values from the earlier studies.
- Switches: `DRY=1` (list cases only), `MESH_ONLY=1` (build and check meshes, no solve), `FORCE=1` (redo finished cases), `WORKDIR=...`. Angles must be integers.
- Behaviour built in: far-field cell size kept constant as the domain changes (`NX`, `NY` scale with size); wing rotated about its own centroid; each step logs stdout **and** stderr to `log.<cmd>`; a case stops at its first failing step and failures are listed at the end; **finished cases are skipped** (finished = `log.pimpleFoam` ends with `End` and `forceCoeffs.dat` exists); `run.sh` is deleted from sweep copies; a `mesh_backup/` copy of the mesh is kept per case; fields are written every 0.5 time units and, after a case finishes, time folders outside **t = 60–120** are deleted (0 and the final time kept). Settings sit at the top of the script.
- Long runs: `AOA="0 15 25" SWEEP=... nohup bash run_domain_sweep.sh > sweep.log 2>&1 &` then `tail -f sweep.log`.
- Post-processing: `python3 post_process.py sweep_domain/AoA* --avg-fraction 0.33 --out upstream_sweep` (averages the last third of the run, about t = 100–150). It prints a table per angle and an "Independence check" comparing the last two cases against a 1% tolerance (also read the `d%max` column).

## 6. What went wrong along the way, and what was done (lessons)

1. **Script syntax error on WSL** (`unexpected token '>'`): fragile backslash line continuations, probably CRLF line endings. Fixed by putting `sed`/`find` commands on single lines.
2. **`locationInMesh` hard-coded to (-5 0 0)** put the seed point outside the mesh for smaller domains (snappyHexMesh fatal error). It is now `(XMIN+1, 0, 0)`, templated per case.
3. **Refinement box1 was fixed at x = -3 while the inlet moved.** It touched the inlet for the smallest domain and sat far from it for the largest. The old `U3_x-7` case showed huge irregular Cl oscillations from t ≈ 100 plus continuity errors; its checkMesh showed higher non-orthogonality (44.8 vs 33.3) and skewness (1.19 vs 0.91) than `U1`. The fix: `refinementBox1` now starts exactly at `XMIN`. This is consistent with the fix working: in the later runs **no case oscillates** (Cl_rms ~ 0 everywhere). But the cause was never isolated by an experiment, so treat it as the most likely explanation, not a proven one. (An intermediate "2-unit buffer" version had no rigorous basis and was dropped.)
4. **A wrong early diagnosis:** the first reading of the flat `U1`/`U2` and wild `U3` signals was "shedding onset takes long to grow". That was wrong; the flat cases are genuinely steady and `U3` was a numerical artifact.
5. **Empty logs / missing results after a run:** only stdout was logged, so failures left empty logs. Now stderr is logged too and success requires `forceCoeffs.dat`.
6. **Mesh disappeared from a finished case** (`constant/polyMesh` held only `cellLevel`, `level0Edge`, `pointLevel`, `surfaceIndex`; ParaView said "contains no meshes"). **The cause was never found.** Likely routes (unconfirmed): re-running the script over a finished case, or the case's own `run.sh` (starts with `foamCleanCase`). Mitigations added: skip finished cases, delete `run.sh` from sweep copies, keep `mesh_backup/`, warn if `points` is missing after meshing. The assistant advised copying the case, rebuilding the mesh (`blockMesh; snappyHexMesh -overwrite; extrudeMesh; createPatch -overwrite`) and checking cells/points against the old log (11858 cells, 24394 points for `AoA0_U1_x-3`) before trusting any picture. The user then reported the issue resolved (the exact steps they took were not shown).
7. **Only the last 10 time folders existed** (`purgeWrite 10`). Now every 0.5 time units, pruned to 60–120 (+ 0 + final). **Note:** the AoA 0–25 upstream runs finished *before* this and have only the last 10 time folders; the skip guard will not redo them, so use `FORCE=1` if ParaView animations of them are wanted.
8. **Post-processing false alarms** (St = 40 from noise, "still transient" on flat signals, % changes computed across different angles). `post_process.py` was rewritten: groups by AoA, shows St only if the flow oscillates, classifies steady / DRIFTING / oscillating, adds Cl_mean and the independence check.

## 7. Results so far (upstream study, AoA 0–25, x = -3 / -5 / -7; XMAX 27, y ±8, endTime 150)

All 18 cases were **steady (Cl_rms < 5e-6)**. Mean Cd read from the user's table (averaged over the last half of the run):

| AoA | Cd x=-3 | Cd x=-5 | Cd x=-7 | change -5 → -7 |
|---|---|---|---|---|
| 0 | 0.40079 | 0.39306 | 0.37232 | -5.3% |
| 5 | 0.41117 | 0.40744 | 0.39392 | -3.3% |
| 10 | 0.44952 | 0.44449 | 0.39885 | -10.3% |
| 15 | 0.51514 | 0.50686 | 0.46918 | -7.4% |
| 20 | 0.59442 | 0.52733 | 0.54984 | +4.3% (not monotonic) |
| 25 | 0.67854 | 0.65784 | 0.60280 | -8.4% |

Cell counts about 11.85k / 12.75k / 13.4k. **Reading:** Cd keeps falling as the inlet moves upstream and is not converged at x = -7. AoA 20 is non-monotonic and its meshes deserve a look in ParaView. `Cl_mean` was not yet seen. A rough channel argument (upstream disturbance fading like exp(-πx/H) with H = 16, so about 55%, 29% and 14% left at 3, 5 and 7 chords) suggests 10–15 chords may be needed; this is an estimate, not a result.
- The script's upstream presets now also include `U4_x-10` and `U5_x-15`, which the assistant recommended. The user later said "we have the upstream data", but **whether those larger cases were run, their results, and the chosen converged upstream distance were not reported in the chat: UP_FIXED = [FILL IN]**.

## 8. Next steps (the user's plan)

1. Confirm the upstream study (check `d%max` and the verdict for each angle) and choose `UP_FIXED` = the largest distance needed across angles.
2. **Downstream study:** `SWEEP=downstream UP_FIXED=<value>`. Suggest 3 representative angles (for example 0, 15, 25) to save time. Caveat: `refinementBox1` ends at x = 25, so for XMAX = 35 or 45 the extra wake lies on the coarse background mesh (cells about 1.3 chords); for XMAX = 15 or 20 the box simply extends past the domain. Consider whether to extend the box with XMAX.
3. **Lateral (half-height) study** is not in the user's stated plan but is recommended: the slip walls at ±8 act like a channel and blockage grows with AoA. It is built in (`SWEEP=lateral`).
4. **Mesh study:** `SWEEP=mesh UP_FIXED=... DOWN_FIXED=... HALF_H=...` (levels coarse to xfine; wing/box levels 3 4 2 2 3, 4 5 3 3 4, 5 6 4 4 5, 6 7 5 5 6). Note the refinement levels were the same in all domain studies (the baseline 4 5 3 3 4).
5. Before long runs: `DRY=1`, then `MESH_ONLY=1` for the chosen angles, and look at one rotated case in ParaView.

## 9. Open questions / things not verified

- Cause of the mesh files disappearing (section 6, item 6).
- Whether the old `U3` instability was really due to the box1 mismatch (section 6, item 3).
- `Cl_mean` and the extended upstream (x = -10, -15) results: not seen by the assistant.
- Any use of an independence tolerance other than 1% (the user has not stated one).
- The scripts are untested against a real OpenFOAM in the assistant's environment (only the user's runs exercised them).

## 10. Working preferences observed

The user wants each change explained, prefers to finish one study before starting the next, works in WSL on a Windows drive, and likes scripts that are safe to re-run. Give exact commands to paste. Keep fixes small and test them (with stubs) before handing them over.

## 11. ParaView quick recipe

`touch case.foam` in the case folder; open it; tick `internalMesh` and the patches (not `Front`/`Back`); Apply; choose the last time; view along +Z. Vorticity: Gradient filter on `U` with "Compute Vorticity". If it says "contains no meshes", `constant/polyMesh` is missing its `points/faces/owner/neighbour/boundary` files: restore from `mesh_backup/` (`cp -r mesh_backup/. constant/polyMesh/`) or rebuild the mesh.
