#!/usr/bin/env python3
"""
post_process.py
================
Post-processes a set of OpenFOAM case directories from a domain- or
mesh-independence sweep (run_domain_sweep.sh / run_mesh_sweep.sh) and reports,
for each case:

  - N_cells        : total cell count, parsed from log.checkMesh
  - Cd_mean        : time-averaged drag coefficient over the last
                      --avg-fraction of the run (default: last 50%)
  - Cl_rms         : RMS of the lift coefficient fluctuation over the same
                      window (a measure of unsteady loading / vortex shedding
                      strength)
  - St             : Strouhal number of the dominant Cl oscillation
                      (from FFT), St = f * lRef / magUInf
  - % change       : relative change of each quantity vs. the previous case
                      in the order given on the command line -- this is the
                      convergence check. Once it drops below your tolerance
                      (commonly 1-2%) that case (or the coarser/smaller one)
                      is considered independent.

USAGE
-----
    python3 post_process.py CASE1 CASE2 CASE3 [CASE4 ...] \
        --out sweep_summary --avg-fraction 0.5 --labels "small,baseline,large,xlarge"

Each CASE must be a path to a case directory containing:
    postProcessing/forces/0/forceCoeffs.dat
    log.checkMesh      (optional, only used to report cell count)

OUTPUT
------
    <out>.csv   summary table
    <out>.png   3-panel convergence plot (Cd_mean, Cl_rms, St vs. case)

Requires: numpy, matplotlib (pandas optional, used if available for a nicer
table print but not required).
"""
import argparse
import os
import re
import sys

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False


def read_force_coeffs(path):
    """Read an OpenFOAM forceCoeffs.dat file.

    Returns (data, meta) where data is an (N,6) array of
    [Time, Cm, Cd, Cl, Clf, Clr] and meta is a dict with any of
    magUInf / lRef / Aref found in the header comments.
    """
    rows = []
    meta = {}
    header_re = re.compile(
        r"#\s*(magUInf|lRef|Aref)\s*:\s*([-+0-9.eE]+)"
    )
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.lstrip().startswith("#"):
                m = header_re.search(line)
                if m:
                    meta[m.group(1)] = float(m.group(2))
                continue
            parts = line.split()
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"No data rows found in {path}")
    return np.array(rows), meta


def read_cell_count(log_checkmesh_path):
    if not os.path.isfile(log_checkmesh_path):
        return None
    with open(log_checkmesh_path) as f:
        for line in f:
            m = re.search(r"\bcells:\s*(\d+)", line)
            if m:
                return int(m.group(1))
    return None


def analyse_case(case_dir, avg_fraction):
    fc_path = os.path.join(case_dir, "postProcessing", "forces", "0", "forceCoeffs.dat")
    if not os.path.isfile(fc_path):
        raise FileNotFoundError(f"Could not find {fc_path}")

    data, meta = read_force_coeffs(fc_path)
    t, Cd, Cl = data[:, 0], data[:, 2], data[:, 3]

    lRef = meta.get("lRef", 1.0)
    magUInf = meta.get("magUInf", 1.0)

    t_start = t[-1] * (1.0 - avg_fraction)
    mask = t >= t_start
    if mask.sum() < 8:
        mask = np.ones_like(t, dtype=bool)

    tm, Cdm, Clm = t[mask], Cd[mask], Cl[mask]

    Cd_mean = float(Cdm.mean())
    Cl_detrend = Clm - Clm.mean()
    Cl_rms = float(np.sqrt(np.mean(Cl_detrend ** 2)))

    # dominant shedding frequency via FFT of the (uniformly sampled) Cl signal
    dt = float(np.median(np.diff(tm)))
    n = len(Cl_detrend)
    St = 0.0
    if n >= 8 and Cl_rms > 1e-9:
        freqs = np.fft.rfftfreq(n, d=dt)
        spec = np.abs(np.fft.rfft(Cl_detrend * np.hanning(n)))
        if len(spec) > 1:
            idx = int(np.argmax(spec[1:]) + 1)
            f_dom = freqs[idx]
            St = f_dom * lRef / magUInf

    # Saturation check: is the oscillation amplitude still growing at the end
    # of the run, or has it plateaued? Compare Cl_rms in the first half vs.
    # second half of the averaging window. A ratio far from 1 means the
    # transient (exponential growth of the shedding instability) hasn't
    # finished yet, so Cd_mean/Cl_rms/St for this case are NOT yet meaningful
    # and shouldn't be compared to a saturated case.
    saturated = True
    growth_ratio = 1.0
    if n >= 16:
        half = n // 2
        rms1 = np.sqrt(np.mean((Clm[:half] - Clm[:half].mean()) ** 2))
        rms2 = np.sqrt(np.mean((Clm[half:] - Clm[half:].mean()) ** 2))
        if rms1 > 1e-9:
            growth_ratio = float(rms2 / rms1)
            saturated = 0.8 <= growth_ratio <= 1.25
        elif rms2 > 1e-9:
            saturated = False  # started oscillating only in the 2nd half -> still growing

    n_cells = read_cell_count(os.path.join(case_dir, "log.checkMesh"))

    return {
        "case": case_dir,
        "N_cells": n_cells,
        "t_window_start": round(t_start, 3),
        "Cd_mean": Cd_mean,
        "Cl_rms": Cl_rms,
        "St": St,
        "saturated": saturated,
        "growth_ratio": growth_ratio,
        "t_full": t,
        "Cd_full": Cd,
        "Cl_full": Cl,
    }


def pct_change(prev, curr):
    if prev is None or prev == 0:
        return None
    return 100.0 * (curr - prev) / abs(prev)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cases", nargs="+", help="Case directories, in sweep order "
                                              "(e.g. small -> large, or coarse -> fine)")
    ap.add_argument("--avg-fraction", type=float, default=0.5,
                    help="Fraction of the run (from the end) used for time-averaging "
                         "(default: 0.5, i.e. last half of the signal)")
    ap.add_argument("--labels", type=str, default=None,
                    help="Comma-separated display labels, one per case, for the plot "
                         "x-axis (defaults to the case directory basenames)")
    ap.add_argument("--out", type=str, default="sweep_summary",
                    help="Output file basename for the .csv and .png (default: sweep_summary)")
    args = ap.parse_args()

    labels = (args.labels.split(",") if args.labels
              else [os.path.basename(os.path.normpath(c)) for c in args.cases])
    if len(labels) != len(args.cases):
        sys.exit("Number of --labels must match number of case directories")

    results = [analyse_case(c, args.avg_fraction) for c in args.cases]

    # ---- print / write table -------------------------------------------
    header = f"{'case':<22}{'N_cells':>10}{'Cd_mean':>12}{'%chg':>8}{'Cl_rms':>12}{'%chg':>8}{'St':>10}{'%chg':>8}  {'saturated?'}"
    print(header)
    print("-" * len(header))

    prev_cd = prev_clrms = prev_st = None
    csv_lines = ["case,label,N_cells,Cd_mean,Cd_pct_change,Cl_rms,Cl_rms_pct_change,St,St_pct_change,saturated,growth_ratio_2nd_vs_1st_half"]
    any_unsaturated = False
    for lab, r in zip(labels, results):
        dcd = pct_change(prev_cd, r["Cd_mean"])
        dcl = pct_change(prev_clrms, r["Cl_rms"])
        dst = pct_change(prev_st, r["St"])
        ncell_str = str(r["N_cells"]) if r["N_cells"] is not None else "n/a"
        sat_str = "yes" if r["saturated"] else "NO -- still transient!"
        if not r["saturated"]:
            any_unsaturated = True
        print(f"{lab:<22}{ncell_str:>10}{r['Cd_mean']:>12.5f}"
              f"{'' if dcd is None else f'{dcd:>7.2f}%'}"
              f"{r['Cl_rms']:>12.5f}"
              f"{'' if dcl is None else f'{dcl:>7.2f}%'}"
              f"{r['St']:>10.5f}"
              f"{'' if dst is None else f'{dst:>7.2f}%'}  {sat_str}")
        csv_lines.append(",".join(str(x) for x in [
            r["case"], lab, r["N_cells"] if r["N_cells"] is not None else "",
            r["Cd_mean"], "" if dcd is None else dcd,
            r["Cl_rms"], "" if dcl is None else dcl,
            r["St"], "" if dst is None else dst,
            r["saturated"], round(r["growth_ratio"], 3),
        ]))
        prev_cd, prev_clrms, prev_st = r["Cd_mean"], r["Cl_rms"], r["St"]

    if any_unsaturated:
        print("\n*** WARNING: at least one case's C_l oscillation had not saturated to a")
        print("    steady periodic amplitude by the end of its run (2nd-half RMS differs")
        print("    from 1st-half RMS by more than ~20-25%). Its Cd_mean/Cl_rms/St are NOT")
        print("    directly comparable to a saturated case -- extend endTime for that case")
        print("    (or all cases) before drawing any domain/mesh-independence conclusion.")
        print("    See the *_timeseries.png plot to check onset/saturation by eye.\n")

    csv_path = args.out + ".csv"
    with open(csv_path, "w") as f:
        f.write("\n".join(csv_lines) + "\n")
    print(f"\nWrote {csv_path}")

    # ---- plot -------------------------------------------------------------
    if HAVE_MPL:
        x = np.arange(len(results))
        cd = [r["Cd_mean"] for r in results]
        clrms = [r["Cl_rms"] for r in results]
        st = [r["St"] for r in results]

        fig, axes = plt.subplots(1, 3, figsize=(13, 4))
        for ax, y, title in zip(axes, [cd, clrms, st],
                                 ["Mean $C_d$", "$C_l$ RMS", "Strouhal number $St$"]):
            ax.plot(x, y, "o-", color="tab:blue")
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=30, ha="right")
            ax.set_title(title)
            ax.grid(True, alpha=0.3)
        fig.suptitle("Convergence check: quantity should plateau across the sweep")
        fig.tight_layout()
        png_path = args.out + ".png"
        fig.savefig(png_path, dpi=150)
        print(f"Wrote {png_path}")

        # Raw time-history plot: one row per case, Cl(t) and Cd(t), with the
        # averaging window shaded. This is what actually tells you whether
        # shedding has saturated or whether the run just ended mid-transient.
        nrows = len(results)
        fig2, axes2 = plt.subplots(nrows, 1, figsize=(9, 2.4 * nrows), sharex=False, squeeze=False)
        for i, (lab, r) in enumerate(zip(labels, results)):
            ax = axes2[i, 0]
            ax.plot(r["t_full"], r["Cl_full"], lw=0.8, color="tab:red", label="$C_l$")
            ax2 = ax.twinx()
            ax2.plot(r["t_full"], r["Cd_full"], lw=0.8, color="tab:blue", alpha=0.6, label="$C_d$")
            ax.axvspan(r["t_window_start"], r["t_full"][-1], color="grey", alpha=0.15)
            tag = "SATURATED" if r["saturated"] else "STILL TRANSIENT"
            ax.set_title(f"{lab}  ({tag}, 2nd/1st-half RMS ratio = {r['growth_ratio']:.2f})",
                         fontsize=10, color=("black" if r["saturated"] else "tab:red"))
            ax.set_ylabel("$C_l$", color="tab:red")
            ax2.set_ylabel("$C_d$", color="tab:blue")
            ax.set_xlabel("t")
            ax.grid(True, alpha=0.3)
        fig2.suptitle("Raw force histories (shaded region = averaging window used above)")
        fig2.tight_layout()
        ts_path = args.out + "_timeseries.png"
        fig2.savefig(ts_path, dpi=150)
        print(f"Wrote {ts_path}")
    else:
        print("matplotlib not available -- skipped plot, CSV written only.")


if __name__ == "__main__":
    main()
