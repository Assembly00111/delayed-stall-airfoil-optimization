#!/usr/bin/env python3
"""
post_process.py
================
Post-processes OpenFOAM sweep cases named  AoA<angle>_<tag>_x-<upstream>
(e.g. AoA10_U2_x-5) and reports, for each AoA separately, how the force
coefficients change as the domain grows. Cases whose names don't match that
pattern are treated as one plain sweep, in the order given.

Per case:
  Cd_mean, Cl_mean : time averages over the last --avg-fraction of the run
  Cl_rms           : RMS fluctuation of Cl in that window
  St               : Strouhal number, ONLY if the flow actually oscillates
                     (Cl_rms > --osc-threshold); otherwise shown as '-'
  state            : steady      Cl barely fluctuates and the means do not drift
                     DRIFTING    steady-looking but the means are still moving
                                 (run longer before trusting the numbers)
                     osc         periodic shedding, amplitude settled
                     osc-GROWING periodic shedding, amplitude still changing

Convergence columns (within one AoA, domains ordered by upstream distance):
  d%prev : change relative to the next-smaller domain
  d%max  : difference from the LARGEST domain (your best estimate)

USAGE
    python3 post_process.py sweep_domain/AoA* --avg-fraction 0.33 --out upstream_sweep

OUTPUT   <out>.csv   <out>.png (convergence vs upstream distance)
         <out>_timeseries.png (Cd(t), Cl(t), one row per AoA)
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

NAME_RE = re.compile(r"^AoA(-?\d+(?:\.\d+)?)_(.+?)_x-(\d+(?:\.\d+)?)$")


def read_force_coeffs(path):
    rows, meta = [], {}
    hdr = re.compile(r"#\s*(magUInf|lRef|Aref)\s*:\s*([-+0-9.eE]+)")
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.lstrip().startswith("#"):
                m = hdr.search(line)
                if m:
                    meta[m.group(1)] = float(m.group(2))
                continue
            try:
                rows.append([float(x) for x in line.split()])
            except ValueError:
                continue
    if not rows:
        raise ValueError("no data rows in " + path)
    return np.array(rows), meta


def read_cell_count(path):
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        for line in f:
            m = re.search(r"\bcells:\s*(\d+)", line)
            if m:
                return int(m.group(1))
    return None


def rms(x):
    return float(np.sqrt(np.mean((x - x.mean()) ** 2)))


def analyse_case(case_dir, avg_fraction, osc_thr, drift_tol):
    fc = os.path.join(case_dir, "postProcessing", "forces", "0", "forceCoeffs.dat")
    if not os.path.isfile(fc):
        raise FileNotFoundError(fc)
    data, meta = read_force_coeffs(fc)
    t, Cd, Cl = data[:, 0], data[:, 2], data[:, 3]
    lRef, U = meta.get("lRef", 1.0), meta.get("magUInf", 1.0)

    t_start = t[-1] * (1.0 - avg_fraction)
    mask = t >= t_start
    if mask.sum() < 16:
        mask = np.ones_like(t, dtype=bool)
        t_start = t[0]
    tm, Cdm, Clm = t[mask], Cd[mask], Cl[mask]
    n = len(tm)
    h = n // 2

    Cd_mean, Cl_mean, Cl_rms = float(Cdm.mean()), float(Clm.mean()), rms(Clm)
    cd_drift = 100 * abs(Cdm[h:].mean() - Cdm[:h].mean()) / max(abs(Cd_mean), 1e-12)
    cl_drift = 100 * abs(Clm[h:].mean() - Clm[:h].mean()) / max(abs(Cl_mean), 1e-2)

    St, growth = None, 1.0
    if Cl_rms > osc_thr:
        dt = float(np.median(np.diff(tm)))
        spec = np.abs(np.fft.rfft((Clm - Clm.mean()) * np.hanning(n)))
        freqs = np.fft.rfftfreq(n, d=dt)
        St = float(freqs[int(np.argmax(spec[1:]) + 1)] * lRef / U)
        r1, r2 = rms(Clm[:h]), rms(Clm[h:])
        growth = r2 / r1 if r1 > 0 else float("inf")
        state = "osc" if 0.8 <= growth <= 1.25 else "osc-GROWING"
    else:
        state = "steady" if max(cd_drift, cl_drift) < drift_tol else "DRIFTING"

    return dict(case=case_dir, N_cells=read_cell_count(os.path.join(case_dir, "log.checkMesh")),
                Cd_mean=Cd_mean, Cl_mean=Cl_mean, Cl_rms=Cl_rms, St=St, state=state,
                Cd_drift=cd_drift, Cl_drift=cl_drift, t_start=t_start, t=t, Cd=Cd, Cl=Cl)


def pct(new, ref):
    if ref is None or abs(ref) < 1e-12:
        return None
    return 100.0 * (new - ref) / abs(ref)


def fmt(v, w, p=5, suffix=""):
    if v is None:
        return f"{'-':>{w}}"
    if suffix == "%" and abs(v) > 999:
        return f"{('>999%' if v > 0 else '<-999%'):>{w}}"
    return f"{v:>{w - len(suffix)}.{p}f}{suffix}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--avg-fraction", type=float, default=0.5,
                    help="fraction of the run (from the end) to average over (default 0.5)")
    ap.add_argument("--osc-threshold", type=float, default=1e-3,
                    help="Cl_rms above this counts as oscillating/shedding (default 1e-3)")
    ap.add_argument("--drift-tol", type=float, default=0.2,
                    help="max %% drift of the means within the window to call it steady (default 0.2)")
    ap.add_argument("--tol", type=float, default=1.0,
                    help="%% change (largest vs next-largest domain) counted as converged (default 1.0)")
    ap.add_argument("--out", default="sweep_summary")
    args = ap.parse_args()

    results = []
    for c in args.cases:
        name = os.path.basename(os.path.normpath(c))
        try:
            r = analyse_case(c, args.avg_fraction, args.osc_threshold, args.drift_tol)
        except (FileNotFoundError, ValueError) as e:
            print(f"SKIPPING {name}: {e}", file=sys.stderr)
            continue
        m = NAME_RE.match(name)
        r["name"] = name
        r["aoa"] = float(m.group(1)) if m else None
        r["tag"] = m.group(2) + "_x-" + m.group(3) if m else name
        r["xup"] = float(m.group(3)) if m else None
        results.append(r)
    if not results:
        sys.exit("No usable cases found.")

    grouped = all(r["aoa"] is not None for r in results)
    if grouped:
        groups = {}
        for r in results:
            groups.setdefault(r["aoa"], []).append(r)
        for a in groups:
            groups[a].sort(key=lambda r: r["xup"])
        groups = dict(sorted(groups.items()))
    else:
        groups = {None: results}

    hdr = (f"{'case':<10}{'cells':>7}{'Cd_mean':>9}{'d%prev':>9}{'d%max':>9}"
           f"{'Cl_mean':>9}{'d%prev':>9}{'d%max':>9}{'Cl_rms':>9}{'St':>7}  state")
    csv = ["case,aoa,x_up,N_cells,Cd_mean,Cd_pct_prev,Cd_pct_vs_largest,Cl_mean,Cl_pct_prev,"
           "Cl_pct_vs_largest,Cl_rms,St,state,Cd_drift_pct,Cl_drift_pct"]
    verdicts = []
    for a, rs in groups.items():
        print()
        print(f"AoA = {a:g} deg" if a is not None else "(cases in the order given)")
        print(hdr)
        print("-" * len(hdr))
        big = rs[-1]
        for i, r in enumerate(rs):
            prev = rs[i - 1] if i > 0 else None
            r["dCd_p"] = pct(r["Cd_mean"], prev["Cd_mean"]) if prev else None
            r["dCl_p"] = pct(r["Cl_mean"], prev["Cl_mean"]) if prev else None
            r["dCd_m"] = pct(r["Cd_mean"], big["Cd_mean"]) if r is not big else None
            r["dCl_m"] = pct(r["Cl_mean"], big["Cl_mean"]) if r is not big else None
            print(f"{r['tag']:<10}{str(r['N_cells'] or 'n/a'):>7}{r['Cd_mean']:>9.5f}"
                  f"{fmt(r['dCd_p'], 9, 2, '%')}{fmt(r['dCd_m'], 9, 2, '%')}"
                  f"{r['Cl_mean']:>9.5f}{fmt(r['dCl_p'], 9, 2, '%')}{fmt(r['dCl_m'], 9, 2, '%')}"
                  f"{r['Cl_rms']:>9.5f}{fmt(r['St'], 7, 4)}  {r['state']}")
            csv.append(",".join("" if v is None else str(v) for v in [
                r["name"], r["aoa"], r["xup"], r["N_cells"], r["Cd_mean"], r["dCd_p"], r["dCd_m"],
                r["Cl_mean"], r["dCl_p"], r["dCl_m"], r["Cl_rms"], r["St"], r["state"],
                r["Cd_drift"], r["Cl_drift"]]))
        if a is not None and len(rs) >= 2:
            last = rs[-1]
            worst = max(abs(last["dCd_p"] or 0), abs(last["dCl_p"] or 0))
            ok = worst < args.tol
            verdicts.append((a, last["dCd_p"], last["dCl_p"], ok, rs[-2]["tag"], last["tag"]))

    if verdicts:
        print(f"\nDomain-independence check: largest vs next-largest domain (tolerance {args.tol:g}%)")
        for a, dcd, dcl, ok, t1, t2 in verdicts:
            print(f"  AoA {a:>4g}: {t1} -> {t2}:  Cd {fmt(dcd, 7, 2, '%').strip():>7}, "
                  f"Cl {fmt(dcl, 7, 2, '%').strip():>7}   "
                  f"{'CONVERGED' if ok else 'NOT converged -- go further upstream'}")
    bad = [r for rs in groups.values() for r in rs if r["state"] in ("DRIFTING", "osc-GROWING")]
    if bad:
        print("\nWARNING: still evolving (numbers not final): " + ", ".join(r["name"] for r in bad))

    with open(args.out + ".csv", "w") as f:
        f.write("\n".join(csv) + "\n")
    print(f"\nWrote {args.out}.csv")

    if not HAVE_MPL:
        return
    keys = list(groups.keys())
    cmap = plt.get_cmap("viridis")
    col = {a: cmap(i / max(len(keys) - 1, 1)) for i, a in enumerate(keys)}

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    for a, rs in groups.items():
        x = [r["xup"] for r in rs] if grouped else list(range(len(rs)))
        lab = f"AoA {a:g}" if a is not None else "sweep"
        ax[0].plot(x, [r["Cd_mean"] for r in rs], "o-", color=col[a], label=lab)
        ax[1].plot(x, [r["Cl_mean"] for r in rs], "o-", color=col[a], label=lab)
        big = rs[-1]
        ax[2].plot(x, [pct(r["Cd_mean"], big["Cd_mean"]) or 0 for r in rs], "o-", color=col[a], label=lab)
    ax[2].axhspan(-args.tol, args.tol, color="green", alpha=0.12)
    for axi, ttl in zip(ax, ["Mean Cd", "Mean Cl", "Cd difference from largest domain (%)"]):
        axi.set_title(ttl)
        axi.set_xlabel("upstream distance (chords)" if grouped else "case #")
        axi.grid(True, alpha=0.3)
    ax[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.out + ".png", dpi=150)
    print(f"Wrote {args.out}.png")

    fig2, axs = plt.subplots(len(keys), 2, figsize=(12, 2.3 * len(keys)), squeeze=False)
    for i, (a, rs) in enumerate(groups.items()):
        for j, (key, ttl) in enumerate([("Cd", "Cd(t)"), ("Cl", "Cl(t)")]):
            axx = axs[i, j]
            lo, hi = np.inf, -np.inf
            for k, r in enumerate(rs):
                axx.plot(r["t"], r[key], lw=0.9, label=r["tag"])
                sel = r["t"] >= 0.3 * r["t"][-1]
                lo, hi = min(lo, r[key][sel].min()), max(hi, r[key][sel].max())
            pad = max(0.1 * (hi - lo), 0.05 * max(abs(hi), abs(lo), 1e-3))
            axx.set_ylim(lo - pad, hi + pad)
            axx.axvline(rs[0]["t_start"], color="grey", ls=":")
            axx.set_title(f"{'AoA ' + format(a, 'g') + ': ' if a is not None else ''}{ttl}", fontsize=9)
            axx.grid(True, alpha=0.3)
            if j == 0:
                axx.legend(fontsize=7, ncol=3)
    fig2.supxlabel("t")
    fig2.tight_layout()
    fig2.savefig(args.out + "_timeseries.png", dpi=130)
    print(f"Wrote {args.out}_timeseries.png")


if __name__ == "__main__":
    main()
