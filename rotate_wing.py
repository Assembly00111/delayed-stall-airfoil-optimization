#!/usr/bin/env python3
"""
rotate_wing.py
==============
Rotate the wing surface (constant/geometry/cylinder_2d.obj) about a chosen
pivot point to set angle of attack, for the airfoil2d_transient case.

WHY A CUSTOM SCRIPT INSTEAD OF surfaceTransformPoints
------------------------------------------------------
OpenFOAM's `surfaceTransformPoints` rotates about the origin (0,0,0) only.
This wing's vertices are NOT centered on the origin (chord runs x=0..1,
y=0..0.085), so rotating about the origin would also translate the body
sideways as a side-effect of "rotating" it -- not what you want for an
angle-of-attack sweep. This script rotates about an arbitrary pivot
(translate -> rotate -> translate back) so the body's position stays fixed
and only its orientation changes.

SIGN CONVENTION
----------------
Free-stream flows in +x (see 0/U). The wing's leading edge is at small x,
trailing edge at x=1 (see area centroid computed at ~(0.417, 0.041)).
A POSITIVE --aoa tips the leading edge UP (clockwise rotation in the
standard x-right/y-up view) -- i.e. increasing --aoa increases the
geometric angle of attack into the oncoming +x flow in the usual sense
(nose up = positive incidence). This was verified visually (see the
--plot option) rather than assumed -- always sanity-check a new geometry
with --plot before trusting the sign.

USAGE
-----
    # find the centroid first (or supply your own pivot):
    python3 rotate_wing.py --find-centroid constant/geometry/cylinder_2d.obj

    # rotate by +10 degrees about that centroid:
    python3 rotate_wing.py constant/geometry/cylinder_2d.obj \
        constant/geometry/cylinder_2d_aoa10.obj \
        --pivot 0.41735 0.04068 --aoa 10 --plot aoa10_check.png

    # generate a full 0/10/20/30 deg family in one go:
    for a in 0 10 20 30; do
      python3 rotate_wing.py constant/geometry/cylinder_2d.obj \
          constant/geometry/cylinder_2d_aoa${a}.obj \
          --pivot 0.41735 0.04068 --aoa $a
    done
"""
import argparse
import math
import sys


def read_obj_vertices(path):
    """Return (lines, vertex_line_indices, vertices) where vertices is a
    list of [x,y,z] for each 'v ...' line, in file order."""
    lines = []
    vidx = []
    verts = []
    with open(path) as f:
        for line in f:
            lines.append(line)
            p = line.split()
            if len(p) >= 4 and p[0] == "v":
                vidx.append(len(lines) - 1)
                verts.append([float(p[1]), float(p[2]), float(p[3])])
    return lines, vidx, verts


def polygon_centroid(xy_points):
    """Area centroid of a closed 2D polygon given as ordered (x,y) points
    (does not need to repeat the first point at the end)."""
    pts = xy_points + [xy_points[0]]
    A = Cx = Cy = 0.0
    for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
        cross = x0 * y1 - x1 * y0
        A += cross
        Cx += (x0 + x1) * cross
        Cy += (y0 + y1) * cross
    A *= 0.5
    if abs(A) < 1e-12:
        raise ValueError("Degenerate polygon (zero area) -- check input points")
    return Cx / (6 * A), Cy / (6 * A), A


def find_centroid(path):
    _, _, verts = read_obj_vertices(path)
    zs = sorted(set(round(v[2], 6) for v in verts))
    if len(zs) < 2:
        raise ValueError(f"Expected two z-layers (thin extrusion), found: {zs}")
    z_pick = zs[-1]
    outline = [(v[0], v[1]) for v in verts if abs(v[2] - z_pick) < 1e-9]
    cx, cy, area = polygon_centroid(outline)
    return cx, cy, area, z_pick, len(outline)


def rotate_obj(in_path, out_path, pivot, aoa_deg, plot_path=None):
    px, py = pivot
    # positive aoa_deg = leading edge up = clockwise in x-right/y-up view
    theta = math.radians(-aoa_deg)
    c, s = math.cos(theta), math.sin(theta)

    lines, vidx, verts = read_obj_vertices(in_path)

    rotated = []
    for (x, y, z) in verts:
        dx, dy = x - px, y - py
        xr = px + dx * c - dy * s
        yr = py + dx * s + dy * c
        rotated.append((xr, yr, z))

    for line_i, (xr, yr, z) in zip(vidx, rotated):
        lines[line_i] = f"v {xr:.6f} {yr:.6f} {z:.6f}\n"

    with open(out_path, "w") as f:
        f.writelines(lines)

    print(f"Wrote {out_path}  (rotated {aoa_deg:+g} deg about pivot ({px:.5f}, {py:.5f}))")

    if plot_path:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        zs = sorted(set(round(v[2], 6) for v in verts))
        z_pick = zs[-1]
        before = [(v[0], v[1]) for v in verts if abs(v[2] - z_pick) < 1e-9]
        after = [(xr, yr) for (xr, yr, z) in rotated if abs(z - z_pick) < 1e-9]

        fig, ax = plt.subplots(figsize=(8, 4))
        bx = before + [before[0]]
        ax_ = after + [after[0]]
        ax.plot(*zip(*bx), "-o", ms=2, color="tab:blue", alpha=0.5, label="original (AoA=0)")
        ax.plot(*zip(*ax_), "-o", ms=2, color="tab:red", label=f"rotated (AoA={aoa_deg:+g} deg)")
        ax.plot(px, py, "k*", ms=14, label="pivot")
        ax.axhline(0, color="grey", lw=0.5)
        ax.annotate("free stream ->", xy=(0.02, 0.95), xycoords="axes fraction")
        ax.set_aspect("equal")
        ax.legend(loc="upper right", fontsize=8)
        ax.set_title(f"rotate_wing.py sanity check: AoA = {aoa_deg:+g} deg about ({px:.3f},{py:.3f})")
        fig.tight_layout()
        fig.savefig(plot_path, dpi=140)
        print(f"Wrote {plot_path} -- CHECK THE SIGN/DIRECTION BEFORE TRUSTING IT")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_obj", nargs="?", help="source .obj (e.g. cylinder_2d.obj)")
    ap.add_argument("output_obj", nargs="?", help="rotated .obj to write")
    ap.add_argument("--pivot", nargs=2, type=float, metavar=("PX", "PY"),
                     help="rotation pivot (x y). Use --find-centroid first if unsure.")
    ap.add_argument("--aoa", type=float, default=0.0,
                     help="angle of attack in degrees (+ = leading edge up)")
    ap.add_argument("--plot", metavar="PNG",
                     help="write a before/after sanity-check plot to this path")
    ap.add_argument("--find-centroid", metavar="OBJ",
                     help="just compute and print the area centroid of OBJ, then exit")
    args = ap.parse_args()

    if args.find_centroid:
        cx, cy, area, z_pick, n = find_centroid(args.find_centroid)
        print(f"Outline read from z={z_pick} layer, {n} points")
        print(f"Enclosed area:  {area:.6f}")
        print(f"Area centroid:  ({cx:.5f}, {cy:.5f})")
        return

    if not args.input_obj or not args.output_obj or not args.pivot:
        ap.error("input_obj, output_obj and --pivot are required unless using --find-centroid")

    rotate_obj(args.input_obj, args.output_obj, tuple(args.pivot), args.aoa, args.plot)


if __name__ == "__main__":
    main()
