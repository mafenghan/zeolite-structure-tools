# -*- coding: utf-8 -*-
"""ring_tables_v2 -- strong window-ring table builder.

The original *ring_tables_v2* part of the ring-table code: symmetry
reduction (spglib), the per-framework pipeline, and the CSV driver.  Window
search/classification lives in ``find_windows``; the supercell gap-fill lives
in ``ring_v2_supercell_fix``; the IZA essential-ring table ships as
``iza_essential_rings.json``.

Usage:
  python ring_tables_v2.py SRC_DIR [-o OUT.csv] [--codes A B C ...]
         [--iza iza.json] [--supercell-root DIR] [--no-fill]

  SRC_DIR           directory with <CODE>/<CODE>.arc subfolders
  -o OUT.csv        output table (default: zeolite_ring_tables.csv)
  --codes           restrict to these framework codes
  --iza             external IZA essential-rings JSON (default: bundled data)
  --supercell-root  reuse existing <CODE>/framework_orig.arc supercells;
                    without it a 2x2x2 supercell is generated on the fly
  --no-fill         skip the supercell gap-fill pass

Output rows: one per inequivalent T site (orbit_size = site multiplicity)
and one per inequivalent ring kind (orbit_size = number of that kind's ring
copies per unit cell), with ring centre, plane normal and O-plane aperture.

Ring recipe (output kind counts follow the IZA essential-ring counts by
construction: candidate kinds are trimmed/padded to those counts, largest
O-plane aperture first):
  * simple T-cycles of the framework's IZA essential-ring sizes (>= 6; odd
    sizes 7/9/11 kept when IZA lists them),
  * closed     (first/last ring atoms bonded in the unwrapped coordinates;
    rejects winding cycles with a net lattice translation),
  * chordless  (chord judged on the ring's unwrapped coordinates),
  * empty-disk (no T projection inside the ring's planar polygon),
  * symmetry-reduced (spglib) to inequivalent kinds, trimmed per size to the
    IZA essential-ring count keeping the largest O-plane apertures,
  * an automatic supercell pass fills kinds that only exist as multi-edge /
    wrap-around rings of the primitive cell (serials still refer to the
    standard file).

Positive control: FAU 12-ring aperture 7.33 A vs IZA free diameter 7.35 A.
"""
import argparse
import csv
import json
import os
import time

import numpy as np

from constants import FALLBACK_SIZES, MIN_SIZE, T_ELEMS
from find_windows import RingFinder, classify, lattice, mic, parse_arc
from ring_v2_supercell_fix import fill_from_supercell

BUNDLED_IZA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "iza_essential_rings.json")


# -------------------------------------------------------------------- spglib
def spglib_data(xyz, L, types):
    import spglib
    frac = np.linalg.solve(L.T, xyz.T).T
    numbers = [{"Si": 14, "Al": 13, "P": 15}.get(t, 6) for t in types]
    ds = spglib.get_symmetry_dataset((L, frac, numbers), symprec=0.5)
    if ds is None:
        return (np.eye(3, dtype=int)[None, ...], np.zeros((1, 3))), \
               np.arange(len(xyz)), ["?"] * len(xyz)
    return ((np.asarray(ds.rotations), np.asarray(ds.translations)),
            np.asarray(ds.equivalent_atoms),
            list(ds.wyckoffs) if hasattr(ds, "wyckoffs") else [])


def build_perms(frac_all, rots, trans):
    lookup = {}
    for i, f in enumerate(frac_all):
        lookup[tuple(np.round(f - np.floor(f + 1e-8), 3))] = i
    perms = []
    for R, t in zip(rots, trans):
        perm = []
        for i, f in enumerate(frac_all):
            g = R @ f + t
            g = g - np.floor(g + 1e-8)
            j = lookup.get(tuple(np.round(g, 3)))
            if j is None:
                d = frac_all - g
                d -= np.round(d)
                j = int(np.argmin(np.abs(d).max(axis=1)))
            perm.append(j)
        perms.append(perm)
    return perms


# ----------------------------------------------------------------- one framework
def process(code, src_dir, super_root, iza, do_fill=True):
    t0 = time.time()
    rec = iza.get(code)
    if rec:
        counts = {int(k): v for k, v in rec.get("counts", {}).items()}
        sizes = [s for s in rec["sizes"] if MIN_SIZE <= s <= rec["max"]]
    else:
        counts = {}
        sizes = list(FALLBACK_SIZES)

    pbc, atoms = parse_arc(os.path.join(src_dir, code, code + ".arc"))
    L = lattice(pbc)
    types = [a[1] for a in atoms]
    xyz_all = np.array([a[2:] for a in atoms])
    t_idx = [i for i, t in enumerate(types) if t in T_ELEMS]
    o_idx = [i for i, t in enumerate(types) if t == "O"]
    T_ids = [atoms[i][0] for i in t_idx]
    T_types = [types[i] for i in t_idx]
    T_xyz = xyz_all[t_idx]
    nT = len(t_idx)
    D = np.linalg.norm(mic(L, T_xyz[:, None, :] - T_xyz[None, :, :]), axis=2)
    adj = [list(np.where((D[i] < 3.4) & (D[i] > 0.1))[0]) for i in range(nT)]

    (rots, trans), eq, wy = spglib_data(T_xyz, L, T_types)
    frac_all = np.linalg.solve(L.T, T_xyz.T).T
    perms = build_perms(frac_all, rots, trans)

    t_sites = {}
    for i, s in enumerate(eq):
        t_sites.setdefault(int(s), []).append(i)     # T-local indices
    t_rows = []
    for j, (_c, loc) in enumerate(sorted(t_sites.items()), 1):
        w = wy[loc[0]] if loc and 0 <= loc[0] < len(wy) else "?"
        t_rows.append((j, len(loc), str(w),
                       sorted(T_ids[k] for k in loc)))

    if not sizes:
        return t_rows, [], 0, round(time.time() - t0, 1)

    lo, hi = min(sizes), max(sizes)
    finder = RingFinder(T_xyz, adj, L, hi)
    finder.run(lo)
    kinds = {}
    inst = {}                 # kind key -> set of per-cell ring centres
    for path, coords in finder.cycles:
        sz = len(path)
        if sz not in sizes:
            continue
        got = classify(path, coords, T_xyz, o_idx, xyz_all, L)
        if got is None:
            continue
        center, nrm, aperture = got
        key = min(tuple(sorted(p[i] for i in path)) for p in perms)
        ckey = tuple(np.round(np.linalg.solve(L.T, center) % 1.0, 3))
        inst.setdefault(key, set()).add(ckey)
        cur = kinds.get(key)
        if cur is None or aperture > cur[3]:
            kinds[key] = (list(path), center, nrm, aperture)

    rows = []
    by = {}
    for key, val in kinds.items():
        by.setdefault(len(val[0]), []).append((key, val))
    for sz in sorted(by):
        group = sorted(by[sz], key=lambda kv: -kv[1][3])
        want = counts.get(sz, min(len(group), 8))
        for key, (path, center, nrm, aperture) in group[:want]:
            rows.append((sz, len(inst.get(key, ())),
                         " ".join(str(T_ids[i]) for i in path),
                         center, nrm, aperture))

    # -------- supercell gap-fill for kind shortfalls --------
    have = {}
    for sz, *_ in rows:
        have[sz] = have.get(sz, 0) + 1
    short = {s: counts.get(s, 8) - have.get(s, 0) for s in sizes}
    short = {s: v for s, v in short.items() if v > 0}
    if do_fill and short:
        rows = fill_from_supercell(code, src_dir, super_root, rows, short,
                                  counts)
    return t_rows, rows, len(finder.cycles), round(time.time() - t0, 1)


def load_iza(path=None):
    """IZA essential-ring table: explicit --iza file or the bundled JSON."""
    p = path or BUNDLED_IZA
    if os.path.isfile(p):
        with open(p, encoding="utf-8") as f:
            return json.load(f)
    return {}


# ----------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description="strong window-ring table")
    ap.add_argument("src", help="directory with <CODE>/<CODE>.arc subfolders")
    ap.add_argument("-o", "--out", default="zeolite_ring_tables.csv")
    ap.add_argument("--codes", nargs="*", default=None)
    ap.add_argument("--iza", default=None, help="IZA essential-rings JSON")
    ap.add_argument("--supercell-root", default=None,
                    help="dir with <CODE>/framework_orig.arc (else 2x2x2 generated)")
    ap.add_argument("--no-fill", action="store_true",
                    help="skip the supercell gap-fill pass")
    a = ap.parse_args()

    iza = load_iza(a.iza)

    codes = sorted(d for d in os.listdir(a.src)
                   if os.path.isfile(os.path.join(a.src, d, d + ".arc")))
    if a.codes:
        codes = [c for c in codes if c in set(a.codes)]

    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("code,record_type,site_id,orbit_size,wyckoff,t_serials,"
                 "cx,cy,cz,nx,ny,nz,aperture\n")
    n_t = n_r = 0
    for code in codes:
        try:
            t_rows, rows, n_raw, secs = process(
                code, a.src, a.supercell_root, iza, do_fill=not a.no_fill)
        except Exception as e:
            print(f"{code}: ERROR {type(e).__name__}: {e}")
            continue
        with open(a.out, "a", encoding="utf-8", newline="\n") as fh:
            w = csv.writer(fh)
            for j, orbit, wy, members in t_rows:
                w.writerow([code, "T", j, orbit, wy,
                            " ".join(map(str, members)),
                            "", "", "", "", "", "", ""])
            for k, (sz, orbit, sers, c, nz, ap) in enumerate(rows, 1):
                w.writerow([code, "RING", k, orbit, "", sers,
                            round(float(c[0]), 6), round(float(c[1]), 6),
                            round(float(c[2]), 6),
                            round(float(nz[0]), 6), round(float(nz[1]), 6),
                            round(float(nz[2]), 6), ap])
        n_t += len(t_rows)
        n_r += len(rows)
        print(f"{code}: T={len(t_rows)} rings={len(rows)} "
              f"sizes={sorted({r[0] for r in rows})} raw={n_raw} {secs}s")
    print(f"DONE codes={len(codes)} T_rows={n_t} RING_rows={n_r} -> {a.out}")


if __name__ == "__main__":
    main()
