# -*- coding: utf-8 -*-
"""ring_v2_supercell_fix -- supercell gap-fill pass.

The original *ring_v2_supercell_fix* part of the ring-table code: build /
load a supercell, canonicalise ring representations, and fill ring kinds that
only exist as multi-edge / wrap-around rings of the primitive cell.
"""
import os

import numpy as np

from constants import T_ELEMS
from find_windows import RingFinder, classify, lattice, mic, parse_arc


def canon_key(ring):
    """ring = [(serial, offset)] canonical over rotation/reversal/translation."""
    n = len(ring)
    best = None
    for seq in (ring, ring[::-1]):
        for k in range(n):
            rot = seq[k:] + seq[:k]
            base = np.array(rot[0][1])
            key = tuple((s, tuple(int(x) for x in (np.array(o) - base)))
                        for s, o in rot)
            if best is None or key < best:
                best = key
    return best


# --------------------------------------------------------------- supercell
def make_supercell(atoms, pbc, repeats=(2, 2, 2)):
    """Replicate the cell; returns (atoms2, pbc2, mapping) with mapping[i] =
    (std_serial, cell_offset) for every atom in the supercell list order."""
    L = lattice(pbc)
    atoms2, mapping = [], []
    a, b, c = repeats
    for i in range(a):
        for j in range(b):
            for k in range(c):
                off = np.array([i, j, k], float) @ L
                for sid, el, x, y, z in atoms:
                    atoms2.append((sid, el, x + off[0], y + off[1], z + off[2]))
                    mapping.append((sid, (i, j, k)))
    pbc2 = [pbc[0] * a, pbc[1] * b, pbc[2] * c, pbc[3], pbc[4], pbc[5]]
    return atoms2, pbc2, mapping


def load_supercell(src_dir, code, super_root):
    """Prefer an existing framework_orig.arc; else generate 2x2x2.

    `mapping` has one entry per ATOM of the supercell (None for non-T), so
    callers can index it with all-atom indices consistently.
    """
    prim = parse_arc(os.path.join(src_dir, code, code + ".arc"))
    if super_root:
        f = os.path.join(super_root, code, "framework_orig.arc")
        if os.path.isfile(f):
            sup = parse_arc(f)
            Lp, Ls = lattice(prim[0]), lattice(sup[0])
            xyz_p = np.array([a[2:] for a in prim[1]])
            t_p = [i for i, a in enumerate(prim[1]) if a[1] in T_ELEMS]
            frac_p = np.linalg.solve(Lp.T, xyz_p[t_p].T).T
            ids_p = [prim[1][i][0] for i in t_p]
            n_cells = round(abs(np.linalg.det(Ls)) / abs(np.linalg.det(Lp)))
            mapping = []
            for a in sup[1]:
                if a[1] not in T_ELEMS:
                    mapping.append(None)   # keep all-atom indexing aligned
                    continue
                f3 = np.linalg.solve(Lp.T, np.array(a[2:]))
                off = np.floor(f3 + 1e-6).astype(int)
                g = f3 - off
                d = frac_p - g
                d -= np.round(d)
                j = int(np.argmin(np.abs(d).max(axis=1)))
                mapping.append((ids_p[j], tuple(int(x) for x in off)))
            return sup[1], sup[0], mapping, n_cells
    atoms2, pbc2, mapping = make_supercell(prim[1], prim[0])
    return atoms2, pbc2, mapping, 8


def fill_from_supercell(code, src_dir, super_root, rows, short, counts):
    """Fill ring kinds missing from the primitive-cell pass using a supercell."""
    atoms2, pbc2, mapping, n_cells = load_supercell(src_dir, code, super_root)
    L = lattice(pbc2)
    types = [a[1] for a in atoms2]
    xyz_all = np.array([a[2:] for a in atoms2])
    t_idx = [i for i, t in enumerate(types) if t in T_ELEMS]
    o_idx = [i for i, t in enumerate(types) if t == "O"]
    T_xyz = xyz_all[t_idx]
    nT = len(t_idx)
    D = np.linalg.norm(mic(L, T_xyz[:, None, :] - T_xyz[None, :, :]), axis=2)
    adj = [list(np.where((D[i] < 3.4) & (D[i] > 0.1))[0]) for i in range(nT)]
    missing = sorted(short)
    finder = RingFinder(T_xyz, adj, L, max(missing))
    finder.run(min(missing))
    Lp = lattice(parse_arc(os.path.join(src_dir, code, code + ".arc"))[0])
    kinds, inst = {}, {}
    for path, coords in finder.cycles:
        sz = len(path)
        if sz not in short:
            continue
        got = classify(path, coords, T_xyz, o_idx, xyz_all, L)
        if got is None:
            continue
        center, nrm, aperture = got
        ring = [mapping[t_idx[i]] for i in path]
        key = canon_key(ring)
        # per-prim-cell centres: the supercell holds n_cells copies
        ckey = tuple(np.round(np.linalg.solve(Lp.T, center) % 1.0, 3))
        inst.setdefault(key, set()).add(ckey)
        cur = kinds.get(key)
        if cur is None or aperture > cur[3]:
            kinds[key] = (ring, center, nrm, aperture)

    merged = [r for r in rows if r[0] not in short]
    by = {}
    for key, val in kinds.items():
        by.setdefault(len(val[0]), []).append((key, val))
    for sz in missing:
        want = counts.get(sz, 8)
        pool = {}
        for r in rows:
            if r[0] == sz:
                pool[tuple(sorted(int(x) for x in r[2].split()))] = \
                    (None, r[3], r[4], r[5], r[1], r)
        for key, (ring, center, nrm, aperture) in by.get(sz, []):
            skey = tuple(sorted(s for s, _o in ring))
            orbit = len(inst.get(key, ()))
            if skey not in pool or aperture > pool[skey][3]:
                pool[skey] = (ring, center, nrm, aperture, orbit, None)
        for (ring, center, nrm, aperture, orbit, src) in \
                sorted(pool.values(), key=lambda v: -v[3])[:want]:
            if src is not None:
                merged.append(src)
            else:
                merged.append((sz, orbit,
                               " ".join(str(s) for s, _o in ring),
                               center, nrm, aperture))
    return sorted(merged, key=lambda r: (r[0], -r[5]))
