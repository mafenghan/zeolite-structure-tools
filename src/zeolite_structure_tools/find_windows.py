# -*- coding: utf-8 -*-
"""find_windows -- arc I/O, lattice geometry, window-ring search.

This is the original *find_windows* part of the ring-table code, unchanged in
function: parse_arc / lattice / mic / RingFinder, plus the window-ring
classification (point_in_poly / classify) that decides which cycles are real
window rings.
"""
import numpy as np

from constants import VDW_O


# ------------------------------------------------------------- arc / lattice
def parse_arc(path):
    pbc, atoms = None, []
    for ln in open(path, encoding="utf-8", errors="replace"):
        p = ln.split()
        if ln.startswith("PBC") and len(p) == 7:
            pbc = [float(v) for v in p[1:7]]
        if len(p) >= 10 and p[4] == "CORE" and p[0].isalpha():
            atoms.append((int(p[-1]), p[0], float(p[1]), float(p[2]), float(p[3])))
    return pbc, atoms


def lattice(pbc):
    a, b, c, al, be, ga = pbc
    al, be, ga = np.deg2rad([al, be, ga])
    av = np.array([a, 0.0, 0.0])
    bv = np.array([b * np.cos(ga), b * np.sin(ga), 0.0])
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    cv = np.array([cx, cy, np.sqrt(max(c * c - cx * cx - cy * cy, 0.0))])
    return np.array([av, bv, cv])


def mic(L, d):
    d = np.asarray(d)
    shape = d.shape
    f = np.linalg.solve(L.T, d.reshape(-1, 3).T).T
    f -= np.round(f)
    return (L.T @ f.T).T.reshape(shape)


class RingFinder:
    """All simple T-cycles of length min_len..max_len (min-start canonical)."""

    def __init__(self, xyz, adj, L, max_len):
        self.xyz, self.adj, self.L, self.max_len = xyz, adj, L, max_len
        self.n = len(xyz)
        self.cycles = []

    def run(self, min_len):
        for start in range(self.n):
            self._dfs(start, [start], [self.xyz[start]], set([start]), min_len)

    def _dfs(self, start, path, coords, visited, min_len):
        last = path[-1]
        depth = len(path)
        if depth >= min_len and start in self.adj[last]:
            self.cycles.append((list(path), np.array(coords)))
        if depth == self.max_len:
            return
        for v in self.adj[last]:
            if v == start or v in visited or v < start:
                continue
            d = mic(self.L, self.xyz[v] - coords[-1])
            coords.append(coords[-1] + d)
            visited.add(v)
            path.append(v)
            self._dfs(start, path, coords, visited, min_len)
            path.pop()
            visited.remove(v)
            coords.pop()


# ------------------------------------------------------------------- geometry
def point_in_poly(p, poly, tol):
    n = len(poly)
    inside = False
    for i in range(n):
        a, b = poly[i], poly[(i + 1) % n]
        ab, ap = b - a, p - a
        cross = ab[0] * ap[1] - ab[1] * ap[0]
        d = np.linalg.norm(ab)
        if d < 1e-9:
            continue
        if abs(cross) / d < tol:
            if -tol <= np.dot(ap, ab) / d <= d + tol:
                return False
        if (a[1] > p[1]) != (b[1] > p[1]):
            if a[0] + (p[1] - a[1]) * ab[0] / ab[1] > p[0]:
                inside = not inside
    return inside


def classify(path, coords, T_xyz, o_idx, xyz_all, L):
    """Strong-ring check + geometry.  None if not a strong window ring."""
    n = len(path)
    # (0) closure: a real ring's first/last atoms are bonded in the unwrapped
    # coords; winding cycles (net lattice translation != 0) fail here.
    if np.linalg.norm(coords[0] - coords[-1]) > 3.4:
        return None
    # (1) chordless, on the ring's UNWRAPPED coords
    for i in range(n):
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue
            if np.linalg.norm(coords[i] - coords[j]) < 3.4:
                return None
    center = coords.mean(axis=0)
    _u, _s, vt = np.linalg.svd(coords - center)
    e1, e2, nrm = vt[0], vt[1], vt[2]
    if nrm[np.argmax(np.abs(nrm))] < 0:
        nrm = -nrm
    # (2) empty disk
    poly = np.array([((c - center) @ e1, (c - center) @ e2) for c in coords])
    rel = mic(L, T_xyz - center)
    proj = np.stack([rel @ e1, rel @ e2], axis=1)
    zdist = np.abs(rel @ nrm)
    rset = set(path)
    for i in range(len(T_xyz)):
        if i in rset or zdist[i] > 2.0:
            continue
        if point_in_poly(proj[i], poly, 0.15):
            return None
    # O-plane aperture
    rel_all = mic(L, xyz_all - center)
    perp = np.abs(rel_all @ nrm)
    radial = np.linalg.norm(rel_all - np.outer(rel_all @ nrm, nrm), axis=1)
    dO = 1e9
    for i in o_idx:
        if perp[i] < 2.5:
            dO = min(dO, radial[i])
    aperture = round(2 * max(dO - VDW_O, 0.0), 3) if dO < 1e8 else 0.0
    return center, nrm, aperture
