# zeolite-ring-table

Strong window-ring table generator for zeolite framework `.arc` files
(BIOSYM archive 2, PBC).  For every framework it writes one CSV row per
inequivalent **T site** and per inequivalent **strong ring** — the ring kind
counts follow the IZA *Essential rings* — with ring centre, plane normal and
O-plane aperture, ready for placing molecules along the standard files.

## Layout (original parts kept)

| file | part | contents |
|---|---|---|
| `find_windows.py` | window search | `parse_arc`, `lattice`, `mic`, `RingFinder`, `point_in_poly`, `classify` |
| `ring_tables_v2.py` | table builder | `spglib_data`, `build_perms`, `process`, `load_iza`, `main` |
| `ring_v2_supercell_fix.py` | supercell gap-fill | `canon_key`, `make_supercell`, `load_supercell`, `fill_from_supercell` |
| `iza_essential_rings.json` | IZA data | essential-ring sizes / kind counts for 226 framework types |
| `constants.py` | shared | `T_ELEMS`, `VDW_O`, `MIN_SIZE`, `FALLBACK_SIZES` |
| `tests/` | self-checks | one file per code part |

## Install

```bash
pip install -r requirements.txt    # numpy, spglib
```

## Usage

```bash
python ring_tables_v2.py SRC_DIR -o ring_tables.csv [--codes A B C]
       [--iza iza.json] [--supercell-root DIR] [--no-fill]
```

| argument | meaning |
|---|---|
| `SRC_DIR` | directory with `<CODE>/<CODE>.arc` subfolders |
| `-o` | output CSV (default `zeolite_ring_tables.csv`) |
| `--codes` | restrict to these framework codes |
| `--iza` | external IZA essential-rings JSON (default: bundled `iza_essential_rings.json`) |
| `--supercell-root` | reuse `<CODE>/framework_orig.arc` supercells; without it a 2×2×2 is generated on the fly |
| `--no-fill` | skip the supercell gap-fill pass |

## Output columns

```
code, record_type, site_id, orbit_size, wyckoff, t_serials,
cx, cy, cz, nx, ny, nz, aperture
```

| column | T rows | RING rows |
|---|---|---|
| `site_id` | row index within the code | row index within the code |
| `orbit_size` | site multiplicity (equivalent T atoms per cell) | copies of that ring kind per unit cell |
| `wyckoff` | Wyckoff letter of the T site | (empty) |
| `t_serials` | serials of the site's atoms | serials of the ring's T atoms |
| `cx,cy,cz` | (empty) | ring centre (cartesian, Å) |
| `nx,ny,nz` | (empty) | ring plane normal |
| `aperture` | (empty) | 2 × (in-plane O distance − 1.35 Å), the O-plane free diameter of the window, Å |

Positive control: FAU's 12-ring aperture comes out 7.33 Å (IZA free
diameter 7.35 Å).

## Ring recipe

A ring kind must be a simple T-cycle that is:

1. **closed** — first/last ring atoms bonded in the unwrapped coordinates
   (rejects winding cycles with a net lattice translation),
2. **chordless** — no T–T bond between non-consecutive ring atoms,
3. **empty-disk** — no T projection inside the ring's planar polygon,

restricted to the framework's IZA essential-ring sizes (>= 6; odd sizes kept
when IZA lists them).  Surviving cycles are reduced to inequivalent kinds with
spglib, then trimmed per size to the IZA essential-kind count, keeping the
largest O-plane apertures.  A supercell pass fills kinds that only exist as
multi-edge / wrap-around rings of the primitive cell; ring serials still refer
to the standard file.

Kind counts follow the embedded IZA counts *by construction* (the output
trims/pads to them); what the recipe decides is **which** rings are kept —
the aperture ranking is the heuristic there.  Verify apertures when in doubt.

## Embedded IZA data

`iza_essential_rings.json` covers 226 framework types (essential-ring sizes
and symmetry-inequivalent kind counts, scraped from the IZA *Database of
Zeolite Structures*, iza-structure.org).  Codes absent from it fall back to
even sizes 6..12, <= 8 kinds each.  The IZA data is subject to the IZA-SC
database terms; redistribute with that in mind.

## Tests

```bash
python tests/test_find_windows.py
python tests/test_find_windows_classify.py
python tests/test_supercell_fix.py
```

Pure-geometry self-checks (no spglib, no structure files needed).
