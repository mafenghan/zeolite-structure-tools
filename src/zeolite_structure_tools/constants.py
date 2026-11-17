# -*- coding: utf-8 -*-
"""Shared constants for the zeolite ring-table tools."""

# T-atom elements recognised in .arc files
T_ELEMS = {"Si", "Al", "P", "Sn", "Ga", "Ge", "B", "Ti", "Zn", "Be"}

# O van der Waals radius used for the O-plane aperture
VDW_O = 1.35

# window rings are >= 6 T atoms; odd sizes 7/9/11 kept when IZA lists them
MIN_SIZE = 6

# practical window range for codes without an IZA record
FALLBACK_SIZES = [6, 8, 10, 12]
