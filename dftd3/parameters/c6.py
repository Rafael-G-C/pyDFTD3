# -*- coding: utf-8 -*-

from ..constants import MAX_CONNECTIVITY, MAX_ELEMENTS
from .pars import PARS


def copyc6(max_elem=MAX_ELEMENTS, maxc=MAX_CONNECTIVITY):
    """Reference systems are read in to compute coordination number dependent dispersion coefficients."""

    c6ab = [[0] * max_elem for _ in range(max_elem)]
    nlines = 32385

    for iat in range(0, max_elem):
        for jat in range(0, max_elem):
            c6ab[iat][jat] = [[0] * maxc for _ in range(maxc)]
    kk = 0

    for nn in range(0, nlines):
        kk = nn * 5
        iadr = 0
        jadr = 0
        iat = int(PARS[kk + 1]) - 1
        jat = int(PARS[kk + 2]) - 1

        while iat > 99:
            iadr = iadr + 1
            iat = iat - 100
        while jat > 99:
            jadr = jadr + 1
            jat = jat - 100

        c6ab[iat][jat][iadr][jadr] = []
        c6ab[iat][jat][iadr][jadr].append(PARS[kk])
        c6ab[iat][jat][iadr][jadr].append(PARS[kk + 3])
        c6ab[iat][jat][iadr][jadr].append(PARS[kk + 4])

        c6ab[jat][iat][jadr][iadr] = []
        c6ab[jat][iat][jadr][iadr].append(PARS[kk])
        c6ab[jat][iat][jadr][iadr].append(PARS[kk + 4])
        c6ab[jat][iat][jadr][iadr].append(PARS[kk + 3])

    return c6ab


C6AB = copyc6(MAX_ELEMENTS, MAX_CONNECTIVITY)
