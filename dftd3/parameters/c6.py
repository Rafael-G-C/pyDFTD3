# -*- coding: utf-8 -*-
#
# pyDFTD3 -- Python implementation of Grimme's D3 dispersion correction.
# Copyright (C) 2020 Rob Paton and contributors.
#
# This file is part of pyDFTD3.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# For information on the complete list of contributors to the
# pyDFTD3, see: <http://github.com/bobbypaton/pyDFTD3/>
#

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
