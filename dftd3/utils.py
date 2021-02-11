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


from math import exp, sqrt

from qcelemental import periodictable

from .parameters import RCOV
from .constants import AU_TO_ANG


def E_to_index(element):
    """Convert element symbol to 0-based index."""
    return periodictable.to_Z(element) - 1


def getMollist(bondmatrix, startatom):
    """From connectivity, establish if there is more than one molecule."""

    # The list of atoms in a molecule
    atomlist = []
    atomlist.append(startatom)
    count = 0

    while count < 100:
        for atom in atomlist:
            for i in range(0, len(bondmatrix[atom])):
                if bondmatrix[atom][i] == 1:
                    alreadyfound = 0
                    for at in atomlist:
                        if i == at:
                            alreadyfound = 1
                    if alreadyfound == 0:
                        atomlist.append(i)
        count = count + 1

    return atomlist


def ncoord(atomtype, coordinates, k1=16, k2=4 / 3):
    """Calculation of atomic coordination numbers.

    Notes
    -----
    The constants ``k1`` and ``k2`` are used to determine fractional connectivities between 2 atoms:

      - ``k1`` is the exponent used in summation
      - ``k2`` is used a fraction of the summed single-bond radii

    These values are copied verbatim from Grimme's code.
    """
    flat_coords = [coordinate * AU_TO_ANG for coordinate in flat_coords]
    cn = []
    natom = len(atomtypes)
    
    if natom != len(coordinates) // 3:
         raise RuntimeError("The size of the coordinates and atom types arrays do not match") 
    for i in range(natom):
        xn = 0.0
        for iat in range(natom):
            if iat != i:
                r = sqrt(
                    (flat_coords[3 * i] - flat_coords[3 * iat]) ** 2
                    + (flat_coords[3 * i + 1] - flat_coords[3 * iat + 1]) ** 2
                    + (flat_coords[3 * i + 2] - flat_coords[3 * iat + 2]) ** 2
                )

                Zi = atomtype[i]
                Ziat = atomtype[iat]

                rco = k2 * (RCOV[Zi] + RCOV[Ziat])
                rr = rco / r
                damp = 1.0 / (1.0 + exp(-k1 * (rr - 1.0)))
                xn = xn + damp
        cn.append(xn)

    return cn


def lin(i1, i2):
    """Linear interpolation."""

    idum1 = max(i1, i2)
    idum2 = min(i1, i2)
    lin = idum2 + idum1 * (idum1 - 1) / 2

    return lin


def getc6(c6ab, mxc, atomtype, cn, a, b, k3=-4.0):
    """Obtain the C6 coefficient for the interaction between atoms A and B.

    Notes
    -----
    The constant ``k3`` is copied verbatim from Grimme's code.
    """

    # atomic charges for atoms A and B, respectively
    iat = atomtype[a]
    jat = atomtype[b]

    c6mem = None
    rsum = 0.0
    csum = 0.0
    c6 = 0.0

    for i in range(0, mxc[iat]):
        for j in range(0, mxc[jat]):
            if isinstance(c6ab[iat][jat][i][j], (list, tuple)):
                c6 = c6ab[iat][jat][i][j][0]
                if c6 > 0:
                    c6mem = c6
                    cn1 = c6ab[iat][jat][i][j][1]
                    cn2 = c6ab[iat][jat][i][j][2]

                    r = (cn1 - cn[a]) ** 2 + (cn2 - cn[b]) ** 2
                    tmp1 = exp(k3 * r)
                    rsum = rsum + tmp1
                    csum = csum + tmp1 * c6

    if rsum > 0:
        c6 = csum / rsum
    else:
        c6 = c6mem

    if c6 is None:
        raise RuntimeError("Computation of C6 failed.")

    return c6
