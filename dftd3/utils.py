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


from dataclasses import InitVar, dataclass, field
from typing import List

import jax.numpy as jnp

from .constants import AU_TO_ANG
from .parameters import BJ_PARMS, RCOV, ZERO_PARMS


def der_order(order: int) -> str:
    """Return correct suffix for order.

    >>> der_order(1)
    ... "1-st"
    >>> der_order(4)
    ... "4-th"
    """
    _suffix = ["st", "nd", "rd"]
    if order - 1 > 3:
        return f"{order}-th"
    else:
        return f"{order}-{_suffix[order - 1]}"


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


def ncoord(charges, coordinates, k1=16, k2=4 / 3):
    """Calculation of atomic coordination numbers.

    Notes
    -----
    The constants ``k1`` and ``k2`` are used to determine fractional connectivities between 2 atoms:

      - ``k1`` is the exponent used in summation
      - ``k2`` is used a fraction of the summed single-bond radii

    These values are copied verbatim from Grimme's code.
    """

    natom = len(charges)

    check_inputs(charges=charges, coordinates=coordinates)

    coordinates = [coordinate * AU_TO_ANG for coordinate in coordinates]
    cn = []

    for i in range(natom):
        xn = 0.0
        for iat in range(natom):
            if iat != i:
                r = jnp.sqrt(
                    (coordinates[3 * i] - coordinates[3 * iat]) ** 2
                    + (coordinates[3 * i + 1] - coordinates[3 * iat + 1]) ** 2
                    + (coordinates[3 * i + 2] - coordinates[3 * iat + 2]) ** 2
                )

                Zi = int(charges[i])
                Ziat = int(charges[iat])

                rco = k2 * (RCOV[Zi] + RCOV[Ziat])
                rr = rco / r
                damp = 1.0 / (1.0 + jnp.exp(-k1 * (rr - 1.0)))
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
    iat = int(atomtype[a])
    jat = int(atomtype[b])

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
                    tmp1 = jnp.exp(k3 * r)
                    rsum = rsum + tmp1
                    csum = csum + tmp1 * c6

    if rsum > 0:
        c6 = csum / rsum
    else:
        c6 = c6mem

    if c6 is None:
        raise RuntimeError("Computation of C6 failed.")

    return c6


def check_inputs(*, charges, coordinates):

    natom = len(charges)

    if natom != len(coordinates) // 3:
        raise RuntimeError(
            f"The size of the coordinates (3*{len(coordinates) // 3}) and atom types ({natom}) arrays do not match!"
        )


@dataclass
class D3Configuration:
    functional: str
    damp: str = "zero"
    nprocs: int = 1
    s6: float = field(init=False)
    rs6: float = field(init=False)
    s8: float = field(init=False)
    a1: float = field(init=False)
    a2: float = field(init=False)
    threebody: bool = False
    bond_index: List[List[int]] = None
    intermolecular: bool = False
    pairwise: bool = False

    # initialization-only variables
    _s6: InitVar[float] = 0.0
    _rs6: InitVar[float] = 0.0
    _s8: InitVar[float] = 0.0
    _a1: InitVar[float] = 0.0
    _a2: InitVar[float] = 0.0

    def __post_init__(self, _s6, _rs6, _s8, _a1, _a2):
        which = self.damp if self.damp == "zero" else "Becke-Johson"
        cfg = f">>> D3-dispersion correction with {which} damping\n"

        if self.damp.casefold() == "zero".casefold():
            if any(map(lambda x: x == 0.0, [_s6, _rs6, _s8])):
                self.s6, self.rs6, self.s8 = ZERO_PARMS[self.functional]
                where = "from database"
            else:
                self.s6 = _s6
                self.rs6 = _rs6
                self.s8 = _s8
                where = "user-defined"
            cfg += f"    - Damping parameters ({where}): s6 = {self.s6}; rs6 = {self.rs6}; s8 = {self.s8}\n"
            # initialize remaining parameters
            self.a1 = _a1
            self.a2 = _a2
        elif self.damp.casefold() == "bj".casefold():
            if any(map(lambda x: x == 0.0, [_s6, _s8, _a1, _a2])):
                self.s6, self.a1, self.s8, self.a2 = BJ_PARMS[self.functional]
                where = "from database"
            else:
                self.s6 = _s6
                self.a1 = _a1
                self.s8 = _s8
                self.a2 = _a2
                where = "user-defined"
            cfg += f"    - Damping parameters: s6 = {self.s6}; s8 = {self.s8}; a1 = {self.a1}; a2 = {self.a2}\n"
            # initialize remaining parameters
            self.rs6 = _rs6
        else:
            raise RuntimeError(f"{self.damp} is an unknown damping scheme.")

        if not self.threebody:
            cfg += "    - 3-body term will not be calculated\n"
        else:
            cfg += "    - Including the Axilrod-Teller-Muto 3-body dispersion term\n"
        if self.intermolecular:
            cfg += "    - Only computing intermolecular dispersion interactions! This is not the total D3-correction\n"

        print(cfg)
