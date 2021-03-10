#!/usr/bin/python
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

"""
This was a little exercise to translate Grimme's D3 Fortran code into Python.
There is no new science as such!  It is possible to implement D3 corrections
directly within most electronic structure packages so the only possible use for
this is pedagogic or to look at individual terms of the total D3-correction
between a pair of atoms or molecules.  This code will read a Gaussian formatted
input/output file and compute the D3-density independent dispersion terms
without modification. Zero and Becke-Johnson damping schemes are both
implemented.

Written by:  Rob Paton and Kelvin Jackson
Last modified:  Mar 20, 2016
"""

from dataclasses import dataclass, field, InitVar
import json
import math
from typing import List

import jax.numpy as jnp
from jax.config import config
from prettytable import PrettyTable

config.update("jax_enable_x64", True)

from .ccParse import *
from .cli import cli
from .constants import (
    ALPHA6,
    ALPHA8,
    AU_TO_ANG,
    MAX_CONNECTIVITY,
    MAX_ELEMENTS,
)
from .parameters import BJ_PARMS, C6AB, R2R4, RAB, ZERO_PARMS
from .utils import check_inputs, getc6, getMollist, lin, ncoord


@dataclass
class D3Configuration:
    functional: str
    damp: str = "zero"
    s6: float = field(init=False)
    rs6: float = field(init=False)
    s8: float = field(init=False)
    a1: float = field(init=False)
    a2: float = field(init=False)
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
        else:
            raise RuntimeError(f"{self.damp} is an unknown damping scheme.")

        if not self.threebody:
            cfg += "    - 3-body term will not be calculated\n"
        else:
            cfg += "    - Including the Axilrod-Teller-Muto 3-body dispersion term\n"
        if self.intermolecular:
            cfg += "    - Only computing intermolecular dispersion interactions! This is not the total D3-correction\n"

        print(cfg)


def d3(
    config: D3Configuration,
    charges: List[float],
    *coordinates: float,
):
    """The code has a faithful implementation of D3-zero and D3-BJ. There are
    also some optional parts that will selectively ignore or scale certain
    interatomic terms. Without these lines our ‘scalefactor’ is set to 1, which
    is equivalent to standard D3 terms.
    """
    # van der Waals attractive R^-6
    attractive_r6_vdw = 0.0
    # van der Waals attractive R^-8
    attractive_r8_vdw = 0.0
    # Axilrod-Teller-Muto 3-body repulsive
    repulsive_abc = 0.0

    # not sure what this is...
    rs8 = 1.0

    natom = len(charges)

    # In case something clever needs to be done wrt inter and intramolecular interactions
    if config.bond_index is not None:
        molAatoms = getMollist(config.bond_index, 0)
        mols = []
        for j in range(natom):
            mols.append(0)
            for atom in molAatoms:
                if atom == j:
                    mols[j] = 1

    mxc = [0]
    for j in range(MAX_ELEMENTS):
        mxc.append(0)
        for k in range(natom):
            if charges[k] > -1:
                for l in range(MAX_CONNECTIVITY):
                    if isinstance(C6AB[j][j][l][l], (list, tuple)):
                        if C6AB[j][j][l][l][0] > 0:
                            mxc[j] = mxc[j] + 1
                break

    # Coordination number based on covalent radii
    cn = ncoord(charges, coordinates)

    # compute C6, C8, and C10 coefficietns from tabulated values (in C6AB) and fractional coordination
    for j in range(natom):
        # C6 coefficient
        C6jj = getc6(C6AB, mxc, charges, cn, j, j)

        z = int(charges[j])

        # C8 coefficient
        C8jj = 3.0 * C6jj * math.pow(R2R4[z], 2.0)

        # C10 coefficient
        C10jj = 49.0 / 40.0 * math.pow(C8jj, 2.0) / C6jj

    icomp = [0] * 100000
    cc6ab = [0] * 100000
    r2ab = [0] * 100000
    dmp = [0] * 100000

    for j in range(natom):
        dist = 0.0
        rr = 0.0
        attractive_r6_term = 0.0
        attractive_r8_term = 0.0
        ## This could be used to 'switch off' dispersion between bonded or geminal atoms ##
        scaling = False
        for k in range(j + 1, natom):
            scalefactor = 1.0

            if config.intermolecular == True:
                if mols[j] == mols[k]:
                    scalefactor = 0
                    print(f"   --- Ignoring interaction between atoms {j+1} and {k+1}")

            if scaling and config.bond_index is not None:
                if config.bond_index[j][k] == 1:
                    scalefactor = 0
                for l in range(natom):
                    if (
                        config.bond_index[j][l] != 0
                        and config.bond_index[k][l] != 0
                        and j != k
                        and config.bond_index[j][k] == 0
                    ):
                        scalefactor = 0
                    for m in range(natom):
                        if (
                            config.bond_index[j][l] != 0
                            and config.bond_index[l][m] != 0
                            and config.bond_index[k][m] != 0
                            and j != m
                            and k != l
                            and config.bond_index[j][m] == 0
                        ):
                            scalefactor = 1 / 1.2

            if k > j:
                ## Pythagoras in 3D to work out distance ##
                totdist = jnp.sqrt(
                    (coordinates[3 * j] - coordinates[3 * k]) ** 2
                    + (coordinates[3 * j + 1] - coordinates[3 * k + 1]) ** 2
                    + (coordinates[3 * j + 2] - coordinates[3 * k + 2]) ** 2
                )

                C6jk = getc6(C6AB, mxc, charges, cn, j, k)

                ## C8 parameters depend on C6 recursively
                atomA = int(charges[j])
                atomB = int(charges[k])

                C8jk = 3.0 * C6jk * R2R4[atomA] * R2R4[atomB]
                C10jk = 49.0 / 40.0 * math.pow(C8jk, 2.0) / C6jk

                # Evaluation of the attractive term dependent on R^-6 and R^-8
                if config.damp.casefold() == "zero".casefold():
                    dist = totdist
                    rr = RAB[atomA][atomB] / dist
                    tmp1 = config.rs6 * rr
                    damp6 = 1 / (1 + 6 * jnp.power(tmp1, ALPHA6))
                    tmp2 = rs8 * rr
                    damp8 = 1 / (1 + 6 * jnp.power(tmp2, ALPHA8))

                    attractive_r6_term = (
                        -config.s6 * C6jk * damp6 / jnp.power(dist, 6) * scalefactor
                    )
                    attractive_r8_term = (
                        -config.s8 * C8jk * damp8 / jnp.power(dist, 8) * scalefactor
                    )
                elif config.damp.casefold() == "bj".casefold():
                    dist = totdist
                    rr = RAB[atomA][atomB] / dist
                    rr = jnp.sqrt(C8jk / C6jk)
                    tmp1 = config.a1 * rr + config.a2
                    damp6 = jnp.power(tmp1, 6)
                    damp8 = jnp.power(tmp1, 8)

                    attractive_r6_term = (
                        -config.s6 * C6jk / (jnp.power(dist, 6) + damp6) * scalefactor
                    )
                    attractive_r8_term = (
                        -config.s8 * C8jk / (jnp.power(dist, 8) + damp8) * scalefactor
                    )
                else:
                    raise RuntimeError(f"{config.damp} is an unknown damping scheme.")

                if config.pairwise and scalefactor != 0:
                    print(
                        f"   --- Pairwise interaction between atoms {j+1} and {k+1}: Edisp = {attractive_r6_term+attractive_r8_term:.6f} kcal/mol",
                    )

                attractive_r6_vdw += attractive_r6_term
                attractive_r8_vdw += attractive_r8_term

                jk = int(lin(k, j))
                icomp[jk] = 1
                cc6ab[jk] = jnp.sqrt(C6jk)
                r2ab[jk] = jnp.power(dist, 2)
                dmp[jk] = jnp.cbrt(1.0 / rr)

    if config.threebody:
        e63 = 0.0
        for iat in range(natom):
            for jat in range(natom):
                ij = int(lin(jat, iat))
                if icomp[ij] == 1:
                    for kat in range(jat, natom):
                        ik = int(lin(kat, iat))
                        jk = int(lin(kat, jat))

                        if (
                            kat > jat
                            and jat > iat
                            and icomp[ik] != 0
                            and icomp[jk] != 0
                        ):
                            rav = (4.0 / 3.0) / (dmp[ik] * dmp[jk] * dmp[ij])
                            tmp = 1.0 / (1.0 + 6.0 * rav ** ALPHA6)

                            c9 = cc6ab[ij] * cc6ab[ik] * cc6ab[jk]
                            d2 = [
                                r2ab[ij],
                                r2ab[jk],
                                r2ab[ik],
                            ]
                            t1 = (d2[0] + d2[1] - d2[2]) / jnp.sqrt(d2[0] * d2[1])
                            t2 = (d2[0] + d2[2] - d2[1]) / jnp.sqrt(d2[0] * d2[2])
                            t3 = (d2[2] + d2[1] - d2[0]) / jnp.sqrt(d2[1] * d2[2])
                            ang = 0.375 * t1 * t2 * t3 + 1.0
                            e63 = e63 + tmp * c9 * ang / (d2[0] * d2[1] * d2[2]) ** 1.50

        repulsive_abc_term = config.s6 * e63
        repulsive_abc += repulsive_abc_term

    return attractive_r6_vdw + attractive_r8_vdw + repulsive_abc


def main():
    # Takes arguments: (1) damping style, (2) s6, (3) rs6, (4) s8, (5) 3-body on/off, (6) input file(s)
    args = cli()

    files = args.infiles

    # prepare table of results
    x = PrettyTable()
    x.field_names = ["Input", "Total (au)"]

    for f in files:
        extension = f.split(".")[1]
        # parse Gaussian input files
        if extension in ["com", "gjf"]:
            data = getinData(f)
        # parse Gaussian output files
        elif extension in ["out", "log"]:
            data = getoutData(f)
        # parse PDB file
        elif extension == "pdb":
            data = getpdbData(f)
        # parse plain text file
        elif extension == "txt":
            data = get_simple_data(f)
        elif extension == "json":
            with open(f, "r") as j:
                data = json.load(j)
        else:
            raise RuntimeError(f"Unrecognized file format {extension}")

        functional = ""
        if isinstance(data, dict):
            charges = data["molecule"]["symbols"]
            # reshape coordinates as (nat, 3) and convert to angstrom
            coordinates = [[0.0, 0.0, 0.0] for _ in range(len(charges))]
            for j in range(3):
                for i in range(len(charges)):
                    coordinates[i][j] = (
                        data["molecule"]["geometry"][3 * i + j] * AU_TO_ANG
                    )
            functional = data["model"]["method"]
        else:
            charges = data.CHARGES
            coordinates = data.CARTESIANS
            functional = data.FUNCTIONAL

        check_inputs(charges=charges, coordinates=coordinates)

        config = D3Configuration(
            functional=functional,
            damp=args.damp,
            _s6=args.s6,
            _rs6=args.rs6,
            _s8=args.s8,
            _a1=args.a1,
            _a2=args.a2,
            threebody=args.three,
            intermolecular=args.inter,
            pairwise=args.pairwise,
        )

        total_vdw = d3(
            config,
            charges,
            *coordinates,
        )

        # convert to atomic units for final printout
        row = [f, total_vdw]
        x.add_row(row)

    print(x)


if __name__ == "__main__":
    main()
