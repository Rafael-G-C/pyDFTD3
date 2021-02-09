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

import math
import json
from prettytable import PrettyTable

from .cli import cli
from .ccParse import *
from .constants import (
    AU_TO_ANG,
    AU_TO_KCAL,
    MAX_CONNECTIVITY,
    MAX_ELEMENTS,
    ALPHA6,
    ALPHA8,
)
from .parameters import BJ_PARMS, R2R4, RAB, ZERO_PARMS, C6AB
from .utils import E_to_index, getc6, getMollist, lin, ncoord

"""
This was a little exercise to translate Grimme's D3
Fortran code into Python. There is no new science as such!
It is possible to implement D3 corrections directly within most
electronic structure packages so the only possible use for this
is pedagogic or to look at individual terms of the total
D3-correction between a pair of atoms or molecules.
This code will read a Gaussian formatted input/output file and
compute the D3-density independent dispersion terms without
modification. Zero and Becke-Johnson damping schemes are both
implemented.

Written by:  Rob Paton and Kelvin Jackson
Last modified:  Mar 20, 2016
"""

## Functional Specific D3 parameters
rs8 = 1.0


def d3(
    atoms,
    coordinates,
    *,
    functional,
    bond_index=None,
    s6=0.0,
    rs6=0.0,
    s8=0.0,
    a1=0.0,
    a2=0.0,
    damp="zero",
    threebody=False,
    intermolecular=False,
    pairwise=False,
    verbose=1,
):

    # van der Waals attractive R^-6
    attractive_r6_vdw = 0.0
    # van der Waals attractive R^-8
    attractive_r8_vdw = 0.0
    # Axilrod-Teller-Muto 3-body repulsive
    repulsive_abc = 0.0

    natom = len(atoms)

    # In case something clever needs to be done wrt inter and intramolecular interactions
    if bond_index is not None:
        molAatoms = getMollist(bond_index, 0)
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
            if E_to_index(atoms[k]) > -1:
                for l in range(MAX_CONNECTIVITY):
                    if isinstance(C6AB[j][j][l][l], (list, tuple)):
                        if C6AB[j][j][l][l][0] > 0:
                            mxc[j] = mxc[j] + 1
                break

    # Coordination number based on covalent radii
    cn = ncoord(natom, atoms, coordinates)

    # compute C6, C8, and C10 coefficietns from tabulated values (in C6AB) and fractional coordination
    for j in range(natom):
        # C6 coefficient
        C6jj = getc6(C6AB, mxc, atoms, cn, j, j)

        z = E_to_index(atoms[j])

        # C8 coefficient
        C8jj = 3.0 * C6jj * math.pow(R2R4[z], 2.0)

        # C10 coefficient
        C10jj = 49.0 / 40.0 * math.pow(C8jj, 2.0) / C6jj

    icomp = [0] * 100000
    cc6ab = [0] * 100000
    r2ab = [0] * 100000
    dmp = [0] * 100000

    ## Compute and output the individual components of the D3 energy correction ##
    if verbose:
        which = damp if damp == "zero" else "Becke-Johson"
        print(f"   D3-dispersion correction with {which} damping:\n")
    # print "\n   Atoms  Types  C6            C8            E6              E8"
    if damp == "zero":
        if s6 == 0.0 or rs6 == 0.0 or s8 == 0.0:
            s6, rs6, s8 = ZERO_PARMS[functional]
        else:
            if verbose:
                print(" manual parameters have been defined")
        if verbose:
            print(f"   Zero-damping parameters: s6 = {s6}; rs6 = {rs6}; s8 = {s8}")

    if damp == "bj":
        if s6 == 0.0 or s8 == 0.0 or a1 == 0.0 or a2 == 0.0:
            s6, a1, s8, a2 = BJ_PARMS[functional]
        else:
            if verbose:
                print(" manual parameters have been defined")
        if verbose:
            print(
                f"   BJ-damping parameters: s6 = {s6}; s8 = {s8}; a1 = {a1}; a2 = {a2}"
            )

    if verbose:
        if threebody == False:
            print("   3-body term will not be calculated")
        else:
            print("   Including the Axilrod-Teller-Muto 3-body dispersion term")
        if intermolecular == True:
            print(
                "   Only computing intermolecular dispersion interactions! This is not the total D3-correction\n"
            )

    for j in range(natom):
        ## This could be used to 'switch off' dispersion between bonded or geminal atoms ##
        scaling = False
        for k in range(j + 1, natom):
            scalefactor = 1.0

            if intermolecular == True:
                if mols[j] == mols[k]:
                    scalefactor = 0
                    print(f"   --- Ignoring interaction between atoms {j+1} and {k+1}")

            if scaling == True and bond_index is not None:
                if bond_index[j][k] == 1:
                    scalefactor = 0
                for l in range(natom):
                    if (
                        bond_index[j][l] != 0
                        and bond_index[k][l] != 0
                        and j != k
                        and bond_index[j][k] == 0
                    ):
                        scalefactor = 0
                    for m in range(natom):
                        if (
                            bond_index[j][l] != 0
                            and bond_index[l][m] != 0
                            and bond_index[k][m] != 0
                            and j != m
                            and k != l
                            and bond_index[j][m] == 0
                        ):
                            scalefactor = 1 / 1.2

            if k > j:
                ## Pythagoras in 3D to work out distance ##
                totdist = math.sqrt(
                    (coordinates[3 * j] - coordinates[3 * k]) ** 2
                    + (coordinates[3 * j + 1] - coordinates[3 * k + 1]) ** 2
                    + (coordinates[3 * j + 2] - coordinates[3 * k + 2]) ** 2
                )

                C6jk = getc6(C6AB, mxc, atoms, cn, j, k)

                ## C8 parameters depend on C6 recursively
                atomA = E_to_index(atoms[j])
                atomB = E_to_index(atoms[k])

                C8jk = 3.0 * C6jk * R2R4[atomA] * R2R4[atomB]
                C10jk = 49.0 / 40.0 * math.pow(C8jk, 2.0) / C6jk

                # Evaluation of the attractive term dependent on R^-6 and R^-8
                if damp == "zero":
                    dist = totdist / AU_TO_ANG
                    rr = RAB[atomA][atomB] / dist
                    tmp1 = rs6 * rr
                    damp6 = 1 / (1 + 6 * math.pow(tmp1, ALPHA6))
                    tmp2 = rs8 * rr
                    damp8 = 1 / (1 + 6 * math.pow(tmp2, ALPHA8))

                    attractive_r6_term = (
                        -s6
                        * C6jk
                        * damp6
                        / math.pow(dist, 6)
                        * AU_TO_KCAL
                        * scalefactor
                    )
                    attractive_r8_term = (
                        -s8
                        * C8jk
                        * damp8
                        / math.pow(dist, 8)
                        * AU_TO_KCAL
                        * scalefactor
                    )

                if damp == "bj":
                    dist = totdist / AU_TO_ANG
                    rr = RAB[atomA][atomB] / dist
                    rr = math.pow((C8jk / C6jk), 0.5)
                    tmp1 = a1 * rr + a2
                    damp6 = math.pow(tmp1, 6)
                    damp8 = math.pow(tmp1, 8)

                    attractive_r6_term = (
                        -s6
                        * C6jk
                        / (math.pow(dist, 6) + damp6)
                        * AU_TO_KCAL
                        * scalefactor
                    )
                    attractive_r8_term = (
                        -s8
                        * C8jk
                        / (math.pow(dist, 8) + damp8)
                        * AU_TO_KCAL
                        * scalefactor
                    )

                if pairwise == True and scalefactor != 0:
                    print(
                        f"   --- Pairwise interaction between atoms {j+1} and {k+1}: Edisp = {attractive_r6_term+attractive_r8_term:.6f} kcal/mol",
                    )

                attractive_r6_vdw += attractive_r6_term
                attractive_r8_vdw += attractive_r8_term

                jk = int(lin(k, j))
                icomp[jk] = 1
                cc6ab[jk] = math.sqrt(C6jk)
                r2ab[jk] = dist ** 2
                dmp[jk] = (1.0 / rr) ** (1.0 / 3.0)

    e63 = 0.0
    for iat in range(natom):
        for jat in range(natom):
            ij = int(lin(jat, iat))
            if icomp[ij] == 1:
                for kat in range(jat, natom):
                    ik = int(lin(kat, iat))
                    jk = int(lin(kat, jat))

                    if kat > jat and jat > iat and icomp[ik] != 0 and icomp[jk] != 0:
                        rav = (4.0 / 3.0) / (dmp[ik] * dmp[jk] * dmp[ij])
                        tmp = 1.0 / (1.0 + 6.0 * rav ** ALPHA6)

                        c9 = cc6ab[ij] * cc6ab[ik] * cc6ab[jk]
                        d2 = [0] * 3
                        d2[0] = r2ab[ij]
                        d2[1] = r2ab[jk]
                        d2[2] = r2ab[ik]
                        t1 = (d2[0] + d2[1] - d2[2]) / math.sqrt(d2[0] * d2[1])
                        t2 = (d2[0] + d2[2] - d2[1]) / math.sqrt(d2[0] * d2[2])
                        t3 = (d2[2] + d2[1] - d2[0]) / math.sqrt(d2[1] * d2[2])
                        ang = 0.375 * t1 * t2 * t3 + 1.0
                        e63 = e63 + tmp * c9 * ang / (d2[0] * d2[1] * d2[2]) ** 1.50

    repulsive_abc_term = s6 * e63 * AU_TO_KCAL
    repulsive_abc += repulsive_abc_term

    return attractive_r6_vdw, attractive_r8_vdw, repulsive_abc


def main():
    # Takes arguments: (1) damping style, (2) s6, (3) rs6, (4) s8, (5) 3-body on/off, (6) input file(s)
    args = cli()

    verbose = args.verbose
    damp = args.damp
    s6 = args.s6
    rs6 = args.rs6
    s8 = args.s8
    bj_a1 = args.a1
    bj_a2 = args.a2
    abc_term = args.three
    intermolecular = args.inter
    pairwise = args.pairwise
    files = args.infiles

    x = PrettyTable()

    fields = ["Input", "D3(R6)", "D3(R8)", "Total (au)"]
    if abc_term:
        fields.insert(3, "D3(3-body)")
    x.field_names = fields

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

        if isinstance(data, dict):
            atoms = data["molecule"]["symbols"]
            # reshape coordinates as (nat, 3) and convert to angstrom
            coordinates = [[0.0, 0.0, 0.0] for _ in range(len(atoms))]
            for j in range(3):
                for i in range(len(atoms)):
                    coordinates[i][j] = (
                        data["molecule"]["geometry"][3 * i + j] * AU_TO_ANG
                    )
            functional = data["model"]["method"]
        else:
            atoms = data.ATOMTYPES
            coordinates = data.CARTESIANS
            functional = data.FUNCTIONAL

        r6, r8, abc = d3(
            atoms,
            coordinates,
            functional=functional,
            s6=s6,
            rs6=rs6,
            s8=s8,
            a1=bj_a1,
            a2=bj_a2,
            damp=damp,
            threebody=abc_term,
            intermolecular=intermolecular,
            pairwise=pairwise,
            verbose=verbose,
        )

        # convert to atomic units for final printout
        attractive_r6_vdw = r6 / AU_TO_KCAL
        attractive_r8_vdw = r8 / AU_TO_KCAL
        repulsive_abc = abc / AU_TO_KCAL

        total_vdw = attractive_r6_vdw + attractive_r8_vdw + repulsive_abc

        row = [f, attractive_r6_vdw, attractive_r8_vdw, total_vdw]
        if abc_term:
            row.insert(3, repulsive_abc)
        x.add_row(row)

    print(x)


if __name__ == "__main__":
    main()
