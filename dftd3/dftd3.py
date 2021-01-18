#!/usr/bin/python
# -*- coding: utf-8 -*-

import math
import sys

# For reading Gaussian formatted input/output files
from .ccParse import *
from .constants import (
    AU_TO_ANG,
    AU_TO_KCAL,
    MAX_CONNECTIVITY,
    MAX_ELEMENTS,
    ALPHA6,
    ALPHA8,
)

# Dependent on parameter file
from .parameters import BJ_PARMS, R2R4, RAB, ZERO_PARMS, copyc6
from .utils import E_to_index, getc6, getMollist, lin, ncoord

# THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Comments and/or additions are welcome (send e-mail to:
# robert.paton@chem.ox.ac.uk

#######################################################################
#                          dftd3.py                                   #
#                                                                     #
#    This was a little exercise to translate Grimme's D3              #
#    Fortran code into Python. There is no new science as such!       #
#    It is possible to implement D3 corrections directly within most  #
#    electronic structure packages so the only possible use for this  #
#    is pedagogic or to look at individual terms of the total         #
#    D3-correction between a pair of atoms or molecules.              #
#    This code will read a Gaussian formatted input/output file and   #
#    compute the D3-density independent dispersion terms without      #
#    modification. Zero and Becke-Johnson damping schemes are both    #
#    implemented.                                                     #
#######################################################################
#######  Written by:  Rob Paton and Kelvin Jackson ####################
#######  Last modified:  Mar 20, 2016 #################################
#######################################################################

## Arrays for attractive and repulsive interactions ##
attractive_vdw = [0]
repulsive_vdw = [0]
total_vdw = [0]

## Functional Specific D3 parameters
rs8 = 1.0

## The computation of the D3 dispersion correction
class calcD3:
    def __init__(
        self,
        data,
        functional,
        s6,
        rs6,
        s8,
        a1,
        a2,
        damp,
        abc,
        intermolecular,
        pairwise,
        verbose,
    ):
        ## Get from pars.py
        c6ab = copyc6(MAX_ELEMENTS, MAX_CONNECTIVITY)

        # verbose = True
        ## Arrays for atoms and Cartesian coordinates ##
        try:
            atomtype = data.ATOMTYPES
            cartesians = data.CARTESIANS
        except:
            atomtype = data.atom_types
            cartesians = data.cartesians
        natom = len(atomtype)

        xco = []
        yco = []
        zco = []
        for at in cartesians:
            xco.append(at[0])
            yco.append(at[1])
            zco.append(at[2])

        ## In case something clever needs to be done wrt inter and intramolecular interactions
        if hasattr(data, "BONDINDEX"):
            molAatoms = getMollist(data.BONDINDEX, 0)
            mols = []
            for j in range(0, natom):
                mols.append(0)
                for atom in molAatoms:
                    if atom == j:
                        mols[j] = 1

        ## Names are pretty obvious...
        self.attractive_r6_vdw = 0.0
        self.attractive_r8_vdw = 0.0
        self.repulsive_abc = 0.0

        mxc = [0]
        for j in range(0, MAX_ELEMENTS):
            mxc.append(0)
            for k in range(0, natom):
                if E_to_index(atomtype[k]) > -1:
                    for l in range(0, MAX_CONNECTIVITY):
                        if isinstance(c6ab[j][j][l][l], (list, tuple)):
                            if c6ab[j][j][l][l][0] > 0:
                                mxc[j] = mxc[j] + 1
                    break

        ## Coordination number based on covalent radii
        cn = ncoord(natom, atomtype, xco, yco, zco)

        ## C6 - Need to calculate these from fractional coordination
        # print "\n           R0(Ang)        CN"
        # print "   #########################"
        for j in range(0, natom):
            C6jj = getc6(c6ab, mxc, atomtype, cn, j, j)

            z = E_to_index(atomtype[j])

            C8jj = 3.0 * C6jj * math.pow(R2R4[z], 2.0)
            # C10 coefficient
            C10jj = 49.0 / 40.0 * math.pow(C8jj, 2.0) / C6jj

        icomp = [0] * 100000
        cc6ab = [0] * 100000
        r2ab = [0] * 100000
        dmp = [0] * 100000

        ## Compute and output the individual components of the D3 energy correction ##
        # print "\n   Atoms  Types  C6            C8            E6              E8"
        if damp == "zero":
            if verbose:
                print("\n   D3-dispersion correction with zero-damping:", end=" ")
            if s6 == 0.0 or rs6 == 0.0 or s8 == 0.0:
                if functional is not None:
                    s6, rs6, s8 = ZERO_PARMS[functional]
                    if verbose:
                        print(
                            f"detected {functional} functional - using default zero-damping parameters",
                        )
                else:
                    if verbose:
                        print(
                            "   WARNING: Damping parameters not specified and no functional could be read!\n"
                        )
                        sys.exit()
            else:
                if verbose:
                    print(" manual parameters have been defined")
            if verbose:
                print(
                    "   Zero-damping parameters:", "s6 =", s6, "rs6 =", rs6, "s8 =", s8
                )

        if damp == "bj":
            if verbose:
                print(
                    "\n   D3-dispersion correction with Becke_Johnson damping:", end=" "
                )
            if s6 == 0.0 or s8 == 0.0 or a1 == 0.0 or a2 == 0.0:
                if functional is not None:
                    s6, a1, s8, a2 = BJ_PARMS[functional]
                    if verbose:
                        print(
                            f"detected {functional} functional - using default BJ-damping parameters",
                        )
                else:
                    if verbose:
                        print(
                            "   WARNING: Damping parameters not specified and no functional could be read!\n"
                        )
                        sys.exit()
            else:
                if verbose:
                    print(" manual parameters have been defined")
            if verbose:
                print(
                    "   BJ-damping parameters:",
                    "s6 =",
                    s6,
                    "s8 =",
                    s8,
                    "a1 =",
                    a1,
                    "a2 =",
                    a2,
                )

        if verbose:
            if abc == False:
                print("   3-body term will not be calculated\n")
            else:
                print("   Including the Axilrod-Teller-Muto 3-body dispersion term\n")
            if intermolecular == True:
                print(
                    "   Only computing intermolecular dispersion interactions! This is not the total D3-correction\n"
                )

        for j in range(0, natom):
            ## This could be used to 'switch off' dispersion between bonded or geminal atoms ##
            scaling = False
            for k in range(j + 1, natom):
                scalefactor = 1.0

                if intermolecular == True:
                    if mols[j] == mols[k]:
                        scalefactor = 0
                        print(
                            "   --- Ignoring interaction between atoms",
                            (j + 1),
                            "and",
                            (k + 1),
                        )

                if scaling == True and hasattr(data, "BONDINDEX"):
                    if data.BONDINDEX[j][k] == 1:
                        scalefactor = 0
                    for l in range(0, natom):
                        if (
                            data.BONDINDEX[j][l] != 0
                            and data.BONDINDEX[k][l] != 0
                            and j != k
                            and data.BONDINDEX[j][k] == 0
                        ):
                            scalefactor = 0
                        for m in range(0, natom):
                            if (
                                data.BONDINDEX[j][l] != 0
                                and data.BONDINDEX[l][m] != 0
                                and data.BONDINDEX[k][m] != 0
                                and j != m
                                and k != l
                                and data.BONDINDEX[j][m] == 0
                            ):
                                scalefactor = 1 / 1.2

                if k > j:
                    ## Pythagoras in 3D to work out distance ##
                    xdist = xco[j] - xco[k]
                    ydist = yco[j] - yco[k]
                    zdist = zco[j] - zco[k]
                    totdist = (
                        math.pow(xdist, 2) + math.pow(ydist, 2) + math.pow(zdist, 2)
                    )
                    totdist = math.sqrt(totdist)

                    C6jk = getc6(c6ab, mxc, atomtype, cn, j, k)

                    ## C8 parameters depend on C6 recursively
                    atomA = E_to_index(atomtype[j])
                    atomB = E_to_index(atomtype[k])

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

                        self.attractive_r6_term = (
                            -s6
                            * C6jk
                            * damp6
                            / math.pow(dist, 6)
                            * AU_TO_KCAL
                            * scalefactor
                        )
                        self.attractive_r8_term = (
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

                        self.attractive_r6_term = (
                            -s6
                            * C6jk
                            / (math.pow(dist, 6) + damp6)
                            * AU_TO_KCAL
                            * scalefactor
                        )
                        self.attractive_r8_term = (
                            -s8
                            * C8jk
                            / (math.pow(dist, 8) + damp8)
                            * AU_TO_KCAL
                            * scalefactor
                        )

                    if pairwise == True and scalefactor != 0:
                        print(
                            "   --- Pairwise interaction between atoms",
                            (j + 1),
                            "and",
                            (k + 1),
                            ": Edisp =",
                            "%.6f"
                            % (self.attractive_r6_term + self.attractive_r8_term),
                            "kcal/mol",
                        )

                    self.attractive_r6_vdw = (
                        self.attractive_r6_vdw + self.attractive_r6_term
                    )
                    self.attractive_r8_vdw = (
                        self.attractive_r8_vdw + self.attractive_r8_term
                    )

                    jk = int(lin(k, j))
                    icomp[jk] = 1
                    cc6ab[jk] = math.sqrt(C6jk)
                    r2ab[jk] = dist ** 2
                    dmp[jk] = (1.0 / rr) ** (1.0 / 3.0)

        e63 = 0.0
        for iat in range(0, natom):
            for jat in range(0, natom):
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
                            d2 = [0] * 3
                            d2[0] = r2ab[ij]
                            d2[1] = r2ab[jk]
                            d2[2] = r2ab[ik]
                            t1 = (d2[0] + d2[1] - d2[2]) / math.sqrt(d2[0] * d2[1])
                            t2 = (d2[0] + d2[2] - d2[1]) / math.sqrt(d2[0] * d2[2])
                            t3 = (d2[2] + d2[1] - d2[0]) / math.sqrt(d2[1] * d2[2])
                            ang = 0.375 * t1 * t2 * t3 + 1.0
                            e63 = e63 + tmp * c9 * ang / (d2[0] * d2[1] * d2[2]) ** 1.50

        self.repulsive_abc_term = s6 * e63 * AU_TO_KCAL
        self.repulsive_abc = self.repulsive_abc + self.repulsive_abc_term


def main():
    # Takes arguments: (1) damping style, (2) s6, (3) rs6, (4) s8, (5) 3-body on/off, (6) input file(s)
    files = []
    verbose = True
    damp = "zero"
    s6 = 0.0
    rs6 = 0.0
    s8 = 0.0
    bj_a1 = 0.0
    bj_a2 = 0.0
    abc_term = False
    intermolecular = False
    pairwise = False
    if len(sys.argv) > 1:
        for i in range(1, len(sys.argv)):
            if sys.argv[i] == "-damp":
                damp = sys.argv[i + 1]
            elif sys.argv[i] == "-s6":
                s6 = float(sys.argv[i + 1])
            elif sys.argv[i] == "-terse":
                verbose = False
            elif sys.argv[i] == "-rs6":
                rs6 = float(sys.argv[i + 1])
            elif sys.argv[i] == "-s8":
                s8 = float(sys.argv[i + 1])
            elif sys.argv[i] == "-a1":
                bj_a1 = float(sys.argv[i + 1])
            elif sys.argv[i] == "-a2":
                bj_a2 = float(sys.argv[i + 1])
            elif sys.argv[i] == "-3body":
                abc_term = True
            elif sys.argv[i] == "-im":
                intermolecular = True
            elif sys.argv[i] == "-pw":
                pairwise = True
            else:
                if len(sys.argv[i].split(".")) > 1:
                    if (
                        sys.argv[i].split(".")[1] == "out"
                        or sys.argv[i].split(".")[1] == "log"
                        or sys.argv[i].split(".")[1] == "com"
                        or sys.argv[i].split(".")[1] == "gjf"
                        or sys.argv[i].split(".")[1] == "pdb"
                        or sys.argv[i].split(".")[1] == "txt"
                    ):
                        files.append(sys.argv[i])

    else:
        print(
            "\nWrong number of arguments used. Correct format: dftd3.py (-damp zero/bj) (-s6 val) (-rs6 val) (-s8 val) (-a1 val) (-a2 val) (-im on/off) (-pw on/off)file(s)\n"
        )
        sys.exit()

    for f in files:
        ## Use ccParse to get the Cartesian coordinates from Gaussian input/output files
        if len(f.split(".com")) > 1 or len(f.split(".gjf")) > 1:
            data = getinData(f)
        if len(f.split(".pdb")) > 1:
            data = getpdbData(f)
        if len(f.split(".out")) > 1 or len(f.split(".log")) > 1:
            data = getoutData(f)
        if len(f.split(".txt")) > 1:
            data = get_simple_data(f)
        fileD3 = calcD3(
            data,
            data.FUNCTIONAL,
            s6,
            rs6,
            s8,
            bj_a1,
            bj_a2,
            damp,
            abc_term,
            intermolecular,
            pairwise,
            verbose,
        )

        attractive_r6_vdw = fileD3.attractive_r6_vdw / AU_TO_KCAL
        attractive_r8_vdw = fileD3.attractive_r8_vdw / AU_TO_KCAL

        # Output includes 3-body term
        if abc_term == True:
            repulsive_abc = fileD3.repulsive_abc / AU_TO_KCAL
            total_vdw = attractive_r6_vdw + attractive_r8_vdw + repulsive_abc
            if verbose:
                print(
                    "\n",
                    " ".rjust(30),
                    "    D3(R6)".rjust(12),
                    "    D3(R8)".rjust(12),
                    "    D3(3-body)".rjust(12),
                    "    Total (au)".rjust(12),
                )
            print(
                "  ",
                f.ljust(30),
                "   %.8f" % attractive_r6_vdw,
                "   %.8f" % attractive_r8_vdw,
                "   %.8f" % repulsive_abc,
                "   %.8f" % total_vdw,
            )

        # Without 3-body term (default)
        else:
            total_vdw = attractive_r6_vdw + attractive_r8_vdw
            if verbose:
                print(
                    "\n",
                    " ".rjust(30),
                    "    D3(R6)".rjust(12),
                    "    D3(R8)".rjust(12),
                    "    Total (au)".rjust(12),
                )
            print(
                "  ",
                f.ljust(30),
                "   %.18f" % attractive_r6_vdw,
                "   %.18f" % attractive_r8_vdw,
                "   %.18f" % total_vdw,
            )


if __name__ == "__main__":
    main()
