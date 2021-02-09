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

###############################################################
#                         ccParse.py                          #
#                Reads compchem job file(s)                   #
###############################################################

import os
import sys

from qcelemental import periodictable


## Check for integer when parsing ##
def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def elementID(massno):
    return periodictable.to_symbol(massno)


# Read Cartesian coordinate data from a PDB file
class getpdbData:
    def __init__(self, file):
        if not os.path.exists(file):
            print("\nFATAL ERROR: Input file [ %s ] does not exist" % file)
            sys.exit()

        def getATOMS(self, inlines):
            self.ATOMTYPES = []
            self.CARTESIANS = []
            for i in range(0, len(inlines)):
                if inlines[i].find("ATOM") > -1 or inlines[i].find("HETATM") > -1:
                    self.ATOMTYPES.append((inlines[i].split()[2]))
                    self.CARTESIANS += [
                        float(coordinate) for coordinate in inlines[i].split()[4:7]
                    ]

        infile = open(file, "r")
        inlines = infile.readlines()
        getATOMS(self, inlines)
        self.FUNCTIONAL = None


# Read Cartesian data from a Gaussian formatted input file (*.com or *.gjf)
class getinData:
    def __init__(self, file):
        if not os.path.exists(file):
            print("\nFATAL ERROR: Input file [ %s ] does not exist" % file)
            sys.exit()

        def getATOMTYPES(self, inlines):
            self.ATOMTYPES = []
            for i in range(0, len(inlines)):
                if inlines[i].find("#") > -1:
                    if len(inlines[i + 1].split()) == 0:
                        start = i + 5
                    if len(inlines[i + 2].split()) == 0:
                        start = i + 6
                    break
            for i in range(start, len(inlines)):
                if len(inlines[i].split()) == 0:
                    break
                else:
                    self.ATOMTYPES.append(inlines[i].split()[0])

        def getCARTESIANS(self, inlines, natoms):
            self.CARTESIANS = []
            for i in range(0, len(inlines)):
                if inlines[i].find("#") > -1:
                    start = i + 5
                    break

            for i in range(start, len(inlines)):
                if len(inlines[i].split()) == 0:
                    break
                elif len(inlines[i].split()) == 4:
                    self.CARTESIANS += [
                        float(coordinate) for coordinate in inlines[i].split()[1:4]
                    ]

        def getMETHOD(self, inlines):
            self.FUNCTIONAL = None
            # looks for a selected group of methods (some of my favourites...)
            for i in range(0, len(inlines)):
                if inlines[i].find("#") > -1:
                    if inlines[i].find("B3LYP") > -1 or inlines[i].find("b3lyp") > -1:
                        self.FUNCTIONAL = "B3LYP"
                    if (
                        inlines[i].find("TPSSTPSS") > -1
                        or inlines[i].find("tpsstpss") > -1
                    ):
                        self.FUNCTIONAL = "TPSSTPSS"

        def getBONDINDEX(self, inlines, natoms):
            conn = []
            connectivity = 0

            for j in range(0, len(inlines)):
                if " 1 " in inlines[j]:
                    startconn = j
                    connectivity = 1
                    break

            if connectivity == 1:
                for j in range(startconn, len(inlines)):
                    conn.append(inlines[j])

            self.BONDINDEX = []
            for j in range(0, natoms):
                self.BONDINDEX.append([0])
                for k in range(0, natoms):
                    self.BONDINDEX[j].append(0)

            for j in range(0, natoms):
                if connectivity == 1:
                    for bonded in conn[j].split():
                        if is_number(bonded) == True:
                            if int(bonded) - 1 != j:
                                self.BONDINDEX[j][int(bonded) - 1] = 1
                                self.BONDINDEX[int(bonded) - 1][j] = 1

        infile = open(file, "r")
        inlines = infile.readlines()
        getATOMTYPES(self, inlines)
        self.NATOMS = len(self.ATOMTYPES)
        getMETHOD(self, inlines)
        getBONDINDEX(self, inlines, self.NATOMS)
        getCARTESIANS(self, inlines, self.NATOMS)


# Read Cartesian data from a Gaussian formatted output file (*.log or *.out)
class getoutData:
    def __init__(self, file):
        if not os.path.exists(file):
            print("\nFATAL ERROR: Output file [ %s ] does not exist" % file)
            sys.exit()

        def getATOMTYPES(self, outlines):
            self.ATOMTYPES = []
            self.CARTESIANS = []
            self.ATOMICTYPES = []

            anharmonic_geom = 0
            for i in range(0, len(outlines)):
                if outlines[i].find("Input orientation") > -1:
                    standor = i
                if outlines[i].find("Standard orientation") > -1:
                    standor = i
                if outlines[i].find("Vib.Av.Geom.") > -1:
                    standor = i
                    anharmonic_geom = 1
                if (
                    outlines[i].find("Distance matrix") > -1
                    or outlines[i].find("Rotational constants") > -1
                ):
                    if outlines[i - 1].find("-------") > -1:
                        self.NATOMS = i - standor - 6

            try:
                standor
            except NameError:
                pass
            else:
                for i in range(standor + 5, standor + 5 + self.NATOMS):
                    self.ATOMTYPES.append(elementID(int(outlines[i].split()[1])))
                    self.ATOMICTYPES.append(int(outlines[i].split()[2]))

                    if anharmonic_geom == 0:
                        if len(outlines[i].split()) > 5:
                            self.CARTESIANS += [
                                float(coordinate)
                                for coordinate in outlines[i].split()[3:6]
                            ]
                        else:
                            self.CARTESIANS += [
                                float(coordinate)
                                for coordinate in outlines[i].split()[2:5]
                            ]
                    if anharmonic_geom == 1:
                        self.CARTESIANS += [
                            float(coordinate) for coordinate in outlines[i].split()[2:5]
                        ]

        def getMETHOD(self, outlines):
            self.FUNCTIONAL = None
            for i in range(0, len(outlines)):
                if outlines[i].strip().find("\GINC") > -1:
                    if len(outlines[i].strip().split("\\")) > 5:
                        func = outlines[i].strip().split("\\")[4]
                        self.FUNCTIONAL = func[1:]

        if os.path.exists(file):
            outfile = open(file, "r")
        outlines = outfile.readlines()
        getMETHOD(self, outlines)
        getATOMTYPES(self, outlines)


class get_simple_data:
    def __init__(self, file):
        self.ATOMTYPES = []
        self.NATOMS = 0
        self.FUNCTIONAL = None
        self.CARTESIANS = []

        def info_getter(lines):
            get_geom = False
            for line in lines:
                if get_geom == True:
                    geometry_atoms = line.split()
                    self.ATOMTYPES.append(geometry_atoms[0])
                    self.CARTESIANS += [
                        float(coordinate) for coordinate in geometry_atoms[1:4]
                    ]
                    self.NATOMS += 1
                elif "geometry" in line:
                    get_geom = True
                    continue

                elif "functional" in line:
                    functional_line = line.split(":")
                    self.FUNCTIONAL = functional_line[1].strip().replace("\n", "")
                else:
                    continue

        with open(file) as simple_data:
            lines = simple_data.readlines()
            info_getter(lines)
