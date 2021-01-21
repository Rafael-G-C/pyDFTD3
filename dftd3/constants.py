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

from qcelemental import constants

"""Global constants in the DFT-D3 model."""

MAX_ELEMENTS = 94
"""int: Maximum number of elements which have D3 parametrization."""

MAX_CONNECTIVITY = 5
"""int: Maximum connectivity."""

AU_TO_ANG = constants.conversion_factor("bohr", "angstrom")

AU_TO_KCAL = constants.conversion_factor("hartree", "kcal per mol")

C6CONV = 1.0 / (constants.conversion_factor("hartree", "J per mol") * AU_TO_ANG ** 6)
"""float: conversion factor for the C6 coefficient (J mol^-1 nm^-6) to atomic units."""

ALPHA6 = 14
"""int: exponent used in distance-dependent damping factor for R^-6 term."""

ALPHA8 = ALPHA6 + 2
"""int: exponent used in distance-dependent damping factor for R^-8 term."""

ALPHA10 = ALPHA8 + 2
"""int: exponent used in distance-dependent damping factor for R^-10 term."""
