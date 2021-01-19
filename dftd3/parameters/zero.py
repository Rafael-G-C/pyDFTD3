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

"""Zero damping parameters."""

ZERO_PARMS = {
    "B1B95": [1.0000, 1.6130, 1.8680],
    "B2GPPLYP": [0.5600, 1.5860, 0.7600],
    "B3LYP": [1.0000, 1.2610, 1.7030],
    "BHLYP": [1.0000, 1.3700, 1.4420],
    "BLYP": [1.0000, 1.0940, 1.6820],
    "BP86": [1.0000, 1.1390, 1.6830],
    "BPBE": [1.0000, 1.0870, 2.0330],
    "mPWLYP": [1.0000, 1.2390, 1.0980],
    "PBE": [1.0000, 1.2170, 0.7220],
    "PBE0": [1.0000, 1.2870, 0.9280],
    "PW6B95": [1.0000, 1.5320, 0.8620],
    "PWB6K": [1.0000, 1.6600, 0.5500],
    "revPBE": [1.0000, 0.9230, 1.0100],
    "TPSS": [1.0000, 1.1660, 1.1050],
    "TPSS0": [1.0000, 1.2520, 1.2420],
    "TPSSh": [1.0000, 1.2230, 1.2190],
    "BOP": [1.0000, 0.9290, 1.9750],
    "MPW1B95": [1.0000, 1.6050, 1.1180],
    "MPWB1K": [1.0000, 1.6710, 1.0610],
    "OLYP": [1.0000, 0.8060, 1.7640],
    "OPBE": [1.0000, 0.8370, 2.0550],
    "oTPSS": [1.0000, 1.1280, 1.4940],
    "PBE38": [1.0000, 1.3330, 0.9980],
    "PBEsol": [1.0000, 1.3450, 0.6120],
    "REVSSB": [1.0000, 1.2210, 0.5600],
    "SSB": [1.0000, 1.2150, 0.6630],
    "B3PW91": [1.0000, 1.1760, 1.7750],
    "BMK": [1.0000, 1.9310, 2.1680],
    "CAMB3LYP": [1.0000, 1.3780, 1.2170],
    "LCwPBE": [1.0000, 1.3550, 1.2790],
    "M052X": [1.0000, 1.4170, 0.0000],
    "M05": [1.0000, 1.3730, 0.5950],
    "M062X": [1.0000, 1.6190, 0.0000],
    "M06HF": [1.0000, 1.4460, 0.0000],
    "M06L": [1.0000, 1.5810, 0.0000],
    "M06": [1.0000, 1.3250, 0.0000],
    "HCTH120": [1.0000, 1.2210, 1.2060],
    "B2PLYP": [0.6400, 1.4270, 1.0220],
    "DSDBLYP": [0.5000, 1.5690, 0.7050],
    "TPSS": [0.7500, 1.5410, 0.8790],
    "PWPB95": [0.8200, 1.5570, 0.7050],
    "revPBE0": [1.0000, 0.9490, 0.7920],
    "revPBE38": [1.0000, 1.0210, 0.8620],
    "rPW86PBE": [1.0000, 1.2240, 0.9010],
}
"""Published parameters (S6, RS6, S8) for zero-damping optimized for different functionals."""
