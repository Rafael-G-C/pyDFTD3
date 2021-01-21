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

"""Becke-Johnson damping parameters."""

BJ_PARMS = {
    "B2PLYP": [0.64, 0.3065, 0.9147, 5.057],
    "B3LYP": [1, 0.3981, 1.9889, 4.4211],
    "B97D": [1, 0.5545, 2.2609, 3.2297],
    "BLYP": [1, 0.4298, 2.6996, 4.2359],
    "BP86": [1, 0.3946, 3.2822, 4.8516],
    "DSDBLYP": [0.5, 0, 0.213, 6.0519],
    "PBE0": [1, 0.4145, 1.2177, 4.8593],
    "PBE": [1, 0.4289, 0.7875, 4.4407],
    "PW6B95": [1, 0.2076, 0.7257, 6.375],
    "PWPB95": [0.82, 0, 0.2904, 7.3141],
    "revPBE0": [1, 0.4679, 1.7588, 3.7619],
    "revPBE38": [1, 0.4309, 1.476, 3.9446],
    "revPBE": [1, 0.5238, 2.355, 3.5016],
    "rPW86PBE": [1, 0.4613, 1.3845, 4.5062],
    "TPSS0": [1, 0.3768, 1.2576, 4.5865],
    "TPSS": [1, 0.4535, 1.9435, 4.4752],
    "B1B95": [1, 0.2092, 1.4507, 5.5545],
    "B2GPPLYP": [0.56, 0, 0.2597, 6.3332],
    "B3PW91": [1, 0.4312, 2.8524, 4.4693],
    "BHLYP": [1, 0.2793, 1.0354, 4.9615],
    "BMK": [1, 0.194, 2.086, 5.9197],
    "BOP": [1, 0.487, 3.295, 3.5043],
    "BPBE": [1, 0.4567, 4.0728, 4.3908],
    "CAMB3LYP": [1, 0.3708, 2.0674, 5.4743],
    "LCwPBE": [1, 0.3919, 1.8541, 5.0897],
    "MPW1B95": [1, 0.1955, 1.0508, 6.4177],
    "MPWB1K": [1, 0.1474, 0.9499, 6.6223],
    "mPWLYP": [1, 0.4831, 2.0077, 4.5323],
    "OLYP": [1, 0.5299, 2.6205, 2.8065],
    "OPBE": [1, 0.5512, 3.3816, 2.9444],
    "oTPSS": [1, 0.4634, 2.7495, 4.3153],
    "PBE38": [1, 0.3995, 1.4623, 5.1405],
    "PBEsol": [1, 0.4466, 2.9491, 6.1742],
    "PTPSS": [0.75, 0, 0.2804, 6.5745],
    "PWB6K": [1, 0.1805, 0.9383, 7.7627],
    "revSSB": [1, 0.472, 0.4389, 4.0986],
    "SSB": [1, -0.0952, -0.1744, 5.217],
    "TPSSh": [1, 0.4529, 2.2382, 4.655],
    "HCTH120": [1, 0.3563, 1.0821, 4.3359],
}
"""Published parameters (S6, S8, a1, a2) for Becke-Johnson-damping optimized for different functionals."""
