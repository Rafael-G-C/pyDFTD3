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

from math import sqrt

"""PBE0/def2-QZVP atomic values for multipole coefficients."""

_R2R4 = [
    8.0589,
    3.4698,
    29.0974,
    14.8517,
    11.8799,
    7.8715,
    5.5588,
    4.7566,
    3.8025,
    3.1036,
    26.1552,
    17.2304,
    17.7210,
    12.7442,
    9.5361,
    8.1652,
    6.7463,
    5.6004,
    29.2012,
    22.3934,
    19.0598,
    16.8590,
    15.4023,
    12.5589,
    13.4788,
    12.2309,
    11.2809,
    10.5569,
    10.1428,
    9.4907,
    13.4606,
    10.8544,
    8.9386,
    8.1350,
    7.1251,
    6.1971,
    30.0162,
    24.4103,
    20.3537,
    17.4780,
    13.5528,
    11.8451,
    11.0355,
    10.1997,
    9.5414,
    9.0061,
    8.6417,
    8.9975,
    14.0834,
    11.8333,
    10.0179,
    9.3844,
    8.4110,
    7.5152,
    32.7622,
    27.5708,
    23.1671,
    21.6003,
    20.9615,
    20.4562,
    20.1010,
    19.7475,
    19.4828,
    15.6013,
    19.2362,
    17.4717,
    17.8321,
    17.4237,
    17.1954,
    17.1631,
    14.5716,
    15.8758,
    13.8989,
    12.4834,
    11.4421,
    10.2671,
    8.3549,
    7.8496,
    7.3278,
    7.4820,
    13.5124,
    11.6554,
    10.0959,
    9.7340,
    8.8584,
    8.0125,
    29.8135,
    26.3157,
    19.1885,
    15.8542,
    16.1305,
    15.6161,
    15.1226,
    16.1576,
]

R2R4 = [sqrt(0.5 * x * sqrt(i + 1)) for i, x in enumerate(_R2R4)]
