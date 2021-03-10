#!/usr/bin/env python
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

import json
from pathlib import Path

import pytest

from dftd3.ccParse import get_simple_data, getinData, getoutData
from dftd3.dftd3 import D3Configuration, d3
from dftd3.jax_diff import D3_derivatives
from dftd3.utils import E_to_index

HERE = Path(__file__).parents[1]


def _from_json(inp):
    with open(inp, "r") as j:
        data = json.load(j)

    charges = [E_to_index(atom) for atom in data["molecule"]["symbols"]]

    # reshape coordinates as (nat, 3) and convert to angstrom
    cartesians = [coordinate for coordinate in data["molecule"]["geometry"]]

    functional = data["model"]["method"]

    return cartesians, charges, functional


def _from_txt(inp):
    data = get_simple_data(inp)
    return data.CARTESIANS, data.CHARGES, data.FUNCTIONAL


def _from_com(inp):
    data = getinData(inp)
    return data.CARTESIANS, data.CHARGES, data.FUNCTIONAL


def _from_log(inp):
    data = getoutData(inp)
    return data.CARTESIANS, data.CHARGES, data.FUNCTIONAL


# reference numbers from canonical repo
@pytest.mark.parametrize(
    "coordinates,charges,functional",
    [
        (_from_txt(HERE / "examples/formic_acid_dimer.txt")),
        (_from_com(HERE / "examples/formic_acid_dimer.com")),
        (_from_log(HERE / "examples/formic_acid_dimer.log")),
        (_from_json(HERE / "examples/formic_acid_dimer.json")),
    ],
    ids=["from_txt", "from_com", "from_log", "from_json"],
)
@pytest.mark.parametrize(
    "damping,ref,order",
    [("zero", -0.005259458983232145, 0), ("bj", -0.009129484886543590, 0)],
    ids=[
        "zero",
        "bj",
    ],
)
def test_CO2H2_2(coordinates, charges, functional, damping, ref, order):
    config = D3Configuration(functional=functional, damp=damping)
    d3_au = d3(config, charges, *coordinates)

    assert d3_au == pytest.approx(ref, rel=1.0e-5)

    d3_dervs = D3_derivatives(order, config, charges, *coordinates)
    assert d3_dervs == pytest.approx(ref, rel=1.0e-5)
