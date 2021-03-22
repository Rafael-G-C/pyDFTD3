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

import numpy as np
import pytest
from qcelemental import periodictable as PT

from dftd3.ccParse import get_simple_data, getinData, getoutData
from dftd3.dftd3 import D3_derivatives, D3Configuration, d3, D3_element_wise
from dftd3.jax_diff import _derv_sequence
from dftd3.utils import der_order

HERE = Path(__file__).parents[1]


def _from_json(inp):
    with open(inp, "r") as j:
        data = json.load(j)

    charges = [PT.to_Z(atom) for atom in data["molecule"]["symbols"]]

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
    "damping,ref",
    [
        # reference numbers from LSDALTON
        ("zero", -0.005259455303),
        ("bj", -0.009129483944),
    ],
    ids=[
        "zero-0",
        "bj-0",
    ],
)
def test_energy(coordinates, charges, functional, damping, ref):
    config = D3Configuration(functional=functional, damp=damping)

    d3_au = d3(config, charges, *coordinates)
    assert d3_au == pytest.approx(ref, rel=1.0e-5)


@pytest.mark.parametrize(
    "damping,ref,order",
    [
        # reference numbers from LSDALTON
        (
            "zero",
            # fmt: off
            np.array(
                [
                    [ 0.000349716853,  0.000062661785, -0.000003311506],
                    [ 0.000439833817,  0.000557690553,  0.000032711342],
                    [ 0.000609561955, -0.000303544148, -0.000003631289],
                    [-0.000501863558,  0.000031792046, -0.000013434572],
                    [-0.000033137414, -0.000168498263, -0.000008039125],
                    [-0.000386522043, -0.000065823761,  0.000002757262],
                    [-0.000310805767, -0.000580921020,  0.000004944208],
                    [-0.000558542699,  0.000283457649, -0.000014030170],
                    [ 0.000506038036,  0.000017788128,  0.000021991826],
                    [-0.000114279180,  0.000165397030, -0.000019957978],
                ],
            ),
            # fmt: on
            1,
        ),
        (
            "bj",
            # fmt: off
            np.array(
                [
                    [ 0.000423636407, -0.000013304710, -0.000001829666],
                    [ 0.000327950989, -0.000203052957, -0.000007080509],
                    [ 0.000288283663,  0.000134917380,  0.000000956567],
                    [ 0.000095919442, -0.000022089153,  0.000000914397],
                    [ 0.000149535193,  0.000158727869, -0.000001062320],
                    [-0.000427595718,  0.000012621136, -0.000001757827],
                    [-0.000341842526,  0.000202703383, -0.000001443181],
                    [-0.000285631580, -0.000137238735,  0.000006350364],
                    [-0.000094892511,  0.000017099502, -0.000003292373],
                    [-0.000135363360, -0.000150383715,  0.000008244549],
                ],
            ),
            # fmt: on
            1,
        ),
    ],
    ids=[
        "zero-1",
        "bj-1",
    ],
)
def test_derivatives(damping, ref, order):
    coordinates, charges, functional = _from_json(
        HERE / "examples/formic_acid_dimer.json"
    )
    config = D3Configuration(functional=functional, damp=damping)
    slices = ((0, 0), (1, 2), (3, 0))

    d3_dervs = D3_derivatives(order, config, charges, *coordinates)
    d3_element_dervs = D3_element_wise(slices, config, charges, *coordinates)

    for element in range(len(d3_element_dervs)):
        assert d3_element_dervs[element] == pytest.approx(
            ref[slices[element]], rel=1.0e-5, abs=1.0e-8
        ), f"matching calculated element {d3_element_dervs[element]} vs {ref[slices[element]]}"

    for i, x in np.ndenumerate(d3_dervs):
        assert x == pytest.approx(
            ref[i], rel=1.0e-5, abs=1.0e-8
        ), f"Element {i} of {der_order(order)} order derivative differs from reference (Delta = {x - ref[i]})"


def test_derv_sequence():
    assert _derv_sequence((3, 2, 1, 0)) == [0, 0, 0, 1, 1, 2]
    assert _derv_sequence((0, 1, 2, 3)) == [1, 2, 2, 3, 3, 3]
    assert _derv_sequence((0, 1, 0, 1)) == [1, 3]
