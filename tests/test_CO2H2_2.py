#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path

import pytest

from dftd3.ccParse import get_simple_data
from dftd3.dftd3 import autokcal, calcD3


def test_CO2H2_2():
    here = Path(__file__).parents[1]
    inp_file = here / "examples/formic_acid_dimer.txt"

    data = get_simple_data(inp_file)
    d3out = calcD3(
        fileData=data,
        functional=data.FUNCTIONAL,
        s6=0.0,
        rs6=0.0,
        s8=0.0,
        a1=0.0,
        a2=0.0,
        damp="zero",
        abc=False,
        intermolecular=False,
        pairwise=False,
        verbose=False,
    )

    d3_au = (d3out.attractive_r6_vdw + d3out.attractive_r8_vdw) / autokcal

    # reference result from Psi4
    # Empirical Dispersion Energy =          -0.0052594600000000
    assert d3_au == pytest.approx(-0.0052594600000000, rel=1.0e-5)
