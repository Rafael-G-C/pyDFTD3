# -*- coding: utf-8 -*-

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
