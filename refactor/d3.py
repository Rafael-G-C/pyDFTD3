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


import jax.numpy as jnp
from jax.config import config

config.update("jax_enable_x64", True)

from qcelemental import constants


def d3(charges, coordinates, rs6, rs8, s6, s8, rab, rcov, r2r4, coefficients):
    bohr_to_angstrom = constants.conversion_factor("bohr", "angstrom")

    coordinates_angstrom = list(map(lambda x: x * bohr_to_angstrom, coordinates))
    coordination_numbers = compute_coordination_numbers(
        coordinates_angstrom, charges, rcov
    )

    result = 0.0
    for j, charge1 in enumerate(charges):
        for k, charge2 in enumerate(charges):
            if k > j:
                dx = coordinates[3 * j + 0] - coordinates[3 * k + 0]
                dy = coordinates[3 * j + 1] - coordinates[3 * k + 1]
                dz = coordinates[3 * j + 2] - coordinates[3 * k + 2]
                dist = jnp.sqrt(dx * dx + dy * dy + dz * dz)

                c6jk = get_c6jk(
                    coefficients[(charge1, charge2)],
                    coordination_numbers[j],
                    coordination_numbers[k],
                )

                rr = rab[(charge1, charge2)] / (dist * bohr_to_angstrom)
                r1 = (
                    -s6
                    * c6jk
                    / (jnp.power(dist, 6) * (1.0 + 6.0 * jnp.power(rs6 * rr, 14)))
                )

                c8jk = 3.0 * c6jk * r2r4[charge1] * r2r4[charge2]
                r2 = (
                    -s8
                    * c8jk
                    / (jnp.power(dist, 8) * (1.0 + 6.0 * jnp.power(rs8 * rr, 16)))
                )
                result += r1 + r2

    return result


def compute_coordination_numbers(coordinates, charges, rcov):
    coordination_numbers = []
    for i, charge1 in enumerate(charges):
        c = 0.0
        for j, charge2 in enumerate(charges):
            if i != j:
                dx = coordinates[3 * i + 0] - coordinates[3 * j + 0]
                dy = coordinates[3 * i + 1] - coordinates[3 * j + 1]
                dz = coordinates[3 * i + 2] - coordinates[3 * j + 2]
                r = jnp.sqrt(dx * dx + dy * dy + dz * dz)
                rco = (4.0 / 3.0) * (rcov[charge1] + rcov[charge2])
                c += 1.0 / (1.0 + jnp.exp(-16.0 * (rco / r - 1.0)))
        coordination_numbers.append(c)
    return coordination_numbers


def get_c6jk(coefficients, cn_j, cn_k):
    csum = 0.0
    rsum = 0.0
    for (c, cn1, cn2) in coefficients:
        r = (cn1 - cn_j) ** 2 + (cn2 - cn_k) ** 2
        t = jnp.exp(-4.0 * r)
        csum += c * t
        rsum += t
    return csum / rsum
