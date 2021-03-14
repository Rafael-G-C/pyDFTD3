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

import argparse
from pathlib import Path


def cli():
    cli = argparse.ArgumentParser(description="Front-end CLI for pyDFTD3")
    cli.add_argument("--verbose", "-v", action="count", default=1)
    cli.add_argument("--version", action="version", version="0.1.0")
    cli.add_argument(
        "--damp",
        action="store",
        default="zero",
        type=str,
        help="damping method: zero (default) or Becke-Johnson",
    )
    cli.add_argument(
        "--s6",
        action="store",
        default="0.0",
        type=float,
        help="document me",
    )
    cli.add_argument(
        "--rs6",
        action="store",
        default="0.0",
        type=float,
        help="document me",
    )
    cli.add_argument(
        "--s8",
        action="store",
        default="0.0",
        type=float,
        help="document me",
    )
    cli.add_argument(
        "--rs8",
        action="store",
        default="0.0",
        type=float,
        help="document me",
    )
    cli.add_argument(
        "--a1",
        action="store",
        default="0.0",
        type=float,
        help="document me",
    )
    cli.add_argument(
        "--a2",
        action="store",
        default="0.0",
        type=float,
        help="document me",
    )
    cli.add_argument(
        "--order",
        action="store",
        default="0",
        type=int,
        help="derivative order to compute",
    )
    cli.add_argument(
        "--three",
        action="store_true",
        help="include Axilrod-Teller-Muto 3-body term",
    )
    cli.add_argument(
        "--inter",
        action="store_true",
        help="document me",
    )
    cli.add_argument(
        "--pairwise",
        action="store_true",
        help="document me",
    )
    cli.add_argument(
        "infiles", nargs=argparse.REMAINDER, type=Path, help="input file(s)"
    )

    return cli.parse_args()
