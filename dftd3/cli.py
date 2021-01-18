# -*- coding: utf-8 -*-

import argparse


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
        "infiles", nargs=argparse.REMAINDER, type=str, help="input file(s)"
    )

    return cli.parse_args()
