# Copied from __main__.py in pip
from __future__ import absolute_import

import os
import sys

# If we are running from a wheel, add the wheel to sys.path
# This allows the usage python pip-*.whl/pip install pip-*.whl

if __package__ == '':
    path = os.path.dirname(os.path.dirname(__file__))
    sys.path.insert(0, path)

import dftd3

if __name__ == '__main__':
    sys.argv.append("/home/kilimanjaro/Documents/Nor/pyDFTD3/examples/formic_acid_dimer.com")
    print(sys.argv)
    sys.exit(dftd3.main())
