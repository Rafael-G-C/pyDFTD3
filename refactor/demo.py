import d3
import pytest
from jax import jacfwd
import time

from data import get_r2r4, get_rcov, get_rab, read_coefficients


coordinates = [
    4.01376695,
    -0.42094027,
    0.00413283,
    2.37659705,
    -2.06988695,
    -0.07118787,
    3.23930049,
    2.02346393,
    0.03394893,
    6.09063455,
    -0.5931548,
    0.05363799,
    1.40721479,
    2.12118734,
    -0.0108338,
    -4.04214497,
    0.34532288,
    -0.03309288,
    -2.55773675,
    2.13450991,
    -0.01370996,
    -3.05560778,
    -2.01961834,
    0.06482139,
    -6.12485749,
    0.3304205,
    -0.12008265,
    -1.22314224,
    -1.95250094,
    0.14005516,
]

charges = [6, 8, 8, 1, 1, 6, 8, 8, 1, 1]


# moved these outside so that the code does not have to know anything about
# functionals and functionals are now the problem of the calling code not the
# problem of the d3 code
rs6 = 1.261
rs8 = 1.0
s6 = 1.0
s8 = 1.703

rab = get_rab()
rcov = get_rcov()
r2r4 = get_r2r4()
coefficients = read_coefficients()

d3_au = d3.d3(charges, coordinates, rs6, rs8, s6, s8, rab, rcov, r2r4, coefficients)
assert d3_au == pytest.approx(-0.0052594552412488435)

f = lambda coordinates: d3.d3(
    charges, coordinates, rs6, rs8, s6, s8, rab, rcov, r2r4, coefficients
)
d = jacfwd  # one could also experiment with jacrev

# compute cubic tensor
t0 = time.perf_counter()
T = d(d(d(f)))(coordinates)
t1 = time.perf_counter()
print(T)
print(f"elapsed: {t1 - t0:0.4f} seconds")
