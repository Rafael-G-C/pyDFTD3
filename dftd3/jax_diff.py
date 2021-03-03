from dftd3.constants import AU_TO_KCAL
from .dftd3 import d3
from itertools import product
from jax.config import config
import numpy as np

config.update("jax_enable_x64", True)

import jax.numpy as jnp
from jax import grad, jacfwd, jacrev


def _derv_sequence(orders):
    sequence = []
    for variable, variable_order in enumerate(orders):
        if variable_order > 0:
            sequence += variable_order * [variable]
    return sequence


def test_derv_sequence():
    assert _derv_sequence((3, 2, 1, 0)) == [0, 0, 0, 1, 1, 2]
    assert _derv_sequence((0, 1, 2, 3)) == [1, 2, 2, 3, 3, 3]
    assert _derv_sequence((0, 1, 0, 1)) == [1, 3]


def derv(fun, variables, orders) -> float:
    """
    fun: function to differentiate which expects a certain number of variables
    variables: list of variables at which to differentiate the function
    orders: [1, 0, 2, 0] means differentate with respect to variable 1 once,
                         and differentiate with respect to variable 3 twice.
    """
    sequence = _derv_sequence(orders)
    functions = [fun]
    for i, order in enumerate(sequence):
        functions.append(grad(functions[i], (order)))
    return functions[-1](*variables)

def distribute(indices, num_variables):
    l = [0 for _ in range(num_variables)]
    for index in indices:
        l[index] += 1
    return l


# this will generate the derivatives layer
def D3_derivatives(charges, functional, damp, order, coordinates):
    """
    1 --> gradient; 2 --> hessian; 3 --> 3rd derivatives; etc...
    """
    dervs = []
    natoms = len(charges)
    num_variables = 3 * natoms

    combo = product(range(num_variables), repeat=order)
    derivative_orders = map(lambda x: distribute(x, num_variables), combo)
    derivative_orders = list(derivative_orders)[0:1]
    for d_order in derivative_orders:
        dervs.append(AU_TO_KCAL *
            derv(d3, [charges, *coordinates], [0] + d_order)
        )

    #dervs = np.array(dervs).reshape((natoms, 3) * order)
    return dervs