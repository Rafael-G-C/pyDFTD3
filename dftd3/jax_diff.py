from jax import grad
from jax.config import config

config.update("jax_enable_x64", True)


def _derv_sequence(orders):
    sequence = []
    for variable, variable_order in enumerate(orders):
        if variable_order > 0:
            sequence += variable_order * [variable]
    return sequence


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
