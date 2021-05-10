import numpy as np
import scipy.optimize as op


def equations(p, s=1e-6, r=6):
    u = 1/r
    E, L = p

    # first:
    first = E ** 2 * (1 - s ** 2 * u ** 2 + 2 * s ** 2 * u ** 3) + E * (
            2 * s * u ** 2 * L - 6 * s * u ** 3 * L) - 1 + 2 * u + 2 * s ** 2 * u ** 3 - 4 * s ** 2 * u ** 4 \
            - s ** 4 * u ** 6 + 2 * s ** 4 * u ** 7 - u ** 2 * L ** 2 + 2 * u ** 3 * L ** 2 + s ** 2 * u ** 6 * L ** 2

    second = E ** 2 * (- 2 * s ** 2 * u + 6 * s ** 2 * u ** 2) + E * (
            4 * s * u * L - 18 * s * u ** 2 * L) + 2 + 6 * s ** 2 * u ** 2 - 16 * s ** 2 * u ** 3 \
             - 6 * s ** 4 * u ** 5 + 14 * s ** 7 * u ** 6 - 2 * u * L ** 2 + 6 * u ** 2 * L ** 2 \
             + 6 * s ** 2 * u ** 5 * L ** 2

    return first, second


def get_com(s, r0, rotation='positive'):
    factor = 1
    if rotation == 'negative':
        factor *= -1
    ig = (np.sqrt(8/9), factor*2*np.sqrt(3))
    res = op.fsolve(equations, ig, args=(s, r0))

    return res[0], r0, res[1]