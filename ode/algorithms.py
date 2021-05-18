import numpy as np
import scipy.optimize as opt


def root(guess, func, args, x0, phi=True):
    """
        (private)
        Support function for the critical angle algorithm
        :param guess: float; (mean) guess of the critical angle
        :param func: function; takes guess and the args as arguments
        :param args: iterable; additional arguments for the func
        :param x0: float; value at which to compare the calculated value to
        :param phi: bool; if True, the calculation will be handled as a phi parameter (periodic in 2pi).
        :return: float; difference of the calculated value and the comparison value.
    """
    # check if the value is pi- or 2pi periodic
    if phi:
        check = 2 * np.pi
    else:
        check = np.pi

    val = func(guess, *args) % check

    if val > check:
        val -= 2 * check

    return val - x0


def find_critical_angle(func, x0, mini, maxi, args, typ='alpha'):
    """
    Wrapper for the critical-angle-finding-algorithm.
    :param func: function; function of which to calculate the critical angle from.
    :param x0: float; comparison value
    :param mini: float; lower bound of the angle
    :param maxi: float; upper bound of the angle
    :param args: iterable; list of additional arguments for the function func
    :param typ: ['alpha', 'beta']; handles if the critical angle should be phi-like or theta-like.
    :return: float; the critical angle / root of the func
    """
    guess = mini + (maxi - mini) / 2

    if typ == 'alpha':
        angle = opt.root_scalar(root, args=(func, args, x0, True), bracket=(mini, maxi), x0=guess)
    elif typ == 'beta':
        angle = opt.root_scalar(root, args=(func, args, x0, False), bracket=(mini, maxi), x0=guess)
    else:
        raise NotImplementedError('False typ specification.')

    return angle.root
