import numpy as np
import scipy.optimize as opt

def root(guess, func, args, x0, phi=True):
    if phi:
        check = 2 * np.pi
    else:
        check = np.pi
    val = func(guess, *args) % (check)
    if val > check:
        val -= 2 * check
    return val - x0

def find_critical_angle(func, x0, mini, maxi, args, typ='alpha'):
    guess = mini + (maxi - mini) / 2

    if typ == 'alpha':
        angle = opt.root_scalar(root, args=(func, args, x0, True), bracket=(mini, maxi), x0=guess)
    elif typ == 'beta':
        angle = opt.root_scalar(root, args=(func, args, x0, False), bracket=(mini, maxi), x0=guess)
    else:
        raise NotImplementedError('False typ specification.')

    return angle.root
