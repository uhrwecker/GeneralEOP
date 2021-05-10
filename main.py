import numpy as np
import matplotlib.pyplot as pl

import time
import os
import configparser as cp

from ode import solver, better_eop
from plotting import plot

from util.coversion import get_rphitheta as rpt

def plotting(s, chi, r0, phi0, alpha, beta, phiobs=np.pi/2, thetaobs=1.25, re=None, ax=None, c='orange'):
    if not ax:
        fig = pl.figure(figsize=(12, 12))
        ax = fig.add_subplot(111, projection='3d')

    sol = solver.ODESolver(s, np.pi/2, phi0, chi, alpha, beta, r0, stop=45)
    if re:
        r, t, p = re
        sol.r0 = r
        sol.theta0 = t
        sol.phi0 = p

    plot.plot_3d(sol, ax, phiobs, thetaobs, r0, np.pi/2, phi0, c)

    return ax

def main():
    chi = 1
    alpha = 2.5 * np.pi / 2
    beta = 2

    robs = 15
    phiobs = np.pi / 2
    thetaobs = 1.25

    r0 = 8
    theta0 = np.pi / 2

    s = 0.
    dir = r'D:/2021_04_26/phi_05_pi/s-1/'

    phi0 = [np.pi/2, np.pi/2, 3.6545057398901672, 3.6545057398901672, 5.770272220879212, 5.770272220879212]
    alpha = [-0.48979591836734687, 3.640330046290727, 2.9217201166180757, 4.398542274052478, 4.075519167693332, 2.3251172015437294]
    beta = [2.082908163265306, 1.0230225456001854, 1.310349854227405, 1.7908892128279885, 1.4096105789254476, 2.0318898670107544]

    c = ['orange'] * 6

    ax = None
    i = 0
    for p, a, b in zip(phi0, alpha, beta):
        ax = plotting(s, chi, r0, p, a, b, ax=ax, c=c[i])
        i += 1
    pl.show()



if __name__ == '__main__':
    main()
