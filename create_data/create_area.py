import numpy as np
import matplotlib.pyplot as pl

import time
import os
import configparser as cp

from ode import solver, better_eop
from plotting import plot

from util.coversion import get_rphitheta as rpt

def plotting(s, chi, r0, phi0, alpha, beta, phiobs=np.pi/2, thetaobs=1.25, re=None):
    r, t, p = re
    fig = pl.figure()
    ax = fig.add_subplot(111, projection='3d')

    sol = solver.ODESolver(s, np.pi/2, phi0, chi, alpha, beta, r0, stop=45)
    sol.r0 = r
    sol.theta0 = t
    sol.phi0 = p

    plot.plot_3d(sol, ax, phiobs, thetaobs, r0, np.pi/2, phi0)
    pl.show()

def main():
    chi = 1
    alpha = 2.5 * np.pi / 2
    beta = 2

    robs = 15
    phiobs = np.pi / 2
    thetaobs = 1.25

    r0 = 8
    theta0 = np.pi / 2
    phi0 = np.pi

    s = -1.
    dir = r'D:/2021_04_26/phi_pi/s-1/'
    dr = 0.2

    res = rpt(r0, phi0, theta0, dr, num=15)

    info = {}

    start = time.time()

    #plotting(s, chi, r0, phi0, alpha, beta, phiobs=np.pi/2, thetaobs=1.25, re=res[0][0])

    for item in res:
        rr, tt, pp = item
        for r, t, p in zip(rr, tt, pp):
            print(r, t, p)
            sol = solver.ODESolver(s, theta0, phi0, chi, alpha, beta, r0, stop=45)

            sol.r0 = r
            sol.theta0 = t
            sol.phi0 = p

            eop = better_eop.EmitterObserverProblem(sol, robs, thetaobs, phiobs)

            if not os.path.isfile(dir + '%.8f_%.8f_%.8f.csv' % (r, t, p)):
                alpha, beta, flag = eop.find_critical_angles(half='right', bmax=1.7, bmin=1.9,
                                                             amax=4.8, amin=4.3)
                info[f'{r}_{p}_{t}'] = flag
                if flag:
                    eop.save_to_csv(alpha, beta, dir)

            else:
                print('Already existing.')

    ti = time.time() - start
    print(f'Took {ti}s ({ti / 60}m)')



if __name__ == '__main__':
    main()
