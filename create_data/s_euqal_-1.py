import numpy as np
import matplotlib.pyplot as pl

import time
import os
import configparser as cp

from ode import solver, better_eop
from plotting import plot

def helping(dir, start, end, num):
    arr = np.linspace(start, end, num=num)
    arr2 = []

    for file in [f for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f))]:
        file = dir + file
        if file.endswith('.ini'):
            config = cp.ConfigParser()
            config.read(file)

            phi0 = float(config['Photon']['phi0'])
            arr2.append(phi0)

    arr2 = np.array(arr2)
    arr3 = []
    for phi in arr:
        if not phi in arr2:
            arr3.append(phi)

    return np.array(arr3)

def compute_4(s, chi, alpha, beta, r0, theta0, robs, thetaobs, phiobs):
    start = time.time()

    info = {}
    for phi0 in np.linspace(3 / 2 * np.pi, 2 * np.pi, num=50):
        sol = solver.ODESolver(s, theta0, phi0, chi, alpha, beta, r0, stop=45)
        eop = better_eop.EmitterObserverProblem(sol, robs, thetaobs, phiobs)  #

        # above:
        #alpha, beta, flag = eop.find_critical_angles(half='right', amin=np.pi / 2, amax=1.5 * np.pi,
        #                                             bmin=np.pi / 2 + 0.3, bmax=np.pi)
        #info[str(phi0)+'above'] = flag
        #if flag:
        #    eop.save_to_csv(alpha, beta, r'D:/2021_04_13/s-1/4/')

        # below:
        alpha, beta, flag = eop.find_critical_angles(half='right', amin=np.pi, amax=4.89,
                                                     bmin=0.85, bmax=1.5)
        info[str(phi0)+'below'] = (alpha, beta, flag)
        if flag:
            eop.save_to_csv(alpha, beta, r'D:/2021_04_13/s-1/4/')

    now = time.time()
    print(f'All those took {now - start}s (or {(now - start) / 60}m).')
    print('Info:')
    print(info)


def compute_3(s, chi, alpha, beta, r0, theta0, robs, thetaobs, phiobs):
    start = time.time()

    info = {}
    for phi0 in np.linspace(np.pi, 3 * np.pi / 2, num=50):
        sol = solver.ODESolver(s, theta0, phi0, chi, alpha, beta, r0, stop=45)
        eop = better_eop.EmitterObserverProblem(sol, robs, thetaobs, phiobs)  #

        # above:
        #alpha, beta, flag = eop.find_critical_angles(half='right', amin=3.6, amax=5,
        #                                             bmin=1.6, bmax=2.4)
        #info[str(phi0)+'above'] = flag
        #if flag:
        #    eop.save_to_csv(alpha, beta, r'D:/2021_04_13/s-1/3/')

        # below:
        alpha, beta, flag = eop.find_critical_angles(half='right', amin=2.8, amax=3.8,
                                                     bmin=0.9, bmax=1.4)
        info[str(phi0)+'below'] = flag
        if flag:
            eop.save_to_csv(alpha, beta, r'D:/2021_04_13/s-1/3/')

    now = time.time()
    print(f'All those took {now - start}s (or {(now - start) / 60}m).')
    print('Info:')
    print(info)
    if False in info.values():
        print('THERE ARE ERRORS!')


def compute_2(s, chi, alpha, beta, r0, theta0, robs, thetaobs, phiobs):
    start = time.time()

    info = {}
    i = 0
    #for phi0 in np.linspace(np.pi / 2, np.pi, num=50):
    #    sol = solver.ODESolver(s, theta0, phi0, chi, alpha, beta, r0, stop=40)
    #    eop = better_eop.EmitterObserverProblem(sol, robs, thetaobs, phiobs)  ###
    #    # above:
    #    alpha, beta, flag = eop.find_critical_angles(half='right', amin=4.4, amax=6,
    #                                                 bmin=1.6, bmax=2.2)
    #    info[str(phi0)+'above'] = flag
    #    if flag:
    #        eop.save_to_csv(alpha, beta, r'D:/2021_04_13/s-1/2/')

    #for phi0 in np.linspace(3*np.pi / 4, np.pi, num=26):
    #    # below:
    #    sol = solver.ODESolver(s, theta0, phi0, chi, alpha, beta, r0, stop=40)
    #    eop = better_eop.EmitterObserverProblem(sol, robs, thetaobs, phiobs)
    ##    alpha, beta, flag = eop.find_critical_angles(half='right', amin=2.95, amax=3.05,
    #                                                 bmin=1.2, bmax=1.4)
    #    info[str(phi0)+'below'] = flag
    #    if flag:
    #        eop.save_to_csv(alpha, beta, r'D:/2021_04_13/s-1/2/')

    phis = helping(r'D:/2021_04_13/s-1/2/', np.pi/2, 3 * np.pi / 4, 25)[::-1]
    print(phis)
    for phi0 in [np.pi/2]:#phis:
        sol = solver.ODESolver(s, theta0, phi0, chi, alpha, beta, r0, stop=40)
        eop = better_eop.EmitterObserverProblem(sol, robs, thetaobs, phiobs)
        alpha, beta, flag = eop.find_critical_angles(half='right', amin=3.615560937897893, amax=3.86542,
                                                     bmin=1.019541796675809, bmax=1.03574537)
        info[str(phi0) + 'below'] = flag
        if flag:
            eop.save_to_csv(alpha, beta, r'D:/2021_04_13/s-1/2/')

        #alpha, beta, flag = eop.find_critical_angles(half='right', amin=3.4439712973205476, amax=3.76542,
        #                                             bmin=0.934534587, bmax=1.0170462174800967)

    now = time.time()
    print(f'All those took {now - start}s (or {(now - start) / 60}m).')
    print('Info:')
    print(info)
    if False in info.values():
        print('THERE ARE ERRORS!')


def compute_1(s, chi, alpha, beta, r0, theta0, robs, thetaobs, phiobs):
    start = time.time()

    info = {}
    #for phi0 in np.linspace(0, np.pi / 2, num=50):
    #    sol = solver.ODESolver(s, theta0, phi0, chi, alpha, beta, r0, stop=35)
    #    eop = better_eop.EmitterObserverProblem(sol, robs, thetaobs, phiobs)  #
    #    # above:
    #    alpha, beta, flag = eop.find_critical_angles(half='right', amin=-0.5, amax=2,
    #                                                 bmin=2, bmax=2.5)
    #    info[str(phi0)+'above'] = flag
    #    if flag:
    #        eop.save_to_csv(alpha, beta, r'D:/2021_04_13/s-1/1/')

    #for phi0 in np.linspace(0, np.pi / 4, num=26):
    #    sol = solver.ODESolver(s, theta0, phi0, chi, alpha, beta, r0, stop=45)
    #    eop = better_eop.EmitterObserverProblem(sol, robs, thetaobs, phiobs)  #
    #    # above:
    #    alpha, beta, flag = eop.find_critical_angles(half='right', amin=4, amax=4.1,
    #                                                 bmin=1.37, bmax=1.45)
    #    info[str(phi0)+'below'] = flag
    #    if flag:
    #        eop.save_to_csv(alpha, beta, r'D:/2021_04_13/s-1/1/')

    phis = helping(r'D:/2021_04_13/s-1/1/', np.pi / 4, np.pi / 2, 25)
    print(phis)
    for phi0 in phis:
        sol = solver.ODESolver(s, theta0, phi0, chi, alpha, beta, r0, stop=45)
        eop = better_eop.EmitterObserverProblem(sol, robs, thetaobs, phiobs)  #
        # above:
        alpha, beta, flag = eop.find_critical_angles(half='right', amax=3.8055248913216775, amin=3.6766988799990963,
                                                     bmax=1.0646566483676932, bmin=1.0321596533292083)
        info[str(phi0)+'below'] = flag
        if flag:
            eop.save_to_csv(alpha, beta, r'D:/2021_04_13/s-1/1/')#

    now = time.time()
    print(f'All those took {now - start}s (or {(now - start) / 60}m).')
    print('Info:')
    print(info)
    if False in info.values():
        print('THERE ARE ERRORS!')


def plotting(s, chi, r0, phi0, alpha, beta, phiobs=np.pi/2, thetaobs=1.25):
    fig = pl.figure()
    ax = fig.add_subplot(111, projection='3d')
    sol = solver.ODESolver(s, np.pi/2, phi0, chi, alpha, beta, r0, stop=45)
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

    s = -1

    compute_1(s, chi, alpha, beta, r0, theta0, robs, thetaobs, phiobs)
    # compute_2(s, chi, alpha, beta, r0, theta0, robs, thetaobs, phiobs)
    # compute_3(s, chi, alpha, beta, r0, theta0, robs, thetaobs, phiobs)
    # compute_4(s, chi, alpha, beta, r0, theta0, robs, thetaobs, phiobs)
    # alpha range: 2.4367274232867313 3.8982937663544246

    # a = (3.7099747776643826 + 3.6100005299999998) / 2
    # b = (1.043015681586237 + 1.0179363187500001) / 2 - 0.00253

    # plotting(s, chi, r0, np.pi/2, a, b)
    #3.7531481481481483, 1.0642469135802468




if __name__ == '__main__':
    main()
