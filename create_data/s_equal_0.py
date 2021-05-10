import numpy as np
import matplotlib.pyplot as pl

import time
import os
import configparser as cp

from ode import solver, better_eop
from plotting import plot

def helping():
    arr = np.linspace(np.pi / 4, np.pi / 2, num=25)
    arr2 = []

    for file in [f for f in os.listdir('D:/2021_04_13/1/') if os.path.isfile(os.path.join('D:/2021_04_13/1/', f))]:
        file = 'D:/2021_04_13/1/' + file
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
        alpha, beta, flag = eop.find_critical_angles(half='right', amin=np.pi / 2, amax=3 / 2 * np.pi,
                                                     bmin=np.pi / 2 + 0.3, bmax=np.pi)
        info[str(phi0)+'above'] = flag
        if flag:
            eop.save_to_csv(alpha, beta, r'D:/2021_04_13/4/')

        # below:
        alpha, beta, flag = eop.find_critical_angles(half='right', amin=np.pi, amax=3 / 2 * np.pi,
                                                     bmin=0.9, bmax=1.4)
        info[str(phi0)+'below'] = (alpha, beta)
        if flag:
            eop.save_to_csv(alpha, beta, r'D:/2021_04_13/4/')

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
        #alpha, beta, flag = eop.find_critical_angles(half='right', amin=3.6, amax=4.8,
        #                                             bmin=1.6, bmax=2.5)
        #info[str(phi0)+'above'] = flag
        #if flag:
        #    eop.save_to_csv(alpha, beta, r'D:/2021_04_13/3/')

        # below:
        alpha, beta, flag = eop.find_critical_angles(half='right', amin=2.8, amax=3.8,
                                                     bmin=0.9, bmax=1.4)
        info[str(phi0)+'below'] = flag
        if flag:
            eop.save_to_csv(alpha, beta, r'D:/2021_04_13/3/')

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
    for phi0 in np.linspace(np.pi / 2, np.pi, num=25):
        sol = solver.ODESolver(s, theta0, phi0, chi, alpha, beta, r0, stop=40)
        eop = better_eop.EmitterObserverProblem(sol, robs, thetaobs, phiobs)  #

        # above:
        alpha, beta, flag = eop.find_critical_angles(half='right', amin=4.4, amax=6,
                                                     bmin=1.6, bmax=2.2)
        info[str(phi0)+'above'] = flag
        if flag:
            eop.save_to_csv(alpha, beta, r'D:/2021_04_13/2/')

    for phi0 in np.linspace(3*np.pi / 4, np.pi, num=25):
        # below:
        sol = solver.ODESolver(s, theta0, phi0, chi, alpha, beta, r0, stop=40)
        eop = better_eop.EmitterObserverProblem(sol, robs, thetaobs, phiobs)
        alpha, beta, flag = eop.find_critical_angles(half='right', amin=2.9, amax=3.1,
                                                     bmin=1.2, bmax=1.4)
        info[str(phi0)+'below'] = flag
        if flag:
            eop.save_to_csv(alpha, beta, r'D:/2021_04_13/2/')


    for phi0 in np.linspace(np.pi / 2, 3 * np.pi / 4, num=25):
        sol = solver.ODESolver(s, theta0, phi0, chi, alpha, beta, r0, stop=40)
        eop = better_eop.EmitterObserverProblem(sol, robs, thetaobs, phiobs)
        alpha, beta, flag = eop.find_critical_angles(half='right', amin=2.95, amax=3.4,
                                                     bmin=1, bmax=1.3)
        info[str(phi0) + 'below'] = flag
        if flag:
            eop.save_to_csv(alpha, beta, r'D:/2021_04_13/2/')

    for phi0 in []: #1.832595714594046 1.8653206380689396 1.8980455615438334 1.930770485018727 1.9962203319685143 2.1598449493429825
        sol = solver.ODESolver(s, theta0, phi0, chi, alpha, beta, r0, stop=40)
        eop = better_eop.EmitterObserverProblem(sol, robs, thetaobs, phiobs)
        alpha, beta, flag = eop.find_critical_angles(half='right', amin=3.02, amax=3.04,
                                                     bmin=1.21, bmax=1.24)
        info[str(phi0) + 'below'] = flag
        if flag:
            eop.save_to_csv(alpha, beta, r'D:/2021_04_13/')

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

        # above:
    #    alpha, beta, flag = eop.find_critical_angles(half='right', amin=0, amax=np.pi,
    #                                                 bmin=1.7, bmax=2.5)
    #    info[str(phi0)+'above'] = flag
    #    if flag:
    #        eop.save_to_csv(alpha, beta, r'D:/2021_04_13/1/')

    #for phi0 in np.linspace(0, np.pi / 4, num=26):
    #    sol = solver.ODESolver(s, theta0, phi0, chi, alpha, beta, r0, stop=45)
    #    eop = better_eop.EmitterObserverProblem(sol, robs, thetaobs, phiobs)  #

        # above:
    #    alpha, beta, flag = eop.find_critical_angles(half='right', amin=4, amax=4.1,
    #                                                 bmin=1.37, bmax=1.45)
    #    info[str(phi0)+'below'] = flag
    #    if flag:
    #        eop.save_to_csv(alpha, beta, r'D:/2021_04_13/1/')

    #for phi0 in np.linspace(np.pi / 4, np.pi/2, num=25):
    #    sol = solver.ODESolver(s, theta0, phi0, chi, alpha, beta, r0, stop=45)
    #    eop = better_eop.EmitterObserverProblem(sol, robs, thetaobs, phiobs)  #
    #    # above:
    #    alpha, beta, flag = eop.find_critical_angles(half='right', amin=3.64, amax=4.05,
    #                                                 bmin=1.02, bmax=1.39)
    #    info[str(phi0)+'below'] = flag
    #    if flag:
    #        eop.save_to_csv(alpha, beta, r'D:/2021_04_13/1/')

    phis = helping()
    print(phis)
    for phi0 in [np.pi/2]:
        sol = solver.ODESolver(s, theta0, phi0, chi, alpha, beta, r0, stop=45)
        eop = better_eop.EmitterObserverProblem(sol, robs, thetaobs, phiobs)  #
        # above:
        alpha, beta, flag = eop.find_critical_angles(half='right', amin=3.65, amax=3.7,
                                                     bmin=1.02, bmax=1.06)
        info[str(phi0)+'below'] = flag
        if flag:
            eop.save_to_csv(alpha, beta, r'D:/2021_04_13/1/')

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

    s = 0

    # compute_1(s, chi, alpha, beta, r0, theta0, robs, thetaobs, phiobs)
    # compute_2(s, chi, alpha, beta, r0, theta0, robs, thetaobs, phiobs)
    # compute_3(s, chi, alpha, beta, r0, theta0, robs, thetaobs, phiobs)
    # compute_4(s, chi, alpha, beta, r0, theta0, robs, thetaobs, phiobs)
    # alpha range: 2.4367274232867313 3.8982937663544246

    #plotting(s, chi, r0, np.pi/4, 3.85, 1.08)




if __name__ == '__main__':
    main()
