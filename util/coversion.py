import numpy as np


def get_phi(x, y, tol):
    sig = []

    for xe, ye in zip(x, y):
        if -tol < xe < tol and -tol < ye < tol:
            sig.append(np.sign(ye) * np.pi / 2)
        elif xe > tol:
            sig.append(np.arctan(ye / xe))
        elif xe < tol:
            sig.append(np.arctan(ye / xe) + np.sign(ye) * np.pi)

    return np.array(sig)


def get_rphitheta(r0, phi0, theta0, dr, tol=1e-5, num=3):
    p = np.linspace(0, 2*np.pi, num=num, endpoint=False)
    t = np.linspace(0, np.pi, num=num)

    res = []

    for t0 in t:
        psi_x = r0 * np.cos(phi0) * np.sin(theta0) + dr * np.cos(p) * np.sin(t0)
        psi_y = r0 * np.sin(phi0) * np.sin(theta0) + dr * np.sin(p) * np.sin(t0)
        psi_z = r0 * np.cos(theta0) + dr * np.cos(t0)

        rho = np.sqrt(psi_x**2 + psi_y**2 + psi_z**2)
        theta = np.arccos(psi_z / rho)
        phi = get_phi(psi_x, psi_y, tol=tol)

        for n in range(len(phi)):
            if -1e-10 < phi[n] < 1e-10:
                phi[n] = 0.

        print(rho, theta, phi)

        res.append((rho, theta, phi % (2 * np.pi)))

    return res