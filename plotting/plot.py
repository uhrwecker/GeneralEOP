import numpy as np

def plot_3d(sol, ax, phi, theta, r0, t0, p0, c):
    sigma, data = sol.solve()

    robs = 15

    r = data[:, 2]
    t = data[:, -4]
    p = data[:, -2]

    #res = -(1 - 2 / r) * data[:, 1]**2 + 1/(1-2/r)*data[:, 3]**2 + r**2 * data[:, 5]**2 + r**2*np.sin(t)**2*data[:, 7]**2


    ax.plot(r * np.cos(p) * np.sin(t),
            r * np.sin(p) * np.sin(t),
            r * np.cos(t), color=c)

    # what I want:
    # r = 8, p = 1, theta = 0.5

    u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    x = 2 * np.cos(u) * np.sin(v)
    y = 2 * np.sin(u) * np.sin(v)
    z = 2 * np.cos(v)
    ax.plot_wireframe(x, y, z, color="r")


    x0 = robs * np.cos(phi) * np.sin(theta)
    y0 = robs * np.sin(phi) * np.sin(theta)
    z0 = robs * np.cos(theta)


    for a in [0]:
        ax.scatter([x0], [y0 * np.cos(a) - z0 * np.sin(a)], [y0 * np.sin(a) + z0 * np.cos(a)], color='blue')

    ax.scatter([r0 * np.cos(p0) * np.sin(t0)],
               [r0 * np.sin(p0) * np.sin(t0)],
               [r0 * np.cos(t0)], color='orange')

    phii = np.linspace(0, 2*np.pi, num=1000)
    ax.plot(r0 * np.cos(phii), r0 * np.sin(phii), 0, color='black', alpha=0.4, lw=1)


    ax.set_xlim(-15, 15)
    ax.set_ylim(-15, 15)
    ax.set_zlim(-10, 10)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')