""" Test the functions from project3.py. """
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

def unittest3_1d():
    """ Here we test all the 1d routines. """
    from project3 import jacobi_step_1d, w_cycle_step_1d, full_mg_1d

    def f_1d(x): return (2*x - 1)**3

    def jacobi_1d(uh, fh, tol=1e-8, kmax=100):
        for k in range(kmax):
            if jacobi_step_1d(uh, fh, 0.5) < tol: return k + 1
        return -1

    def mg_w_cycle_1d(uh, fh, alpha1, alpha2, tol=1e-8, kmax=100):
        for k in range(kmax):
            if w_cycle_step_1d(uh, fh, 2/3, alpha1, alpha2) < tol: return k + 1
        return -1

    def w_i(i, x): return np.sin(i * x * np.pi)

    kmax, n, graphics = 50, 1024, False
    h = 2/n
    x = np.linspace(0, 1, n+1)
    fh = np.zeros_like(x)
    uh = 0.1 * w_i(1, x) + 0.5 * w_i(n//2, x) + w_i(n-3, x)
    if graphics: 
        plt.plot(x, uh, color = "blue", label = "u_h^0")

    jacobi_1d(uh, fh, 0, kmax)
    if graphics: 
        plt.plot(x, uh, color = "red", label = f"u_h^{kmax}")
        plt.legend()
        plt.show()
    error = la.norm(uh, np.inf)
    print(f"The error after {kmax} Jacobi steps is {error:6.4e}.")

    uh = w_i(n-1, x)
    k = jacobi_1d(uh, fh, 1e-6, kmax)
    error = la.norm(uh, np.inf)
    print(f"Needed {k} Jacobi steps to reduce the error to {error:6.4e}.")

    uh = 0.1 * w_i(1, x) + 0.5 * w_i(n//2, x) + w_i(n-3, x)
    k = mg_w_cycle_1d(uh, fh, 2, 1, 1e-8, kmax)
    error = la.norm(uh, np.inf)
    print(f"Needed {k} W-cycle steps to reduce the error to {error:6.4e}.")

    fh[:] = f_1d(x)
    fh[0], fh[-1] = 0, 0
    uh[0], uh[-1] = fh[0], fh[-1]
    uh[1:n] = 0
    k = mg_w_cycle_1d(uh, fh, 1, 1, 1e-8, kmax)
    if graphics: 
        plt.plot(x, uh, color = "blue")
        plt.show()
    print(f"Needed {k} W-cycle steps to converge.")

    uh[1:n] = 0
    smax = full_mg_1d(uh, fh, 2/3, 1, 1, 2)
    if graphics: 
        plt.plot(x, uh, color = "blue")
        plt.show()
    print(f"The pseudo-residual after a full MG step is {smax:6.4e}.")


def unittest3_2d():
    #Here we test all the 2d routines.
    from project3 import jacobi_step_2d, w_cycle_step_2d, full_mg_2d
    from mpl_toolkits.mplot3d import Axes3D

    def f_2d(x, y): return 3*x**4 - y**2

    def jacobi_2d(uh, fh, tol=1e-8, kmax=100):
        for k in range(kmax):
            if jacobi_step_2d(uh, fh, 0.5) < tol: return k + 1
        return -1
    
    def mg_w_cycle_2d(uh, fh, alpha1, alpha2, tol=1e-8, kmax=100):
        for k in range(kmax):
            if w_cycle_step_2d(uh, fh, 2/3, alpha1, alpha2) < tol: return k + 1
        return -1
    
    def w_ij(i, j, x, y): 
        return np.sin(i * x * np.pi) * np.sin(j * y * np.pi)

    kmax, n, graphics = 50, 64, False
    h = 4/n
    x = np.linspace(0, 1, n+1)
    y = np.linspace(0, 1, n+1)
    x, y = np.meshgrid(x, y)
    fh = np.zeros_like(x)
    uh = (0.1 * w_ij(1, 1, x, y) + 0.5 * w_ij(n//2, n//2, x, y) 
          + w_ij(n-3, n-3, x, y))
    if graphics: 
        fig = plt.figure()
        ax = fig.add_subplot(211, projection="3d")
        ax.set_title("u_h^0")
        ax.plot_surface(x, y, uh)
    jacobi_2d(uh, fh, 0, kmax)
    if graphics: 
        ax = fig.add_subplot(212, projection="3d")
        ax.set_title(f"u_h^{kmax}")
        ax.plot_surface(x, y, uh)
        plt.show()
    error = la.norm(uh.flat, np.inf)
    print(f"The error after {kmax} Jacobi steps is {error:6.4e}.")

    uh = w_ij(n-1, n-1, x, y)
    k = jacobi_2d(uh, fh, 1e-6, kmax)
    error = la.norm(uh.flat, np.inf)
    print(f"Needed {k} Jacobi steps to reduce the error to {error:6.4e}.")
    
    uh = (0.1 * w_ij(1, 1, x, y) + 0.5 * w_ij(n//2, n//2, x, y) 
          + w_ij(n-3, n-3, x, y))
    k = mg_w_cycle_2d(uh, fh, 2, 1, 1e-8, kmax)
    error = la.norm(uh.flat, np.inf)
    print(f"Needed {k} W-cycle steps to reduce the error to {error:6.4e}.")

    fh = np.zeros((n+1,n+1))
    fh[1:n,1:n] = f_2d(x[1:n,1:n], y[1:n,1:n])
    uh = np.zeros((n+1,n+1))
    uh[1:n,1:n] = 0
    k = mg_w_cycle_2d(uh, fh, 1, 1, 1e-8, kmax)
    if graphics: 
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.plot_surface(x, y, uh)
        plt.show()
    print(f"Needed {k} W-cycle steps to converge.")
    
    uh[1:n,1:n] = 0
    smax = full_mg_2d(uh, fh, 2/3, 1, 1, 1)
    if graphics: 
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.plot_surface(x, y, uh)
        plt.show()
    print(f"The pseudo-residual after a full MG step is {smax:6.4e}.")

    
print("****** 1d tests ********")
unittest3_1d()

print("****** 2d tests ********")
unittest3_2d()
