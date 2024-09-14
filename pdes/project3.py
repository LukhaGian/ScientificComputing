"""
    Project3
    Name: Gianmaria
    Surname: Lucca
    MAT: 241440
"""

import numpy as np
"""
We'll solve the equation
 -Δu + (2 - |x|)u = f            (*)
 with boundary conditions u = 0
"""

########################
"""
In 1D, f(x) = (2x - 1)^3
"""
def jacobi_step_1d(uh, fh, omega):
    """
    function that performs one step of the weighted
    Jacobi method with weight omega for the linear system
    arising from the discretization of the equation (*)
    It outputs the new iterate uh^(k+1) (in place of uh) and the 
    pseudo-residual |uh^(k+1) - uh^(k)|
    """
    n = len(uh) - 1
    x = np.linspace(0, 1, n+1)
    h = 1/n
    den = 2 + (h**2)*(2 - np.abs(x))
    v = uh.copy()
    smax = 0
    for i in range(1,n):
        # Std Jacobi step
        v[i] = (((h**2)*fh[i]) + uh[i-1] + uh[i+1])/den[i]
        # convex comb, weighted step
        v[i] = (1-omega)*uh[i] + omega*v[i]
        smax = max(smax, abs(v[i] - uh[i]))
    uh[:] = v    
    return smax


def two_grid_jacobi(uh, fh, omega, nu):
    """
    Function that implements the 2 grid correction scheme
    Steps:
    1) Pre smoothing A_h u_h = f_h ==> Obtain u_h^k (By weighted Jacobi Method)
    2) Compute the residual r_h = f_h - A_h u_h^k
    3) Cast the residual in the coarser grid
    4) Solve A_2h e_2h = r_2h with nu weighted Jacobi
    5) Cast the error in the fine grid and u_h^(k+1) = u_h(k) + casted error
    6) Post smoothing A_h u_h = f_h with inital guess u_h^(k+1)
    """
    #1)
    jacobi_step_1d(uh, fh, omega)
    #2)
    n = len(uh) - 1
    x = np.linspace(0, 1, n+1)
    h = 1/n
    rh = np.zeros_like(x)
    for i in range(1, n):
        rh[i] = fh[i] + (uh[i-1] - 2*uh[i] + uh[i+1])/(h**2) - (2 - abs(x[i]))*uh[i]
    #3)
    n_half = int(n/2)
    r2h = np.zeros(n_half + 1)
    e2h = np.zeros(n_half + 1) # inital guess
    for j in range(1,n_half):
        r2h[j] = (0.25)*(rh[2*j-1] + 2*rh[2*j] + rh[2*j+1])
    #4)
    for _ in range(nu):
        jacobi_step_1d(e2h, r2h, omega)
    #5)
    eh = np.zeros_like(x)
    for p in range(0,int(n/2 + 1)):
        eh[2*p] = e2h[p]
    for p in range(0,int(n/2)):
        eh[2*p+1] = (0.5)*(e2h[p] + e2h[p+1])
    # last step
    v = uh + eh
    uh = v.copy()
    #6) 
    smax = jacobi_step_1d(uh, fh ,omega)
    return smax


def w_cycle_step_1d(uh, fh, omega, alpha1, alpha2):
    """
    Function that performs a W-cycle, using alpha1 pre smoothing steps with weighted Jacobi
    and alpha2 post smoothing with weighted Jacobi
    """
    mu = 2 # W cycle
    n = len(uh) -1
    h = 1/n
    if n == 2: # base case, only 1 interior node ==> Exact solve
        uh[1] = fh[1]/((2/(h**2))+ 1.5)
        smax = 0
    else:
        # relax alpha1 times the system
        for _ in range(alpha1):
            jacobi_step_1d(uh, fh, omega)
        # cast into the coarser grid (f_h - A_h u_h)
        # let v = f_h - A_h u_h
        x = np.linspace(0, 1, n+1)
        v = np.zeros_like(x)
        for i in range(1, n):
            v[i] = fh[i] + (uh[i-1] - 2*uh[i] + uh[i+1])/(h**2) - (2 - abs(x[i]))*uh[i]
        n_half = int(n/2)
        f2h = np.zeros(n_half+1)    
        for i in range(1, n_half):
            f2h[i] = (0.25)*(v[2*i-1] + 2*v[2*i] + v[2*i+1])
        # initialize u2h
        u2h = np.zeros(n_half+1)
        # now we recursively call ourmethod mu times
        for _ in range(mu):
            w_cycle_step_1d(u2h, f2h, omega, alpha1, alpha2)
        # Correction of uh
        # cast u2h into the finer grid
        u2h_cast = np.zeros_like(x)
        for i in range(0, n_half):
            u2h_cast[2*i] = u2h[i]
            u2h_cast[2*i+1] = (0.5)*(u2h[i]+ u2h[i+1])
        # last copy
        u2h_cast[2*n_half] = u2h[n_half]   
        w = uh.copy()
        w = uh + u2h_cast
        uh[:] = w.copy()
        # post smoothing
        for _ in range(alpha2-1):
            jacobi_step_1d(uh, fh, omega)  
        smax = jacobi_step_1d(uh, fh, omega)
    return smax


def full_mg_1d(uh, fh, omega, alpha1, alpha2, nu):
    """
    Function that performs a full multigrid step, performing nu W-cycle steps as defined before
    """
    n = len(fh)-1
    h = 1/n
    if n == 2: # Base case, exact solve
        uh[1] = fh[1]/((2/(h**2))+ 1.5)
        smax = 0
    else:
        # cast fh into the coarser grid
        n_half = int(n/2)
        f2h = np.zeros(n_half+1)
        for i in range(1,n_half):
            f2h[i] = (0.25)*(fh[2*i-1] + 2*fh[2*i] + fh[2*i+1])
        #### we initalize u2h
        u2h = np.zeros(n_half+1)    
        # recursive call
        full_mg_1d(u2h, f2h, omega, alpha1, alpha2, nu)
        # cast u2h in the finer grid into uh
        u2h_cast = np.zeros(n+1)    
        for i in range(0, n_half):
            u2h_cast[2*i] = u2h[i]
            u2h_cast[2*i+1] = (0.5)*(u2h[i]+ u2h[i+1])
        # last copy
        u2h_cast[2*n_half] = u2h[n_half]
        uh[:] = u2h_cast.copy()
        # we now call nu times the W-cycle
        for _ in  range(nu-1):
            w_cycle_step_1d(uh, fh, omega, alpha1, alpha2)   
        # last call
        smax = w_cycle_step_1d(uh, fh, omega, alpha1, alpha2)
    return smax

########################
"""
In 2D, f(x, y) = 3x^4 - y^2
"""
def jacobi_step_2d(uh, fh, omega):
    """
    Function that performs a weighted Jacobi step
    Outputs: the pseudo residual |uh^(k+1) - uh^(k)| and update the new iterate uh^(k+1)
    """
    n = len(uh)-1 # here we return the number of rows of the discretization, which is the same as the number of columns
    h = 1/n
    x = np.linspace(0,1,n+1)
    y = np.linspace(0,1,n+1)
    x,y = np.meshgrid(x,y)
    normx = normcalc(x,y)
    den = 4 + (h**2)*(2 - normx)
    smax = 0
    v = uh.copy()
    for i in range(1, n):
        for j in range(1, n):
            # Std Jacobi step
            v[i,j] = (((h**2)*fh[i,j]) + uh[i-1,j] + uh[i+1,j] + uh[i,j-1] + uh[i,j+1])/den[i,j]
            # convex comb, weigthed step
            v[i,j] = (1-omega)*uh[i,j] + omega*v[i,j]
            smax = max(smax, abs(v[i,j] - uh[i,j]))       
    uh[:] = v
    return smax


def normcalc(x,y):
    """
    Retruns the Euclidean norm element wise of the matrix
    """
    return np.sqrt(x**2 + y**2)


def w_cycle_step_2d(uh, fh, omega, alpha1, alpha2):
    """
    Function that does a full W-cycle using alpha1 post smoothing steps with
    weighted Jacobi and alpha2 post smoothing steps with weighted Jacobi
    """
    mu = 2
    n = len(uh) - 1
    h = 1/n
    if n == 2: # base case, only 1 interior node ==> Exact solve
        uh[1,1] = fh[1,1]/((4/h**2) + (2 - 1/pow(2,0.5)))
        smax = 0
    else:
        # relax alpha1 times the system
        for _ in range(alpha1):
            jacobi_step_2d(uh, fh, omega)
        # cast into the coarser grid (f_h - A_h u_h)
        # let v = f_h - A_h u_h
        x = np.linspace(0, 1, n+1)
        y = np.linspace(0, 1, n+1)
        v = np.zeros((n+1,n+1))
        x,y = np.meshgrid(x,y)
        normx = normcalc(x,y)
        for i in range(1, n):
            for j in range(1, n):
                v[i,j] = fh[i,j] + (uh[i-1,j] + uh[i+1,j] - 4*uh[i,j] + uh[i,j-1] + uh[i,j+1])/(h**2) - (2 - normx[i,j])*uh[i,j]
        n_half = int(n/2)
        f2h = np.zeros((n_half+1,n_half+1))
        # cast into f2h
        for i in range(1, n_half):
            for j in range(1, n_half):
                f2h[i,j] = (1/16)*(v[2*i-1,2*j-1] + v[2*i-1,2*j+1] + v[2*i+1,2*j-1] + v[2*i+1,2*j+1] + 2*(v[2*i,2*j-1] + v[2*i,2*j+1] + v[2*i-1,2*j] + v[2*i+1,2*j]) + 4*v[2*i,2*j])
        u2h = np.zeros((n_half+1,n_half+1))
        # recursive call mu times
        for _ in range(mu):
            w_cycle_step_2d(u2h, f2h, omega, alpha1, alpha2)
        # correction of uh
        u2h_cast = np.zeros((n+1,n+1))
        for i in range(0, n_half+1):
            for j in range(0, n_half+1):
                u2h_cast[2*i,2*j] = u2h[i,j]
                if i < n_half:
                    u2h_cast[2*i+1,2*j] = (0.5)*(u2h[i,j] + u2h[i+1,j])
                if j < n_half:    
                    u2h_cast[2*i,2*j+1] = (0.5)*(u2h[i,j] + u2h[i,j+1])
                if (i < n_half) and (j < n_half):    
                    u2h_cast[2*i+1,2*j+1] = (0.25)*(u2h[i,j] + u2h[i+1,j] + u2h[i,j+1] + u2h[i+1,j+1])
        w = uh.copy()
        w = uh + u2h_cast
        uh[:] = w.copy()
        # post smoothing
        for _ in range(alpha2-1):
            jacobi_step_2d(uh, fh, omega)
        smax = jacobi_step_2d(uh, fh, omega)    
    return smax


def full_mg_2d(uh, fh, omega, alpha1, alpha2, nu):
    """
    Function that performs a full multigrid step, performing nu W-cycle steps as defined before
    """
    n = len(fh) -1
    h = 1/n
    if n == 2: # base case, exact solve
        uh[1,1] = fh[1,1]/((4/h**2) + (2 - 1/pow(2,0.5)))
        smax = 0
    else:
        # cast fh into the coarser grid
        n_half = int(n/2)
        f2h = np.zeros((n_half+1,n_half+1))
        for i in range(1,n_half):
            for j in range(1,n_half):
                f2h[i,j] = (1/16)*(fh[2*i-1,2*j-1] + fh[2*i-1,2*j+1] + fh[2*i+1,2*j-1] + fh[2*i+1,2*j+1] + 2*(fh[2*i,2*j-1] + fh[2*i,2*j+1] + fh[2*i-1,2*j] + fh[2*i+1,2*j]) + 4*fh[2*i,2*j])
        u2h = np.zeros((n_half+1,n_half+1))
        # recursive call
        full_mg_2d(u2h, f2h, omega, alpha1, alpha2, nu)
        # cast u2h in the finer grid and update uh
        u2h_cast = np.zeros((n+1,n+1))
        for i in range(0, n_half+1):
            for j in range(0, n_half+1):
                u2h_cast[2*i,2*j] = u2h[i,j]
                if i < n_half:
                    u2h_cast[2*i+1,2*j] = (0.5)*(u2h[i,j] + u2h[i+1,j])
                if j < n_half:    
                    u2h_cast[2*i,2*j+1] = (0.5)*(u2h[i,j] + u2h[i,j+1])
                if (i < n_half) and (j < n_half):    
                    u2h_cast[2*i+1,2*j+1] = (0.25)*(u2h[i,j] + u2h[i+1,j] + u2h[i,j+1] + u2h[i+1,j+1])
        uh[:] = u2h_cast.copy()
        # now nu times W-cycle
        for _ in range(nu-1):
            w_cycle_step_2d(uh, fh, omega, alpha1, alpha2)
        # last call
        smax = w_cycle_step_2d(uh, fh, omega, alpha1, alpha2)
    return smax                   


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from time import time

    def f_1d(x):
        return (2*x - 1)**3
    
    def f_2d(x, y):
        return 3*x**4 - y**2
    
    # Tests
    """
    #2)#########
    tol = 1e-8
    for l in range(3,8):
        n = 2**l
        x = np.linspace(0,1,n+1)
        uh = np.zeros_like(x)
        fh = np.zeros_like(x)
        fh[1:n-1] = f_1d(x[1:n-1])
        omega = 2/3
        k = 1
        flag = True # flagger
        tic = time()
        while flag:
            q = jacobi_step_1d(uh, fh, omega)
            k += 1
            if q < tol:
                flag = False
        tac = time() - tic
        #plt.plot(x,uh)    
        #plt.show()    
        print(f"For l = {l}: {k} iterations, time = {tac:.5f} seconds")        
    ########
        
    #3) 2 grid correction scheme
    tol = 1e-8
    for l in range(3,8):
        n = 2**l
        x = np.linspace(0,1,n+1)
        uh = np.zeros_like(x)
        fh = np.zeros_like(x)
        fh[1:n-1] = f_1d(x[1:n-1])
        omega = 2/3
        nu = 10
        k = 1
        flag = True # flagger
        tic = time()
        while flag:
            q = two_grid_jacobi(uh, fh, omega, nu)
            k += 1
            if q < tol:
                flag = False
        tac = time() - tic
        #plt.plot(x,uh)    
        #plt.show()    
        print(f"For l = {l}: {k} iterations, time = {tac:.5f} seconds")    
    # 4) W-cycle
    for l in range(3,15):
        n = 2**l
        x = np.linspace(0,1,n+1)
        uh = np.zeros_like(x)
        fh = np.zeros_like(x)
        fh[1:n-1] = f_1d(x[1:n-1])
        omega = 2/3
        k = 1
        flag = True # flagger
        tic = time()
        while flag:
           q = w_cycle_step_1d(uh, fh, omega, 2,1)
           k += 1
           if q < 1e-8:
               flag = False
        tac = time() - tic      
        print(f"\item For $l = {l}:$ {k} iterations, time = {tac:.5f} seconds ")   

    #6) full multigrid
    for l in range(14,15):
        n = 2**l
        x = np.linspace(0,1,n+1)
        uh = np.zeros_like(x)
        fh = np.zeros_like(x)
        fh[1:n-1] = f_1d(x[1:n-1])
        omega = 2/3
        alpha1 = 1
        alpha2 = 1
        nu = 1
        tic = time()
        smax = full_mg_1d(uh, fh, omega, alpha1, alpha2, nu)
        tac = time() - tic     
        # compute the residual
        h = 1/n 
        res = np.zeros_like(x)
        residual = 0
        for i in range(1,n):
            res[i] = fh[i] + (uh[i-1]-2*uh[i]+uh[i+1])/(h**2) - (2 - abs(x[i]))*uh[i]
            residual = max(residual, abs(res[i]))
            if residual == abs(res[i]):
                ind = i   
        print(f"\item For $l = {l}:$ pseudo-residual = {smax},\\\ residual = {residual}, time = {tac:.5f} seconds")
        #print(x[ind], ind)
    plt.plot(x,uh, markersize='22', color='red')
    plt.title(r'Approximation $\tilde{u}$ to the solution $u$ over the domain $\overline{\Omega}$',fontsize="15")
    plt.xlabel(r'x', fontsize="14")
    plt.ylabel(r'$\tilde{u}$', fontsize="14")   
    plt.show()
    print(min(uh))
    
    #2D ###################################
    # Jacobi
    tol = 1e-8
    for l in range(2,7):
        n = 2**l
        x = np.linspace(0,1,n+1)
        y = np.linspace(0,1,n+1)
        x,y = np.meshgrid(x,y)
        uh = np.zeros((n+1,n+1))
        fh = np.zeros((n+1,n+1))
        #fh = f_2d(x,y)
        fh[1:-1,1:-1] = f_2d(x[1:-1,1:-1],y[1:-1,1:-1])
        #for i in range(0,n+1):
        #    fh[i,0] = 0
        #    fh[i,-1] = 0
        #    fh[0,i] = 0  
        #    fh[-1,i] = 0        
        omega = 2/3
        k = 1
        flag = True # flagger
        tic = time()
        while flag:
            q = jacobi_step_2d(uh, fh, omega)
            k += 1
            if q < tol:
                flag = False
        tac = time() - tic    
        print(f"\item For $l = {l}:$ {k} iterations, time = {tac:.5f} seconds ")
    
    # W-cycle
    for l in range(2,9):
        n = 2**l
        x = np.linspace(0,1,n+1)
        y = np.linspace(0,1,n+1)
        x,y = np.meshgrid(x,y)
        uh = np.zeros((n+1,n+1))
        fh = np.zeros((n+1,n+1))
        #fh = f_2d(x,y)
        fh[1:-1,1:-1] = f_2d(x[1:-1,1:-1],y[1:-1,1:-1])

        omega = 2/3
        k = 1
        flag = True # flagger
        tic = time()
        while flag:
           q = w_cycle_step_2d(uh, fh, omega, 1,2)
           k += 1
           if q < 1e-8:
               flag = False
        tac = time() - tic      
        print(f"\item For $l = {l}:$ {k} iterations, time = {tac:.5f} seconds ")
    
    #5) Full multigrid
    for l in range(2,9):
        n = 2**l
        x = np.linspace(0,1,n+1)
        y = np.linspace(0,1,n+1)
        x,y = np.meshgrid(x,y)
        normx = normcalc(x,y)
        uh = np.zeros((n+1,n+1))
        fh = np.zeros((n+1,n+1))
        fh[1:-1,1:-1] = f_2d(x[1:-1,1:-1],y[1:-1,1:-1])

        omega = 2/3
        alpha1 = 1
        alpha2 = 1
        nu = 1
        tic = time()
        smax = full_mg_2d(uh, fh, omega, alpha1, alpha2, nu)
        tac = time() - tic    
        h = 1/n 
        res = np.zeros((n+1,n+1))
        residual = 0
        for i in range(1,n):
            for j in range(1,n):
                res[i,j] = fh[i,j] + (uh[i-1,j]-4*uh[i,j]+uh[i+1,j]+uh[i,j-1]+uh[i,j+1])/(h**2) - (2 - normx[i,j])*uh[i,j]
                residual = max(residual, abs(res[i,j]))
        print(f"\item For $l = {l}:$ pseudo-residual = {smax},\\\ residual = {residual}, time = {tac:.5f} seconds")
    
    print(min((uh.flatten())))    
    from matplotlib import cm 
    from matplotlib.ticker import LinearLocator
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(x, y, uh, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter('{x:.02f}')
    plt.title(r'Approximation $\tilde{u}$ to the solution $u$ over the domain $\overline{\Omega}$',fontsize="15")
    plt.xlabel(r'$x$', fontsize="14")
    plt.ylabel(r'$y$', fontsize="14")  
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()
    """