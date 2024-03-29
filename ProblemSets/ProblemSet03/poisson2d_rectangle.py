# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 20:16:41 2014
Script to solve the poisson equation with manufactured solution:
u = exp(x+y/2)
over the [0,1]x[0,2] with uniform mesh and dirichlet boundary conditions.

@author: tfuser
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import pylab
from mpl_toolkits.mplot3d import Axes3D
import time

def u_exact(x,y):
    return np.exp(x+y/2)
    
def f(x,y):
    return 1.25*np.exp(x+y/2)
    
def grid_norm2(f,h):
    return np.sqrt(h)*np.linalg.norm(f, 2)

def calcSolution(h,show_matri,show_result):
    ax = 0.0
    bx = 2.0
    ay = 0.0
    by = 1.0
    mx = 2.0/h-1
    my = 1.0/h-1
    x = np.linspace(ax,bx,mx+2)   # grid points x including boundaries
    y = np.linspace(ay,by,my+2)   # grid points y including boundaries

    X,Y = np.meshgrid(x,y)     # 2d arrays of x,y values
    X = X.T                    # transpose so that X(i,j),Y(i,j) are
    Y = Y.T                    # coordinates of (i,j) point

    Xint = X[1:-1,1:-1]        # interior points
    Yint = Y[1:-1,1:-1]
    rhs = f(Xint,Yint)         # evaluate f at interior points for right hand side
                           # rhs is modified below for boundary conditions.

    # set boundary conditions around edges of usoln array:

    usoln = np.zeros(X.shape)
    usoln[:,0] = u_exact(x,ay)
    usoln[:,-1] = u_exact(x,by)
    usoln[0,:] = u_exact(ax,y)
    usoln[-1,:] = u_exact(bx,y)

    # adjust the rhs to include boundary terms: 
    rhs[:,0] -= usoln[1:-1,0] / h**2
    rhs[:,-1] -= usoln[1:-1,-1] / h**2
    rhs[0,:] -= usoln[0,1:-1] / h**2
    rhs[-1,:] -= usoln[-1,1:-1] / h**2


    # convert the 2d grid function rhs into a column vector for rhs of system:
    F = rhs.reshape((mx*my,1))
    
    # form matrix A:
    Ix = sp.eye(mx,mx)
    Iy = sp.eye(my,my)
    ex = np.ones(mx)
    ey = np.ones(my)
    T = sp.spdiags([ey,-4.*ey,ey],[-1,0,1],my,my)
    S = sp.spdiags([ex,ex],[-1,1],mx,mx)
    A = (sp.kron(Ix,T) + sp.kron(S,Iy)) / h**2    
    A = A.tocsr()
    
    show_matrix = True
    if (show_matrix):
        pylab.figure()
        pylab.spy(A,marker='.')
        
    # Solve the linear system:
    tic = time.time()
    uvec = spsolve(A, F)
    toc = time.time()
    
    # reshape vector solution uvec as a grid function and
    # insert this interior solution into usoln for plotting purposes:
    # (recall boundary conditions in usoln are already set)
    
    usoln[1:-1, 1:-1] = uvec.reshape( (mx,my) )
    
    show_result = True
    if show_result:
        # plot results:
        pylab.figure()
        ax = Axes3D(pylab.gcf())
        ax.plot_surface(X,Y,usoln, rstride=1, cstride=1, cmap=pylab.cm.jet)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('u')
        #pylab.axis([a, b, a, b])
        #pylab.daspect([1 1 1])
        pylab.title('Surface plot of computed solution')
        
        pylab.show(block=False)
        
def main():
    calcSolution(1.0/16.0,True,True)
            
if __name__ == "__main__":
    main()
