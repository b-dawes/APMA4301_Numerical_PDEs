# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 18:42:43 2014
Script to solve the poisson equation with manufactured solution:
u = exp(x+y/2)
over the unit square with uniform mesh and dirichlet boundary conditions.

Also plots relative and absolute errors
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

def calcSolution(m,show_matrix,show_result):
    a = 0.0
    b = 1.0
    m = 16                   # number of interior points in each direction
    h = (b-a)/(m+1)
    x = np.linspace(a,b,m+2)   # grid points x including boundaries
    y = np.linspace(a,b,m+2)   # grid points y including boundaries


    X,Y = np.meshgrid(x,y)     # 2d arrays of x,y values
    X = X.T                    # transpose so that X(i,j),Y(i,j) are
    Y = Y.T                    # coordinates of (i,j) point

    Xint = X[1:-1,1:-1]        # interior points
    Yint = Y[1:-1,1:-1]
    rhs = f(Xint,Yint)         # evaluate f at interior points for right hand side
                           # rhs is modified below for boundary conditions.

    # set boundary conditions around edges of usoln array:
    usoln = np.zeros(X.shape)
    usoln[:,0] = u_exact(x,a)
    usoln[:,-1] = u_exact(x,b)
    usoln[0,:] = u_exact(a,y)
    usoln[-1,:] = u_exact(b,y)

    # adjust the rhs to include boundary terms: 
    rhs[:,0] -= usoln[1:-1,0] / h**2
    rhs[:,-1] -= usoln[1:-1,-1] / h**2
    rhs[0,:] -= usoln[0,1:-1] / h**2
    rhs[-1,:] -= usoln[-1,1:-1] / h**2


    # convert the 2d grid function rhs into a column vector for rhs of system:
    F = rhs.reshape((m*m,1))
    
    # form matrix A:
    I = sp.eye(m,m)
    e = np.ones(m)
    T = sp.spdiags([e,-4.*e,e],[-1,0,1],m,m)
    S = sp.spdiags([e,e],[-1,1],m,m)
    A = (sp.kron(I,T) + sp.kron(S,I)) / h**2
    A = A.tocsr()
    
    if show_matrix:
        pylab.figure()
        pylab.spy(A,marker='.')
        
    # Solve the linear system:
    tic = time.time()
    uvec = spsolve(A, F)
    toc = time.time()
    
    # reshape vector solution uvec as a grid function and
    # insert this interior solution into usoln for plotting purposes:
    # (recall boundary conditions in usoln are already set)
        
    usoln[1:-1, 1:-1] = uvec.reshape( (m,m) )
    
    # Find errors
    abs_err = grid_norm2(usoln-u_exact(x,y),h)
    rel_err = abs_err/grid_norm2(usoln,h)
    print "m = {0}".format(m)
    print "Absolute error = {0:10.3e}, relative error = {1:10.3e}".format(abs_err,rel_err)
    print 'Elapsed Time = {0} s'.format(toc-tic)
    
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
    
    return [abs_err,rel_err]
            
def main():
    m_arr = [16,32,64,128,256]
    abs_err = []
    rel_err = []
    h = []
    for m in m_arr:
        i,j = calcSolution(m,False,False)
        abs_err.append(i)
        rel_err.append(j)
        h.append(1.0/m)
    calcSolution(32,True,True)
    pylab.figure()
    pylab.loglog(h,abs_err,'-o')
    pylab.loglog(h,rel_err,'-d')

if __name__ == "__main__":
    main()
