# -*- coding: utf-8 -*-
"""
Created on Sun Sep 21 17:43:51 2014

@author: tfuser
"""

import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import numpy as np
import pylab

from diffMatrix import setD

def u_exact(x):
    """ 
    exact solution for this problem
    -x^2/2+3x/2
    """
    
    # return exact solution
    return -x**2/2.0+3*x/2.0
    
def f_exact(x):
    """
    return rhs for manufactured solution given by u_exact
    """
    return 1
    
def grid_norm2(f,h):
    """calculate grid L2 norm given discrete function f with uniform spacing h
    """
    
    return np.sqrt(h)*np.linalg.norm(f, 2)
    
# plot solution and errors
def plotpoisson1d(x,u_solve,title=None):
    """ Quick test routine to display solution of 
        input: x: numpy array of mesh-points
               f: function pointer to function to differentiate
    """  

    # compare solution to analytic
    pylab.figure()
    pylab.plot(x,u_solve,'rx-',x,u_exact(x),'b+-')
    pylab.title("compare solutions N={0}".format(len(x)))
    pylab.legend(["u_h", "u_exact"],loc="best")
    pylab.xlabel("x")
    pylab.ylabel("u")
    pylab.grid()
        
    # plot errors
    err = u_solve - u_exact(x)    
    pylab.figure()
    pylab.plot(x,err)
    pylab.xlabel("x")
    pylab.ylabel("err")
    pylab.title("errors, N={0}".format(len(x)))
    if title:
        pylab.title(title)
    pylab.show(block=False)

def plotconvergence(N_ar,abs_err_ar,rel_err_ar):
    """ make pretty convergence plots
    """
    # plot absolute and relative errors against h
    # convert from lists to numpy arrays for plotting
    h = 1./np.array(N_ar)
    abs_err = np.array(abs_err_ar)
    rel_err = np.array(rel_err_ar)
    
    # calculate best-fit polynomial to log(h), log(abs_err)
    p = np.polyfit(np.log(h),np.log(abs_err),1)
    pylab.figure()
    pylab.loglog(h,abs_err,'o-',h,rel_err,'d-')
    pylab.xlabel("Mesh width")
    pylab.ylabel("Error")
    pylab.title("Errors in numeric solution of the Poisson equation")
    pylab.legend(["Absolute Error","Relative error"],loc="best")
    pylab.grid()
    pylab.show(block=False)
    

def main():
    # set number of mesh intervals (mesh points is N+1)
    N_ar = [8,16,32,64,128,256,512,1024]
    # initialize lists for storing output
    abs_err_ar = []
    rel_err_ar = []
    for N in N_ar:
        # set numpy grid array to be evenly spaced    
        x = np.linspace(0,1,N+1)
        
        u_true = u_exact(x)
    
        # get full second order second derivative operator
        A = -setD(2,x)
        
        # set rhs vector
        f = f_exact(x)
        
        # set solution vector with boundary values 
        u_bound = np.zeros(x.shape)        
        u_bound[0] = 0
        u_bound[-1] = 1
        u_solve = u_bound
    
        # lift rhs
        f = f - A*u_bound
    
        # solve interior points
        u_solve[1:-1] = spsolve(A[1:-1,1:-1],f[1:-1])
        
        # calculate error and errornorm
        h = 1./N
        err = u_solve - u_true
        abs_err = grid_norm2(err,h)
        rel_err = abs_err/grid_norm2(u_solve,h)
        
       
        print 'N=', N, 'abs_err=', abs_err, 'rel_err=',rel_err
        # collect results for later plotting 
        abs_err_ar.append(abs_err)
        rel_err_ar.append(rel_err)
        
    # plot it out
    plotpoisson1d(x,u_solve)
    plotconvergence(N_ar,abs_err_ar,rel_err_ar)
    

if __name__ == "__main__":
    main()
