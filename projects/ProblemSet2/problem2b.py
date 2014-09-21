# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 21:35:32 2014

Script to solve problem 2b on problem set 2.
Calculates and plots the errors in the 1st and 2nd derivatives of f=x^2+sin4pix
over the interval [0,1] as a function of the mesh width.

@author: tfuser
"""

import scipy.sparse as sp
import numpy as np
import pylab
from fdcoeffF import fdcoeffF

def setD(k,x):
    """ 
    example function for setting k'th order sparse differentiation matrix over 
    arbitrary mesh x with a 3 point stencil
    
    input:
        k = degree of derivative <= 2
        x = numpy array of coordinates >=3 in length
    returns:
        D sparse differention matrix
    """
    
    assert(k < 5) # check to make sure k < 3
    assert(len(x) > 4)
    
    N = len(x)
    # initialize a sparse NxN matrix in "lil" (linked list) format
    D = sp.lil_matrix((N,N))
    print x[0:5]
    # assign the k'th derivative at x[0] and x[1]
    D[0,0:5] = fdcoeffF(k,x[0],x[0:5])
    D[1,0:5] = fdcoeffF(k,x[1],x[0:5])
    # assign centered k'th ordered derivatives in the interior
    for i in xrange(2,N-2):
        D[i,i-2:i+3] = fdcoeffF(k,x[i],x[i-2:i+3])
    # assign one sided k'th derivative at end point x[-1]
    D[N-1,-5:] = fdcoeffF(k,x[N-1],x[-5:])
    D[N-2,-5:] = fdcoeffF(k,x[N-2],x[-5:])
    
    # convert to csr (compressed row storage) and return
    return D.tocsr()
    
def calcError(f,N,k,DkF):
    '''
    Method to calculate the absolute mesh error in the kth derivative over 
    a uniform mesh from 0 to 1 containing N points.
    inputs:
        f - input function
        N - number of points in mesh
        k - order of derivative
        DkF - kth derivative of f
    '''
    x = np.linspace(0,1,N+1)
    h = 1.0/N
    
    Dk = setD(k,x)
    
    return np.sqrt(h*np.sum((Dk*f(x)-DkF(x))**2))   

def main():
    # Test function and first and second derivatives
    def f(x):
        return x**2 + np.sin(4*np.pi*x)
    
    def d1f(x):
        return 2*x+4*np.pi*np.cos(4*np.pi*x)
        
    def d2f(x):
        return 2-16*np.pi**2*np.sin(4*np.pi*x)
    
    # Number of points in mesh    
    N = [8,16,32,64,128,256,512,1024]

    # Error in first and second derivatives for each N
    eh1 = []
    eh2 = []  
    # Mesh width for each N
    h = []
    
    # Calc errors and h
    for n in N:
        eh1.append(calcError(f,n,1,d1f))
        eh2.append(calcError(f,n,2,d2f))     
        h.append(1.0/n)
    
    print eh1
    print eh2
    
    # Plot error in derivatives
    pylab.figure(1)
    pylab.axis([0.0008,.2,0.000000005,500])
    pylab.loglog(h,eh1,'-o',h,eh2,'-d')
    pylab.xlabel('Mesh width')
    pylab.ylabel('Absolute mesh error')
    pylab.title('Error in derivatives of $f=x^2+\sin(4\pi x)$ (5pt stencil)')
    pylab.legend(['First derivative','Second derivative'],loc='best')
    pylab.grid()
    pylab.draw()
    pylab.show()


if __name__ == "__main__":
    main()
