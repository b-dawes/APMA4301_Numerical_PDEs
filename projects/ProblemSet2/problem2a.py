# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 22:29:42 2014

@author: tfuser
"""

import scipy.sparse as sp
import numpy as np
import pylab
import math
import matplotlib.pyplot as plt

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
    
    assert(k < 3) # check to make sure k < 3
    assert(len(x) > 2)
    
    N = len(x)
    # initialize a sparse NxN matrix in "lil" (linked list) format
    D = sp.lil_matrix((N,N))
    # assign the one-sided k'th derivative at x[0]
    D[0,0:3] = fdcoeffF(k,x[0],x[0:3])
    # assign centered k'th ordered derivatives in the interior
    for i in xrange(1,N-1):
        D[i,i-1:i+2] = fdcoeffF(k,x[i],x[i-1:i+2])
    # assign one sided k'th derivative at end point x[-1]
    D[N-1,-3:] = fdcoeffF(k,x[N-1],x[-3:])
    
    # convert to csr (compressed row storage) and return
    return D.tocsr()

# quicky test program
def plotDf(x,f,title=None):
    """ Quick test routine to display derivative matrices and plot
        derivatives of an arbitrary function f(x)
        input: x: numpy array of mesh-points
               f: function pointer to function to differentiate
    """  
    # calculate first and second derivative matrices
    D1 = setD(1,x) 
    D2 = setD(2,x)
    
    print D1
    print D2
    
    # show sparsity pattern of D1
    pylab.figure()
    pylab.spy(D1,precision=1.e-6)
    pylab.title("Sparsity pattern of D1")
        
    # plot a function and it's derivatives
    y = f(x)    
    pylab.figure()
    pylab.plot(x,y,x,D1*y,x,D2*y,x,2*x+4*math.pi*cos(4*math.pi*x))
    pylab.legend(['f','D1*f','D2*f'],loc="best")
    if title:
        pylab.title(title)
    pylab.show(block=False)
    
def calcError(f,N,k,DkF):
    x = np.linspace(0,1,N+1)
    h = 1.0/N
    
    Dk = setD(k,x)
    return sqrt(h)*normL2(Dk*f(x)-DkF(x))
    
def normL2(x):
    norm = 0    
    for i in range(len(x)):
        norm += x[i]**2
    return sqrt(norm)        

def main():
    def f(x):
        return x**2 + sin(4*pi*x)
    
    def d1f(x):
        return 2*x+4*pi*cos(4*pi*x)
        
    def d2f(x):
        return 2+16*pi**2*sin(4*pi*x)
        
    #plotDf(x,f,"f=x^2 + sin(4\pi x)")

    N = [8,16,32,64,128,256,512,1024]

    eh1 = []
    eh2 = []   
    h = []
    for n in N:
        eh1.append(calcError(f,n,1,d1f))
        eh2.append(calcError(f,n,2,d2f))     
        h.append(1.0/n)
    
    print h
    print eh2
    
    pylab.figure()
    pylab.axis([0.0001,1,0.0001,10])
    pylab.loglog(h,eh1,'-o')
    pylab.legend(['$e_{h1}$'],loc="best")
    
    pylab.figure()
    pylab.axis([0.0001,1,210,225])    
    pylab.loglog(h,eh2,'-o')
    pylab.legend(['$e_{h2}$'],loc='best')

if __name__ == "__main__":
    main()
