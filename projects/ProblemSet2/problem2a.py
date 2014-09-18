# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 22:29:42 2014
Script to solve problem 2a on problem set 2.
Calculates and plots the errors in the 1st and 2nd derivatives of f=x^2+sin4pix
over the interval [0,1] as a function of the mesh width.
@author: tfuser
"""

import scipy.sparse as sp
import numpy as np
import pylab

from diffMatrix import setD
    
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
    
    # Plot error in derivatives
    pylab.figure(1)
    pylab.axis([0.0008,.2,0.0001,100])
    pylab.loglog(h,eh1,'-o',h,eh2,'-d')
    pylab.xlabel('Mesh width')
    pylab.ylabel('Absolute mesh error')
    pylab.title('Error in derivatives of $f=x^2+\sin(4\pi x)$')
    pylab.legend(['First derivative','Second derivative'],loc='best')
    pylab.grid()
    pylab.draw()
    pylab.show()


if __name__ == "__main__":
    main()
