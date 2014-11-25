# -*- coding: utf-8 -*-
"""
Stability diagrams:  draws stability diagrams for various ODE schemes
Created on Tue Nov 13 10:13:08 2012

@author: mspieg
"""
import numpy as np
import pylab as pl

def plotR(x,y,R,title=None):
    pl.contourf(x,y,np.abs(R),[0.,1.])
    pl.colorbar()
    pl.axis('tight')
    pl.grid()
    if title:
        pl.title(title)
    return
    
bnd=13.    
x=np.linspace(-bnd,bnd,100)
y=np.linspace(-bnd,bnd,100)
X,Y = np.meshgrid(x,y)
Z = X + Y*1j
V=np.array([1.])

# Forward Euler
R = 1.+Z
pl.figure()
plotR(x,y,R,'Forward Euler')
pl.savefig('LATEX/Euler.pdf')
pl.show()

# Backward Euler
R = 1./(1.-Z)
pl.figure()
plotR(x,y,R,'Backward Euler')
pl.savefig('LATEX/BackwardEuler.pdf')
pl.show()

# Midpt
R = 1.+Z+Z**2/2.
pl.figure()
plotR(x,y,R,'Midpoint')
pl.savefig('LATEX/Midpoint.pdf')
pl.show()

# Improved Euler
R = 1.+Z+Z**2/2.
pl.figure()
plotR(x,y,R,'Imporved Euler')
pl.savefig('LATEX/ImprovedEuler.pdf')
pl.show()

# Trap
R = (1.+Z/2.)/(1.-Z/2.)
pl.figure()
plotR(x,y,R,'Trapezoidal')
pl.savefig('LATEX/Trapezoidal.pdf')
pl.show()

# RK4
R = 1.+Z+Z**2/2.+Z**3/6.+Z**4/24
pl.figure()
plotR(x,y,R,'RK4')
pl.savefig('LATEX/RK4.pdf')
pl.show()

# TR-BDF2
R = (1.+5.*Z/12.)/(1-7.*Z/12.+Z**2/12.)
pl.figure()
plotR(x,y,R,'TR-BDF2')
pl.savefig('LATEX/TRBDF2.pdf')
pl.show()