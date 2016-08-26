from __future__ import division,print_function

import os, sys, inspect
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a)

from classes import *
#from ElectronInMagneticBottle import *

from math import pi
from scipy import integrate

def calcBcoils(px, py, pz):
    eta = float(N*I*10**(-5))          #Constant "muo * I * 10^2 / 4 * pi" [T*cm]
    
    print(N, I, R, d)
    
    #Equations to Integrate
    lfdBx = lambda a: eta*(-R*pz*sin(a) - R*py*cos(a) + R**2) / (((px + d)**2 + (py - R*cos(a))**2 + (pz - R*sin(a))**2)**(3/2))
    rtdBx = lambda a: eta*(-R*pz*sin(a) - R*py*cos(a)+R**2) / (((px - d)**2 + (py - R*cos(a))**2 + (pz - R*sin(a))**2)**(3/2))
    lfdBy = lambda a: eta*(R*cos(a)*(px + d)) / (((px + d)**2 + (py - R*cos(a))**2 + (pz - R*sin(a))**2)**(3/2))
    rtdBy = lambda a: eta*(R*cos(a)*(px - d)) / (((px - d)**2 + (py - R*cos(a))**2 + (pz - R*sin(a))**2)**(3/2))
    lfdBz = lambda a: eta*(R*sin(a)*(px + d)) / (((px + d)**2 + (py - R*cos(a))**2 + (pz - R*sin(a))**2)**(3/2))
    rtdBz = lambda a: eta*(R*sin(a)*(px - d)) / (((px - d)**2 + (py - R*cos(a))**2 + (pz - R*sin(a))**2)**(3/2))
    
    # Integrate Functions Iteratively
    lfBx = integrate.quad(lfdBx, 0, 2*pi)
    lfBy = integrate.quad(lfdBy, 0, 2*pi)
    lfBz = integrate.quad(lfdBz, 0, 2*pi)
    rtBx = integrate.quad(rtdBx, 0, 2*pi)
    rtBy = integrate.quad(rtdBy, 0, 2*pi)
    rtBz = integrate.quad(rtdBz, 0, 2*pi)

    # Find total B
    bx = lfBx[0] + rtBx[0]
    by = lfBy[0] + rtBy[0]
    bz = lfBz[0] + rtBz[0]
    
    return bx, by, bz

Px = -4.75
Py = 0
Pz = 0

ver = 4

if ver == 4:
    wireCoils = WireCoilPair(None, 1, 1, 5, 5)
    print(wireCoils.N, wireCoils.I, wireCoils.R, wireCoils.d)
    Bx, By, Bz = wireCoils.calcBatP(Px, Py, Pz)
elif ver == 3:
    d = 5
    R = 5
    N = 1                              #Number of turns of wire in coil [None]
    I = 1                              #Current [A]
    Bx, By, Bz = calcBcoils(Px, Py, Pz)
    
print(Bx, By, Bz)