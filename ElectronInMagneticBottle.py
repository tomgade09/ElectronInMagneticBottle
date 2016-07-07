# Tom Gade
# Electron in a Magnetic Bottle
# v3.2.0 06Jul16
# Works on Python 2.7 ONLY (as of now), and VPython 6
# Written as part of a Senior research project at the Citadel, 2008-2009

from __future__ import division
from math import pi
from scipy import integrate

import os, sys, inspect
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a)

# Only use one of these at a time to avoid namespace conflict
from VPyDraw import *
#from pygletDraw import *
#from OpenGLDraw import *

muo = float(4*pi*10**(-7))         #Permeability of Free Space [T*m/A]
eom = 1.76*10**11                  #Charge to mass ratio [C/kg]
dt = 5*10**-9                      #Time slice [s]
d = 5                              #Half distance between coils [cm]
R = 5                              #Radius of coil [cm]

def calcBcoils(px, py, pz):
    N = 1                              #Number of turns of wire in coil [None]
    I = 1                              #Current [A]
    eta = float(N*I*10**(-5))          #Constant "muo * I * 10^2 / 4 * pi" [T*cm]
    
    #Equations to Integrate
    lfdBx = lambda a: eta*(-R*pz*sin(a)-R*py*cos(a)+R**2)/(((px+d)**2+
                                                       (py-R*cos(a))**2+
                                                       (pz-R*sin(a))**2)**(3/2))
    rtdBx = lambda a: eta*(-R*pz*sin(a)-R*py*cos(a)+R**2)/(((px-d)**2+
                                                       (py-R*cos(a))**2+
                                                       (pz-R*sin(a))**2)**(3/2))
    lfdBy = lambda a: eta*(R*cos(a)*(px+d))/(((px+d)**2+
                                        (py-R*cos(a))**2+
                                        (pz-R*sin(a))**2)**(3/2))
    rtdBy = lambda a: eta*(R*cos(a)*(px-d))/(((px-d)**2+
                                        (py-R*cos(a))**2+
                                        (pz-R*sin(a))**2)**(3/2))
    lfdBz = lambda a: eta*(R*sin(a)*(px+d))/(((px+d)**2+
                                        (py-R*cos(a))**2+
                                        (pz-R*sin(a))**2)**(3/2))
    rtdBz = lambda a: eta*(R*sin(a)*(px-d))/(((px-d)**2+
                                        (py-R*cos(a))**2+
                                        (pz-R*sin(a))**2)**(3/2))
    
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
    
def calcV(bx, by, bz, vx, vy, vz):
    # Calculating B's Effect on the Electron's Velocity
    vXB_x = vy*bz - vz*by     #cmpts of v cross B [T*cm/s]
    vXB_y = vz*bx - vx*bz
    vXB_z = vx*by - vy*bx
    ax = -eom * vXB_x         #cmpts of acceleration [cm/s2]
    ay = -eom * vXB_y
    az = -eom * vXB_z
    vx += ax*dt               #cmpts of new velocity [cm/s]
    vy += ay*dt
    vz += az*dt
    
    return vx, vy, vz
    
def calcP(px, py, pz, vx, vy, vz):
    # Iterate Time Forward dt, Calculate New Position
    px += vx*dt              #cmpts of new position [cm]
    py += vy*dt
    pz += vz*dt
    
    return px, py, pz

def drawBlines(Po_x, Po_y, Po_z, linedist):
    px = Po_x
    py = Po_y
    pz = Po_z
    
    while px < 5:
        bx, by, bz = calcBcoils(px, py, pz)
        #normfact = linedist / sqrt(bx**2 + by**2 + bz**2)
        normfact = 10000
        
        bxn = bx * normfact
        byn = by * normfact
        bzn = bz * normfact
        
        drawLines(px, py, pz, bxn, byn, bzn)
        
        px += bxn
        py += byn
        pz += bzn

def main():
    ind = 0                            #Index (Calculation Counter)
    t = 0                              #Initial time [s]
    Vx = 1000                          #Initial Electron Velocity [cm/s]
    Vy = 1000
    Vz = 1000
    Px = -4.75                         #Initial Electron Position [cm]
    Py = 0
    Pz = 0

    # Draw some things
    scene1, relTclock = drawWindow(1920, 1080, Px, Py, Pz, t, d, R)    
    drawWireLoopPair(d, R)
    electron1 = drawParticlePic(Px,Py,Pz,10,100)
    drawBlines(-5, 0, -3.5, 0.1)
    drawBlines(-5, -3.5, 0, 0.1)
    drawBlines(-5, 0, 3.5, 0.1)
    drawBlines(-5, 3.5, 0, 0.1)
        
    while t <= .001:
        FPSrate(10000)
        Bx, By, Bz = calcBcoils(Px, Py, Pz)
        Vx, Vy, Vz = calcV(Bx, By, Bz, Vx, Vy, Vz)
        Px, Py, Pz = calcP(Px, Py, Pz, Vx, Vy, Vz)
        t += dt                  #time increase [s]
        ind += 1                 #index increase
        updatePic(electron1,relTclock,Px,Py,Pz,t)

if __name__ == "__main__":
    main()