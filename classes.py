from __future__ import division

# Import visual library.  Only use one of these at a time to avoid namespace conflicts.
#from VPyDraw import *
from OpenGLDraw import *

from scipy import integrate
from math import *
import numpy as np
from vectortools import *
import version
import ctypes, os, sys, inspect

classpath = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))

class Particle(object):
    """Define a particle to be placed in the specified magnetic field.
    
    Attributes of the particle:
    - Charge (q) [C]
    - Mass [kg]
    - Position Vector [px,py,pz] [cm]
    - Velocity Vector [vx,vy,vz] [m/s]
    
    Objects associated with the particle:
    - wind - The window object where the particle is being drawn
        - Note: if not using visualization, set this to 'None'
    - pic - The visual point object representing the position of the particle"""
    def __init__(self, windObj, charge, mass, po, vo, name=None):
        """Initiate a particle object."""
        self.wind = windObj
        self.q = charge
        self.mass = mass
        self.p = po
        self.v = vo
        self.name = name
        self.eom = self.q / self.mass
        self.pic = None

    def initDraw(self, intrvl, traillng, Dcolor=(0,1,0)):
        """Draw a point object representing the particle in the object specified by self.wind.  intrvl represents how often to 'draw' a point.  traillng represents how long of a 'trail' to leave behind the current position of the particle."""
        if self.pic is not None:
            print "Pic has already been initialized.  Use updDraw to change the position."
            return
        self.pic = drawParticlePic(self.wind, self.p, intrvl, traillng, Dcolor)
        return self.pic
        
    def updDraw(self):
        """Update the location of the point object drawn with initDraw.  Obviously, it can't be updated if it hasn't been initialized.  Use the initDraw function first."""
        if self.pic is None:
            print "Pic has not been initialized.  Use initDraw to create a picture of the particle first."
            return
        updateParticlePic(self.wind, self.pic, self.p)
    
    def calcBatP(self, pB):
        """Calculate B at a point pB due to this particle."""
        pB = pB[:] #If particle.p is passed in, need to change pB to a list vs a pointer
        pB -= np.array(self.p)
        if abs(pB[0]) <= 10e-15 and abs(pB[1]) <= 10e-15 and abs(pB[2]) <=10e-15:
            return 0, 0, 0
        c = 10e-5 * self.q / sqrt(pB[0]**2 + pB[1]**2 + pB[2]**2)**3
        b = np.cross(self.v, pB)
        b = np.array(b) * c
        return b[0], b[1], b[2]
    
    def __updV(self, b, dt):
        """Calculate the new velocity of the particle based on the specified B field."""
        dv = np.cross(self.v, b) * self.eom * dt
        self.v += dv
        
    def updP(self, b, dt):
        """Calculate the new position based on the particle's velocity."""
        self.__updV(b, dt)
        #self.p += self.v * dt #For some reason, doesn't work, but would be quicker.
        for i in range(3):
            self.p[i] += self.v[i] * dt
    
    def foRKvCrossB(self, BFieldObj, h): #Highly experimental!  Not sure if algorithm is implemented right.
        k1 = self.eom * np.cross(self.v,BFieldObj.totalBatP(self.p)) * h
        k2 = self.eom * np.cross(self.v + k1 / 2, BFieldObj.totalBatP(self.p +
            np.array(self.v) * h / 2)) * h
        k3 = self.eom * np.cross(self.v + k2 / 2, BFieldObj.totalBatP(self.p +
            np.array(self.v) * h / 2)) * h
        k4 = self.eom * np.cross(self.v + k3, BFieldObj.totalBatP(self.p +
            np.array(self.v) * h)) * h
    
        self.v += (k1 + 2 * (k2 + k3) + k4) / 6
        self.p += self.v * h
            
class Electron(Particle):
    """Define an electron as a specific type of 'Particle'"""
    def __init__(self, windObj, po, vo):
        """Values pulled from https://en.wikipedia.org/wiki/Electron, July 15, 2016."""
        Particle.__init__(self, windObj, -1.60217657e-19, 9.10938356e-31, po, vo)
        
class Positron(Particle):
    """Define a positron as a specific type of 'Particle'"""
    def __init__(self, windObj, po, vo):
        """Values pulled from https://en.wikipedia.org/wiki/Positron, July 17, 2016."""
        Particle.__init__(self, windObj, 1.60217657e-19, 9.10938291e-31, po, vo)

class Proton(Particle):
    """Define a proton as a specific type of 'Particle'"""
    def __init__(self, windObj, po, vo):
        """Values pulled from https://en.wikipedia.org/wiki/Positron, July 17, 2016."""
        Particle.__init__(self, windObj, 1.60217662e-19, 1.67262190e-27, po, vo)
        
class WireCoilPair(object):
    """Define a pair of wire coils to create a B field.
    
    Attributes of the wire coils:
    - Number of wire turns: N
    - Current running through coils: I [A]
    - Radius of coils: R [cm]
    - Distance of center of loop from origin: d [cm]
    
    Note: This builds a pair of coils, parallel to one another, offset from the origin by a distance d.  There is no way to build a single wire loop, or to have the loops offset from one another.  Hence the name 'Wire Coil Pair'"""
    def __init__(self, windObj, Cpair, axis, N, I, R, d, name=None):
        """Initiate a WireCoilPair object."""
        self.wind = windObj
        self.Cpair = np.array(Cpair)
        self.axis = np.array(axis)
        self.N = N; self.I = I; self.R = R; self.d = d
        self.name = name
        self.cst = float(self.N * self.I * 10**(-5))
        self.pic = None
        
        if self.axis[0] != 0 and self.axis[1] == 0 and self.axis[2] == 0:
            #X axis - no rotation, don't run code that's not necessary
            self.axiscf_theta = self.axiscf_phi = self.axis_phi = 0
            self.axis_theta = pi / 2
        elif self.axis[0] == 0 and self.axis[1] == 0 and self.axis[2] != 0:
            #Z axis - causes problems in spherical coordinate system
            self.axiscf_theta = - pi / 2; self.axiscf_phi = 0
            self.axis_theta = self.axis_phi = 0
        else:
            self.axis_rho, self.axis_theta, self.axis_phi = cartesianToSpherical(
                self.axis)
            self.axiscf_theta = self.axis_theta - (pi / 2)
            self.axiscf_phi = self.axis_phi

        if self.axis_theta == 0 and self.axis_phi == 0: #Z Axis calculate center points
            self.Cright = [self.Cpair[0], self.Cpair[1], self.Cpair[2] + self.d]
            self.Cleft = [self.Cpair[0], self.Cpair[1], self.Cpair[2] - self.d]
        else: #All other cases calculate center points for loops
            self.Cright = sphericalToCartesian(self.d, self.axis_theta, self.axis_phi)
            if self.axis_phi > pi:
                self.Cleft = sphericalToCartesian(self.d, pi-self.axis_theta,
                    self.axis_phi-pi)
            else:
                self.Cleft = sphericalToCartesian(self.d, pi-self.axis_theta,
                    self.axis_phi+pi)
            self.Cright += self.Cpair; self.Cleft += self.Cpair

    def initDraw(self):
        """Draw the pair of Wire Coils."""
        self.pic = drawWireCoilPair(self.wind, self.Cpair, self.axis, self.Cleft,
            self.Cright, self.R)
    
    def calcBatPinC(self, p):
        """Calculate the B field as a result of the wire coils at a position P."""
        pcnt = [p[0] - self.Cpair[0], p[1] - self.Cpair[1], p[2] - self.Cpair[2]]
        ppr = rotateVector(pcnt, -self.axiscf_theta, -self.axiscf_phi)
        lib = ctypes.CDLL(classpath + '/c/WireCoilPairB.so') # For Linux
        lib.dBx.restype = ctypes.c_double
        lib.dBy.restype = ctypes.c_double
        lib.dBz.restype = ctypes.c_double
        lib.dBx.argtypes = (ctypes.c_int, ctypes.c_double)
        lib.dBy.argtypes = (ctypes.c_int, ctypes.c_double)
        lib.dBz.argtypes = (ctypes.c_int, ctypes.c_double)
        
        c1 = -self.cst * self.R * ppr[2]
        c2 = -self.cst * self.R * ppr[1]
        c3 = self.cst * self.R**2
        c4lf = ppr[0]**2+ppr[1]**2+ppr[2]**2 + 2*ppr[0]*self.d + self.d**2 + self.R**2
        c4rt = ppr[0]**2+ppr[1]**2+ppr[2]**2 - 2*ppr[0]*self.d + self.d**2 + self.R**2
        c5 = -2 * ppr[1] * self.R
        c6 = -2 * ppr[2] * self.R
        c7lf = self.cst * self.R * (ppr[0] + self.d)
        c7rt = self.cst * self.R * (ppr[0] - self.d)
        
        lfBx = integrate.quad(lib.dBx, 0, 2*pi, 
            args=(c1, c2, c3, c4lf, c5, c6))
        lfBy = integrate.quad(lib.dBy, 0, 2*pi, 
            args=(c7lf, c4lf, c5, c6))
        lfBz = integrate.quad(lib.dBz, 0, 2*pi, 
            args=(c7lf, c4lf, c5, c6))
        rtBx = integrate.quad(lib.dBx, 0, 2*pi, 
            args=(c1, c2, c3, c4rt, c5, c6))
        rtBy = integrate.quad(lib.dBy, 0, 2*pi, 
            args=(c7rt, c4rt, c5, c6))
        rtBz = integrate.quad(lib.dBz, 0, 2*pi, 
            args=(c7rt, c4rt, c5, c6))
        
        bx = lfBx[0] + rtBx[0]
        by = lfBy[0] + rtBy[0]
        bz = lfBz[0] + rtBz[0]
        
        if self.axiscf_theta == 0 and self.axiscf_phi == 0:
            return bx, by, bz
            
        return rotateVector([bx,by,bz], self.axiscf_theta, self.axiscf_phi)
    
    def calcBatP(self, p):
        """Calculate the B field as a result of the wire coils at a position P."""
        pcnt = [p[0] - self.Cpair[0], p[1] - self.Cpair[1], p[2] - self.Cpair[2]]
        ppr = rotateVector(pcnt, -self.axiscf_theta, -self.axiscf_phi)
        #ToDo Add X, Z Axis conditions
        #ToDo Separate loop pair and define loops one-by-one
        #Equations to Integrate
        lfdBx = lambda a: self.cst*(-self.R*ppr[2]*sin(a) - self.R*ppr[1]*cos(a) + 
            self.R**2) / (((ppr[0] + self.d)**2 + (ppr[1] - self.R*cos(a))**2 + (ppr[2] - 
            self.R*sin(a))**2)**(3/2))
        rtdBx = lambda a: self.cst*(-self.R*ppr[2]*sin(a) - self.R*ppr[1]*cos(a) + 
            self.R**2) / (((ppr[0] - self.d)**2 + (ppr[1] - self.R*cos(a))**2 + (ppr[2] - 
            self.R*sin(a))**2)**(3/2))
        lfdBy = lambda a: self.cst*(self.R*cos(a)*(ppr[0] + self.d)) / (((ppr[0] + 
            self.d)**2 + (ppr[1] - self.R*cos(a))**2 + (ppr[2] -
            self.R*sin(a))**2)**(3/2))
        rtdBy = lambda a: self.cst*(self.R*cos(a)*(ppr[0] - self.d)) / (((ppr[0] - 
            self.d)**2 + (ppr[1] - self.R*cos(a))**2 + (ppr[2] -
            self.R*sin(a))**2)**(3/2))
        lfdBz = lambda a: self.cst*(self.R*sin(a)*(ppr[0] + self.d)) / (((ppr[0] +
            self.d)**2 + (ppr[1] - self.R*cos(a))**2 + (ppr[2] -
            self.R*sin(a))**2)**(3/2))
        rtdBz = lambda a: self.cst*(self.R*sin(a)*(ppr[0] - self.d)) / (((ppr[0] -
            self.d)**2 + (ppr[1] - self.R*cos(a))**2 + (ppr[2] -
            self.R*sin(a))**2)**(3/2))
    
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
        
        if self.axiscf_theta == 0 and self.axiscf_phi == 0:
            return bx, by, bz
            
        return rotateVector([bx,by,bz], self.axiscf_theta, self.axiscf_phi)

class BField(object):
    """Define a B Field object containing the elements in BObjList."""
    def __init__(self, windObj, BObjList=[], name=None):
        """Initialize a B Field object."""
        self.windObj = windObj
        self.BObjList = BObjList[:]
        self.name = name
        #ToDo Define Particle list separate from BObjList
        
    def totalBatP(self, p):
        """Calculate total B at p due to all objects in BObjList.  Calculates them one at a time and adds them together."""
        Bx = By = Bz = 0
        for BObj in self.BObjList:
            bx, by, bz = BObj.calcBatP(p)
            Bx += bx; By += by; Bz += bz
        
        return [Bx, By, Bz]
    
    def drawBlines(self, windObj, p, pupbound=[None,None,None],
        plobound=[None,None,None], numiter=None, linelength=None, multlng=None):
        """Draw B field lines starting at po and ending at ####."""
        loopind = 0
        BoundBool = True
        while BoundBool:
            for i in [pupbound, plobound]:
                #ToDo Write code to check lower bound - maybe negative both sides
                if i[0] is not None:
                    BoundBool = BoundBool and p[0] <= i[0]
                if i[1] is not None:
                    BoundBool = BoundBool and p[1] <= i[1]
                if i[2] is not None:
                    BoundBool = BoundBool and p[2] <= i[2]
            if numiter is not None:
                loopind += 1
                BoundBool = BoundBool and loopind < numiter
            
            bx, by, bz = self.totalBatP(p)
            
            if linelength is not None:
                Blen, Bth, Bphi = cartesianToSpherical([bx,by,bz])
                bx, by, bz = sphericalToCartesian(linelength, Bth, Bphi)
            elif multlng is not None:
                bx *= multlng; by *= multlng; bz *= multlng
            
            drawLine(windObj, p, [bx,by,bz])
            p[0] += bx; p[1] += by; p[2] += bz