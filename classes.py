from __future__ import division,print_function

# Import visual library.  Only use one of these at a time to avoid namespace conflicts.
#from VPyDraw import *
#from OpenGLDraw import * #Only use for MinVR - Virtual Reality Display

from scipy import integrate
from math import *
import numpy as np
from vectortools import *
import version
import ctypes, os, sys, inspect

classpath = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))

from sys import platform as _platform
libName = 'WireCoilB'
if _platform == "linux" or _platform == "linux2":
    libFilePath = classpath + '/lib/' + libName + '.so'
elif _platform == "win32":
    libFilePath = classpath + '/lib/' + libName + '.dll'
elif _platform == "darwin": #All of MagBottlePy is untested on Darwin
    libFilePath = classpath + '/lib/' + libName + '.dylib'

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
    
    def calcBatP(self, pB):
        """Calculate B at a point pB due to this particle."""
        #pB = pB[:] #If particle.p is passed in, need to change pB to a list vs a pointer
        #pB -= np.array(self.p)
        pB = [pB[0] - self.p[0], pB[1] - self.p[1], pB[2] - self.p[2]]
        if abs(pB[0]) <= 10e-15 and abs(pB[1]) <= 10e-15 and abs(pB[2]) <= 10e-15:
            return 0, 0, 0 #B due to point charge is 0 at the location of the pt charge
        c = 10e-5 * self.q / sqrt(pB[0]**2 + pB[1]**2 + pB[2]**2)**3
        #b = np.cross(self.v, pB)
        #b = np.array(b) * c
        #b = [(self.v[1] * pB[2] - self.v[2] * pB[1]) * c,
            #(self.v[2] * pB[0] - self.v[0] * pB[2]) * c,
            #(self.v[0] * pB[1] - self.v[1] * pB[0]) * c]
        B = cross3DandMult(self.v, pB, c)
        return B[0], B[1], B[2]
    
    #@profile #For running with line/memory profiler
    def __updV(self, b, dt):
        """Calculate the new velocity of the particle based on the specified B field."""
        #dv = np.cross(self.v, b) * self.eom * dt #Apparently not faster
        #self.v += dv
        #dv = [self.v[1]*b[2]-self.v[2]*b[1], #Below cuts 20% time off
            #self.v[2]*b[0]-self.v[0]*b[2],
            #self.v[0]*b[1]-self.v[1]*b[0]]
        #for i in range(3):
            #self.v[i] += (dv[i] * self.eom * dt)
        #self.v = [self.v[0] + (self.v[1] * b[2] - self.v[2] * b[1]) * self.eom * dt,
            #self.v[1] + (self.v[2] * b[0] - self.v[0] * b[2]) * self.eom * dt,
            #self.v[2] + (self.v[0] * b[1] - self.v[1] * b[0]) * self.eom * dt]
        vpr = cross3DandMult(self.v, b, self.eom * dt)
        self.v = [self.v[0] + vpr[0], self.v[1] + vpr[1], self.v[2] + vpr[2]]
        
    #@profile #For running with line/memory profiler
    def updP(self, b, dt):
        """Calculate the new position based on the particle's velocity."""
        self.__updV(b, dt)
        #self.p += self.v * dt #For some reason, doesn't work.
        #for i in range(3): #Slower code than below
            #self.p[i] += self.v[i] * dt
        self.p = [self.p[0] + self.v[0] * dt, self.p[1] + self.v[1] * dt,
            self.p[2] + self.v[2] * dt]
    
    #@profile #For running with line/memory profiler
    def foRKvCrossB(self, BFieldObj, h): #Highly experimental!  Not sure if algorithm is implemented right.
    #ToDo verify code/algorithm
    #ToDo spin out totalBatP into separate processes
    # Two methods execute at similar speeds - in Py2.7, top executes about 5% faster
    # In Py 3.5, top executes about 10% slower
    # Not exactly sure what to pick
    # Even with MKL on a Core2Duo, numpy(bottom) is slower - not sure about newer procs
        B1 = BFieldObj.totalBatP(self.p[:])
        P23 = [self.p[0] + self.v[0] * h / 2, self.p[1] + self.v[1] * h / 2, self.p[2] + 
            self.v[2] * h / 2]
        B23 = BFieldObj.totalBatP(P23)
        P4 = [self.p[0] + self.v[0] * h, self.p[1] + self.v[1] * h, self.p[2] + 
            self.v[2] * h]
        B4 = BFieldObj.totalBatP(P4)
        V1 = self.v[:]
        k1 = cross3DandMult(V1,B1,self.eom * h)
        V2 = [self.v[0] + k1[0] / 2, self.v[1] + k1[1] / 2, self.v[2] + k1[2] / 2]
        k2 = cross3DandMult(V2,B23,self.eom * h)
        V3 = [self.v[0] + k2[0] / 2, self.v[1] + k2[1] / 2, self.v[2] + k2[2] / 2]
        k3 = cross3DandMult(V3,B23,self.eom * h)
        V4 = [self.v[0] + k3[0], self.v[1] + k3[1], self.v[2] + k3[2]]
        k4 = cross3DandMult(V4,B4,self.eom * h)
        
        k = [(k1[0] + 2 * (k2[0] + k3[0]) + k4[0]) / 6,(k1[1] + 2 * (k2[1] + k3[1]) + 
            k4[1]) / 6, (k1[2] + 2 * (k2[2] + k3[2]) + k4[2]) / 6]
        
        self.v = [self.v[0] + k[0], self.v[1] + k[1], self.v[2] + k[2]]
        self.p = [self.p[0] + self.v[0] * h, self.p[1] + self.v[1] * h, self.p[2] + 
            self.v[2] * h]
        
        #B2 = BFieldObj.totalBatP(self.p + np.array(self.v) * h / 2)
        #k11 = self.eom * np.cross(self.v,BFieldObj.totalBatP(self.p)) * h
        #k22 = self.eom * np.cross(self.v + k11 / 2, B2) * h
        #k33 = self.eom * np.cross(self.v + k22 / 2, B2) * h
        #k44 = self.eom * np.cross(self.v + k33, BFieldObj.totalBatP(self.p +
            #np.array(self.v) * h)) * h
            
        #self.v += (k11 + 2 * (k22 + k33) + k44) / 6
        #self.p += self.v * h
        
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
    def __init__(self, windObj, Cpair, axis, N, I, R, d, useC=False, name=None):
        """Initiate a WireCoilPair object."""
        self.wind = windObj
        self.Cpair = Cpair
        self.axis = axis
        self.N = N; self.I = I; self.R = R; self.d = d
        self.name = name
        self.cst = float(self.N * self.I * 10**(-5))
        self.pic = None
        self.useC = useC
                
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
                self.Cleft = sphericalToCartesian(self.d, pi - self.axis_theta,
                    self.axis_phi - pi)
            else:
                self.Cleft = sphericalToCartesian(self.d, pi - self.axis_theta,
                    self.axis_phi + pi)
            #self.Cright += self.Cpair; self.Cleft += self.Cpair
            self.Cright = [self.Cright[0] + self.Cpair[0], self.Cright[1] + self.Cpair[1],
                self.Cright[2] + self.Cpair[2]]
            self.Cleft = [self.Cleft[0] + self.Cpair[0], self.Cleft[1] + self.Cpair[1],
                self.Cleft[2] + self.Cpair[2]]
    
    def calcBatPinC(self, p):
        """Calculate the B field as a result of the wire coils at a position P."""
        pcnt = [p[0] - self.Cpair[0], p[1] - self.Cpair[1], p[2] - self.Cpair[2]]
        ppr = rotateVector(pcnt, -self.axiscf_theta, -self.axiscf_phi)
        
        lib = ctypes.CDLL(libFilePath)        
        lib.dBx.restype = ctypes.c_double
        lib.dBy.restype = ctypes.c_double
        lib.dBz.restype = ctypes.c_double
        lib.dBx.argtypes = (ctypes.c_int, ctypes.c_double)
        lib.dBy.argtypes = (ctypes.c_int, ctypes.c_double)
        lib.dBz.argtypes = (ctypes.c_int, ctypes.c_double)
        
        c1 = -self.cst * self.R * ppr[2]
        c2 = -self.cst * self.R * ppr[1]
        c3 = self.cst * self.R**2
        d2 = ppr[0]**2 + ppr[1]**2 + ppr[2]**2
        c4lf = d2 + 2*ppr[0]*self.d + self.d**2 + self.R**2
        c4rt = d2 - 2*ppr[0]*self.d + self.d**2 + self.R**2
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
    
    def calcBatPinPy(self, p):
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
    
    def calcBatP(self, p):
        if self.useC == True:
            return self.calcBatPinC(p)
        return self.calcBatPinPy(p)

class BField(object):
    """Define a B Field object containing the elements in BObjList."""
    def __init__(self, windObj, BObjList=[], PartList=[], name=None):
        """Initialize a B Field object."""
        self.windObj = windObj
        self.BObjList = BObjList[:]
        self.particleList = PartList[:]
        self.name = name
        
    def totalBatP(self, p):
        """Calculate total B at p due to all objects in BObjList.  Calculates them one at a time and adds them together."""
        Bx = By = Bz = 0
        for BObj in self.BObjList:
            bx, by, bz = BObj.calcBatP(p)
            Bx += bx; By += by; Bz += bz
        for BObj in self.particleList:
            bx, by, bz = BObj.calcBatP(p)
            Bx += bx; By += by; Bz += bz
        
        return [Bx, By, Bz]