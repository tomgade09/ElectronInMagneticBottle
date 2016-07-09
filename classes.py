from __future__ import division

from VPyDraw import *
from scipy import integrate
from math import pi

class Particle(object):
    """Define a particle to be placed in the specified magnetic field.
    
    Attributes of the particle:
    - Charge
    - Mass
    - Position Vector (x,y,z)
    - Velocity Vector (vx,vy,vz)"""
    
    def __init__(self, windObj, charge, mass, po, vo):
        """Initiate a particle object."""
        self.charge = charge
        self.mass = mass
        self.px = po[0]
        self.py = po[1]
        self.pz = po[2]
        self.vx = vo[0]
        self.vy = vo[1]
        self.vz = vo[2]
        self.eom = self.charge / self.mass
        self.pic = None
        self.wind = windObj

    def initDraw(self, intrvl, traillng):
        if self.pic != None:
            print "Pic has already been initialized.  Use updDraw to change the position."
            return
        self.pic = drawParticlePic(self.wind, (self.px, self.py, self.pz), intrvl, traillng)
        
        return self.pic
        
    def updDraw(self):
        if self.pic == None:
            print "Pic has not been initialized.  Use initDraw to create a picture of the particle first."
            return
        updateParticlePic(self.wind, self.pic, (self.px, self.py, self.pz))
    
    def __updV(self, bx, by, bz, dt):
        """Calculate the new velocity of the particle based on the specified B field."""
        vXB_x = self.vy * bz - self.vz * by
        vXB_y = self.vz * bx - self.vx * bz
        vXB_z = self.vx * by - self.vy * bx
        ax = self.eom * vXB_x
        ay = self.eom * vXB_y
        az = self.eom * vXB_z
        self.vx += ax*dt
        self.vy += ay*dt
        self.vz += az*dt
        
    def updP(self, bx, by, bz, dt):
        """Calculate the new position based on the particle's velocity."""
        self.__updV(bx, by, bz, dt)
        self.px += self.vx * dt
        self.py += self.vy * dt
        self.pz += self.vz * dt
        
class WireCoilPair(object):
    """Define a pair of wire coils to create a B field.
    
    Attributes of the wire coils:
    - Number of wire turns: N
    - Current running through coils: I [A]
    - Radius of coils: R [cm]
    - Distance of center of loop from origin: d [cm]
    
    Note: This builds a pair of coils, parallel to one another, offset from the origin by a distance d.  There is no way to build a single wire loop, or to have the loops offset from one another.  Hence the name 'Wire Coil Pair'"""
    
    def __init__(self, windObj, C, axis, N, I, R, d):
        self.Cx = C[0]
        self.Cy = C[1]
        self.Cz = C[2]
        self.axis_x = axis[0]
        self.axis_y = axis[1]
        self.axis_z = axis[2]
        #self.axis_theta = 
        #self.axis_phi = 
        self.N = N
        self.I = I
        self.R = R
        self.d = d
        self.wind = windObj
        self.cst = float(self.N * self.I * 10**(-5))

    def initDraw(self):
        drawWireCoilPair(self.wind, self.d, self.R)
        
    def calcBatP(self, p):
        """Calculate the B field as a result of the wire coils at a position P."""    
        #Equations to Integrate
        lfdBx = lambda a: self.cst*(-self.R*p[2]*sin(a) - self.R*p[1]*cos(a) + self.R**2) / (((p[0] + self.d)**2 + (p[1] - self.R*cos(a))**2 + (p[2] - self.R*sin(a))**2)**(3/2))
        rtdBx = lambda a: self.cst*(-self.R*p[2]*sin(a) - self.R*p[1]*cos(a) + self.R**2) / (((p[0] - self.d)**2 + (p[1] - self.R*cos(a))**2 + (p[2] - self.R*sin(a))**2)**(3/2))
        lfdBy = lambda a: self.cst*(self.R*cos(a)*(p[0] + self.d)) / (((p[0] + self.d)**2 + (p[1] - self.R*cos(a))**2 + (p[2] - self.R*sin(a))**2)**(3/2))
        rtdBy = lambda a: self.cst*(self.R*cos(a)*(p[0] - self.d)) / (((p[0] - self.d)**2 + (p[1] - self.R*cos(a))**2 + (p[2] - self.R*sin(a))**2)**(3/2))
        lfdBz = lambda a: self.cst*(self.R*sin(a)*(p[0] + self.d)) / (((p[0] + self.d)**2 + (p[1] - self.R*cos(a))**2 + (p[2] - self.R*sin(a))**2)**(3/2))
        rtdBz = lambda a: self.cst*(self.R*sin(a)*(p[0] - self.d)) / (((p[0] - self.d)**2 + (p[1] - self.R*cos(a))**2 + (p[2] - self.R*sin(a))**2)**(3/2))
    
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
        
class BField(object):
    """Calculate a number of B Field related things based on all B field producing elements."""
    
    def __init__(self, windObj, BObjList):
        self.BObjList = BObjList
        self.windObj = windObj
        
    def totalBatP(self, p):
        Bx = 0
        By = 0
        Bz = 0
        for BObj in self.BObjList:
            bx, by, bz = BObj.calcBatP((p[0], p[1], p[2]))
            Bx += bx
            By += by
            Bz += bz
        
        return Bx, By, Bz
    
    def drawBlines(self, windObj, Po, linedist):
        px = Po[0]
        py = Po[1]
        pz = Po[2]
    
        while px < 5:
            bx, by, bz = self.totalBatP((px, py, pz))
            #normfact = linedist / sqrt(bx**2 + by**2 + bz**2)
            normfact = 10000
        
            bxn = bx * normfact
            byn = by * normfact
            bzn = bz * normfact
        
            drawLine(windObj, (px, py, pz), (bxn, byn, bzn))
        
            px += bxn
            py += byn
            pz += bzn