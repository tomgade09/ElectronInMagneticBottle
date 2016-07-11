from __future__ import division

from VPyDraw import *
from scipy import integrate
from math import *

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
        self.N = N
        self.I = I
        self.R = R
        self.d = d
        self.wind = windObj
        self.cst = float(self.N * self.I * 10**(-5))
        
        if self.axis_x != 0 and self.axis_y == 100 and self.axis_z == 0:
            self.axiscf_theta = self.axiscf_phi = 0
        else:
            self.axis_rho, self.axis_theta, self.axis_phi = cartesianToSpherical(
                self.axis_x, self.axis_y, self.axis_z)
        
            self.axiscf_theta = (pi / 2) - self.axis_theta
            self.axiscf_phi = -self.axis_phi

    def initDraw(self):
        cntrt = sphericalToCartesian(self.d, self.axis_theta, self.axis_phi)
        if self.axis_phi > pi:
            cntlf = sphericalToCartesian(self.d, pi - self.axis_theta, self.axis_phi - pi)
        else:
            cntlf = sphericalToCartesian(self.d, pi - self.axis_theta, self.axis_phi + pi)
        cntrt = (cntrt[0] + self.Cx, cntrt[1] + self.Cy, cntrt[2] + self.Cz)
        cntlf = (cntlf[0] + self.Cx, cntlf[1] + self.Cy, cntlf[2] + self.Cz)
        drawWireCoilPair(self.wind, (self.Cx, self.Cy, self.Cz), 
            (self.axis_x, self.axis_y, self.axis_z), cntlf, cntrt, self.R)
        
    def calcBatP(self, p):
        """Calculate the B field as a result of the wire coils at a position P."""
        ppr = (p[0] - self.Cx, p[1] - self.Cy, p[2] - self.Cz)
        
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
        
        #if self.axiscf_theta == 0 and self.axiscf_phi == 0:
            #return bx, by, bz
        
        brho, btheta, bphi = cartesianToSpherical(bx, by, bz)
        
        bxprime, byprime, bzprime = sphericalToCartesian(
            brho, btheta - self.axiscf_theta, bphi - self.axiscf_phi)
        
        return bxprime, byprime, bzprime

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

def cartesianToSpherical(x, y, z):
    rho = sqrt(x**2 + y**2 + z**2)
    theta = acos(z / rho)
    phi = atan2(y, x)
    
    return rho, theta, phi

def sphericalToCartesian(rho, theta, phi):
    x = rho * sin(theta) * cos(phi)
    y = rho * sin(theta) * sin(phi)
    z = rho * cos(theta)
    
    return x, y, z