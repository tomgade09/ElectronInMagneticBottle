from __future__ import division

from VPyDraw import *
from scipy import integrate
from math import *
import numpy

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
    def __init__(self, windObj, charge, mass, po, vo):
        """Initiate a particle object."""
        self.wind = windObj
        self.q = charge
        self.mass = mass
        self.p = po
        self.v = vo
        self.eom = self.q / self.mass
        self.pic = None

    def initDraw(self, intrvl, traillng):
        """Draw a point object representing the particle in the object specified by self.wind.  intrvl represents how often to 'draw' a point.  traillng represents how long of a 'trail' to leave behind the current position of the particle."""
        if self.pic != None:
            print "Pic has already been initialized.  Use updDraw to change the position."
            return
        self.pic = drawParticlePic(self.wind, self.p, intrvl, traillng)
        return self.pic
        
    def updDraw(self):
        """Update the location of the point object drawn with initDraw.  Obviously, it can't be updated if it hasn't been initialized.  Use the initDraw function first."""        
        if self.pic == None:
            print "Pic has not been initialized.  Use initDraw to create a picture of the particle first."
            return
        updateParticlePic(self.wind, self.pic, self.p)
    
    def calcBatP(self, pB):
        """Calculate B at a point pB due to this particle."""
        pB = pB[:] #If particle.p is passed in, need to change pB to a list vs a pointer
        for i in range(len(pB)):
            pB[i] -= self.p[i]
        if abs(pB[0]) <= 10e-15 and abs(pB[1]) <= 10e-15 and abs(pB[2]) <=10e-15:
            return 0, 0, 0
        plen = sqrt(pB[0]**2 + pB[1]**2 + pB[2]**2)
        c = 10e-7 * self.q / plen**3 #Should be in m?  Need to get units right.
        b = numpy.cross(self.v, pB)
        for i in range(len(d)):
            b[i] *= c
        
        return b[0], b[1], b[2]
    
    def __updV(self, b, dt):
        """Calculate the new velocity of the particle based on the specified B field."""
        a = numpy.cross(self.v, b)
        for i in range(len(self.v)):
            self.v[i] += self.eom * a[i] * dt
        
    def updP(self, b, dt):
        """Calculate the new position based on the particle's velocity."""
        self.__updV(b, dt)
        for i in range(len(self.p)):
            self.p[i] += self.v[i] * dt
            
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
    def __init__(self, windObj, C, axis, N, I, R, d):
        """Initiate a WireCoilPair object."""
        self.C = C
        self.axis = axis
        self.N = N; self.I = I; self.R = R; self.d = d
        self.wind = windObj
        self.cst = float(self.N * self.I * 10**(-5))
        
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
                self.axis[0], self.axis[1], self.axis[2])
            self.axiscf_theta = self.axis_theta - (pi / 2)
            self.axiscf_phi = self.axis_phi

    def initDraw(self):
        """Draw the pair of Wire Coils."""
        if self.axis_theta == 0 and self.axis_phi == 0: #Z Axis calculate center points
            cntrt = [self.C[0], self.C[1], self.C[2] + self.d]
            cntlf = [self.C[0], self.C[1], self.C[2] - self.d]
        else: #All other cases calculate center points for loops
            cntrt = sphericalToCartesian(self.d, self.axis_theta, self.axis_phi)
            if self.axis_phi > pi:
                cntlf = sphericalToCartesian(self.d, pi-self.axis_theta, self.axis_phi-pi)
            else:
                cntlf = sphericalToCartesian(self.d, pi-self.axis_theta, self.axis_phi+pi)
            #cntrt = [cntrt[0] + self.C[0], cntrt[1] + self.C[1], cntrt[2] + self.C[2]]
            #cntlf = [cntlf[0] + self.C[0], cntlf[1] + self.C[1], cntlf[2] + self.C[2]]
            for i in range(len(cntrt)):
                cntrt[i] += self.C[i]; cntlf[i] += self.C[i]
        drawWireCoilPair(self.wind, self.C, self.axis, cntlf, cntrt, self.R)
        
    def calcBatP(self, p):
        """Calculate the B field as a result of the wire coils at a position P."""
        pcnt = [p[0] - self.C[0], p[1] - self.C[1], p[2] - self.C[2]]
        ppr = rotateVector(pcnt, -self.axiscf_theta, -self.axiscf_phi)
########Add X, Z Axis conditions
        
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
            
        bxprime, byprime, bzprime = rotateVector(bx, by, bz, self.axiscf_theta, 
            self.axiscf_phi)
        
        return bxprime, byprime, bzprime

class BField(object):
    """Define a B Field object containing the elements in BObjList."""
    def __init__(self, windObj, BObjList=[]):
        """Initialize a B Field object."""
        self.BObjList = BObjList
        self.windObj = windObj
        
    def totalBatP(self, p):
        """Calculate total B at p due to all objects in BObjList.  Calculates them one at a time and adds them together."""
        Bx = By = Bz = 0
        for BObj in self.BObjList:
            bx, by, bz = BObj.calcBatP(p)
            Bx += bx; By += by; Bz += bz
        
        return [Bx, By, Bz]
    
    def drawBlines(self, windObj, po, linedist):
        """Draw B field lines starting at po and ending at ####."""
        px = po[0]; py = po[1]; pz = po[2]
#########Need better boundary conditions, but not sure how to define at this time
        #Maybe pass in an argument for bounds
        while px < 5:
            bx, by, bz = self.totalBatP([px, py, pz])
            #normfact = linedist / sqrt(bx**2 + by**2 + bz**2)
            normfact = 10000
            bx *= normfact; by *= normfact; bz *= normfact        
            drawLine(windObj, [px, py, pz], [bx, by, bz])
            px += bx; py += by; pz += bz

def cartesianToSpherical(a):
    """Convert cartesian coords to spherical."""
    rho = sqrt(a[0]**2 + a[1]**2 + a[2]**2)
    if a[0] == 0 and a[1] == 0 and a[2] != 0: #Z axis case
        theta = phi = 0
        return rho, theta, phi
    
    theta = acos(a[2] / rho)
    phi = atan2(a[1], a[0])
    
    return rho, theta, phi

def sphericalToCartesian(rho, theta, phi):
    """Convert spherical coords to cartesian."""
    if theta == 0: #Z axis case
        x = y = 0
        z = rho
    else:
        x = rho * sin(theta) * cos(phi)
        y = rho * sin(theta) * sin(phi)
        z = rho * cos(theta)
    
    return [x, y, z]
    
def rotateVector(v, rot_theta, rot_phi):
    """Rotate a vector, v by rot_theta and rot_phi.  You can guess which corresponds to theta and which corresponds to phi."""
    if rot_theta == -pi / 2: #Z axis case
        xprm = - v[2]; yprm = v[1]; zprm = v[0]
        return xprm, yprm, zprm
    elif rot_theta == rot_phi == 0: #X axis case
        return v
    tmp_rho, tmp_theta, tmp_phi = cartesianToSpherical(v)
    vprm = sphericalToCartesian(tmp_rho, tmp_theta + rot_theta, tmp_phi +
        rot_phi)
    
    return vprm