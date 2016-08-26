from __future__ import division

from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import math
import numpy as np

#####Code within comments by Dan Orban.  Some modifications by Tom Gade.
###Class for sphere and code to draw
def drawSphere(r,color=(0,1,0),ndiv=100):
    quadGLObj = gluNewQuadric()
    gluQuadricNormals(quadGLObj, GLU_SMOOTH)
    glColor(color)
    gluSphere(quadGLObj,r,ndiv,ndiv)
    
class sphere(object):
    """Define a sphere object to be drawn by OpenGL."""
    def __init__(self, pos=[0,0,0], radius=1.0, color=[0,1,0], make_trail=False,
        trail_type="points", interval=1, retain=0):
        
        self.trail = []
        self.trail.append(pos)
        self.radius = radius
        self.color = color
        self.retain = retain * interval
        self.interval = interval
        self.make_trail = make_trail
        self.trail_type = trail_type
        
    def draw(self):
        for p in self.trail[::self.interval]:
            glMatrixMode(GL_MODELVIEW)
            glPushMatrix()
            glTranslatef(p[0], p[1], p[2])
            drawSphere(self.radius, self.color)
            glPopMatrix()
            
    def setPos(self, pos):
        if (not(self.make_trail)):
            self.trail = []
        if (len(self.trail) > self.retain):
            self.trail.pop(0)
        self.trail.append((pos[0],pos[1],pos[2])) #ensure appended object is a tuple

###Ring Class and Draw code
def Torus(numc, numt, radius=1, thickness=0.01,color=(0.75,0.75,0.75)):
    twopi = 2 * math.pi;
    for i in range(0, numc):
        glBegin(GL_QUAD_STRIP)
        glColor3f(color[0], color[1], color[2]);
        for j in range(0, numt+1):
            for k in range(1,-1, -1):
                s = (i + k) % numc + 0.5
                t = j % numt
                x = (radius+thickness*math.cos(s*twopi/numc))*math.cos(t*twopi/numt)
                y = (radius+thickness*math.cos(s*twopi/numc))*math.sin(t*twopi/numt)
                z = thickness * math.sin(s * twopi / numc)
                glVertex3f(x, y, z)
        glEnd()
        
class ring(object):
    """Create a ring object."""
    def __init__(self, pos=(0,0,0), axis=(1,0,0), radius=1, thickness=0.01):
        self.pos = pos
        self.axis = axis
        self.radius = radius
        self.thickness = thickness
        
    def draw(self):
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        xyz = np.array([0,0,-1])
        c = np.cross(xyz,self.axis)
        glTranslatef(self.pos[0], self.pos[1], self.pos[2])
        glRotatef(90, c[0], c[1], c[2])
        Torus(5, 30, self.radius, self.thickness)
        glPopMatrix()

###Line draw code (set of two vertices)
def Line(verts): #Original code by Danny Orban.
    glBegin(GL_LINE_STRIP)
    for vertex in verts:
        glVertex3fv(vertex)
    glEnd()
#####End Dan's code.

###Functions called on by classes.py
#Core functions to be defined:

#def drawWindow(wd, ht, cent, axistype="points"): Need to decide how I'm drawing window - GLUT or PyGame?
    #Could I sneak the code in for lighting, depth, perspective, etc into here?

#def drawParticlePic(windObj, po, intrvl, traillng, col):
    ##somehow select the window (don't know how to specify in OpenGL)
    #particle = sphere(pos=po, interval=intrvl, retain=traillng, color=col)
    #particle.draw()
    
    #return particle

#def updateParticlePic(windObj, partObj, p):
    ##somehow select the window (don't know how to specify in OpenGL)
    #partObj.setPos(p)
    #partObj.draw

#def drawTimeClock(windObj, po, to):

#def updateTimeClock(windObj, relclockObj, relt):

#def drawLine(windObj, p, laxis, headwd=0.005, headln=0.001, shaftwd=0.01,
    #col=[1,0,0]):
    ##somehow select the window (don't know how to specify in OpenGL)
    ##set color to col
    #vert = []
    #vert.append(p)
    #vert.append((p[0] + laxis[0], p[1] + laxis[1], p[2] + laxis[2]))
    #Line(vert)
    
#def drawWireCoilPair(windObj, C, norm, cntlf, cntrt, R):
    ##somehow select the window (don't know how to specify in OpenGL)
    ##set color to gray (maybe (0.85,0.85,0.85) but play with it)
    #rleft = ring(pos=cntlf, axis=norm, radius=R)
    #rright = ring(pos=cntrt, axis=norm, radius=R)
    
    #return [rleft, rright]

#def FPSrate(fps):

#def displayFlip():