#!/usr/bin/python

# Use file's current directory and append path to MinVR Module
import sys, os, inspect
fileDir = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
os.chdir(fileDir)
#sys.path.append("../../plugins/Python/src/python")

# Import MinVR and OpenGL Dependencies
#from MinVR import *
from OpenGL.GL import *
from OpenGL.GLU import *
import math

def Line(verts):
    glBegin(GL_LINE_STRIP)
    for vertex in verts:
        glVertex3fv(vertex)
        glEnd()

def drawSphere(s,color=(1.0,0,0),ndiv=100):
    quad = gluNewQuadric()
    gluQuadricNormals(quad, GLU_SMOOTH)
    glColor(color)
    gluSphere(quad,s,ndiv,ndiv)
    
class sphere(object):
    def __init__(self, pos=[0,0,0], radius = 1.0, color=[0,1,0], make_trail = False,
        trail_type="points", interval = 1, retain=0):
        
        self.trail = []
        self.trail.append(pos)
        self.radius = radius
        self.color = color
        self.retain = retain
        self.interval = interval
        self.make_trail = make_trail
        
    def draw(self):
        for p in self.trail[::self.interval]:
            glPushMatrix()
            glTranslatef(p[0], p[1], p[2])
            drawSphere(self.radius, self.color)
            glPopMatrix()
            
    def setPos(self, pos):
        if (not(self.make_trail)):
            self.trail = []
        if (len(self.trail) > self.retain):
            self.trail.pop(0)
        self.trail.append(pos)

def Torus(numc, numt, radius=1, thickness=0.01):
    twopi = 2 * math.pi;
    for i in range(0, numc):
        glBegin(GL_QUAD_STRIP)
        for j in range(0, numt+1):
            for k in range(1,-1, -1):
                s = (i + k) % numc + 0.5
                t = j % numt
                x = (radius+thickness*math.cos(s*twopi/numc))*math.cos(t*twopi/numt);
                y = (radius+thickness*math.cos(s*twopi/numc))*math.sin(t*twopi/numt);
                z = thickness * math.sin(s * twopi / numc);
                glVertex3f(x, y, z)
        glEnd()
        
class ring(object):
    def __init__(self, pos=(0,0,0), axis=(0,0,1), radius=1, thickness=1):
        self.pos = pos
        self.axis = axis
        self.radius = radius
        self.thickness = thickness
        
    def draw(self):
        glPushMatrix()
        x = 0
        y = 0
        z = -1
        c = (y * self.axis[2] - z * self.axis[1], z * self.axis[0] - x * self.axis[2], x * self.axis[1] - y * self.axis[0])
        glRotatef(90.0, c[0], c[1], c[2]);
        glTranslatef(self.pos[0], self.pos[1], self.pos[2])
        Torus(5, 30, self.radius, self.thickness)
        glPopMatrix()