from __future__ import division

import os, sys, inspect
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '/../')
sys.path.append(a + '/../vis/')
os.chdir(a)

from classes import *
#Only use one of these at a time to avoid namespace conflicts
#from VPyDraw import *
from OpenGLDraw import *
import pygame
from pygame.locals import *

size = width, height = 1280, 720

pygame.init()
scene1 = pygame.display.set_mode(size, DOUBLEBUF|OPENGL)

electron1 = sphere(pos=(0,0,0), radius = 0.01, interval=2, retain=20, make_trail=True)
r1 = ring(pos=(-5,0,0), axis=(0,0,1), radius=5, thickness=0.01)
r2 = ring(pos=(5,0,0), axis=(1,0,0), radius=5, thickness=0.01)

while True:
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glEnable(GL_DEPTH_TEST) #Depth test, checks whether things are in front of another
    glDepthFunc(GL_LEQUAL) #Won't display things behind
    glClearDepth(1.0)
    glEnable(GL_LIGHTING)
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)
    glEnable(GL_COLOR_MATERIAL)

    lightKa = (.2, .2, .2, 1.0) #Lighting
    lightKd = (.7, .7, .7, 1.0)
    lightKs = (1, 1, 1, 1)
    glLightfv(GL_LIGHT0, GL_AMBIENT, lightKa)
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightKd)
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightKs)
    lightPos = (0.5, 0, 3, 1)
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos)
    glEnable(GL_LIGHT0)

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) #Clears screen
    width = 1280
    height = 720
    ratio = width / height;
    glViewport(0, 0, width, height); #Set to same as window ht, wd - set viewable area
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION); #Setting perspective of view
    glLoadIdentity(); #Load identity matrix
    gluPerspective(45, ((1.0*width)/(1.0*height)), 0.1, 50.0) #

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0,0.0,-10.0) #Moving eye out 5 units
    
    electron1.draw()
    r1.draw()
    r2.draw()
    pygame.display.flip()