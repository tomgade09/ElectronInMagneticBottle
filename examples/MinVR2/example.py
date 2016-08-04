#!/usr/bin/python

from __future__ import division

# Use file's current directory and append path to MinVR Module
import sys, os, inspect
fileDir = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
os.chdir(fileDir)
sys.path.append("../..")
sys.path.append("../../vis")

import numpy as np

# --------- MinVR Implementation -----------------

# Import MinVR and OpenGL Dependencies
from MinVR2 import *
from OpenGL.GL import *
from OpenGLDraw import *
from classes import *

def DrawPic(ObjList):
	for obj in ObjList:
		obj.draw()

# Application class which handles events and rendering
class App(VREventHandler, VRRenderHandler):
	def __init__(self, PicList, pos):
		self.loop = True
		self.rotateAngle = 0.0
		self.PicList = PicList
		self.pos = pos
		self.direction = 1.0

	# Called when an event is passed from MinVR
	def onVREvent(self, eventName):
		print eventName
		if eventName == "/KbdEsc_Down":
			self.loop = False
		elif eventName == "/KbdRight_Down" or eventName == "/KbdRight_Repeat":
			self.rotateAngle += 0.1
		elif eventName == "/KbdLeft_Down" or eventName == "/KbdLeft_Repeat":
			self.rotateAngle -= 0.1

	# Renders the scene
	def onVRRenderScene(self, renderState):
		b=[0,0,0]
		b[0],b[1],b[2] = W1.calcBatP([Electron1.p[0],Electron1.p[1],Electron1.p[2]])
		Electron1.updP(b,5e-8)
		self.PicList[0].setPos((Electron1.p[0],Electron1.p[1],Electron1.p[2]))
		glClear(GL_COLOR_BUFFER_BIT)
		width = renderState.getValue("WindowWidth","/")
		height = renderState.getValue("WindowHeight","/")
		ratio = width / height
		glViewport(0, 0, width, height)
		glClear(GL_COLOR_BUFFER_BIT)
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		gluPerspective(60, width/height, 0.1, 50.0)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		glTranslatef(0.0,0.0, -15)
		DrawPic(self.PicList)

# ----------- Main program ------------------

Electron1 = Electron(None,[-4.75,0,0],[1000,1000,1000])
W1 = WireCoilPair(None,[0,0,0],[1,0,0],1,1,5,5)

# Define Objects
e1 = sphere(pos=[0,0,0],radius=0.1)
r1 = ring(pos=[-5,0,0],radius=5)
r2 = ring(pos=[5,0,0],radius=5)

VisObjList = [e1, r1, r2]
pos = np.array([0.0,0.0,0.0])
direction = 1

# Create application
app = App(VisObjList, pos)

# Create VRMain instance passing in vrsetup configuration
#config = sys.argv[1]
config = 'desktop.xml'
vrmain = VRMain(config)

# Add event handler and render handler
vrmain.addEventHandler(app)
vrmain.addRenderHandler(app)

# Main loop
while app.loop:
	vrmain.mainloop()

# Shutdown MinVR
vrmain.shutdown()