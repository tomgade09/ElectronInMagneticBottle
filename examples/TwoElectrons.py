# Tom Gade
# Electron in a Magnetic Bottle
# v4.1.2 11Jul16
# Works on Python 2.7 ONLY (as of now), and VPython 6
# Written as part of a Senior research project at the Citadel, 2008-2009

from __future__ import division

import os, sys, inspect
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '/../')
sys.path.append(a + '/../vis/')

from classes import *
# Only use one of these at a time to avoid namespace conflicts
from VPyDraw import *
#from pygletDraw import *
#from OpenGLDraw import *

def main():
    ind = 0                            #Index (Calculation Counter)
    t = 0                              #Initial time [s]
    BList = []
    dt = 5*10**-9
    
    e1center = (5, 0.25, 5)
    e1vel = (-1000, 1000, 1000)
    
    # Draw some things
    windObj1 = drawWindow(1920, 1080, e1center)
    relclockObj1 = drawTimeClock(windObj1, t)
    
    wireCoils = WireCoilPair(windObj1, (5, 5, 5), (0, 1, 0), 1, 1, 5, 5)
    wireCoils.initDraw()
    BList.append(wireCoils)
    
    electron1 = Particle(windObj1, -1.76*10**11, 1, e1center, e1vel)
    electron1.initDraw(10, 50)
    
    #electron2 = Particle(windObj1, -1.76*10**11, 1, (4.75, 0, 0), (-1000, -1000, -1000))
    #electron2.initDraw(10, 50)
    #electron2.color=color.yellow
    
    B = BField(windObj1, BList)
    #Use rotateVector to rotate to appropriate start point
    #B.drawBlines(windObj1, (-5, 0, -3.5), 0.1)
    #B.drawBlines(windObj1, (-5, -3.5, 0), 0.1)
    #B.drawBlines(windObj1, (-5, 0, 3.5), 0.1)
    #B.drawBlines(windObj1, (-5, 3.5, 0), 0.1)
    
    #while (-10 <= x <= 10) and (-10 <= y <= 10) and (-10 <= z <= 10):
    while t < 0.001:
        FPSrate(10000)
        Bx, By, Bz = B.totalBatP((electron1.px, electron1.py, electron1.pz))
        #Bx, By, Bz = wireCoils.calcBatP((electron1.px, electron1.py, electron1.pz))
        #Bxs, Bys, Bzs = wireCoils.calcBatP((electron2.px, electron2.py, electron2.pz))
        
        electron1.updP(Bx, By, Bz, dt)
        #electron2.updP(Bxs, Bys, Bzs, dt)
        t += dt                  #time increase [s]
        ind += 1                 #index increase
        electron1.updDraw()
        #electron2.updDraw()
        updateTimeClock(windObj1, relclockObj1, t)

if __name__ == "__main__":
    main()