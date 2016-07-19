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
    dt = 5*10**-9
    
    e1center = [-4.75,0,0]
    e1vel = [1000,1000,1000]
    wccenter = [0,0,0]
    loopaxis = [1,0,0]
    
    # Draw some things
    windObj1 = drawWindow(1920, 1080, e1center)
    relclockObj1 = drawTimeClock(windObj1, t)
    
    wireCoils = WireCoilPair(windObj1, wccenter, loopaxis, 1, 1, 5, 5)
    wireCoils.initDraw()
    
    electron1 = Electron(windObj1, e1center, e1vel)
    electron1.initDraw(10, 50)
    
    B = BField(windObj1)
    B.BObjList.append(wireCoils)
    B.BObjList.append(electron1)
    #Need to update drawBlines for appropriate boundary conditions, then uncomment below
    #for i in [[-5,0,-3.5],[-5,-3.5,0],[-5,0,3.5],[-5,3.5,0]]:
        #j = rotateVector(i,wireCoils.axis_theta,wireCoils.axis_phi)
        #j += wireCoils.C
        #B.drawBlines(windObj1, j, 0.1)
    
    while ((-10 + wccenter[0]) <= electron1.p[0] <= (10 + wccenter[0])) and ((-10 + 
            wccenter[1]) <= electron1.p[1] <= (10 + wccenter[1])) and ((-10 + 
            wccenter[2]) <= electron1.p[2] <= (10 + wccenter[2])):      
        FPSrate(10000)
        Barray = B.totalBatP(electron1.p)
        electron1.updP(Barray, dt)
        electron1.updDraw()

        t += dt                  #time increase [s]
        ind += 1                 #index increase
        updateTimeClock(windObj1, relclockObj1, t)
        #print electron1.p
    
    while True:
        FPSrate(10000)

if __name__ == "__main__":
    main()