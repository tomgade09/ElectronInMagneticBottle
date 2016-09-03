from __future__ import division

import os, sys, inspect, time
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '/../')
sys.path.append(a + '/../vis/')
os.chdir(a)

from classes import *
# Only use one of these at a time to avoid namespace conflicts
#from VPyDraw import *
from OpenGLDraw import *

def main():
    ind = 0                            #Index (Calculation Counter)
    t = 0                              #Initial time [s]
    dt = 1*10**-7
    
    e1center = [-4.75,0,0]
    e1vel = [1000,1000,1000]
    wccenter = [0,0,0]
    loopaxis = [1,0,0]

    # Draw some things
    #windObj1 = drawWindow(1920, 1080, e1center)
    #relclockObj1 = drawTimeClock(windObj1, [-6.5,0,0], t)
    windObj1 = None
    wireCoils = WireCoilPair(windObj1, wccenter, loopaxis, 1, 1, 5, 5, useC=True)
    #wireCoils.initDraw()
    
    electron1 = Electron(windObj1, e1center, e1vel)
    #electron1.initDraw(10, 50)
    
    B = BField(None)
    B.BObjList.append(wireCoils)
    B.particleList.append(electron1)
    #####Below needs a test
    #for i in [[-5,0,-3.5],[-5,-3.5,0],[-5,0,3.5],[-5,3.5,0]]:
        #j = rotateVector(i,wireCoils.axiscf_theta,wireCoils.axiscf_phi) + wireCoils.Cpair
        #B.drawBlines(windObj1, j, pupbound=[5,None,None], multlng=10000)
    
    start = time.time()
    while ind < 2500000: #Need to specify a number of iterations or loop would never end
    #while ((-10 + wccenter[0]) <= electron1.p[0] <= (10 + wccenter[0])) and ((-10 + 
            #wccenter[1]) <= electron1.p[1] <= (10 + wccenter[1])) and ((-10 + 
            #wccenter[2]) <= electron1.p[2] <= (10 + wccenter[2])):
        FPSrate(10000)
        
        electron1.foRKvCrossB(B,dt)
        #electron1.updDraw()
        t += dt
        ind += 1
        #updateTimeClock(windObj1, relclockObj1, t)
        if ind % 1000 == 0:
            #windObj1.center = (electron1.p[0], electron1.p[1], electron1.p[2])
            print(ind, time.time() - start)
    
    while True:
        FPSrate(30)

if __name__ == "__main__":
    main()