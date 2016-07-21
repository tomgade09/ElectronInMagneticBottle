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
    dt = 5*10**-8
    
    evel = [1000,1000,1000]
    e1center = [-10.75,0,0]
    e2center = [1.25,0,0]
    wc1center = [-6,0,0]
    wc2center = [6,0,0]
    loopaxis = [1,0,0]
    
    # Draw some things
    windObj1 = drawWindow(1920, 1080, [0,0,0])
    relclockObj1 = drawTimeClock(windObj1, [0,6.5,0], t)
    
    wireCoils1 = WireCoilPair(windObj1, wc1center, loopaxis, 1, 1, 5, 5)
    wireCoils1.initDraw()
    
    wireCoils2 = WireCoilPair(windObj1, wc2center, loopaxis, 1, 1, 5, 5)
    wireCoils2.initDraw()
    
    electron1 = Electron(windObj1, e1center, evel)
    electron1.initDraw(10, 50, color.red)
    
    electron2 = Electron(windObj1, e2center, evel)
    electron2.initDraw(10, 50, color.yellow)
    
    #Define 2 separate B Fields that don't interact with each other.
    B1 = BField(windObj1)
    B1.BObjList.append(wireCoils1)
    B1.BObjList.append(electron1)
    
    B2 = BField(windObj1)
    B2.BObjList.append(wireCoils2)
    B2.BObjList.append(electron2)
    #####Below needs a test
    for i in [[-5,0,-3.5],[-5,-3.5,0],[-5,0,3.5],[-5,3.5,0]]:
        j = rotateVector(i,wireCoils.axis_theta,wireCoils.axis_phi)
        j += wireCoils.Cpair
        k = sphericalToCartesian(2*wireCoils.d,wireCoils.axis_theta,wireCoils.axis_phi)
        B.drawBlines(windObj1, j, pupbound=k, multlng=10000)
    
    while ((-10 + wc1center[0]) <= electron1.p[0] <= (10 + wc1center[0])) and ((-10 + 
            wc1center[1]) <= electron1.p[1] <= (10 + wc1center[1])) and ((-10 + 
            wc1center[2]) <= electron1.p[2] <= (10 + wc1center[2])):
        FPSrate(10000)
        electron1.updP(B1.totalBatP(electron1.p), dt)
        electron1.updDraw()
        electron2.updP(B2.totalBatP(electron2.p), dt)
        electron2.updDraw()
        
        t += dt                  #time increase [s]
        ind += 1                 #index increase
        updateTimeClock(windObj1, relclockObj1, t)
        #print electron1.p
        #print electron2.p
        
    while True:
        FPSrate(10000)

if __name__ == "__main__":
    main()