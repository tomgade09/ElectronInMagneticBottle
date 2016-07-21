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

e1center = [-4.75,0,0]
e1vel = [10000,10000,10000]
wccenter = [0,0,0]
loopaxis = [1,0,0]

ind = 0
t = 0
dt = 5*10**-9

windObj1 = drawWindow(1920, 1080, e1center)
relclockObj1 = drawTimeClock(windObj1, [-6.5,0,0], t)

wireCoils1 = WireCoilPair(windObj1, wccenter, loopaxis, 1, 1, 5, 5)
wireCoils1.initDraw()

wireCoils2 = WireCoilPair(windObj1, wccenter, [0,1,0], 1, 1, 5, 5)
wireCoils2.initDraw()

electron1 = Electron(windObj1, e1center, e1vel)
electron1.initDraw(10, 50)

B = BField(windObj1)
B.BObjList.append(wireCoils1)
B.BObjList.append(wireCoils2)

for i in [[-5,0,-3.5],[-5,-3.5,0],[-5,0,3.5],[-5,3.5,0]]:
    j = rotateVector(i,wireCoils1.axiscf_theta,wireCoils1.axiscf_phi) + wireCoils1.Cpair
    B.drawBlines(windObj1, j, pupbound=[5,5,5], multlng=10000)
    
while ((-10 + wccenter[0]) <= electron1.p[0] <= (10 + wccenter[0])) and ((-10 + 
        wccenter[1]) <= electron1.p[1] <= (10 + wccenter[1])) and ((-10 + 
        wccenter[2]) <= electron1.p[2] <= (10 + wccenter[2])):      
    FPSrate(10000)
    Barray = B.totalBatP(electron1.p)
    electron1.updP(Barray, dt)
    electron1.updDraw()

    t += dt
    ind += 1
    updateTimeClock(windObj1, relclockObj1, t)
    
while True:
    FPSrate(30)