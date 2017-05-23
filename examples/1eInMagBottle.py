from __future__ import division

import os, sys, inspect, time
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '/../')
sys.path.append(a + '/../vis/')
os.chdir(a)

from VPyDraw import *

def main():
    ind = 0                            #Index (Calculation Counter)
    dt = 1*10**-7                      #Loses sufficient resolution at any faster (larger dt) than 1e-7
                                       #Set smaller for greater accuracy/larger for quicker simulation
    
    e1center = [-4,3,0]
    e1vel = [500,4000,4000]            #Need high vparallel to vperpendicular (to B field) ratio to confine
    #e2center = [-4,0,0]
    wccenter = [0,0,0]
    loopaxis = [1,0,0]

    # Draw some things
    windObj1 = drawWindow(1920, 1080, [0,0,0])
    relclockObj1 = drawTimeClock(windObj1, [-6.5,0,0], 0)
    
    wireCoils = WireCoilPair(windObj1, wccenter, loopaxis, 1, 1, 5, 5, useC=True)
    wireCoils.initDraw()
    
    electron1 = Electron(windObj1, e1center, e1vel)
    electron1.initDraw(1, 5000)
    points(pos=e1center, size=5, color=color.red) #draw a point at electron1 start point - for reference

    #electron2 = Electron(windObj1, e2center, e1vel)
    #electron2.initDraw(1, 5000)
    #points(pos=e2center, size=5, color=color.red)

    Field = Fields(windObj1, dt, [wireCoils], [electron1])#, electron2])
    
    for i in [[-5,3.5,0],[-5,0,3.5],[-5,-3.5,0],[-5,0,-3.5]]:
        j = rotateVector(i,wireCoils.axiscf_theta,wireCoils.axiscf_phi) + wireCoils.Cpair
        Field.drawBlines(j, numiter=5750, multlng=500)
    
    windObj1.center = (electron1.p[0], electron1.p[1], electron1.p[2])
    
    wireCoils.drawBfromGC(electron1, Field)

    start = time.time()
    while ind < 2500000: #Need to specify a number of iterations or loop would never end
    #while ((-10 + wccenter[0]) <= electron1.p[0] <= (10 + wccenter[0])) and ((-10 + 
            #wccenter[1]) <= electron1.p[1] <= (10 + wccenter[1])) and ((-10 + 
            #wccenter[2]) <= electron1.p[2] <= (10 + wccenter[2])):
        FPSrate(10000)
        
        Field.updateParticleP_V()
        electron1.updDraw()
        #electron2.updDraw()
        
        ind += 1 #Increment Index
        
        updateTimeClock(windObj1, relclockObj1, Field.t)
        pauseOnKey('\n', windObj1) #Pause on enter key
        if ind % 1000 == 0:
            windObj1.center = (electron1.p[0], electron1.p[1], electron1.p[2])
            #print(ind, time.time() - start)
    
    print("Done!")
    while True: #Done with calculations - don't close window, leave it open to continue looking at it
        FPSrate(30)

if __name__ == "__main__":
    main()