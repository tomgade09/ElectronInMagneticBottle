from __future__ import division

import os, sys, inspect, time
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '/../')
sys.path.append(a + '/../vis/')
os.chdir(a)

from VPyDraw import *

def main():
    ind = 0                            #Index (Calculation Counter)
    t = 0                              #Initial time [s]
    dt = 1*10**-7                      #Loses sufficient resolution at any faster (larger dt) than 1e-7
    
    e1center = [-4,3.5,0]
    e1vel = [1000,4000,4000]         #Need high vparallel to vperpendicular (to B field) ratio to confine
    wccenter = [0,0,0]
    loopaxis = [1,0,0]

    # Draw some things
    windObj1 = drawWindow(1920, 1080, [0,0,0])
    relclockObj1 = drawTimeClock(windObj1, [-6.5,0,0], t)
    
    wireCoils = WireCoilPair(windObj1, wccenter, loopaxis, 1, 1, 5, 5, useC=True)
    wireCoils.initDraw()
    
    electron1 = Electron(windObj1, e1center, e1vel)
    electron1.initDraw(1, 5000)
    
    #points(pos=(0,5.78,0), size=5, color=color.cyan)

    B = BField(windObj1)
    B.BObjList.append(wireCoils)
    #B.BObjList.append(electron1) #Messes with mag field lines.  Don't need for now.
    
    for i in [[-5,3.5,0],[-5,0,3.5],[-5,-3.5,0],[-5,0,-3.5]]:#,e1center]:#, [-5,-2.5*np.sqrt(2),2.5*np.sqrt(2)],[-5,-2.5*np.sqrt(2),-2.5*np.sqrt(2)],[-5,2.5*np.sqrt(2),-2.5*np.sqrt(2)]]:
        j = rotateVector(i,wireCoils.axiscf_theta,wireCoils.axiscf_phi) + wireCoils.Cpair
        B.drawBlines(j, numiter=5750, multlng=500)#pupbound=[5,None,None]
    
    windObj1.center = (electron1.p[0], electron1.p[1], electron1.p[2])
    
    start = time.time()
    while ind < 2500000: #Need to specify a number of iterations or loop would never end
    #while ((-10 + wccenter[0]) <= electron1.p[0] <= (10 + wccenter[0])) and ((-10 + 
            #wccenter[1]) <= electron1.p[1] <= (10 + wccenter[1])) and ((-10 + 
            #wccenter[2]) <= electron1.p[2] <= (10 + wccenter[2])):
        FPSrate(10000)
        
        vdrift = np.array(B.totalDriftsF(B, wireCoils, electron1, dt)) * dt / electron1.mass
        electron1.foRKvCrossB(B, dt)
        electron1.v = [electron1.v[0] + vdrift[0], electron1.v[1] + vdrift[1], electron1.v[2] + vdrift[2]]
        electron1.updDraw()
        t += dt
        ind += 1
        updateTimeClock(windObj1, relclockObj1, t)
        #pauseOnKey('\n')
        if windObj1.kb.keys:
            key = windObj1.kb.getkey()
            if key == '\n':
                windObj1.waitfor('\n keydown')
                windObj1.kb.queue = []
        if ind % 1000 == 0:
            windObj1.center = (electron1.p[0], electron1.p[1], electron1.p[2])
            #print(ind, time.time() - start)
    
    print("Done!")
    while True:
        FPSrate(30)

if __name__ == "__main__":
    main()