from __future__ import division

import os, sys, inspect, time
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '/../')
sys.path.append(a + '/../vis/')
os.chdir(a)

from classes import *

def main():
    ind = 0                            #Index (Calculation Counter)
    t = 0                              #Initial time [s]
    dt = 5*10**-9
    
    e1center = [-4.75,0,0]
    e1vel = [1000,1000,1000]
    wccenter = [0,0,0]
    loopaxis = [1,0,0]
    
    wireCoils = WireCoilPair(None, wccenter, loopaxis, 1, 1, 5, 5, useC=True)
    
    electron1 = Electron(None, e1center, e1vel)
    
    B = BField(None)
    B.BObjList.append(wireCoils)
    
    start = time.time()
    while ((-10 + wccenter[0]) <= electron1.p[0] <= (10 + wccenter[0])) and ((-10 + 
            wccenter[1]) <= electron1.p[1] <= (10 + wccenter[1])) and ((-10 + 
            wccenter[2]) <= electron1.p[2] <= (10 + wccenter[2])):
        Barray = B.totalBatP(electron1.p)
        electron1.updP(Barray, dt)
        
        t += dt
        ind += 1
        if ind % 1000 == 0:
            print(ind, electron1.p, time.time() - start)

if __name__ == "__main__":
    main()