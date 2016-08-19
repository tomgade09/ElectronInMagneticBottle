from __future__ import division

import os, sys, inspect
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '/../')
sys.path.append(a + '/../vis/')
os.chdir(a)

from classes import *

import csv
import time

def main(deltat,foldstring,usec):
    ind = 0                            #Index (Calculation Counter)
    t = 0                              #Initial time [s]
    dt = deltat
    
    e1center = [-4.75,0,0]
    e1vel = [1000,1000,1000]
    wccenter = [0,0,0]
    loopaxis = [1,0,0]
    pos = []
    
    print deltat0
    
    windObj1 = None    
    wireCoils = WireCoilPair(windObj1, wccenter, loopaxis, 1, 1, 5, 5)    
    electron1 = Electron(windObj1, e1center, e1vel)
    
    B = BField(windObj1)
    B.BObjList.append(wireCoils)
    #B.BObjList.append(electron1)
    
    start = time.time()
    
    while t<=0.000001:
        if usec != True:
            Barray = wireCoils.calcBatP(electron1.p)
        elif usec == True:
            Barray = wireCoils.calcBatPinC(electron1.p)
        
        electron1.updP(Barray, dt)

        t += dt
        ind += 1
        if ind % int(1e-8 / deltat) == 0:
            tmp = [ind, t, time.time()-start]
            print tmp
            pos.append(np.array([t,electron1.p[0],electron1.p[1],electron1.p[2]]))

    resultfile = open(foldstring + str(deltat) + str(usec) + '.csv','w+')    
    wr = csv.writer(resultfile)
    wr.writerows(pos)
    
    end = time.time()
    tottime = end - start
    print tottime
    print deltat, usec
    
    return [deltat, tottime, float(usec)]

if __name__ == "__main__":
    fstr = './TestC/'
    timed = []
    timed.append(main(1e-8,fstr,False))
    timed.append(main(1e-8,fstr,True))
    for j in range(-9,-12,-1):
        for i in [5,2,1]:
            for k in [False, True]:
                timed.append(main(i*(10**j),fstr,k))
    
    timefile = open(fstr + 'time.csv','w+')
    wrt = csv.writer(timefile)
    wrt.writerows(timed)