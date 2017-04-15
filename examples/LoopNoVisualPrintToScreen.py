from __future__ import division,print_function

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
    dt = 1e-10
    
    e1center = [-4.75,1,1]
    e1vel = [1000,1000,1000]
    wccenter = [0,0,0]
    loopaxis = [1,0,0]
    pos = []
    
    print(deltat)
    
    windObj1 = None
    wireCoils = WireCoilPair(windObj1, wccenter, loopaxis, 1, 1, 5, 5)    
    electron1 = Electron(windObj1, e1center, e1vel)
    
    B = BField(windObj1)
    B.BObjList.append(wireCoils)
    
    start = time.time()
    
    while t<=5.0001e-6:
        wireCoils.useC=True
        Barrayc = wireCoils.calcBatP(electron1.p)
        #print("B(c++): ",Barrayc[0],Barrayc[1],Barrayc[2])
        #wireCoils.useC=False
        #Barraypy = wireCoils.calcBatP(electron1.p)
        #print("B(py ): ",Barraypy[0],Barraypy[1],Barraypy[2])
        electron1.updP(Barrayc, dt)
        #if (Barraypy[1] != 0):
            #Ratio = [Barrayc[0]/Barraypy[0],Barrayc[1]/Barraypy[1],Barrayc[2]/Barraypy[2]]
            #print("Ratio: ",Ratio[0],Ratio[1],Ratio[2])
        #print("-----")
        
        t += dt
        if ind % 1000 == 0:
            tmp = [ind, t, time.time()-start]
            #print(tmp)
            print("Position: ",electron1.p[0],electron1.p[1],electron1.p[2]," | Ind: ",ind," | dT: ", dt)
            print("Velocity: ",electron1.v[0], electron1.v[1],electron1.v[2]," | Time: ", t)
            print("------")
            pos.append(np.array([t,electron1.p[0],electron1.p[1],electron1.p[2]]))
        ind += 1
    
    end = time.time()
    tottime = end - start
    print(tottime)
    print(deltat, usec)
    
    return [deltat, tottime, float(usec)]

if __name__ == "__main__":
    fstr = './../__TestResults__/TestNewC/'
    timed = []
    #timed.append(main(1e-8,fstr,False))
    timed.append(main(1e-10,fstr,True))
    #timed.append(main(1e-10,fstr,False))