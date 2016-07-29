from __future__ import division

import os, sys, inspect
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '/../')
sys.path.append(a + '/../vis/')
os.chdir(a)

from classes import *

import cProfile, pstats, io

ind = 0                            #Index (Calculation Counter)
t = 0                              #Initial time [s]
dt = 5*10**-8
    
e1center = [-4.75,0,0]
e1vel = [1000,1000,1000]
wccenter = [0,0,0]
loopaxis = [1,0,0]
windObj1 = None

wireCoils = WireCoilPair(windObj1, wccenter, loopaxis, 1, 1, 5, 5)
electron1 = Electron(windObj1, e1center, e1vel)
B = BField(windObj1)
B.BObjList.append(wireCoils)
B.BObjList.append(electron1)

pr = cProfile.Profile()
pr.enable()
Barray = B.totalBatP(electron1.p)
electron1.updP(Barray, dt)
pr.disable()

s = io.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print s.getvalue()