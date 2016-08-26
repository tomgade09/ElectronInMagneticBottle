from __future__ import division,print_function

import os, sys, inspect
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '/../')

from classes import *

e1center = (0, -4.75, 0)
e1vel = (-1000, 1000, 1000)
BList = []
dt = 5*10**-9
ind = 0
axis = (0, 1, 0)

wireCoils = WireCoilPair(None, (0, 0, 0), axis, 1, 1, 5, 5)
BList.append(wireCoils)
electron1 = Particle(None, -1.76*10**11, 1, e1center, e1vel)
B = BField(None, BList)

print("Axis:")
print(axis)

for i in range(0,5,1):
    Bx, By, Bz = wireCoils.calcBatP((electron1.px, electron1.py, electron1.pz))
    electron1.updP(Bx, By, Bz, dt)
    ind += 1
    
    print(ind)
    print("B:")
    print(Bx, By, Bz)
    print("V:")
    print(electron1.vx, electron1.vy, electron1.vz)
    print("P:")
    print(electron1.px, electron1.py, electron1.pz)
    print("==============")