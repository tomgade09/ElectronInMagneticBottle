from __future__ import division

import os, sys, inspect
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '/../')
sys.path.append(a + '/../vis/')

from classes import *
from math import pi

dt = 5e-9

parray = []
varray = []
Barray = []
ind = 0

po_x = (-4.75,0,0)
po_y = rotateVector(po_x[0], po_x[1], po_x[2], 0, pi/2)
#po_y = (po_x[0] + 5, po_x[1] + 5, po_x[2] + 5)
#po_z = rotateVector(po_x[0], po_x[1], po_x[2], -pi/2, 0)
po_z = ((po_x[0] + 5, po_x[1] + 5, po_x[2] + 5))
po_rot = rotateVector(po_x[0], po_x[1], po_x[2], -pi/4, pi / 4)

v_x = (1000,1000,1000)
v_y = rotateVector(v_x[0], v_x[1], v_x[2], 0, pi/2)
#v_z = rotateVector(v_x[0], v_x[1], v_x[2], -pi/2, 0)
v_z = v_x
v_rot = rotateVector(v_x[0], v_x[1], v_x[2], -pi/4, pi / 4)

wx = WireCoilPair(None, (0,0,0), (1,0,0), 1, 1, 5, 5)
wy = WireCoilPair(None, (0,0,0), (0,1,0), 1, 1, 5, 5)
wz = WireCoilPair(None, (-5,-5,-5), (1,0,0), 1, 1, 5, 5)
wrot = WireCoilPair(None, (0,0,0), (1,1,0), 1, 1, 5, 5)

elec_x = Particle(None, -1.602e-19, 9.109e-31, po_x, v_x)
elec_y = Particle(None, -1.602e-19, 9.109e-31, po_y, v_y)
elec_z = Particle(None, -1.602e-19, 9.109e-31, po_z, v_z)
elec_rot = Particle(None, -1.602e-19, 9.109e-31, po_rot, v_rot)

B_x = BField(None, [wx])
B_y = BField(None, [wy])
B_z = BField(None, [wz])
B_rot = BField(None, [wrot])

for i in range(0,9999,1):
    bx_x, by_x, bz_x = B_x.totalBatP((elec_x.px, elec_x.py, elec_x.pz))
    bx_y, by_y, bz_y = B_y.totalBatP((elec_y.px, elec_y.py, elec_y.pz))
    bx_z, by_z, bz_z = B_z.totalBatP((elec_z.px, elec_z.py, elec_z.pz))
    #bx_rot, by_rot, bz_rot = B_rot.totalBatP((elec_rot.px, elec_rot.py, elec_rot.pz))
    
    Barray.append([ind, [bx_x, by_x, bz_x], [bx_y, by_y, bz_y], [bx_z, by_z, bz_z]])
    
    elec_x.updP(bx_x, by_x, bz_x, dt)
    elec_y.updP(bx_y, by_y, bz_y, dt)
    elec_z.updP(bx_z, by_z, bz_z, dt)
    #elec_rot.updP(bx_rot, by_rot, bz_rot, dt)
    
    parray.append([ind, [elec_x.px, elec_x.py, elec_x.pz], [elec_y.px, elec_y.py, elec_y.pz], [elec_z.px, elec_z.py, elec_z.pz]])
    varray.append([ind, [elec_x.vx, elec_x.vy, elec_x.vz], [elec_y.vx, elec_y.vy, elec_y.vz], [elec_z.vx, elec_z.vy, elec_z.vz]])
    
    ind += 1

print elec_x.px - elec_y.py, elec_x.py + elec_y.px, elec_x.pz - elec_x.pz
print elec_x.px - elec_z.pz, elec_x.py - elec_z.py, elec_x.pz - elec_z.px

#if __name__ == "__main__":
    #main()