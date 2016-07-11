from __future__ import division

import os, sys, inspect
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '/../')

from classes import *

px = 10
py = 100
pz = 1000

thetacf = 0
phicf = pi / 2

#rho, theta, phi = cartesianToSpherical(px, py, pz)
#print rho, theta, phi
#x, y, z = sphericalToCartesian(rho, theta + thetacf, phi + phicf)

#print x, y, z

xx, yy, zz = rotateVector(px, py, pz, thetacf, phicf)

print xx, yy, zz