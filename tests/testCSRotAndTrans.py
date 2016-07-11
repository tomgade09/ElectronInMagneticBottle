from __future__ import division

import os, sys, inspect
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '\\..\\')

from classes import *

px = 1000
py = 1000
pz = 1000

thetacf = -pi / 2
phicf = 0

rho, theta, phi = cartesianToSpherical(px, py, pz)
print rho, theta, phi
x, y, z = sphericalToCartesian(rho, theta + thetacf, phi + phicf)

print x, y, z