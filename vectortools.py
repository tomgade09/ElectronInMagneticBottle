from __future__ import division

from math import *

def cartesianToSpherical(a):
    """Convert cartesian coords to spherical."""
    rho = sqrt(a[0]**2 + a[1]**2 + a[2]**2)
    if a[0] == 0 and a[1] == 0 and a[2] != 0: #Z axis case
        return rho, 0, 0
    theta = acos(a[2] / rho)
    phi = atan2(a[1], a[0])
    
    return rho, theta, phi

def sphericalToCartesian(rho, theta, phi):
    """Convert spherical coords to cartesian."""
    if theta == 0: #Z axis case
        return 0, 0, rho
    else:
        x = rho * sin(theta) * cos(phi)
        y = rho * sin(theta) * sin(phi)
        z = rho * cos(theta)
    
    return [x, y, z]
    
def rotateVector(v, rot_theta, rot_phi):
    """Rotate a vector, v by rot_theta and rot_phi.  You can guess which corresponds to theta and which corresponds to phi."""
    if rot_theta == -pi / 2: #Z axis case
        xprm = - v[2]; yprm = v[1]; zprm = v[0]
        return xprm, yprm, zprm
    elif rot_theta == rot_phi == 0: #X axis case
        return v
    tmp_rho, tmp_theta, tmp_phi = cartesianToSpherical(v)
    
    return sphericalToCartesian(tmp_rho, tmp_theta + rot_theta, tmp_phi + rot_phi)