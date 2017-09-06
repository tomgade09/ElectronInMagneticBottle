from __future__ import division

from math import *

def fourthOrderRungeKutta_dy(callback, funcArgs): #requires [h, t_0 = 0, y_0x, y_0y, y_0z, ... whatever else you want here]
    #calculates a differential value for y - time starts at 0, return value is dy NOT y total

    if (funcArgs[1] != 0.0):
        print("Warning: t_0 is not set to 0.\nfourthOrderRungeKutta_dy returns differential value (dy) NOT y. Setting t_0 equal to 0.")
        funcArgs[1] = 0.0

    y_0x = funcArgs[2]                            #[dt, t_0, vx, vy, vz, q, m, px_0, py_0, pz_0] in this implementation
    y_0y = funcArgs[3]                            # h is dt, t_0 starts at 0, y is velocity, dy / dt is acceleration
    y_0z = funcArgs[4]                            # callback must return dy / dt - f(t, y) = dy / dt

    k1 = callback(funcArgs)                       # k1 = f(t_n, y_n), units of dy / dt

    funcArgs[1] = funcArgs[0] / 2                 # k2 = f(t_n + h/2, y_n + h/2 * k1)
    funcArgs[2] = y_0x + k1[0] * funcArgs[1]
    funcArgs[3] = y_0y + k1[1] * funcArgs[1]
    funcArgs[4] = y_0z + k1[2] * funcArgs[1]
    k2 = callback(funcArgs)

    funcArgs[2] = y_0x + k2[0] * funcArgs[1]      # k3 = f(t_n + h/2, y_n + h/2 * k2)
    funcArgs[3] = y_0y + k2[1] * funcArgs[1]
    funcArgs[4] = y_0z + k2[2] * funcArgs[1]
    k3 = callback(funcArgs)

    funcArgs[1] = funcArgs[0]                     # k4 = f(t_n + h, y_n + h * k3)
    funcArgs[2] = y_0x + k3[0] * funcArgs[1]
    funcArgs[3] = y_0y + k3[1] * funcArgs[1]
    funcArgs[4] = y_0z + k3[2] * funcArgs[1]
    k4 = callback(funcArgs)

    #normally, FORK is "y_0 + dy (the below value)", but this code just returns dy
    return (k1 + 2 * k2 + 2 * k3 + k4) * funcArgs[0] / 6

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
    x = rho * sin(theta) * cos(phi)
    y = rho * sin(theta) * sin(phi)
    z = rho * cos(theta)
    
    return [x, y, z]
    
def rotateVector(v, rot_theta, rot_phi):
    """Rotate a cartesian vector, v by rot_theta and rot_phi.  You can guess which corresponds to theta and which corresponds to phi."""
    if rot_theta == -pi / 2: #Z axis case
        xprm = - v[2]; yprm = v[1]; zprm = v[0]
        return xprm, yprm, zprm
    elif rot_theta == rot_phi == 0: #X axis case
        return v
    tmp_rho, tmp_theta, tmp_phi = cartesianToSpherical(v)
    
    return sphericalToCartesian(tmp_rho, tmp_theta + rot_theta, tmp_phi + rot_phi)
    
def cross3DandMult(a,b,c):
    return [(a[1] * b[2] - a[2] * b[1]) * c, (a[2] * b[0] - a[0] * b[2]) * c,
        (a[0] * b[1] - a[1] * b[0]) * c]