# Tom Gade
# Electron in a Magnetic Bottle
# v2.0.21Nov13
# Works on Python 2.7 ONLY (as of now), and VPython 6
# Written as part of a research project at the Citadel, 2008-2009

# Modules
from __future__ import division
from math import pi
from visual import *
from scipy import integrate

# Constants and index
d = 5                              #Half distance between coils [cm]
R = 5                              #Radius of coil [cm]
N = 1                              #Number of turns of wire in coil [None]
I = 1                              #Current [A]
muo = float(4*pi*10**(-7))         #Permeability of Free Space [T*m/A]
eta = float(N*I*10**(-5))          #Constant "muo * I * 10^2 / 4 * pi" [T*cm]
eom = 1.76*10**11                  #Charge to mass ratio [C/kg]
Px = -4.75                         #Initial Electron Position [cm]
Py = 0
Pz = 0
vx = 100                           #Initial Electron Velocity [cm/s]
vy = 100
vz = 100
dt = 5*10**-9                      #time slice [s]
t = 0                              #initial time [s]
Bx = 0                             #Comps of Total B at P [T]
By = 0
Bz = 0
ind = 0                              #Index (Calculation Counter)

# Graphical Necessities
scene1 = display(title='Electron in Magnetic Bottle', autocenter=0, width=1920,
                 height=1200, center=(Px,Py,Pz), exit=0, range = (15,15,15))
xaxpt=[0,1,2,3,4,5,6,7,8,9,10]
yaxpt=[0,1,2,3,4,5,6,7,8,9,10]
zaxpt=[0,1,2,3,4,5,6,7,8,9,10]
ring(pos=(-d,0,0), axis=(1,0,0), radius=R, thickness=0.01)
ring(pos=(d,0,0), axis=(1,0,0), radius=R, thickness=0.01)
xlbl = label(pos=(10,1,0), text='x')
ylbl =  label(pos=(1,10,0), text='y')
zlbl =  label(pos=(0,1,10), text='z')
for i in range(0,10,1):
    points(pos=(xaxpt[i],0,0), size=5, color=color.cyan)
    points(pos=(0,yaxpt[i],0), size=5, color=color.cyan)
    points(pos=(0,0,zaxpt[i]), size=5, color=color.cyan)
#for i in range(0,100,1):
    #points(pos=(Px+i*10**-4,0,0), size=1, color=color.cyan)
#for i in range(0,250,5):
    #points(pos=(Px+i*10**-4,0,0), size=3, color=color.blue)
#for i in range(0,500,10):
    #points(pos=(Px+i*10**-4,0,0), size=5, color=color.orange)
#label(pos=(-5.5,0.2,0),height=0.1,text='cyan =\n1*10^-4 cm',
      #color=color.cyan,box=false)
#label(pos=(-5.5,0,0), height=0.1, text='blue =\n5*10^-4 cm',
      #color=color.blue,box=false)
#label(pos=(-5.5,-0.2,0), height=0.1, text='orange =\n1*10^-3 cm',
      #color=color.orange,box=false)

electron=sphere(pos=(Px,Py,Pz), radius=0.0000001, color=color.green,
    make_trail=True, trail_type="points", interval=1, retain=50000)

# Approximation for Magnetic Field Lines
flxpt = []
flypt = []
flxpt = range(-50,50,1)  #Field Line Approximation
for i in range(len(flxpt)):
    flxpt.append(flxpt.pop(0)/10)
for i in range(len(flxpt)):
    flypt.append(4/cosh(.3*flxpt[i]))
for i in range(len(flxpt)):
    points(pos=(flxpt[i],flypt[i],0), size=5, color=color.red)
    points(pos=(flxpt[i],-flypt[i],0), size=5, color=color.red)
    points(pos=(flxpt[i],0,flypt[i]), size=5, color=color.red)
    points(pos=(flxpt[i],0,-flypt[i]), size=5, color=color.red)

# Motion Calculator and Graphing
while t <= .001:
    rate(10000)
    # Define Functions Which Will Be Integrated Later in the Script
    lfdBx = lambda a: eta*(-R*Pz*sin(a)-R*Py*cos(a)+R**2)/(((Px+d)**2+
                                                       (Py-R*cos(a))**2+
                                                       (Pz-R*sin(a))**2)**(3/2))
    rtdBx = lambda a: eta*(-R*Pz*sin(a)-R*Py*cos(a)+R**2)/(((Px-d)**2+
                                                       (Py-R*cos(a))**2+
                                                       (Pz-R*sin(a))**2)**(3/2))
    lfdBy = lambda a: eta*(R*cos(a)*(Px+d))/(((Px+d)**2+
                                        (Py-R*cos(a))**2+
                                        (Pz-R*sin(a))**2)**(3/2))
    rtdBy = lambda a: eta*(R*cos(a)*(Px-d))/(((Px-d)**2+
                                        (Py-R*cos(a))**2+
                                        (Pz-R*sin(a))**2)**(3/2))
    lfdBz = lambda a: eta*(R*sin(a)*(Px+d))/(((Px+d)**2+
                                        (Py-R*cos(a))**2+
                                        (Pz-R*sin(a))**2)**(3/2))
    rtdBz = lambda a: eta*(R*sin(a)*(Px-d))/(((Px-d)**2+
                                        (Py-R*cos(a))**2+
                                        (Pz-R*sin(a))**2)**(3/2))

    # Integrate Functions Iteratively
    lfBx = integrate.quad(lfdBx, 0, 2*pi)
    lfBy = integrate.quad(lfdBy, 0, 2*pi)
    lfBz = integrate.quad(lfdBz, 0, 2*pi)
    rtBx = integrate.quad(rtdBx, 0, 2*pi)
    rtBy = integrate.quad(rtdBy, 0, 2*pi)
    rtBz = integrate.quad(rtdBz, 0, 2*pi)

    # Add B from Left Loop to Right Loop (Using Superposition Principle)
    Bx = lfBx[0] + rtBx[0]
    By = lfBy[0] + rtBy[0]
    Bz = lfBz[0] + rtBz[0]

    # Calculating B's Effect on the Electron's Velocity
    vxBx = vy*Bz - vz*By     #cmpts of v cross B [T*cm/s]
    vxBy = vz*Bx - vx*Bz
    vxBz = vx*By - vy*Bx
    ax = -eom * vxBx         #cmpts of acceleration [cm/s2]
    ay = -eom * vxBy
    az = -eom * vxBz
    vx += ax*dt              #cmpts of new velocity [cm/s]
    vy += ay*dt
    vz += az*dt

    # Iterate Time Forward dt, Calculate New Position, Index
    Px += vx*dt              #cmpts of new position [cm]
    Py += vy*dt
    Pz += vz*dt
    t += dt                  #time increase [s]
    ind += 1                 #index increase
    electron.pos=(Px,Py,Pz)
