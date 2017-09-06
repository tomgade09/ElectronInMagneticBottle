from __future__ import division

import os, sys, inspect, time
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '/../')
sys.path.append(a + '/../vis/')
sys.path.append(a + '/../cExtension/')
os.chdir(a)

useCLib = True

from VPyDraw import *

#def initSimulationPy(): #BField* initSimulation(double dt, double WC1x, double WC1y, double WC1z, double WC2x, double WC2y, double WC2z, int N, double I, double R)
#def createAnElectronPy(): #Particle** createAnElectron()
#def updateParticlePosition(particleArray, particleIndex, BFieldptr): #void updateParticlePosition(Particle** particleArray, int particleIndex, BField* B)
#def returnSimTimePy(BFieldPtr): #double returnSimTime(BField* B)
#def returnParticlePositionPy(BFieldPtr): #double* returnParticlePosition(Particle** particleArray, int particleIndex)
#def returnParticleVelocityPy(BFieldPtr): #double* returnParticleVelocity(Particle** particleArray, int particleIndex)

def cmain():
    ind = 0
    dt = 1*10**-7

    BFieldPtr = initSimulationPy(dt, -5.0, 0.0, 0.0, 5.0, 0.0, 0.0, 1, 1.0, 5.0)
    ParticlePtr = createAnElectronPy()

    windObj1 = drawWindow(1920, 1080, returnParticlePositionPy(ParticlePtr, 0))
    relclockObj1 = drawTimeClock(windObj1, [-6.5,0,0], 0)

    electron1 = Electron(windObj1, returnParticlePositionPy(ParticlePtr, 0), returnParticleVelocityPy(ParticlePtr, 0))
    electron1.initDraw(1, 250)
    points(pos=returnParticlePositionPy(ParticlePtr, 0), size=5, color=color.red)

    #for i in [[-5,3.5,0],[-5,0,3.5],[-5,-3.5,0],[-5,0,-3.5]]:
        #Field.drawBlines(i, numiter=5750, multlng=500)

    #wireCoils.drawBfromGC(electron1, Field, numiter=2250)

    while ind < 2500000:
        FPSrate(10000)
        
        updateParticlePositionPy(ParticlePtr, 0, BFieldPtr)
        electron1.p = returnParticlePositionPy(ParticlePtr, 0)
        electron1.updDraw()
        
        ind += 1
        
        updateTimeClock(windObj1, relclockObj1, returnSimTimePy(BFieldPtr))
        pauseOnKey('\n', windObj1)
        if ind % 1000 == 0:
            windObj1.center = (electron1.p[0], electron1.p[1], electron1.p[2])
    
    print("Done!")
    while True:
        FPSrate(30)

def main():
    ind = 0
    dt = 1*10**-7
    
    e1center = [-4,3,0]
    e1vel = [500,4000,4000]
    wccenter = [0,0,0]
    loopaxis = [1,0,0]

    windObj1 = drawWindow(1920, 1080, e1center)
    relclockObj1 = drawTimeClock(windObj1, [-6.5,0,0], 0)

    wireCoils = WireCoilPair(windObj1, wccenter, loopaxis, 1, 1, 5, 5, useC=False)
    wireCoils.initDraw()
    
    electron1 = Electron(windObj1, e1center, e1vel)
    #electron1.initDraw(1, 250)
    points(pos=e1center, size=5, color=color.red)

    electron2 = Electron(windObj1, e1center, e1vel)
    electron2.initDraw(1, 250, color.yellow)

    Field = Fields(windObj1, dt, [wireCoils], [electron1])
    
    for i in [[-5,3.5,0],[-5,0,3.5],[-5,-3.5,0],[-5,0,-3.5]]:
        j = rotateVector(i,wireCoils.axiscf_theta,wireCoils.axiscf_phi) + wireCoils.Cpair
        Field.drawBlines(j, numiter=5750, multlng=500)
    
    wireCoils.drawBfromGC(electron1, Field, numiter=2250)

    #C++ stuff
    BFieldPtr = initSimulationPy(dt, -5.0, 0.0, 0.0, 5.0, 0.0, 0.0, 1, 1.0, 5.0)
    ParticlePtr = createAnElectronPy()

    start = time.time()
    while ind < 2500000:
        FPSrate(10000)
        
        #C++ stuff
        updateParticlePositionPy(ParticlePtr, 0, BFieldPtr)
        electron2.p = returnParticlePositionPy(ParticlePtr, 0)
        electron2.updDraw()

        #Field.updateParticleP_V()
        #electron1.updDraw()
        
        ind += 1
        
        updateTimeClock(windObj1, relclockObj1, Field.t)
        pauseOnKey('\n', windObj1)
        if ind % 1000 == 0:
            windObj1.center = (electron1.p[0], electron1.p[1], electron1.p[2])
    
    print("Done!")
    while True:
        FPSrate(30)

if __name__ == "__main__":
    main()
    #cmain()