import ctypes

from sys import platform as _platform
import os, inspect
libName = 'MagBottleC'
apipath = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
if _platform == "linux" or _platform == "linux2":
    libFilePath = apipath + '/../lib/' + libName + '.so'
elif _platform == "win32":
    libFilePath = apipath + '/../lib/' + libName + '.dll'
elif _platform == "darwin": #All of MagBottlePy is untested on Darwin, but should work
    libFilePath = apipath + '/../lib/' + libName + '.dylib'

MBCapi = ctypes.CDLL(libFilePath)

def initSimulationC(dt, WC1x, WC1y, WC1z, WC2x, WC2y, WC2z, N, I, R):
    #DLLEXPORT BField* initSimulation(double dt, double WC1x, double WC1y, double WC1z, double WC2x, double WC2y, double WC2z, int N, double I, double R)
    MBCapi.initSimulation.argtypes = (ctypes.c_double, 
                                      ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                      ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                      ctypes.c_int, ctypes.c_double, ctypes.c_double)
    MBCapi.initSimulation.restype = ctypes.c_void_p #hopefully fake out Python - don't need to access
    
    return MBCapi.initSimulation(dt, WC1x, WC1y, WC1z, WC2x, WC2y, WC2z, N, I, R)

def createAnElectronC(p_o, v_o): #DLLEXPORT Particle** createAnElectron(double Px, double Py, double Pz, double Vx, double Vy, double Vz)
    MBCapi.createAnElectron.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double)
    MBCapi.createAnElectron.restype = ctypes.c_void_p #hopefully fake out Python - don't need to access

    return MBCapi.createAnElectron(p_o[0], p_o[1], p_o[2], v_o[0], v_o[1], v_o[2])

def updateParticlePositionC(particle, BFieldptr): #DLLEXPORT void updateParticlePosition(Particle* particle, BField* B)
    MBCapi.updateParticlePosition.argtypes = (ctypes.c_void_p, ctypes.c_void_p)
    MBCapi.updateParticlePosition.restype = None
    
    MBCapi.updateParticlePosition(particle, BFieldptr)

def returnSimTimeC(BFieldPtr): #DLLEXPORT double returnSimTime(BField* B)
    MBCapi.returnSimTime.argtypes = (ctypes.c_void_p,) #hopefully fake out Python - don't need to access
    MBCapi.returnSimTime.restype = ctypes.c_double

    return MBCapi.returnSimTime(BFieldPtr)

def returnParticleMassC(particle): #DLLEXPORT double returnParticleMass(Particle* particle)
    MBCapi.returnParticleMass.argtypes = (ctypes.c_void_p,) #hopefully fake out Python - don't need to access
    MBCapi.returnParticleMass.restype = ctypes.c_double
    
    return MBCapi.returnParticleMass(particle)

def returnParticleChargeC(particle): #DLLEXPORT double returnParticleCharge(Particle** particleArray, int particleIndex)
    MBCapi.returnParticleCharge.argtypes = (ctypes.c_void_p,) #hopefully fake out Python - don't need to access
    MBCapi.returnParticleCharge.restype = ctypes.c_double
    
    return MBCapi.returnParticleCharge(particle)

def returnParticlePositionC(particle): #DLLEXPORT double* returnParticlePosition(Particle** particleArray, int particleIndex)
    MBCapi.returnParticlePosition.argtypes = (ctypes.c_void_p,) #hopefully fake out Python - don't need to access
    MBCapi.returnParticlePosition.restype = ctypes.POINTER(ctypes.c_double)

    a = MBCapi.returnParticlePosition(particle)
    ret = [a[0], a[1], a[2]]

    deleteDblArray(a)

    return ret

def returnParticleVelocityC(particle): #DLLEXPORT double* returnParticleVelocity(Particle** particleArray, int particleIndex)
    MBCapi.returnParticleVelocity.argtypes = (ctypes.c_void_p,) #hopefully fake out Python - don't need to access
    MBCapi.returnParticleVelocity.restype = ctypes.POINTER(ctypes.c_double)

    a = MBCapi.returnParticleVelocity(particle)
    ret = [a[0], a[1], a[2]]

    deleteDblArray(a)

    return ret

def returnTotalBatPC(BFieldPtr, P):
    MBCapi.returnTotalBatP.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double)
    MBCapi.returnTotalBatP.restype = ctypes.POINTER(ctypes.c_double)

    a = MBCapi.returnTotalBatP(BFieldPtr, P[0], P[1], P[2])
    ret = [a[0], a[1], a[2]]

    deleteDblArray(a)

    return ret

def incTimeC(BFieldPtr):
    MBCapi.increaseTimeBydt.argtypes = (ctypes.c_void_p,)
    MBCapi.increaseTimeBydt.restype = None

    MBCapi.increaseTimeBydt(BFieldPtr)

def deleteDblArray(arrayPtr):
    MBCapi.deleteDblArray.argtypes = (ctypes.c_void_p,)
    MBCapi.deleteDblArray.restype = None

    MBCapi.deleteDblArray(arrayPtr)

#DLLEXPORT void deleteObjects(BObject* B, Particle** particleArray, int numberOfParticles)

def cInit(self):
    #check wirecoils for C speedup code
    for BObj in self.BObjList:
        if BObj.useC:
            print("Warning: Both Wire Coil Pair and overall simulation is set to use C speedup.  These two are not compatible.  Setting WireCoil variable useC to False.")
            BObj.useC = False
    self.Cobject = ctypes.c_void_p #Store the pointer to the C Fields object here

    if (self.BObjList == []):
        print("Error: add wire coils to list at instantiation of Fields.  Fields is improperly formed and will not function correctly.")
        return

    #Assumes Wire Coil Pair is the first (maybe only) object in BObjList (Python)
    wcp = self.BObjList[0]
    self.Cobject = initSimulationC(self.dt, wcp.Cleft[0], wcp.Cleft[1], wcp.Cleft[2], wcp.Cright[0], wcp.Cright[1], wcp.Cright[2], wcp.N, wcp.I, wcp.R)

    for Part in self.particleList: #For now, particles need to be electrons
        if (not((Part.q == -1.60217657e-19) and (Part.mass == 9.10938356e-31))):
            print("Error: For right now, only electrons are supported with C code.  Results will be inaccurate.")
        Part.Cobject = ctypes.c_void_p #Store the pointer to a single C particle object here (1 C object per Python Object)
        Part.Cobject = createAnElectronC(Part.p, Part.v)