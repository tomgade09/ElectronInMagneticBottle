#ifndef MBC_API_H
#define MBC_API_H

#include "bfield/WireCoil_c.h"
#include "particle/Particle_c.h"
#include "bfield/BField_c.h"
#include "tools/vectortools.h"

#define DLLEXPORT extern "C" __declspec(dllexport)

DLLEXPORT BField* initSimulation(double dt, double WC1x, double WC1y, double WC1z, double WC2x, double WC2y, double WC2z, int N, double I, double R);
DLLEXPORT Particle* createAnElectron(double Px, double Py, double Pz, double Vx, double Vy, double Vz);
DLLEXPORT void    updateParticlePosition(Particle* particle, BField* B);
DLLEXPORT double  returnSimTime(BField* B);
DLLEXPORT double  returnParticleMass(Particle* particle);
DLLEXPORT double  returnParticleCharge(Particle* particle);
DLLEXPORT double* returnParticlePosition(Particle* particle);
DLLEXPORT double* returnParticleVelocity(Particle* particle);
DLLEXPORT double* returnTotalBatP(BField* B, double Px, double Py, double Pz);
DLLEXPORT void    increaseTimeBydt(BField* B);
DLLEXPORT void    deleteDblArray(double* arrayPtr);
//DLLEXPORT void    deleteObjects(BObject* B, Particle** particleArray, int numberOfParticles);

#endif