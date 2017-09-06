#include "api.h"

///initialization functions
DLLEXPORT BField* initSimulation(double dt, double WC1x, double WC1y, double WC1z, double WC2x, double WC2y, double WC2z, int N, double I, double R)
{
	//create B Field object
	BField* B = new BField(dt);

	//creating Wire Coil objects
	dblArray3_t WC1cent{ WC1x, WC1y, WC1z };
	dblArray3_t WC2cent{ WC2x, WC2y, WC2z };
	dblArray3_t WCaxis{ 1.0, 0.0, 0.0 };
	WireCoil* WC1 = new WireCoil(WC1cent, WCaxis, N, I, R, 0, "Left");
	WireCoil* WC2 = new WireCoil(WC2cent, WCaxis, N, I, R, 1, "Right");

	B->addBObj(WC1);
	B->addBObj(WC2);

	return B;
}

DLLEXPORT Particle* createAnElectron(double Px, double Py, double Pz, double Vx, double Vy, double Vz)
{//for right now, I'm only doing electrons
	dblArray3_t e1posinit{ Px, Py, Pz };
	dblArray3_t e1velinit{ Vx, Vy, Vz };
	Particle* Electron = new Particle(-1.60217657e-19, 9.10938356e-31, e1posinit, e1velinit);

	return Electron;
}

///functions to run in loop
DLLEXPORT void updateParticlePosition(Particle* particle, BField* B)
{
	B->updateParticleP_V(particle, 37);
}

///access functions
DLLEXPORT double returnSimTime(BField* B)
{
	return B->getTime();
}

DLLEXPORT double returnParticleMass(Particle* particle)
{
	return particle->getM();
}

DLLEXPORT double returnParticleCharge(Particle* particle)
{
	return particle->getq();
}

DLLEXPORT double* returnParticlePosition(Particle* particle)
{
	dblArray3_t tmp = particle->getP();

	double* ret = new double[3];
	ret[0] = tmp[0];
	ret[1] = tmp[1];
	ret[2] = tmp[2];

	return ret;
}

DLLEXPORT double* returnParticleVelocity(Particle* particle)
{
	dblArray3_t tmp = particle->getV();

	double* ret = new double[3];
	ret[0] = tmp[0];
	ret[1] = tmp[1];
	ret[2] = tmp[2];
	
	return ret;
}

DLLEXPORT double* returnTotalBatP(BField* B, double Px, double Py, double Pz)
{
	dblArray3_t P{ nulldA3_t };
	dblArray3_t tmp{ nulldA3_t };

	double* ret = new double[3];

	P = { Px, Py, Pz };
	tmp = B->totalBatP(P);

	ret[0] = tmp[0];
	ret[1] = tmp[1];
	ret[2] = tmp[2];

	return ret;
}

DLLEXPORT void increaseTimeBydt(BField* B)
{
	B->forceIncTime();
}

///delete allocated memory
DLLEXPORT void deleteDblArray(double* arrayPtr)
{
	delete[] arrayPtr;
}

/*DLLEXPORT void deleteObjects(BObject* B, Particle** particleArray, int numberOfParticles)
{
	for (int iii = 0; iii < numberOfParticles; iii++)
	{
		delete particleArray[iii];
	}
	delete[] particleArray;
	delete B;
}*/