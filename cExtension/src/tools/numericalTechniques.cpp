#include "bfield/BField_c.h"

//forward decls for these functions are in "BField/BField_c.h"

dblArray3_t fourthOrderRungeKutta_dy(FuncPtr3D_t funcptr, BField* B, double* funcArgs, int argsLen)
{//requires [h, t_0 = 0, y_0x, y_0y, y_0z, ... whatever else you want here]
	
	//check that t_0 is 0.0
	
	double y_0x{ funcArgs[2] };
	double y_0y{ funcArgs[3] };
	double y_0z{ funcArgs[4] };

	dblArray3_t k1{ nulldA3_t };
	dblArray3_t k2{ nulldA3_t };
	dblArray3_t k3{ nulldA3_t };
	dblArray3_t k4{ nulldA3_t };

	k1 = funcptr(B, funcArgs, argsLen);

	funcArgs[1] = funcArgs[0] / 2;
	funcArgs[2] = y_0x + k1[0] * funcArgs[1];
	funcArgs[3] = y_0y + k1[1] * funcArgs[1];
	funcArgs[4] = y_0z + k1[2] * funcArgs[1];
	k2 = funcptr(B, funcArgs, argsLen);

	funcArgs[2] = y_0x + k2[0] * funcArgs[1];
	funcArgs[3] = y_0y + k2[1] * funcArgs[1];
	funcArgs[4] = y_0z + k2[2] * funcArgs[1];
	k3 = funcptr(B, funcArgs, argsLen);

	funcArgs[1] = funcArgs[0];
	funcArgs[2] = y_0x + k3[0] * funcArgs[1];
	funcArgs[3] = y_0y + k3[1] * funcArgs[1];
	funcArgs[4] = y_0z + k3[2] * funcArgs[1];
	k4 = funcptr(B, funcArgs, argsLen);

	return (k1 + (k2 + k3) * 2 + k4) * funcArgs[0] / 6;
}

dblArray3_t Fwrapper(BField* B, double* funcArgs, int argsLen)
{//[dt, t_0, vx, vy, vz, q, m, px_0, py_0, pz_0]
	dblArray3_t Fl{ nulldA3_t };
	dblArray3_t Fm{ nulldA3_t };

	Fl = B->calcLorentzF(funcArgs, argsLen);
	Fm = B->calcMirrorF(funcArgs, argsLen);
	
	return (Fl + Fm) / funcArgs[6]; //converts to acceleration
}