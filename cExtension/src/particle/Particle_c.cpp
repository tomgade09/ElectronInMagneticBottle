#include "particle/Particle_c.h"
#include "tools/vectortools.h"
#include <cmath>

//private functions
void Particle::updV(const dblArray3_t& B, double dt)
{
	velocity_m += cross3dVectorsAndMultByConst(velocity_m, B, qOVERm_m * dt);
}

void Particle::updV(const dblArray3_t& dv)
{ velocity_m += dv; }

//public functions
void Particle::updP(const dblArray3_t& B, double dt)
{
	updV(B, dt);
	position_m += (velocity_m * dt);
}

void Particle::updP(const dblArray3_t& dV, double dt, bool newmethod)
{
	updV(dV);
	position_m += (velocity_m * dt);
}

dblArray3_t Particle::calcBatP(const dblArray3_t& P, int norder)
{
	dblArray3_t r(P - position_m);
	bool boo = 1;
	for (unsigned int iii = 0; iii < 3; iii++)
		boo &= (r[iii] < 1e-14);
	if (boo)
		return{ 0.0, 0.0, 0.0 }; //no mag field "inside/at" point charge
	return cross3dVectorsAndMultByConst(velocity_m, r, (SCALE * charge_m / pow(dA3len(r), 3)));
}