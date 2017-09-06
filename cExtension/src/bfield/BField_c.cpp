#include "bfield/BField_c.h"
#include "tools/vectortools.h"
#include "particle/Particle_c.h"
#include <iostream>

dblArray3_t BField::totalBatP(const dblArray3_t& P, int norder, bool WarnFlag)
{
	dblArray3_t B{ 0.0, 0.0, 0.0 };
	
	if (!bObjPtrList_m[0])
	{
		if (WarnFlag) { std::cout << "Warning (totalBatP): No B Objects or charged particles added to this BField instance.  Returning B = {0,0,0}." << std::endl; }
		return nulldA3_t;
	}

	for (unsigned int iii = 0; iii < bObjPtrList_m.size(); iii++)
	{
		dblArray3_t c{nulldA3_t};
		if (!(bObjPtrList_m[iii]))
		{
			if (WarnFlag) { std::cout << "Warning (totalBatP): nullptr in bObjPtrList_m at index " << std::to_string(iii) << ".  Skipping this index." << std::endl; }
			continue;
		}
		c = bObjPtrList_m[iii]->calcBatP(P,norder);

		B += c;
	}

	incTime();

	return B;
}

void BField::addBObj(BObject* bobjptr)
{
	if (!bObjPtrList_m[0] && (bObjPtrList_m.size() == 1))
	{
		bObjPtrList_m.pop_back();
	}
	bObjPtrList_m.push_back(bobjptr);
	std::cout << "Added BObj.  Size: " << std::to_string(bObjPtrList_m.size()) << std::endl;
}

dblArray3_t BField::totaldVatP(Particle* partobj, double dt, int norder)
{//[dt, t_0, vx, vy, vz, q, m, px_0, py_0, pz_0]
	dblArray3_t v{ nulldA3_t };
	dblArray3_t p{ nulldA3_t };
	v = partobj->getV();
	p = partobj->getP();

	double* funcArgs = new double[10];
	funcArgs[0] = dt;
	funcArgs[1] = 0.0;
	funcArgs[2] = v[0];
	funcArgs[3] = v[1];
	funcArgs[4] = v[2];
	funcArgs[5] = partobj->getq();
	funcArgs[6] = partobj->getM();
	funcArgs[7] = p[0];
	funcArgs[8] = p[1];
	funcArgs[9] = p[2];

	dblArray3_t res;
	res = fourthOrderRungeKutta_dy(&Fwrapper, this, funcArgs, 10);

	delete[] funcArgs;
	return res;
}

void BField::updateParticleP_V(Particle* partobj, int norder)
{
	//Need to add updating particle P_V of all particles in a list, then increasing time - at some point
	partobj->updP(totaldVatP(partobj, dt_m, norder), dt_m, true);
}

dblArray3_t BField::calcMirrorF(double* funcArgs, int argsLen)
{//[dt, t_0, vx, vy, vz, q, m, px_0, py_0, pz_0]
	dblArray3_t vtmp{ funcArgs[2], funcArgs[3], funcArgs[4] };
	dblArray3_t ptmp{ funcArgs[7], funcArgs[8], funcArgs[9] };
	ptmp = ptmp + (vtmp * funcArgs[1]);

	dblArray3_t Bp = totalBatP(ptmp);
	double Bplen = dA3len(Bp);
	double BpdotV = Bp[0] * vtmp[0] + Bp[1] * vtmp[1] + Bp[2] * vtmp[2];

	double scaledLength = 1e-8; //total length of vector in the direction of B that is used to calculate B(p + ds) and B(p - ds)
	dblArray3_t halfds = Bp * (scaledLength / Bplen);

	dblArray3_t Bpminus = totalBatP(ptmp - halfds);
	dblArray3_t Bpplus = totalBatP(ptmp + halfds);

	double vperp2 = pow(vtmp[0], 2) + pow(vtmp[1], 2) + pow(vtmp[2], 2) - pow(BpdotV / Bplen, 2);
	double mu = funcArgs[6] * vperp2 / (2 * Bplen);
	dblArray3_t FgradB = (Bpplus - Bpminus) * -mu;

	for (unsigned int jjj = 0; jjj < 3; jjj++)
	{
		if (abs(Bp[jjj]) < 1e-20)
		{
			FgradB[jjj] = 0;
			continue;
		}
		FgradB[jjj] /= (2 * scaledLength);
	}

	return FgradB;
}

dblArray3_t BField::calcLorentzF(double* funcArgs, int argsLen)
{
	dblArray3_t vtmp{ funcArgs[2], funcArgs[3], funcArgs[4] };
	dblArray3_t ptmp{ funcArgs[7], funcArgs[8], funcArgs[9] };
	ptmp = ptmp + (vtmp * funcArgs[1]);

	dblArray3_t B{ nulldA3_t };
	B = totalBatP(ptmp);
	dblArray3_t vxB{ nulldA3_t };
	vxB[0] = vtmp[1] * B[2] - vtmp[2] * B[1];
	vxB[1] = vtmp[2] * B[0] - vtmp[0] * B[2];
	vxB[2] = vtmp[0] * B[1] - vtmp[1] * B[0];

	return vxB * funcArgs[5];
}