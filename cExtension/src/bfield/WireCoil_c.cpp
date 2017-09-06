#include "bfield/WireCoil_c.h"
#include "integrate/gauss_legendre.h"
#include "tools/vectortools.h"
#include <cmath>
#include <iostream>

//private functions

void WireCoil::Init()
{
	dblArray3_t a{ { 1.0, M_PI_2, 0.0 } };
	dblArray3_t b(cartesianToSpherical(coilAxis_m));
	axisRot_m[0] = 1.0;

	for (unsigned int iii = 1; iii < 3; iii++)
		axisRot_m[iii] = (a[iii] - b[iii]);
	//somehow account for loop pair centers not being at {0,0,0}
	//not sure if I should allow a dblArray3_t to be passed, and let the user worry about figuring out the right value
	//or perhaps passing in a pointer to the paired loop and a T/F to allow the user to pair or not with another loop
	//but for now, just going to leave it and ensure symmetry in the examples
	d_m = dA3len(coilCenter_m); //coilCenter represents a point at the center of the loop, but here, is also the vector from origin to the center (cx-0, cy-0, cz-0)
	if (!rightLoop_m)			//need to be careful here...BatP point P simply adds/subtracts d to Px...need to translate, then rotate before calcing B
		d_m *= -1;
}

void WireCoil::setIntegConst(const dblArray3_t& P)
{
	c1_m = -constant_m * R_m * P[2];
	c2_m = -constant_m * R_m * P[1];
	c3_m =  constant_m * pow(R_m,2);
	a4_m =  pow(P[0],2) + pow(P[1],2) + pow(P[2],2);
	c4_m =  a4_m + 2 * P[0] * d_m + pow(d_m, 2) + pow(R_m, 2);
	c5_m = -2 * P[1] * R_m;
	c6_m = -2 * P[2] * R_m;
	c7_m =  constant_m * R_m * (P[0] + d_m);
}

//public functions

double WireCoil::dBx(double x, void* data, double var[6]) //for gauss_legendre
{
	//return ((c1_m * sin(x) + c2_m * cos(x) + c3_m) / pow((c4_m + c5_m * cos(x) + c6_m * sin(x)), 1.5));
	return ((var[0] * sin(x) + var[1] * cos(x) + var[2]) / pow((var[3] + var[4] * cos(x) + var[5] * sin(x)), 1.5));
}

double WireCoil::dBy(double x, void* data, double var[6]) //for gauss_legendre
{
	//return (c7_m * cos(x) / pow((c4_m + c5_m * cos(x) + c6_m * sin(x)), 1.5));
	return (var[0] * cos(x) / pow((var[1] + var[2] * cos(x) + var[3] * sin(x)), 1.5));
}

double WireCoil::dBz(double x, void* data, double var[6]) //for gauss_legendre
{
	//return (c7_m * sin(x) / pow((c4_m + c5_m * cos(x) + c6_m * sin(x)), 1.5));
	return (var[0] * sin(x) / pow((var[1] + var[2] * cos(x) + var[3] * sin(x)), 1.5));
}

dblArray3_t WireCoil::calcBatP(const dblArray3_t& P, int norder)
{
	setIntegConst(P);

	dblArray3_t B{ 0.0, 0.0, 0.0 };

	typedef double(*SFP)(double, void*, double[]); //for gauss_legendre
	SFP FP[3] = { dBx, dBy, dBz };
	double dBxCst[6]{ c1_m, c2_m, c3_m, c4_m, c5_m, c6_m };
	double dByCst[6]{ c7_m, c4_m, c5_m, c6_m, 0.0, 0.0 };
	double *dBCst[] { dBxCst, dByCst, dByCst };
	
	for (unsigned int jjj = 0; jjj < 3; jjj++)
		B[jjj] = gauss_legendre(norder, FP[jjj], NULL, 0, 2 * M_PI, dBCst[jjj]); //for gauss_legendre, n=5 to match default SciPy quadpack value
	//32-37 yields some pretty close values...around this range, the numbers don't change up to 14 digits and match python almost perfectly
	//around the 40s the numbers start fluctuating again

	return B;
}