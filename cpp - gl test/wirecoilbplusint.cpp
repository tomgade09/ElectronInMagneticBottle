#include <array>
#include <cmath>
#include "gauss_legendre.h"
#define M_PI       3.14159265358979323846   // pi
using dblArray3_t = std::array<double, 3>;

#if defined(WIN32)
#define SHARED_API extern "C" __declspec(dllexport)
#else
#define SHARED_API
#endif

double dBx(double x, void* data, double var[6]) //for gauss_legendre
{
	//return ((c1_m * sin(x) + c2_m * cos(x) + c3_m) / pow((c4_m + c5_m * cos(x) + c6_m * sin(x)), 1.5));
	return ((var[0] * sin(x) + var[1] * cos(x) + var[2]) / pow((var[3] + var[4] * cos(x) + var[5] * sin(x)), 1.5));
}

double dBy(double x, void* data, double var[6]) //for gauss_legendre
{
	//return (c7_m * cos(x) / pow((c4_m + c5_m * cos(x) + c6_m * sin(x)), 1.5));
	return (var[0] * cos(x) / pow((var[1] + var[2] * cos(x) + var[3] * sin(x)), 1.5));
}

double dBz(double x, void* data, double var[6]) //for gauss_legendre
{
	//return (c7_m * sin(x) / pow((c4_m + c5_m * cos(x) + c6_m * sin(x)), 1.5));
	return (var[0] * sin(x) / pow((var[1] + var[2] * cos(x) + var[3] * sin(x)), 1.5));
}
/*
dblArray3_t calcBatP(const double Px, const double Py, const double Pz, const int N_m, const double I_m, const double R_m, const double d_m)
{

	dblArray3_t B{ 0.0, 0.0, 0.0 };
	dblArray3_t P{ Px, Py, Pz };

	double constant_m = N_m * I_m * 1e-5;

	double c1_m = { -constant_m * R_m * P[2] };
	double c2_m = { -constant_m * R_m * P[1] };
	double c3_m = { constant_m * pow(R_m,2) };
	double a4_m = { pow(P[0],2) + pow(P[1],2) + pow(P[2],2) };
	double c4_m = { a4_m + 2 * P[0] * d_m + pow(d_m, 2) + pow(R_m, 2) };
	double c5_m = { -2 * P[1] * R_m };
	double c6_m = { -2 * P[2] * R_m };
	double c7_m = { constant_m * R_m * (P[0] + d_m) };

	typedef double(*SFP)(double, void*, double[]); //for gauss_legendre
	
	SFP FP[3] = { dBx, dBy, dBz };
	double dBxCst[6]{ c1_m, c2_m, c3_m, c4_m, c5_m, c6_m };
	double dByCst[6]{ c7_m, c4_m, c5_m, c6_m, 0.0, 0.0 };
	double *dBCst[]{ dBxCst, dByCst, dByCst };

	for (int jjj = 0; jjj < 3; jjj++)
	{
		for (int iii = 2; iii <= 128; iii++)
		{
			B[jjj] = gauss_legendre(iii, FP[jjj], NULL, 0, 2 * M_PI, dBCst[jjj]); //for gauss_legendre
		}
	}

	return B;
}
*/
SHARED_API double calcBxatP(const double Px, const double Py, const double Pz, const int N_m, const double I_m, const double R_m, const double d_m)
{

	double B{ 0.0 };
	dblArray3_t P{ Px, Py, Pz };

	double constant_m = N_m * I_m * 1e-5;

	double c1_m = { -constant_m * R_m * P[2] };
	double c2_m = { -constant_m * R_m * P[1] };
	double c3_m = { constant_m * pow(R_m,2) };
	double a4_m = { pow(P[0],2) + pow(P[1],2) + pow(P[2],2) };
	double c4_m = { a4_m + 2 * P[0] * d_m + pow(d_m, 2) + pow(R_m, 2) };
	double c5_m = { -2 * P[1] * R_m };
	double c6_m = { -2 * P[2] * R_m };
	double c7_m = { constant_m * R_m * (P[0] + d_m) };

	double dBxCst[6]{ c1_m, c2_m, c3_m, c4_m, c5_m, c6_m };

	for (int iii = 2; iii <= 128; iii++)
	{
		B = gauss_legendre(iii, dBx, NULL, 0, 2 * M_PI, dBxCst); //for gauss_legendre
	}

	return B;
}

SHARED_API double calcByatP(const double Px, const double Py, const double Pz, const int N_m, const double I_m, const double R_m, const double d_m)
{

	double B{ 0.0 };
	dblArray3_t P{ Px, Py, Pz };

	double constant_m = N_m * I_m * 1e-5;

	double c1_m = { -constant_m * R_m * P[2] };
	double c2_m = { -constant_m * R_m * P[1] };
	double c3_m = { constant_m * pow(R_m,2) };
	double a4_m = { pow(P[0],2) + pow(P[1],2) + pow(P[2],2) };
	double c4_m = { a4_m + 2 * P[0] * d_m + pow(d_m, 2) + pow(R_m, 2) };
	double c5_m = { -2 * P[1] * R_m };
	double c6_m = { -2 * P[2] * R_m };
	double c7_m = { constant_m * R_m * (P[0] + d_m) };

	double dByCst[6]{ c7_m, c4_m, c5_m, c6_m, 0.0, 0.0 };

	for (int iii = 2; iii <= 128; iii++)
	{
		B = gauss_legendre(iii, dBy, NULL, 0, 2 * M_PI, dByCst); //for gauss_legendre
	}

	return B;
}

SHARED_API double calcBzatP(const double Px, const double Py, const double Pz, const int N_m, const double I_m, const double R_m, const double d_m)
{

	double B{ 0.0 };
	dblArray3_t P{ Px, Py, Pz };

	double constant_m = N_m * I_m * 1e-5;

	double c1_m = { -constant_m * R_m * P[2] };
	double c2_m = { -constant_m * R_m * P[1] };
	double c3_m = { constant_m * pow(R_m,2) };
	double a4_m = { pow(P[0],2) + pow(P[1],2) + pow(P[2],2) };
	double c4_m = { a4_m + 2 * P[0] * d_m + pow(d_m, 2) + pow(R_m, 2) };
	double c5_m = { -2 * P[1] * R_m };
	double c6_m = { -2 * P[2] * R_m };
	double c7_m = { constant_m * R_m * (P[0] + d_m) };

	double dByCst[6]{ c7_m, c4_m, c5_m, c6_m, 0.0, 0.0 };

	for (int iii = 2; iii <= 128; iii++)
	{
		B = gauss_legendre(iii, dBz, NULL, 0, 2 * M_PI, dByCst); //for gauss_legendre
	}

	return B;
}