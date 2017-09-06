#ifndef BFIELD_WIRECOIL_H
#define BFIELD_WIRECOIL_H
#include "tools/constsandtypes.h"
#include "bfield/BObject_c.h"
#include "particle/Particle_c.h"
#include "bfield/BField_c.h"

class WireCoil : public BObject
{
private:
	std::string name_m;
	int* window_m; //change the type once you know how and have one defined
	int* pic_m; //change the type once you know how and have one defined
	dblArray3_t coilCenter_m; //center point of wire loop
	dblArray3_t coilAxis_m; //normal vector from center of loop perpendicular to loop plane along positive B (use rhr to relate to I flow direction)
	dblArray3_t axisRot_m; //theta, phi to rotate coord system, abbreviated CS: (CS)' - (CS) = axisRot_m
	int N_m; //number of wire wraps
	double I_m; //current flowing through loop
	double R_m; //loop radius
	double d_m; //distance along x axis...define this way?  is needed if each loop defined separately?
	double constant_m; //constant that is easier to calculate once
	bool rightLoop_m; //is this the left or right loop? right loop = 1, left loop = 0 - of course left/right are relative - B flows from left to right loop

	double c1_m{ 0.0 }; //variables for numerical integration
	double c2_m{ 0.0 };
	double c3_m{ 0.0 };
	double a4_m{ 0.0 };
	double c4_m{ 0.0 };
	double c5_m{ 0.0 };
	double c6_m{ 0.0 };
	double c7_m{ 0.0 };

	

public:
	WireCoil(const dblArray3_t& coilCenter, const dblArray3_t& coilAxis, int N, double I, double R, bool rightLoop, std::string name = "", int* window = nullptr, int* pic = nullptr) 
		: coilCenter_m(coilCenter), coilAxis_m(coilAxis), N_m{ N }, I_m{ I }, R_m{ R }, rightLoop_m{ rightLoop }, name_m { name	}, window_m{ window }, pic_m{ pic }, 
		  constant_m{ N * I * SCALE }, d_m{0.0} { Init(); }

	void Init();

	void setIntegConst(const dblArray3_t& P); //set back to private after compile
	static double dBx(double x, void* data, double var[6]); //for gauss_legendre
	static double dBy(double x, void* data, double var[6]);
	static double dBz(double x, void* data, double var[6]);

	double getc1_m() { return c1_m; };
	double getc2_m() { return c2_m; };
	double getc3_m() { return c3_m; };
	double geta4_m() { return a4_m; };
	double getc4_m() { return c4_m; };
	double getc5_m() { return c5_m; };
	double getc6_m() { return c6_m; };
	double getc7_m() { return c7_m; };

	dblArray3_t calcBatP(const dblArray3_t& P, int norder);
};

#endif