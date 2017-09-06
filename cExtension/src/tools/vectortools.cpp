#include <cmath>
#include <iostream>
#include "tools/vectortools.h"

//Physics (ISO) convention for spherical used: phi = angle from x axis, theta = angle from z axis
//Range: 0 <= phi < 2*PI; 0 <= theta < PI
//All angles in radians unless otherwise noted
//From here out, the word "vector" refers to a 1 x 3 array (represented by std::array<double, 3>), NOT std::vector
//        | z
//        |
//        |   @ (rho, theta, phi)
//        |__/-
//  theta>| / -      << 0 to Pi
//        |/  -
//--------------------- y
//       / -  -
//  phi>/___- -      << 0 to 2*Pi
//     /     --
//    /       -
//   /x

dblArray3_t cartesianToSpherical(const dblArray3_t& cart)
{ //convert a cartesian array {x, y, z} to a spherical array {rho, theta, phi}
	if ((cart[0] == 0.0) && (cart[1] == 0.0) && (cart[2] > 0.0)) //pos z-axis vector case
		return{ cart[2], 0.0, 0.0 };
	else if ((cart[0] == 0.0) && (cart[1] == 0.0) && (cart[2] < 0.0)) //neg z-axis vector case
		return{ -cart[2], M_PI, 0.0 };
	else if ((cart[0] == 0.0) && (cart[1] == 0.0) && (cart[2] == 0.0)) //zero vector case
		return nulldA3_t;

	double rho = sqrt(pow(cart[0], 2.0) + pow(cart[1], 2.0) + pow(cart[2], 2.0));
	double theta = acos((cart[2] / rho));
	double phi = atan2(cart[1], cart[0]);

	if (phi < 0.0)
		phi += (2 * M_PI);

	return{ rho, theta, phi };
}

dblArray3_t sphericalToCartesian(const dblArray3_t& sphere, int warnFlag) //warnFlag default = 0 (in fwd declaration)
{ //convert a spherical vector {rho, theta, phi} to a cartesian vector {x, y, z}
	if (sphere[1] == 0.0 && sphere[2] == 0.0) //z-axis case
		return{ 0.0, 0.0, sphere[0] };
	else if (sphere[2] == 0.0 && (abs(1 - (sphere[1] / M_PI)) < 1e-14)) //neg z-axis case, case accounts for float rounding error
		return{ 0.0, 0.0, -sphere[0] };
	else if (sphere[0] < 0.0)
	{
		if (warnFlag == 1) { std::cout << "WARNING: Negative rho passed to sphericalToCartesian.  You can't do that!  Returning 0 vector." << std::endl; }
		return nulldA3_t; //could also do a Pi/2 +/- Phi (depending on quadrant) and Pi + Theta...but shouldn't be passing in a negative rho
	}

	double x = (sphere[0] * sin(sphere[1]) * cos(sphere[2]));
	double y = (sphere[0] * sin(sphere[1]) * sin(sphere[2]));
	double z = (sphere[0] * cos(sphere[1]));

	return{ x, y, z };
}

dblArray3_t rotateCartesianVector(const dblArray3_t& cart, double rotTheta, double rotPhi, int warnFlag) //warnFlag default = 0 (in fwd declaration)
{ //rotate a cartesian vector cart by Theta and Phi in spherical, return 3d cartesian vector
	if ((rotTheta == 0.0) && (rotPhi == 0.0)) //no rotation
		return cart;

	dblArray3_t sphere = cartesianToSpherical(cart);
	sphere[1] += rotTheta;
	sphere[2] += rotPhi;

	if (sphere[1] < 0.0) //negative theta value - flip to positive, rotate phi PI radians
	{
		if (warnFlag == 1) { std::cout << "Warning: Negative theta.  Converting to equivalent vector with 0 <= theta < PI." << std::endl; }
		sphere[1] = -sphere[1];
		sphere[2] += M_PI;
	}
	else if (sphere[1] > M_PI) //theta greater than PI - subtract from 2 * PI to get offset from PI radians, rotate phi PI radians
	{
		if (warnFlag == 1) { std::cout << "Warning: theta > PI.  Converting to equivalent vector with 0 <= theta < PI." << std::endl; }
		sphere[1] = (2 * M_PI - sphere[1]);
		sphere[2] += M_PI;
	}

	if (sphere[2] < 0.0) //negative phi value - add 2*PI to put within 0 <= phi < 2*PI
	{
		if (warnFlag == 1) { std::cout << "Warning: Negative phi.  Converting to equivalent vector with 0 <= phi < 2 * PI." << std::endl; }
		sphere[2] += (2 * M_PI);
	}
	else if (sphere[2] > (2 * M_PI))
	{
		if (warnFlag == 1) { std::cout << "Warning: phi > 2 * PI.  Converting to equivalent vector with 0 <= phi < 2 * PI." << std::endl; }
		sphere[2] -= (2 * M_PI);
	}

	return sphericalToCartesian(sphere);
} 

dblArray3_t cross3dVectorsAndMultByConst(const dblArray3_t& a, const dblArray3_t& b, double c)
{ //return 3D vector v = a(3D vector) x b(3D vector) * c (double)
	return{ (a[1] * b[2] - a[2] * b[1]) * c, 
		    (a[2] * b[0] - a[0] * b[2]) * c,
		    (a[0] * b[1] - a[1] * b[0]) * c };
}

double dA3len(const dblArray3_t& a)
	{ return sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2)); }

void dA3coutSci(const dblArray3_t& var)
	{ std::cout << "{ " << std::scientific << var[0] << ", " << var[1] << ", " << var[2] << " }" << std::endl; }