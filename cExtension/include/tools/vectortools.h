#ifndef VECTORTOOLS_H
#define VECTORTOOLS_H
#include "tools/constsandtypes.h"

dblArray3_t cartesianToSpherical(const dblArray3_t& cart);
dblArray3_t sphericalToCartesian(const dblArray3_t& sphere, int warnFlag=0);
dblArray3_t rotateCartesianVector(const dblArray3_t& cart, double rotTheta, double rotPhi, int warnFlag=0);
dblArray3_t cross3dVectorsAndMultByConst(const dblArray3_t& a, const dblArray3_t& b, double c = 1.0);
double dA3len(const dblArray3_t& a);
void dA3coutSci(const dblArray3_t& var);

#endif