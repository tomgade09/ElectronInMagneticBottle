#ifndef CONSTSANDTYPES_H
#define CONSTSANDTYPES_H
#include <array>

#define SCALE 1e-5
#define M_PI       3.14159265358979323846   // pi
#define M_PI_2     1.57079632679489661923   // pi/2
#define M_PI_4     0.785398163397448309616  // pi/4 ^^ from <math.h>
using dblArray3_t = std::array<double, 3>;
#define nulldA3_t { 0.0, 0.0, 0.0 }

dblArray3_t operator+(const dblArray3_t& a, const dblArray3_t& b);
void operator+=(dblArray3_t& a, const dblArray3_t& b);

dblArray3_t operator-(const dblArray3_t& a, const dblArray3_t& b);
void operator-=(dblArray3_t& a, const dblArray3_t& b);

dblArray3_t operator*(const dblArray3_t& a, int b);
dblArray3_t operator*(const dblArray3_t& a, double b);
void operator*=(dblArray3_t& a, int b);
void operator*=(dblArray3_t& a, double b);

dblArray3_t operator/(const dblArray3_t& a, int b);
dblArray3_t operator/(const dblArray3_t& a, double b);
void operator/=(dblArray3_t& a, int b);
void operator/=(dblArray3_t& a, double b);

std::ostream& operator<<(std::ostream& out, const dblArray3_t& a);

#endif