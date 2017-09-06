#include "tools/constsandtypes.h"

dblArray3_t operator+(const dblArray3_t& a, const dblArray3_t& b) { return{ a[0] + b[0], a[1] + b[1], a[2] + b[2] }; }
void operator+=(dblArray3_t& a, const dblArray3_t& b) { a = (a + b); }

dblArray3_t operator-(const dblArray3_t& a, const dblArray3_t& b) { return{ a[0] - b[0], a[1] - b[1], a[2] - b[2] }; }
void operator-=(dblArray3_t& a, const dblArray3_t& b) { a = (a - b); }

dblArray3_t operator*(const dblArray3_t& a, int b) { return{ a[0] * b, a[1] * b, a[2] * b }; }
dblArray3_t operator*(const dblArray3_t& a, double b) { return{ a[0] * b, a[1] * b, a[2] * b }; }
void operator*=(dblArray3_t& a, int b) { a = (a * b); }
void operator*=(dblArray3_t& a, double b) { a = (a * b); }

dblArray3_t operator/(const dblArray3_t& a, int b) { return{ a[0] / b, a[1] / b, a[2] / b }; }
dblArray3_t operator/(const dblArray3_t& a, double b) { return{ a[0] / b, a[1] / b, a[2] / b }; }
void operator/=(dblArray3_t& a, int b) { a = (a / b); }
void operator/=(dblArray3_t& a, double b) { a = (a / b); }

std::ostream& operator<<(std::ostream& out, const dblArray3_t& a)
{
	out << "{ " << a[0] << ", " << a[1] << ", " << a[2] << " }";
	return out;
}