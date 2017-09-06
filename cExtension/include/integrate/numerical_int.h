#ifndef NUMERICAL_INT_H
#define NUMERICAL_INT_H

//pulled from https://rosettacode.org/wiki/Numerical_integration#C.2B.2B on Dec 10, 2016 at 1:00 AM CST
// the integration routine
template<typename Method, typename F, typename Float>
double integrate(F f, Float a, Float b, int steps, Method m, double var[])
{
	double s = 0;
	double h = (b - a) / steps;
	for (unsigned int i = 0; i < steps; ++i)
		s += m(f, a + h*i, h, var);
	return h*s;
}

// methods
class rectangular
{
public:
	enum position_type { left, middle, right };
	rectangular(position_type pos) : position(pos) {}
	template<typename F, typename Float>
	double operator()(F f, Float x, Float h, double var[]) const
	{
		switch (position)
		{
		case left:
			return f(x,var);
		case middle:
			return f(x + h / 2,var);
		case right:
			return f(x + h,var);
		}
	}
private:
	const position_type position;
};

class trapezium
{
public:
	template<typename F, typename Float>
	double operator()(F f, Float x, Float h, double var[]) const
	{
		return (f(x,var) + f(x + h,var)) / 2;
	}
};

class simpson
{
public:
	template<typename F, typename Float>
	double operator()(F f, Float x, Float h, double var[]) const
	{
		return (f(x,var) + 4 * f(x + h / 2,var) + f(x + h,var)) / 6;
	}
};
#endif