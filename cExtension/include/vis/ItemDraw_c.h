#ifndef ITEMDRAW_C_H
#define ITEMDRAW_C_H
#include <list>
#include <iostream>
#include "tools/constsandtypes.h"

class Pic
{
public:
	virtual ~Pic() {}
	virtual void draw(const dblArray3_t& P) = 0;
};

class ParticlePic : public Pic
{
private:
	double radius_m;
	std::list<dblArray3_t> trail_m;
	dblArray3_t color_m;
	int retain_m;
	int interval_m;
	bool maketrail_m;

	void setPos(const dblArray3_t& P);
	void incrementInd();

public:
	ParticlePic(int interval = 1, int retain = 0, double radius = 0.05, dblArray3_t color = { 0.0, 1.0, 0.0 }) :
		radius_m{ radius }, color_m( color ), retain_m{ retain }, interval_m{ interval }, maketrail_m{ (retain > 1) }
	{ 
		if (retain > 1)
		{
			for (int iii = 0; iii < retain; iii++)
				trail_m.push_back(nulldA3_t);
			std::cout << "trail_m.size() " << trail_m.size() << std::endl;
		}
		else
			trail_m.push_back(nulldA3_t);
		if (interval < 1)
		{
			std::cout << "Can't set interval to less than 1.  Setting to 1." << std::endl;
			interval_m = 1;
		}
	}

	void draw(const dblArray3_t& P);
};

class WireCoilPic : public Pic
{
private:
	dblArray3_t cent_m;
	dblArray3_t axis_m;
	double radius_m;
	double thickness_m;

public:
	WireCoilPic(dblArray3_t center, dblArray3_t axis, double radius, double thickness = 0.01) :
		cent_m( center ), axis_m( axis ), radius_m{ radius }, thickness_m{ thickness } {}
	
	void draw(const dblArray3_t& dummy = nulldA3_t);
};

#endif