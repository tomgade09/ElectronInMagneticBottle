#ifndef BFIELD_BOBJECT_H
#define BFIELD_BOBJECT_H
#include "tools/constsandtypes.h"

class BObject
{
public:
	virtual ~BObject() {}
	virtual dblArray3_t calcBatP(const dblArray3_t& P, int norder) = 0;
};

#endif