#pragma once
#include "solution.h"
class IVP
{
public:
	virtual solution solve(const function<valarray<db>(const valarray<db> &,db)> &f,db L,db R,const valarray<db> &xl,int N) const=0;
};