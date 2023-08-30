#pragma once
#include "name.h"
#include "bits/stdc++.h"
using namespace std;
class solution
{
	vector<valarray<db>> a;
	db l,r,k;
public:
	solution(const vector<valarray<db>> &_a,const db &_l,const db &_r)
		:a(_a),l(_l),r(_r),k((db)(a.size()-1)/(r-l)) {}
	/*solution(vector<valarray<db>> &&_a,const db &_l,const db &_r)
		:a(_a),l(_l),r(_r),k((db)(a.size()-1)/(r-l)) {}*/
	const valarray<db> &operator()(db x) const
	{
		return a[clamp<int>(round((x-l)*k),0,a.size()-1)];
	}
};