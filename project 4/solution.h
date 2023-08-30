#pragma once
#include "bits/stdc++.h"
using namespace std;
class solution
{
	vector<valarray<db>> y;
	db xl,tl,ik,ih;
public:
	solution(const vector<valarray<db>> &y,db xl,db xr,db tl,db tr):y(y),xl(xl),tl(tl),ik((y.size()-1)/(tr-tl)),ih((y[0].size()-1)/(xr-xl)){}
	db operator()(db x,db t) const
	{
		return y[clamp<int>(round((t-tl)*ik),0,y.size()-1)][clamp<int>(round((x-xl)*ih),0,y[0].size()-1)];
	}
	void debug() const
	{
		for (const auto &v:y)
		{
			for (const db &x:v) cout<<x<<' ';
			cout<<endl;
		}
	}
};