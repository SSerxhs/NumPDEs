#pragma once
#include "bits/stdc++.h"
using namespace std;
using db=double;
vector<vector<db>> tran(const vector<vector<pair<int,db>>> &A)
{
	int n=A.size(),i;
	vector r(n,vector<db>(n));
	for (i=0; i<n; i++) for (auto [y,v]:A[i]) r[i][y]+=v;
	return r;
}
valarray<db> relax(const vector<vector<pair<int,db>>> &A,valarray<db> x0,valarray<db> b,int T)
{
	int n=x0.size(),i;
	if (n==0) throw invalid_argument("Incorrect length of x0.");
	if (b.size()!=n) throw invalid_argument("Incorrect length of b.");
	if (A.size()!=n) throw invalid_argument("Incorrect size of A.");
	valarray<db> D(n);
	for (i=0; i<n; i++)
	{
		for (auto [y,v]:A[i])
		{
			if (y<0||y>=n) throw invalid_argument("Incorrect coordinate of A.");
			if (i==y) D[i]+=v;
		}
	}
	const db w=2./3;
	D=w/D;
	b*=D;
	while (T--)
	{
		valarray<db> r(n);
		for (i=0; i<n; i++) for (auto [y,v]:A[i]) if (i!=y) r[i]-=v*x0[y];
		x0=r*D+(1-w)*x0+b;
	}
	return x0;
}
valarray<db> multiply(const vector<vector<pair<int,db>>> &A,const valarray<db> &x)
{
	int n=A.size(),i;
	if (x.size()!=n) throw invalid_argument("Incorrect length of x.");
	valarray<db> r(n);
	for (i=0; i<n; i++) for (auto [y,v]:A[i])
	{
		if (y<0||y>=n) throw invalid_argument("Incorrect coordinate of A.");
		r[i]+=x[y]*v;
	}
	return r;
}
class restriction
{
public:
	virtual valarray<db> operator() (const valarray<db> &a) const=0;
};
class interpolation
{
public:
	virtual valarray<db> operator() (const valarray<db> &a) const=0;
};
class boundary_value
{
protected:
	db v;
public:
	boundary_value(db x) :v(x) {}
	virtual db value() const=0;
};
class function_value :public boundary_value
{
public:
	function_value(db x) :boundary_value(x) {}
	db value() const { return v; }
};
class derivative_value :public boundary_value
{
public:
	derivative_value(db x) :boundary_value(x) {}
	db value() const { return v; }
};