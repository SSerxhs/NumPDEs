#pragma once
#include "bits/stdc++.h"
#include "matrix.h"
using namespace std;
using db=double;
class full_weighting :public restriction
{
public:
	valarray<db> operator()(const valarray<db> &a) const
	{
		int n=a.size(),i;
		if (n<5||n!=(1<<__lg(n)|1)) throw invalid_argument("Incorrect length.");
		n=n/2+1;
		valarray<db> r(n);
		r[0]=a[0]; r[n-1]=a[n*2-2];
		for (i=1; i+1<n; i++) r[i]=(a[i*2-1]+a[i*2]*2+a[i*2+1])/4;
		return r;
	}
};
class injection :public restriction
{
public:
	valarray<db> operator()(const valarray<db> &a) const
	{
		int n=a.size(),i;
		if (n<5||n!=(1<<__lg(n)|1)) throw invalid_argument("Incorrect length.");
		n=n/2+1;
		valarray<db> r(n);
		for (i=0; i<n; i++) r[i]=a[i*2];
		return r;
	}
};
class linear :public interpolation
{
public:
	valarray<db> operator()(const valarray<db> &a) const
	{
		int n=a.size(),i;
		if (n<3||n!=(1<<__lg(n)|1)) throw invalid_argument("Incorrect length.");
		valarray<db> r(n*2-1);
		for (i=0; i<n; i++) r[i*2]=a[i];
		for (i=1; i<n; i++) r[i*2-1]=(a[i]+a[i-1])/2;
		return r;
	}
};
class quadratic :public interpolation//?
{
public:
	valarray<db> operator()(const valarray<db> &a) const
	{
		int n=a.size(),i;
		if (n<3||n!=(1<<__lg(n)|1)) throw invalid_argument("Incorrect length.");
		valarray<db> r(n*2-1);
		for (i=0; i<n; i++) r[i*2]=a[i];
		for (i=1; i*2<n; i++) r[i*2-1]=(3*a[i-1]+6*a[i]-a[i+1])/8;
		for (; i<n; i++) r[i*2-1]=(-a[i-2]+6*a[i-1]+3*a[i])/8;
		return r;
	}
};
function<db(db)> ZERO_FUNC=[&](db x) ->db { return 0; };
vector<vector<pair<int,db>>> get_A(int n,const boundary_value &vl,const boundary_value &vr)
{
	int i;
	vector<vector<pair<int,db>>> A(n);
	db h=1.0/(n-1);
	const function_value &tmp=function_value(0);
	if (typeid(vl)==typeid(tmp)) A[0]={{0,1}};
	else A[0]={{0,-1.5/h},{1,2/h},{2,-0.5/h}};
	for (i=1; i<n-1; i++) A[i]={{i-1,1/h/h},{i,-2/h/h},{i+1,1/h/h}};
	if (typeid(vr)==typeid(tmp)) A[n-1]={{n-1,1}};
	else A[n-1]={{n-3,0.5/h},{n-2,-2/h},{n-1,1.5/h}};
	return A;
}
class cycles
{
public:
	virtual void operator()(valarray<db> &v,valarray<db> f,const boundary_value &vl,const boundary_value &vr,const restriction &rs,const interpolation &inte,int v1=2,int v2=2) const=0;
};
class V_cycle :public cycles
{
public:
	// template<typename func> V_cycle(const func &rhs) :cycles(rhs) {}
	void operator()(valarray<db> &v,valarray<db> f,const boundary_value &vl,const boundary_value &vr,const restriction &rs,const interpolation &inte,int v1=2,int v2=2) const
	{
		int n=v.size(),i;
		if (n<3||n!=(1<<__lg(n)|1)) throw invalid_argument("Incorrect length of v.");
		if (n!=f.size()) throw invalid_argument("Incorrect length of f.");
		const auto &A=get_A(n,vl,vr);
		v=relax(A,v,f,v1);
		if (n>3)
		{
			valarray f2=rs(f-multiply(A,v));
			valarray<db> vv(n/2+1);
			V_cycle()(vv,f2,vl,vr,rs,inte,v1,v2);
			v=v+inte(vv);
		}
		v=relax(A,v,f,v2);
	}
};
class FMG :public cycles
{
public:
	void operator()(valarray<db> &v,valarray<db> f,const boundary_value &vl,const boundary_value &vr,const restriction &rs,const interpolation &inte,int v1=2,int v2=2) const
	{
		int n=v.size(),i;
		if (n<3||n!=(1<<__lg(n)|1)) throw invalid_argument("Incorrect length of v.");
		if (n!=f.size()) throw invalid_argument("Incorrect length of f.");
		if (n==3) v=0;
		else
		{
			valarray f2=rs(f);
			v.resize(n/2+1);
			FMG()(v,f2,vl,vr,rs,inte,v1,v2);
			v=inte(v);
		}
		V_cycle()(v,f,vl,vr,rs,inte,v1,v2);
	}
};

class solution
{
	valarray<db> v;
	int times;
	db residual;
public:
	template<typename func1,typename func2,typename func3>
	solution(const func1 &_f,const func2 &boundary,int n,int T,db eps,const restriction &rs,const interpolation &inte,const cycles &cyc,const func3 &guess,int v1=2,int v2=2) :v(n)
	{
		if (n<3||n!=(1<<__lg(n)|1)) throw invalid_argument("Incorrect n.");
		int i;
		for (i=0; i<n; i++) v[i]=guess(i/(n-1.));
		valarray<db> f(n);
		function_value tmp(0);
		unique_ptr<boundary_value> vl,vr;
		vl=boundary(0); vr=boundary(1);
		if (typeid(*vr)!=typeid(tmp)&&typeid(*vl)!=typeid(tmp)) vl=unique_ptr<boundary_value>(&tmp);
		f[0]=vl->value();
		f[n-1]=vr->value();
		for (i=1; i+1<n; i++) f[i]=_f(i/(n-1.));
		const auto &A=get_A(n,*vl,*vr);
		f-=multiply(A,v);
		for (times=1; (times<=T||abs(f).sum()>=eps)&&times<=100; times++)
		{
			valarray<db> r(n);
			cyc(r,f,*vl,*vr,rs,inte,v1,v2);
			v+=r;
			f-=multiply(A,r);
		}
		residual=abs(f).sum();
	}
	template<typename func1,typename func2>
	solution(const func1 &_f,const func2 &boundary,int n,int T,db eps,const restriction &rs,const interpolation &inte,const cycles &cyc,int v1=2,int v2=2) :solution(_f,boundary,n,T,eps,rs,inte,cyc,ZERO_FUNC,v1,v2) {}
	db operator()(db x) const { return v[(int)round(x*(v.size()-1))]; }
	int get_times() const { return times; }
	db get_residual() const { return residual; }
};