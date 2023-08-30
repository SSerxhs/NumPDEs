#pragma once
#include "bits/stdc++.h"
#include "matrix.h"
#include "debug.h"
using namespace std;
using db=double;
class full_weighting2 :public restriction
{
public:
	valarray<db> operator()(const valarray<db> &a) const
	{
		int m=a.size(),n=round(sqrt(m)),i,j;
		if (n*n!=m||n<5||n!=(1<<__lg(n)|1)) throw invalid_argument("Incorrect length.");
		m=n; n=n/2+1;
		valarray<db> r(n*n);
		for (i=0; i<n; i++)
		{
			r[i]=a[i*2];
			r[(n-1)*n+i]=a[(m-1)*m+i*2];
			r[i*n]=a[i*2*m];
			r[i*n+n-1]=a[i*2*m+m-1];
		}
		for (i=1; i+1<n; i++) for (j=1; j+1<n; j++) r[i*n+j]=a[i*2*m+j*2];
		return r;
	}
};
class injection2 :public restriction
{
public:
	valarray<db> operator()(const valarray<db> &a) const
	{
		int m=a.size(),n=round(sqrt(m)),i,j;
		if (n*n!=m||n<5||n!=(1<<__lg(n)|1)) throw invalid_argument("Incorrect length.");
		m=n; n=n/2+1;
		valarray<db> r(n*n);
		for (i=0; i<n; i++) for (j=0; j<n; j++) r[i*n+j]=a[i*2*m+j*2];
		return r;
	}
};
class linear2 :public interpolation
{
public:
	valarray<db> operator()(const valarray<db> &a) const
	{
		int m=a.size(),n=round(sqrt(m)),i,j;
		if (n*n!=m||n<3||n!=(1<<__lg(n)|1)) throw invalid_argument("Incorrect length.");
		m=n; n=n*2-1;
		valarray<db> r(n*n);
		for (i=0; i<m; i++) for (j=0; j<m; j++) r[i*2*n+j*2]=a[i*m+j];
		for (i=0; i<m; i++) for (j=1; j<m; j++) r[i*2*n+j*2-1]=(a[i*m+j-1]+a[i*m+j])/2;
		for (i=1; i<m; i++) for (j=0; j<m; j++) r[(i*2-1)*n+j*2]=(a[(i-1)*m+j]+a[i*m+j])/2;
		for (i=1; i<m; i++) for (j=1; j<m; j++) r[(i*2-1)*n+j*2-1]=(a[(i-1)*m+j-1]+a[(i-1)*m+j]+a[i*m+j-1]+a[i*m+j])/4;
		return r;
	}
};
vector<vector<pair<int,db>>> get_A2(int n)
{
	int i,j;
	vector<vector<pair<int,db>>> A(n*n);
	db h=1.0/(n-1);
	const function_value &tmp=function_value(0);
	A[0]={{0,1}}; A[n-1]={{n-1,1}};
	A[(n-1)*n]={{(n-1)*n,1}}; A[n*n-1]={{n*n-1,1}};
	for (i=1; i<n-1; i++) for (j=1; j<n-1; j++) A[i*n+j]={{(i-1)*n+j,1/h/h},{i*n+j-1,1/h/h},{i*n+j,-4/h/h},{i*n+j+1,1/h/h},{(i+1)*n+j,1/h/h}};
	for (i=1; i<n-1; i++)
	{
		A[i]={{i,1}};
		A[i*n]={{i*n,1}};
		A[i*n+n-1]={{i*n+n-1,1}};
		A[(n-1)*n+i]={{(n-1)*n+i,1}};
	}
	return A;
}
class cycles2
{
public:
	virtual void operator()(valarray<db> &v,valarray<db> f,const restriction &rs,const interpolation &inte,int v1=2,int v2=2) const=0;
};
class V_cycle2 :public cycles2
{
public:
	// template<typename func> V_cycle(const func &rhs) :cycles(rhs) {}
	void operator()(valarray<db> &v,valarray<db> f,const restriction &rs,const interpolation &inte,int v1=2,int v2=2) const
	{
		int n=round(sqrt(v.size())),i;
		const auto &A=get_A2(n);
		v=relax(A,v,f,v1);
		if (n>3)
		{
			valarray f2=rs(f-multiply(A,v));
			valarray<db> vv((n/2+1)*(n/2+1));
			V_cycle2()(vv,f2,rs,inte,v1,v2);
			v=v+inte(vv);
		}
		v=relax(A,v,f,v2);
	}
};
class FMG2 :public cycles2
{
public:
	void operator()(valarray<db> &v,valarray<db> f,const restriction &rs,const interpolation &inte,int v1=2,int v2=2) const
	{
		int n=v.size(),i;
		if (n==3) v=0;
		else
		{
			valarray f2=rs(f);
			v.resize((n/2+1)*(n/2+1));
			FMG2()(v,f2,rs,inte,v1,v2);
			v=inte(v);
		}
		V_cycle2()(v,f,rs,inte,v1,v2);
	}
};
function<db(db,db)> ZERO_FUNC2=[&](db x,db y) ->db { return 0; };
class solution2
{
	valarray<db> v;
	int times;
	db residual;
public:
	template<typename func1,typename func2,typename func3>
	solution2(const func1 &_f,const func2 &boundary,int n,int T,db eps,const restriction &rs,const interpolation &inte,const cycles2 &cyc,const func3 &guess,int v1=2,int v2=2) :v(n *n)
	{
		if (n<3||n!=(1<<__lg(n)|1)) throw invalid_argument("Incorrect n.");
		int i,j;
		for (i=0; i<n; i++) for (j=0; j<n; j++) v[i*n+j]=guess(i/(n-1.),j/(n-1.));
		valarray<db> f(n*n);
		f[0]=0; f[n-1]=0; f[(n-1)*n]=0; f[n*n-1]=0;
		for (i=1; i+1<n; i++) for (j=1; j+1<n; j++) f[i*n+j]=_f(i/(n-1.),j/(n-1.));
		for (i=0; i<n; i++) for (j=0; j<n; j++) if (i==0||j==0||i==n-1||j==n-1) f[i*n+j]=boundary(i/(n-1.),j/(n-1.));
		const auto &A=get_A2(n);
		f-=multiply(A,v);
		for (times=1; (times<=T||abs(f).sum()>=eps)&&times<=1000; times++)
		{
			valarray<db> r(n*n);
			cyc(r,f,rs,inte,v1,v2);
			v+=r;
			f-=multiply(A,r);
		}
		residual=abs(f).sum();
	}
	template<typename func1,typename func2>
	solution2(const func1 &_f,const func2 &boundary,int n,int T,db eps,const restriction &rs,const interpolation &inte,const cycles2 &cyc,int v1=2,int v2=2) :solution2(_f,boundary,n,T,eps,rs,inte,cyc,ZERO_FUNC2,v1,v2) {}
	db operator()(db x,db y) const { int n=round(sqrt(v.size())); return v[(int)round(x*(n-1))*n+(int)round(y*(n-1))]; }
	int get_times() const { return times; }
	db get_residual() const { return residual; }
};