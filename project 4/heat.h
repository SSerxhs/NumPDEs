#pragma once
#include "linear.h"
#include "solution.h"
#include "bits/stdc++.h"
using namespace std;
class theta_method
{
	db theta;
public:
	theta_method(db theta):theta(theta){}
	template<typename func> solution solve(const func &boundary,int m,int n,db nu,db r) const
	{
		db h=1.0/m;
		db k=r*h*h/nu;
		vector y(n+1,valarray<db>(m+1));
		{
			auto &y0=y[0];
			for (int i=0;i<=m;i++)
			{
				y0[i]=boundary((db)i/m,0);
			}
		}
		for (int j=1;j<=n;j++)
		{
			const auto &y0=y[j-1];
			valarray<db> l(m),d(m+1),u(m),b(m+1);
			l=-theta*r;
			d=1+2*theta*r;
			u=-theta*r;
			db t=j*k;
			d[0]=1;b[0]=boundary(0,t);
			d[m]=1;b[m]=boundary(1,t);
			l[m-1]=0;u[0]=0;
			for (int i=1;i<m;i++)
			{
				b[i]=(1-theta)*r*y0[i-1]+(1-2*(1-theta)*r)*y0[i]+(1-theta)*r*y0[i+1];
			}
			Tridiagonal_Matrix(l,d,u,b).swap(y[j]);
		}
		return solution(y,0,1,0,n*k);
	}
};
class FTCS:public theta_method
{
public:
	FTCS():theta_method(0){}
};
class Crank_Nicolson:public theta_method
{
public:
	Crank_Nicolson():theta_method(0.5){}
};
class BTCS:public theta_method
{
public:
	BTCS():theta_method(1){}
};
class RKM
{
	vector<valarray<db>> a;
	valarray<db> b,c;
public:
	RKM(const vector<valarray<db>> &a,const valarray<db> &b,const valarray<db> &c)
	:a(a),b(b),c(c){}
	template<typename func> solution solve(const func &boundary,int m,int n,db nu,db r) const
	{
		int s=a.size();
		db h=1.0/m;
		db k=r*h*h/nu;
		vector U(n+1,valarray<db>(m+1));
		{
			auto &U0=U[0];
			for (int i=0;i<=m;i++)
			{
				U0[i]=boundary((db)i/m,0);
			}
		}
		for (int index=1;index<=n;index++)
		{
			const auto &Un=U[index-1];
			vector coef((m+1)*s,valarray<db>((m+1)*s));
			valarray<db> rhs((m+1)*s);
			for (int i=0;i<s;i++)
			{
				coef[i*(m+1)][i*(m+1)]=1;
				coef[i*(m+1)+m][i*(m+1)+m]=1;
				for (int j=1;j<m;j++)
				{
					int x=i*(m+1)+j;
					auto &co=coef[x];
					rhs[x]=Un[j+1]-2*Un[j]+Un[j-1];
					co[x]+=k/r;
					for (int l=0;l<s;l++)
					{
						x=l*(m+1)+j;
						co[x+1]-=k*a[i][l];
						co[x]+=2*k*a[i][l];
						co[x-1]-=k*a[i][l];
					}
				}
			}
			const auto &y=Gauss(coef,rhs);
			for (int i=0;i<=m;i++)
			{
				for (int j=0;j<s;j++) U[index][i]+=b[j]*y[j*(m+1)+i];
			}
			U[index]=Un+k*U[index];
		}
		return solution(U,0,1,0,n*k);
	}
};
class RKM2:public RKM
{
public:
	RKM2():RKM(vector{
		valarray<db>{5/12.,-1/12.},
		valarray<db>{3/4.,1/4.}},
		valarray<db>{3/4.,1/4.},
		valarray<db>{1/3.,1}){}
};
class GLRKM1:public RKM
{
public:
	GLRKM1():RKM(vector{
		valarray<db>{1/2.}},
		valarray<db>{1},
		valarray<db>{1/2.}){}
};