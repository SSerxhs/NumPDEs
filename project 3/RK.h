#pragma once
#include "IVP.h"
#include "name.h"
#include "bits/stdc++.h"
using namespace std;
class RK:IVP
{
protected:
	vector<valarray<db>> a;
	valarray<db> b,c;
	RK(const vector<valarray<db>> &_a,const valarray<db> &_b,const valarray<db> &_c):a(_a),b(_b),c(_c) {}
};
class explicit_RK:public RK
{
public:
	explicit_RK(const vector<valarray<db>> &a,const valarray<db> &b,const valarray<db> &c):RK(a,b,c) {}
	solution solve(const function<valarray<db>(const valarray<db> &,db)> &f,db L,db R,const valarray<db> &xl,int N) const
	{
		vector<valarray<db>> x(N+1,xl);
		db step=(R-L)/N;
		int s=a.size(),m=xl.size();
		vector<valarray<db>> y(s,valarray<db>(m));
		valarray<db> cur(m);
		for (int i=0; i<N; i++)
		{
			x[i+1]=x[i];
			for (int j=0; j<s; j++)
			{
				cur=0;
				for (int k=0; k<j; k++) cur+=a[j][k]*y[k];
				y[j]=f(x[i]+step*cur,L+(i+c[j])*step);
				// cerr<<b<<' '<<y<<endl;
				x[i+1]+=b[j]*step*y[j];
			}
		}
		return solution(x,L,R);
	}
};
class implicit_RK:public RK
{
public:
	implicit_RK(const vector<valarray<db>> &a,const valarray<db> &b,const valarray<db> &c):RK(a,b,c) {}
	solution solve(const function<valarray<db>(const valarray<db> &,db)> &f,db L,db R,const valarray<db> &xl,int N) const
	{
		vector<valarray<db>> x(N+1,xl);
		db step=(R-L)/N;
		int s=a.size(),m=xl.size();
		vector<valarray<db>> y(s,valarray<db>(m));
		auto tmp=y;
		valarray<db> cur(m);
		for (int i=0; i<N; i++)
		{
			db err;
			int cnt=0;;
			do
			{
				err=0;
				for (int j=0; j<s; j++)
				{
					cur=0;
					for (int k=0; k<s; k++) cur+=a[j][k]*y[k];
					tmp[j]=f(x[i]+step*cur,L+(i+c[j])*step);
					cmax(err,abs(y[j]-tmp[j]).max());
				}
				swap(y,tmp);
			} while (err>1e-14&&++cnt<=1000);
			x[i+1]=x[i];
			for (int j=0; j<s; j++) x[i+1]+=b[j]*step*y[j];
		}
		return solution(x,L,R);
	}
};
class Gauss_Legendre_RK:public implicit_RK
{
public:
	const static vector<vector<valarray<db>>> A;
	const static vector<valarray<db>> B;
	const static vector<valarray<db>> C;
	Gauss_Legendre_RK(int s):implicit_RK(A[s-1],B[s-1],C[s-1]) {}
};
class embedded_RK:public explicit_RK
{
protected:
	valarray<db> d;
public:
	embedded_RK(const vector<valarray<db>> &_a,const valarray<db> &_b,const valarray<db> &_d,const valarray<db> &_c):d(_d),explicit_RK(_a,_b,_c) {}
};
class classical_RK:public explicit_RK
{
public:classical_RK():explicit_RK(vector<valarray<db>>{
	valarray<db>{0,0,0,0},
		valarray<db>{0.5,0,0,0},
		valarray<db>{0,0.5,0,0},
		valarray<db>{0,0,1,0}},
	valarray<db>{1.0/6,1.0/3,1.0/3,1.0/6},
	valarray<db>{0,0.5,0.5,1})
{
}
};
class ESDIRK:public implicit_RK
{
public:
	ESDIRK():implicit_RK(vector<valarray<db>>{
		valarray<db>{0,0,0,0,0,0},
			valarray<db>{1/4.,1/4.,0,0,0,0},
			valarray<db>{8611/62500.,-1743/31250.,1/4.,0,0,0},
			valarray<db>{5012029/34652500.,-654441/2922500.,174375/388108.,1/4.,0,0},
			valarray<db>{15267082809/155376265600.,-71443401/120774400.,730878875/902184768.,2285395/8070912.,1/4.,0},
			valarray<db>{82889/524892.,0,15625/83664.,69875/102672.,-2260/8211.,1/4.}},
		valarray<db>{82889/524892.,0,15625/83664.,69875/102672.,-2260/8211.,1/4.},
		valarray<db>{0,1/2.,89/250.,31/50.,17/20.,1})
	{
	}
};
const vector<vector<valarray<db>>> Gauss_Legendre_RK::A={
	vector<valarray<db>>{
	valarray<db>{1/2.}
},vector<valarray<db>>{
	valarray<db>{1/4.,(3-2*sqrtl(3))/12},
		valarray<db>{(3+2*sqrtl(3))/12,1/4.}
},vector<valarray<db>>{
	valarray<db>{5/36.,2/9.-sqrtl(15)/15,5/36.-sqrtl(15)/30},
		valarray<db>{5/36.+sqrtl(15)/24,2/9.,5/36.-sqrtl(15)/24},
		valarray<db>{5/36.+sqrtl(15)/30,2/9.+sqrtl(15)/15,5/36.},
}
};
const vector<valarray<db>> Gauss_Legendre_RK::B={
	valarray<db>{1},
	valarray<db>{1/2.,1/2.},
	valarray<db>{5/18.,4/9.,5/18.}
};
const vector<valarray<db>> Gauss_Legendre_RK::C={
	valarray<db>{1/2.},
	valarray<db>{(3-sqrtl(3))/6,(3+sqrtl(3))/6},
	valarray<db>{(5-sqrtl(15)/10,1/2.,(5+sqrtl(15))/10)}
};
class Fehlberg_embedded_RK:public embedded_RK
{
public:
	Fehlberg_embedded_RK():embedded_RK(vector<valarray<db>>{
		valarray<db>{0,0,0,0,0,0},
			valarray<db>{1/4.,0,0,0,0,0},
			valarray<db>{3/32.,9/32.,0,0,0,0},
			valarray<db>{1932/2197.,-7200/2197.,7296/2197.,0,0,0},
			valarray<db>{439/216.,-8,3680/513.,-845/4104.,0,0},
			valarray<db>{-8/27.,2,-3544/2565.,1859/4104.,-11/40.,0}
	},
		valarray<db>{25/216.,0,1408/2565.,2197/4104.,-1/5.,0},
		valarray<db>{16/135.,0,6656/12825.,28561/56340.,-9/50.,2/55.},
		valarray<db>{0,1/4.,3/8.,12/13.,1,1/.2})
	{
	}
};
class Dormand_Prince_embedded_RK:public embedded_RK
{
public:Dormand_Prince_embedded_RK():embedded_RK(vector<valarray<db>>{
	valarray<db>{0,0,0,0,0,0,0},
		valarray<db>{1/5.,0,0,0,0,0,0},
		valarray<db>{3/40.,9/40.,0,0,0,0,0},
		valarray<db>{44/45.,-56/15.,32/9.,0,0,0,0},
		valarray<db>{19372/6561.,-25360/2187.,64448/6561.,-212/729.,0,0,0},
		valarray<db>{9017/3168.,-355/33.,46732/5247.,49/176.,-5103/18656.,0,0},
		valarray<db>{35/384.,0,500/1113.,125/192.,-2187/6784.,11/84.,0}
},
valarray<db>{35/384.,0,500/1113.,125/192.,-2187/6784.,11/84.,0},
valarray<db>{5179/57600.,0,7571/16695.,393/640.,-92091/339200.,187/2100.,1/40.},
valarray<db>{0,1/5.,3/10.,4/5.,8/9.,1,1})
{
}
};