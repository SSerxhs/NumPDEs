#pragma once
#include "name.h"
#include "IVP.h"
#include "RK.h"
#include "bits/stdc++.h"
#include "Gauss.h"
using namespace std;
class LMM:IVP
{
protected:
	vector<db> alpha,beta;
public:
	LMM(int s):alpha(s),beta(s) {}
};
class explicit_LMM:public LMM
{
protected:
	explicit_LMM(int s):LMM(s) {}
	template<typename func1,typename func2> solution iteration(const func1 &f,db L,db R,const valarray<db> &xl,int N,const func2 &method) const
	{
		db step=(R-L)/N;
		int s=alpha.size()-1,m=xl.size();
		vector<valarray<db>> x(N+1),kfx(N+1);
		if (s==1)
		{
			x[0]=xl;
			kfx[0]=f(x[0],L)*step;
		}
		else
		{
			solution g=method.solve(f,L,L+s*step,xl,N*s);
			for (int i=0; i<s; i++)
			{
				x[i]=g(L+i*step);
				kfx[i]=f(x[i],L+i*step)*step;
			}
		}
		for (int i=0; i+s<=N; i++)
		{
			valarray<db> &tmp=x[i+s];
			tmp.resize(m);
			for (int j=0; j<s; j++) tmp+=kfx[i+j]*beta[j]-x[i+j]*alpha[j];
			kfx[i+s]=f(tmp,L+i*step)*step;
		}
		return solution(x,L,R);
	}
};
class implicit_LMM:public LMM
{
protected:
	implicit_LMM(int s):LMM(s) {}
	template<typename func1,typename func2> solution iteration(const func1 &f,db L,db R,const valarray<db> &xl,int N,const func2 &method) const
	{
		db step=(R-L)/N;
		int s=alpha.size()-1;
		vector<valarray<db>> x(N+1),kfx(N+1);
		if (s==1)
		{
			x[0]=xl;
			kfx[0]=f(x[0],L)*step;
		}
		else
		{
			solution g=method.solve(f,L,L+s*step,xl,N*s);
			for (int i=0; i<s; i++)
			{
				x[i]=g(L+i*step);
				kfx[i]=f(x[i],L+i*step)*step;
			}
		}
		for (int i=0; i+s<=N; i++)
		{
			x[i+s]=x[i+s-1];
			kfx[i+s]=f(x[i+s],L+i*step)*step;
			valarray<db> tmp;
			int cnt=0;
			do
			{
				tmp=kfx[i+s]*beta[s];
				for (int j=0; j<s; j++) tmp+=kfx[i+j]*beta[j]-x[i+j]*alpha[j];
				swap(x[i+s],tmp);
				kfx[i+s]=f(x[i+s],L+i*step)*step;
			} while (abs(x[i+s]-tmp).max()>1e-14&&++cnt<=1000);
		}
		return solution(x,L,R);
	}
};
class Adams_Bashforth:public explicit_LMM
{
public:
	Adams_Bashforth(int p):Adams_Bashforth(p,p) {}
	Adams_Bashforth(int s,int p):explicit_LMM(s+1)
	{
		++s; ++p;
		vector<vector<db>> A;
		vector<db> b,coefa(s,1),coefb(s,0);
		int i,j;
		for (i=0; i<p; i++)
		{
			vector<db> tmp(s*2);
			for (j=0; j<s; j++)
			{
				tmp[j]=coefa[j];
				coefa[j]*=j/(i+1.);
				tmp[j+s]=coefb[j];
				if (i==0) coefb[j]=-1; else coefb[j]*=(db)j/i;
			}
			A.push_back(tmp);
			b.push_back(0);
		}
		{
			vector<db> tmp(s*2);
			tmp[s-1]=1;
			A.push_back(tmp);
			b.push_back(1);
			for (i=0; i<s-2; i++)
			{
				vector<db> tmp(s*2);
				tmp[i]=1;
				A.push_back(tmp);
				b.push_back(0);
			}
		}
		{
			vector<db> tmp(s*2);
			tmp[s*2-1]=1;
			A.push_back(tmp);
			b.push_back(0);
		}
		assert(A.size()<=s*2);
		int lim=s*2-A.size();
		for (i=0; i<lim; i++)
		{
			vector<db> tmp(s*2);
			tmp[s+i]=1;
			A.push_back(tmp);
			b.push_back(0);
		}
		auto x=Gauss(A,b);
		copy(x.begin(),x.begin()+s,alpha.begin());
		copy(x.begin()+s,x.end(),beta.begin());
	}
	solution solve(const function<valarray<db>(const valarray<db> &,db)> &f,db L,db R,const valarray<db> &xl,int N) const { return iteration(f,L,R,xl,N,Gauss_Legendre_RK(3)); }
};
class Adams_Moulton:public implicit_LMM
{
public:
	Adams_Moulton(int p):Adams_Moulton(max(1,p-1),p) {}
	Adams_Moulton(int s,int p):implicit_LMM(s+1)
	{
		++s; ++p;
		vector<vector<db>> A;
		vector<db> b,coefa(s,1),coefb(s,0);
		int i,j;
		for (i=0; i<p; i++)
		{
			vector<db> tmp(s*2);
			for (j=0; j<s; j++)
			{
				tmp[j]=coefa[j];
				coefa[j]*=j/(i+1.);
				tmp[j+s]=coefb[j];
				if (i==0) coefb[j]=-1; else coefb[j]*=(db)j/i;
			}
			A.push_back(tmp);
			b.push_back(0);
		}
		{
			vector<db> tmp(s*2);
			tmp[s-1]=1;
			A.push_back(tmp);
			b.push_back(1);
			for (i=0; i<s-2; i++)
			{
				vector<db> tmp(s*2);
				tmp[i]=1;
				A.push_back(tmp);
				b.push_back(0);
			}
		}
		assert(A.size()<=s*2);
		int lim=s*2-A.size();
		for (i=0; i<lim; i++)
		{
			vector<db> tmp(s*2);
			tmp[s+i]=1;
			A.push_back(tmp);
			b.push_back(0);
		}
		auto x=Gauss(A,b);
		copy(x.begin(),x.begin()+s,alpha.begin());
		copy(x.begin()+s,x.end(),beta.begin());
	}
	solution solve(const function<valarray<db>(const valarray<db> &,db)> &f,db L,db R,const valarray<db> &xl,int N) const { return iteration(f,L,R,xl,N,Gauss_Legendre_RK(3)); }
};
class BDF:public implicit_LMM
{
public:
	BDF(int p):BDF(p,p) {}
	BDF(int s,int p):implicit_LMM(s+1)
	{
		++s; ++p;
		vector<vector<db>> A;
		vector<db> b,coefa(s,1),coefb(s,0);
		int i,j;
		for (i=0; i<p; i++)
		{
			vector<db> tmp(s*2);
			for (j=0; j<s; j++)
			{
				tmp[j]=coefa[j];
				coefa[j]*=j/(i+1.);
				tmp[j+s]=coefb[j];
				if (i==0) coefb[j]=-1; else coefb[j]*=(db)j/i;
			}
			A.push_back(tmp);
			b.push_back(0);
		}
		{
			vector<db> tmp(s*2);
			tmp[s-1]=1;
			A.push_back(tmp);
			b.push_back(1);
			for (i=0; i<s-1; i++)
			{
				vector<db> tmp(s*2);
				tmp[i+s]=1;
				A.push_back(tmp);
				b.push_back(0);
			}
		}
		assert(A.size()<=s*2);
		int lim=s*2-A.size();
		for (i=0; i<lim; i++)
		{
			vector<db> tmp(s*2);
			tmp[s+i]=1;
			A.push_back(tmp);
			b.push_back(0);
		}
		auto x=Gauss(A,b);
		copy(x.begin(),x.begin()+s,alpha.begin());
		copy(x.begin()+s,x.end(),beta.begin());
	}
	solution solve(const function<valarray<db>(const valarray<db> &,db)> &f,db L,db R,const valarray<db> &xl,int N) const { return iteration(f,L,R,xl,N,Gauss_Legendre_RK(3)); }
};
class Euler:public Adams_Bashforth
{
public:
	Euler():Adams_Bashforth(1) {}
};