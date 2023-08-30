#pragma once
#include "1D.h"
#include "2D.h"
template<int dim> class BPV;
template<> class BPV<1> :public solution
{
public:
	template<typename func1,typename func2,typename func3>
	BPV<1>(const func1 &_f,const func2 &boundary,int n,int T,db eps,const restriction &rs,const interpolation &inte,const cycles &cyc,const func3 &guess,int v1=2,int v2=2) :solution(_f,boundary,n,T,eps,rs,inte,cyc,guess,v1,v2) {}
	template<typename func1,typename func2>
	BPV<1>(const func1 &_f,const func2 &boundary,int n,int T,db eps,const restriction &rs,const interpolation &inte,const cycles &cyc,int v1=2,int v2=2):solution(_f,boundary,n,T,eps,rs,inte,cyc,v1,v2) {}
};
template<> class BPV<2> :public solution2
{
public:
	template<typename func1,typename func2,typename func3>
	BPV<2>(const func1 &_f,const func2 &boundary,int n,int T,db eps,const restriction &rs,const interpolation &inte,const cycles2 &cyc,const func3 &guess,int v1=2,int v2=2) :solution2(_f,boundary,n,T,eps,rs,inte,cyc,guess,v1,v2) {}
	template<typename func1,typename func2>
	BPV<2>(const func1 &_f,const func2 &boundary,int n,int T,db eps,const restriction &rs,const interpolation &inte,const cycles2 &cyc,int v1=2,int v2=2):solution2(_f,boundary,n,T,eps,rs,inte,cyc,v1,v2) {}

};