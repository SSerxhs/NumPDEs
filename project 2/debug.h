#pragma once
#include "bits/stdc++.h"
using namespace std;
template<typename typC,typename typD> istream &operator>>(istream &cin,pair<typC,typD> &a) { return cin>>a.first>>a.second; }
template<typename typC> istream &operator>>(istream &cin,vector<typC> &a) { for (auto &x:a) cin>>x; return cin; }
template<typename typC> ostream &operator<<(ostream &cout,const vector<typC> &a)
{
	cout<<'[';
	int n=a.size();
	if (n)
	{
		cout<<a[0]; for (int i=1; i<n; i++) cout<<", "<<a[i];

	}
	cout<<"]";
	return cout;
}
template<typename Tuple,size_t N>
struct tuple_print
{
	static void print(const Tuple &t,std::ostream &os)
	{
		tuple_print<Tuple,N-1>::print(t,os);
		os<<", "<<std::get<N-1>(t);
	}
};
template<typename Tuple>
struct tuple_print<Tuple,1>
{
	static void print(const Tuple &t,std::ostream &os)
	{
		os<<std::get<0>(t);
	}
};
template<typename... Args>
ostream &operator<<(std::ostream &os,const std::tuple<Args...> &t)
{
	os<<"{";
	tuple_print<decltype(t),sizeof...(Args)>::print(t,os);
	os<<"}";
	return os;
}
template<typename T1,typename T2> ostream &operator<<(ostream &cout,const pair<T1,T2> &x)
{
	return cout<<"{"<<x.first<<", "<<x.second<<"}";
}
template<typename typC,typename typD> bool cmin(typC &x,const typD &y) { if (y<x) { x=y; return 1; } return 0; }
template<typename typC,typename typD> bool cmax(typC &x,const typD &y) { if (x<y) { x=y; return 1; } return 0; }
template<typename typC> vector<typC> range(typC l,typC r,typC step=1) { assert(step>0); int n=(r-l+step-1)/step,i; vector<typC> res(n); for (i=0; i<n; i++) res[i]=l+step*i; return res; }
template<typename typC> ostream &operator<<(ostream &cout,const valarray<typC> &a)
{
	cout<<'[';
	int n=a.size();
	if (n)
	{
		cout<<a[0]; for (int i=1; i<n; i++) cout<<", "<<a[i];

	}
	cout<<"]";
	return cout;
}
