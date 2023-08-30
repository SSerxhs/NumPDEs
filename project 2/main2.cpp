#include "bits/stdc++.h"
#include "debug.h"
#include "BPV.h"
using namespace std;
using db=double;
int main(int argc,char *argv[])
{
	auto t1=chrono::steady_clock::now().time_since_epoch().count();
	cout<<"Function: exp(sin(x))\n";
	cout<<" n  32 64 128 256\n";
	auto F=[&](db x,db y) { return exp(y+sin(x)); };
	auto DxF=[&](db x,db y) { assert(x==0||x==1); return exp(y+sin(x))*cos(x); };
	auto DyF=[&](db x,db y) { assert(y==0||y==1); return exp(y+sin(x)); };
	auto D2F=[&](db x,db y) { x=sin(x); return exp(y+x)*(x+2)*(1-x); };
	int n=64,i,j;
	BPV<2> s(D2F,F,n+1,0,1e-8,injection2(),linear2(),V_cycle2(),2,2);
	// for (i=0; i<=n; i++) for (j=0; j<=n; j++) cout<<'('<<i*1./n<<", "<<j*1./n<<") = "<<s(i*1./n,j*1./n)<<' '<<F(i*1./n,j*1./n)<<'\n';
	cout<<s.get_residual()<<endl;
	cout<<chrono::steady_clock::now().time_since_epoch().count()-t1<<endl;
}