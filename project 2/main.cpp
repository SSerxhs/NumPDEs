#include "bits/stdc++.h"
#include "debug.h"
#include "BPV.h"
using namespace std;
using db=double;
int main(int argc,char *argv[])
{
	{
		const db pi=acos(-1);
		cout<<"Function: exp(sin(x))\n";
		cout<<" n  32 64 128 256\n";
		auto f=[&](db x) { return exp(sin(x)); };
		auto Df=[&](db x) { return exp(sin(x))*cos(x); };
		auto DDf=[&](db x) { return exp(sin(x))*(cos(x)*cos(x)-sin(x)); };
		auto boundary_f=[&](db x)
		{
			return unique_ptr<boundary_value>(new function_value(f(x)));
			// if (x==0) return unique_ptr<boundary_value>(new function_value(f(x)));
			// return unique_ptr<boundary_value>(new derivative_value(Df(x)));
		};
		unique_ptr<restriction> pt1;
		unique_ptr<interpolation> pt2;
		unique_ptr<cycles> pt3;
		for (int b1:{0,1}) for (int b2:{0,1}) for (int b3:{0,1})
		{
			if (b1) pt1.reset(new injection); else pt1.reset(new full_weighting);
			if (b2) pt2.reset(new quadratic); else pt2.reset(new linear);
			if (b3) pt3.reset(new FMG); else pt3.reset(new V_cycle);
			cout<<b1<<b2<<b3;
			for (int n:{32,64,128,256})
			{
				BPV<1> s(DDf,boundary_f,n+1,0,1e-8,*pt1,*pt2,*pt3,2,2);
				cout<<' '<<s.get_times();
			}
			cout<<'\n';
		}
		cout<<"after 15\n";
		for (int b1:{0,1}) for (int b2:{0,1}) for (int b3:{0,1})
		{
			if (b1) pt1.reset(new injection); else pt1.reset(new full_weighting);
			if (b2) pt2.reset(new quadratic); else pt2.reset(new linear);
			if (b3) pt3.reset(new FMG); else pt3.reset(new V_cycle);
			cout<<b1<<b2<<b3;
			for (int n:{32,64,128,256})
			{
				BPV<1> s(DDf,boundary_f,n+1,15,1e100,*pt1,*pt2,*pt3,2,2);
				cout<<' '<<s.get_residual();
			}
			cout<<'\n';
		}
		cout<<"error:";
		for (int n:{32,64,128,256})
		{
			int i;
			db err=0;
			BPV<1> s(DDf,boundary_f,n+1,15,1e100,full_weighting(),quadratic(),FMG(),2,2);
			for (i=1; i<n; i++) err+=abs(s(i*1./n)-f(i*1./n));
			cout<<" "<<err;
		}
		cout<<'\n';
	}
	{
		const db pi=acos(-1);
		cout<<"Function: exp(sin(x))\n";
		cout<<" n  32 64 128 256\n";
		auto f=[&](db x) { return exp(sin(x)); };
		auto Df=[&](db x) { return exp(sin(x))*cos(x); };
		auto DDf=[&](db x) { return exp(sin(x))*(cos(x)*cos(x)-sin(x)); };
		auto boundary_f=[&](db x)
		{
			// return unique_ptr<boundary_value>(new function_value(f(x)));
			if (x==0) return unique_ptr<boundary_value>(new function_value(f(x)));
			return unique_ptr<boundary_value>(new derivative_value(Df(x)));
		};
		unique_ptr<restriction> pt1;
		unique_ptr<interpolation> pt2;
		unique_ptr<cycles> pt3;
		for (int b1:{0,1}) for (int b2:{0,1}) for (int b3:{0,1})
		{
			if (b1) pt1.reset(new injection); else pt1.reset(new full_weighting);
			if (b2) pt2.reset(new quadratic); else pt2.reset(new linear);
			if (b3) pt3.reset(new FMG); else pt3.reset(new V_cycle);
			cout<<b1<<b2<<b3;
			for (int n:{32,64,128,256})
			{
				BPV<1> s(DDf,boundary_f,n+1,0,1e-8,*pt1,*pt2,*pt3,2,2);
				cout<<' '<<s.get_times();
			}
			cout<<'\n';
		}
		cout<<"after 30\n";
		for (int b1:{0,1}) for (int b2:{0,1}) for (int b3:{0,1})
		{
			if (b1) pt1.reset(new injection); else pt1.reset(new full_weighting);
			if (b2) pt2.reset(new quadratic); else pt2.reset(new linear);
			if (b3) pt3.reset(new FMG); else pt3.reset(new V_cycle);
			cout<<b1<<b2<<b3;
			for (int n:{32,64,128,256})
			{
				BPV<1> s(DDf,boundary_f,n+1,30,1e100,*pt1,*pt2,*pt3,2,2);
				cout<<' '<<s.get_residual();
			}
			cout<<'\n';
		}
		cout<<"error:";
		for (int n:{32,64,128,256})
		{
			int i;
			db err=0;
			BPV<1> s(DDf,boundary_f,n+1,30,1e100,full_weighting(),quadratic(),FMG(),2,2);
			for (i=1; i<n; i++) err+=abs(s(i*1./n)-f(i*1./n));
			cout<<" "<<err;
		}
		cout<<'\n';
	}
	cout<<"-----------------------------------------------\n";
	{
		const db pi=acos(-1);
		cout<<"Function: sin(pi*x)\n";
		cout<<" n  32 64 128 256\n";
		auto f=[&](db x) { return sin(pi*x); };
		auto Df=[&](db x) { return pi*cos(pi*x); };
		auto DDf=[&](db x) { return -pi*pi*sin(pi*x); };
		auto boundary_f=[&](db x)
		{
			return unique_ptr<boundary_value>(new function_value(f(x)));
			// if (x==0) return unique_ptr<boundary_value>(new function_value(f(x)));
			// return unique_ptr<boundary_value>(new derivative_value(Df(x)));
		};
		unique_ptr<restriction> pt1;
		unique_ptr<interpolation> pt2;
		unique_ptr<cycles> pt3;
		for (int b1:{0,1}) for (int b2:{0,1}) for (int b3:{0,1})
		{
			if (b1) pt1.reset(new injection); else pt1.reset(new full_weighting);
			if (b2) pt2.reset(new quadratic); else pt2.reset(new linear);
			if (b3) pt3.reset(new FMG); else pt3.reset(new V_cycle);
			cout<<b1<<b2<<b3;
			for (int n:{32,64,128,256})
			{
				BPV<1> s(DDf,boundary_f,n+1,0,1e-8,*pt1,*pt2,*pt3,2,2);
				cout<<' '<<s.get_times();
			}
			cout<<'\n';
		}
		cout<<"after 15\n";
		for (int b1:{0,1}) for (int b2:{0,1}) for (int b3:{0,1})
		{
			if (b1) pt1.reset(new injection); else pt1.reset(new full_weighting);
			if (b2) pt2.reset(new quadratic); else pt2.reset(new linear);
			if (b3) pt3.reset(new FMG); else pt3.reset(new V_cycle);
			cout<<b1<<b2<<b3;
			for (int n:{32,64,128,256})
			{
				BPV<1> s(DDf,boundary_f,n+1,15,1e100,*pt1,*pt2,*pt3,2,2);
				cout<<' '<<s.get_residual();
			}
			cout<<'\n';
		}
		cout<<"error:";
		for (int n:{32,64,128,256})
		{
			int i;
			db err=0;
			BPV<1> s(DDf,boundary_f,n+1,15,1e100,full_weighting(),quadratic(),FMG(),2,2);
			for (i=1; i<n; i++) err+=abs(s(i*1./n)-f(i*1./n));
			cout<<" "<<err;
		}
		cout<<'\n';
	}
	{
		const db pi=acos(-1);
		cout<<"Function: sin(pi*x)\n";
		cout<<" n  32 64 128 256\n";
		auto f=[&](db x) { return sin(pi*x); };
		auto Df=[&](db x) { return pi*cos(pi*x); };
		auto DDf=[&](db x) { return -pi*pi*sin(pi*x); };
		auto boundary_f=[&](db x)
		{
			// return unique_ptr<boundary_value>(new function_value(f(x)));
			if (x==0) return unique_ptr<boundary_value>(new function_value(f(x)));
			return unique_ptr<boundary_value>(new derivative_value(Df(x)));
		};
		unique_ptr<restriction> pt1;
		unique_ptr<interpolation> pt2;
		unique_ptr<cycles> pt3;
		for (int b1:{0,1}) for (int b2:{0,1}) for (int b3:{0,1})
		{
			if (b1) pt1.reset(new injection); else pt1.reset(new full_weighting);
			if (b2) pt2.reset(new quadratic); else pt2.reset(new linear);
			if (b3) pt3.reset(new FMG); else pt3.reset(new V_cycle);
			cout<<b1<<b2<<b3;
			for (int n:{32,64,128,256})
			{
				BPV<1> s(DDf,boundary_f,n+1,0,1e-8,*pt1,*pt2,*pt3,2,2);
				cout<<' '<<s.get_times();
			}
			cout<<'\n';
		}
		cout<<"after 30\n";
		for (int b1:{0,1}) for (int b2:{0,1}) for (int b3:{0,1})
		{
			if (b1) pt1.reset(new injection); else pt1.reset(new full_weighting);
			if (b2) pt2.reset(new quadratic); else pt2.reset(new linear);
			if (b3) pt3.reset(new FMG); else pt3.reset(new V_cycle);
			cout<<b1<<b2<<b3;
			for (int n:{32,64,128,256})
			{
				BPV<1> s(DDf,boundary_f,n+1,30,1e100,*pt1,*pt2,*pt3,2,2);
				cout<<' '<<s.get_residual();
			}
			cout<<'\n';
		}
		cout<<"error:";
		for (int n:{32,64,128,256})
		{
			int i;
			db err=0;
			BPV<1> s(DDf,boundary_f,n+1,30,1e100,full_weighting(),quadratic(),FMG(),2,2);
			for (i=1; i<n; i++) err+=abs(s(i*1./n)-f(i*1./n));
			cout<<" "<<err;
		}
		cout<<'\n';
	}
	cout<<"-----------------------------------------------\n";
	{
		const db pi=acos(-1);
		cout<<"Function: 1/(x+1)\n";
		cout<<" n  32 64 128 256\n";
		auto f=[&](db x) { return 1/(x+1); };
		auto Df=[&](db x) { return -1/pow(x+1,2); };
		auto DDf=[&](db x) { return 2/pow(x+1,3); };
		auto boundary_f=[&](db x)
		{
			return unique_ptr<boundary_value>(new function_value(f(x)));
			// if (x==0) return unique_ptr<boundary_value>(new function_value(f(x)));
			// return unique_ptr<boundary_value>(new derivative_value(Df(x)));
		};
		unique_ptr<restriction> pt1;
		unique_ptr<interpolation> pt2;
		unique_ptr<cycles> pt3;
		for (int b1:{0,1}) for (int b2:{0,1}) for (int b3:{0,1})
		{
			if (b1) pt1.reset(new injection); else pt1.reset(new full_weighting);
			if (b2) pt2.reset(new quadratic); else pt2.reset(new linear);
			if (b3) pt3.reset(new FMG); else pt3.reset(new V_cycle);
			cout<<b1<<b2<<b3;
			for (int n:{32,64,128,256})
			{
				BPV<1> s(DDf,boundary_f,n+1,0,1e-8,*pt1,*pt2,*pt3,2,2);
				cout<<' '<<s.get_times();
			}
			cout<<'\n';
		}
		cout<<"after 15\n";
		for (int b1:{0,1}) for (int b2:{0,1}) for (int b3:{0,1})
		{
			if (b1) pt1.reset(new injection); else pt1.reset(new full_weighting);
			if (b2) pt2.reset(new quadratic); else pt2.reset(new linear);
			if (b3) pt3.reset(new FMG); else pt3.reset(new V_cycle);
			cout<<b1<<b2<<b3;
			for (int n:{32,64,128,256})
			{
				BPV<1> s(DDf,boundary_f,n+1,15,1e100,*pt1,*pt2,*pt3,2,2);
				cout<<' '<<s.get_residual();
			}
			cout<<'\n';
		}
		cout<<"error:";
		for (int n:{32,64,128,256})
		{
			int i;
			db err=0;
			BPV<1> s(DDf,boundary_f,n+1,15,1e100,full_weighting(),quadratic(),FMG(),2,2);
			for (i=1; i<n; i++) err+=abs(s(i*1./n)-f(i*1./n));
			cout<<" "<<err;
		}
		cout<<'\n';
	}
	{
		const db pi=acos(-1);
		cout<<"Function: 1/(x+1)\n";
		cout<<" n  32 64 128 256\n";
		auto f=[&](db x) { return 1/(x+1); };
		auto Df=[&](db x) { return -1/pow(x+1,2); };
		auto DDf=[&](db x) { return 2/pow(x+1,3); };
		auto boundary_f=[&](db x)
		{
			// return unique_ptr<boundary_value>(new function_value(f(x)));
			if (x==0) return unique_ptr<boundary_value>(new function_value(f(x)));
			return unique_ptr<boundary_value>(new derivative_value(Df(x)));
		};
		unique_ptr<restriction> pt1;
		unique_ptr<interpolation> pt2;
		unique_ptr<cycles> pt3;
		for (int b1:{0,1}) for (int b2:{0,1}) for (int b3:{0,1})
		{
			if (b1) pt1.reset(new injection); else pt1.reset(new full_weighting);
			if (b2) pt2.reset(new quadratic); else pt2.reset(new linear);
			if (b3) pt3.reset(new FMG); else pt3.reset(new V_cycle);
			cout<<b1<<b2<<b3;
			for (int n:{32,64,128,256})
			{
				BPV<1> s(DDf,boundary_f,n+1,0,1e-8,*pt1,*pt2,*pt3,2,2);
				cout<<' '<<s.get_times();
			}
			cout<<'\n';
		}
	}
	{
		const db pi=acos(-1);
		auto f=[&](db x) { return exp(sin(x)); };
		auto Df=[&](db x) { return exp(sin(x))*cos(x); };
		auto DDf=[&](db x) { return exp(sin(x))*(cos(x)*cos(x)-sin(x)); };
		auto boundary_f=[&](db x)
		{
			return unique_ptr<boundary_value>(new function_value(f(x)));
			// if (x==0) return unique_ptr<boundary_value>(new function_value(f(x)));
			// return unique_ptr<boundary_value>(new derivative_value(Df(x)));
		};
		unique_ptr<restriction> pt1;
		unique_ptr<interpolation> pt2;
		unique_ptr<cycles> pt3;
		for (int b1:{0,1}) for (int b2:{0,1}) for (int b3:{0,1})
		{
			if (b1) pt1.reset(new injection); else pt1.reset(new full_weighting);
			if (b2) pt2.reset(new quadratic); else pt2.reset(new linear);
			if (b3) pt3.reset(new FMG); else pt3.reset(new V_cycle);
			cout<<b1<<b2<<b3;
			for (int n:{32,64,128,256})
			{
				for (db eps=1e-8; eps>=2.2e-16; eps*=0.1)
				{
					BPV<1> s(DDf,boundary_f,n+1,0,eps,*pt1,*pt2,*pt3,2,2);
					if (s.get_residual()>eps)
					{
						cout<<' '<<eps;
						break;
					}
				}
			}
			cout<<'\n';
		}
	}
}