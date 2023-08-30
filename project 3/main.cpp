#include "bits/stdc++.h"
#include "name.h"
#include "debug.h"
#include "Gauss.h"
#include "LMM.h"
#include "RK.h"
#include "factory.h"
using namespace std;
const db mu=0.012277471;
factory<IVP> &initialize_factory()
{
	auto &s=factory<IVP>::create_factory();
	for (int p:{1,2,3,4})
	{
		s.insert("Adams-Bashforth("s+to_string(p)+")"s,shared_ptr<IVP>((IVP *)(new Adams_Bashforth(p))));
	}
	for (int p:{2,3,4,5})
	{
		s.insert("Adams-Moulton("s+to_string(p)+")"s,shared_ptr<IVP>((IVP *)(new Adams_Moulton(p))));
	}
	for (int p:{1,2,3,4})
	{
		s.insert("BDF("s+to_string(p)+")"s,shared_ptr<IVP>((IVP *)(new BDF(p))));
	}
	s.insert("classical_RK"s,shared_ptr<IVP>((IVP *)(new classical_RK())));
	s.insert("ESDIRK"s,shared_ptr<IVP>((IVP *)(new ESDIRK())));
	for (int p:{1,2,3})
	{
		s.insert("Gauss-Legendre RK("s+to_string(p)+")"s,shared_ptr<IVP>((IVP *)(new Gauss_Legendre_RK(p))));
	}
	s.insert("classical RK"s,shared_ptr<IVP>((IVP *)(new classical_RK())));
	s.insert("ESDIRK"s,shared_ptr<IVP>((IVP *)(new ESDIRK())));
	s.insert("Fehlberg embedded RK"s,shared_ptr<IVP>((IVP *)(new Fehlberg_embedded_RK())));
	s.insert("Dormand-Prince embedded RK"s,shared_ptr<IVP>((IVP *)(new Dormand_Prince_embedded_RK())));
	s.insert("Euler"s,shared_ptr<IVP>((IVP *)(new Euler())));
	return s;
}
class model
{
public:
	valarray<db> x0;
	db T;
	model(const valarray<db> &_x0,const db &_T):x0(_x0),T(_T) {}
};
class question
{
public:
	string name;
	int id,n,t;
};
istream &operator>>(istream &in,question &o)
{
	string s;
	if (getline(in,s))
	{
		s=s.substr(s.find('"')+1);
		o.name=s.substr(0,s.find('"'));
		s=s.substr(s.find(',')+1);
		o.id=stoi(s);
		s=s.substr(s.find(',')+1);
		o.n=stoi(s);
		s=s.substr(s.find(',')+1);
		o.t=stoi(s);
	}
	return in;
}
int main()
{
	cout<<scientific<<setprecision(3);
	const auto f=[](const valarray<db> &x,db t)
	{
		db tmp0=x[0]+mu-1;
		db tmp1=x[1]*x[1]+x[2]*x[2]+tmp0*tmp0;
		db tmp2=x[1]*x[1]+x[2]*x[2]+(tmp0+1)*(tmp0+1);
		tmp1=mu/sqrt(tmp1*tmp1*tmp1);
		tmp2=(1-mu)/sqrt(tmp2*tmp2*tmp2);
		return valarray<db>{
			x[3],
				x[4],
				x[5],
				2*x[4]+x[0]-(x[0]+mu-1)*tmp1-(x[0]+mu)*tmp2,
				-2*x[3]+x[1]-x[1]*(tmp1+tmp2),
				-x[2]*(tmp1+tmp2)
		};
	};
	vector<model> Q;
	Q.emplace_back(valarray<db>{0.994,0,0,0,-2.0015851063790825224,0},17.06521656015796);
	Q.emplace_back(valarray<db>{0.879779227778,0,0,0,-0.379677780949,0},19.140540691377);
	auto &ID=initialize_factory();
	auto solve=[&](const model &q,int n,int t,shared_ptr<IVP> p)
	{
		const auto &[x0,T]=q;
		n=n*t;
		return p->solve(f,0,T*t,x0,n);
	};
	auto plot=[&]()
	{
		ifstream in("plot.in"s);
		question tmp;
		while (in>>tmp)
		{
			auto [name,id,n,t]=tmp;
			ofstream out(name+".in");
			out<<fixed<<setprecision(9);
			auto res=solve(Q[id],n,t,ID[name]);
			// n*=Q[id].T;
			db step=Q[id].T/n;
			int k=1;
			out<<n*t/k+1<<endl;
			for (int i=0; i<=n*t; i+=k)
			{
				auto tmp=res(step*i);
				out<<tmp[0]<<' '<<tmp[1]<<'\n';
			}
		}
	};
	auto error=[&](int id)
	{
		ifstream in("error"+to_string(id)+".in"s);
		ofstream out("error"+to_string(id)+".out"s);
		question tmp;
		while (in>>tmp)
		{
			const auto &[name,id,n,t]=tmp;
			auto begin=chrono::steady_clock::now();
			auto res=solve(Q[id],n,t,ID[name]);
			auto end=chrono::steady_clock::now();
			out<<scientific<<setprecision(2)<<name<<": step = "<<n<<", error norm = "<<abs(res(Q[id].T)-Q[id].x0).max()<<fixed<<setprecision(3)<<", run time = "<<((chrono::duration<double>)(end-begin)).count()<<" seconds, "<<"error = "<<scientific<<setprecision(2)<<valarray(res(Q[id].T)-Q[id].x0)<<"."<<endl;
		}
	};
	auto error_latex=[&](int id)
	{
		ifstream in("error"+to_string(id)+".in"s);
		ofstream out("error"+to_string(id)+".tex"s);
		question tmp;
		out<<"\\begin{tabular}{c|c|c|c|c|c|c}\n";
		out<<"名字&$N$&误差 $1$&用时 $1$&误差 $2$&用时 $2$&收敛阶\\\\\n";
		while (in>>tmp)
		{
			db er;
			{
				const auto &[name,id,n,t]=tmp;
				auto begin=chrono::steady_clock::now();
				auto res=solve(Q[id],n,t,ID[name]);
				auto end=chrono::steady_clock::now();
				out<<scientific<<setprecision(2)<<"\\texttt{"<<name<<"}&$"<<n<<"$&$"<<(er=abs(res(Q[id].T)-Q[id].x0).max())<<"$&$"<<fixed<<setprecision(3)<<((chrono::duration<double>)(end-begin)).count()<<"$&$";
			}
			in>>tmp;
			{
				const auto &[name,id,n,t]=tmp;
				auto begin=chrono::steady_clock::now();
				auto res=solve(Q[id],n,t,ID[name]);
				auto end=chrono::steady_clock::now();
				out<<scientific<<setprecision(2)<<abs(res(Q[id].T)-Q[id].x0).max()<<"$&$"<<fixed<<setprecision(3)<<((chrono::duration<double>)(end-begin)).count()<<"$&$"<<log2(er/abs(res(Q[id].T)-Q[id].x0).max())<<"$\\\\"<<endl;
			}
		}
		out<<"\\end{tabular}\n";
	};
	auto bsearch=[&]()
	{
		ifstream in("bsearch.in"s);
		ofstream out("bsearch.out"s);
		question tmp;
		while (in>>tmp)
		{
			auto [name,id,n,t]=tmp;
			n=800;
			while (n<=10'000'000&&abs(solve(Q[id],n,t,ID[name])(Q[id].T)-Q[id].x0).max()>1e-3)
			{
				n*=2;
			}
			int l=n/2,r=n;
			while (r-l>1000)
			{
				n=(l+r)/2;
				if (abs(solve(Q[id],n,t,ID[name])(Q[id].T)-Q[id].x0).max()>1e-3) l=n+1; else r=n;
			}
			n=(l+r)/2;
			auto begin=chrono::steady_clock::now();
			auto res=solve(Q[id],l,t,ID[name]);
			auto end=chrono::steady_clock::now();
			out<<scientific<<setprecision(5)<<name<<": step = "<<n<<", run time in ["<<((chrono::duration<double>)(end-begin)).count();
			begin=chrono::steady_clock::now();
			res=solve(Q[id],r,t,ID[name]);
			end=chrono::steady_clock::now();
			out<<", "<<((chrono::duration<double>)(end-begin)).count()<<"] seconds."<<endl;
		}
	};
	auto bsearch2=[&]()
	{
		ifstream in("bsearch2.in"s);
		ofstream out("bsearch2.out"s);
		question tmp;
		while (in>>tmp)
		{
			auto [name,id,n,t]=tmp;
			n=100;
			while (n<=10'000'000&&abs(solve(Q[id],n,t,ID[name])(Q[id].T)-Q[id].x0).max()>1e-3)
			{
				n*=2;
			}
			int l=1,r=n;
			while (l<r)
			{
				n=(l+r)/2;
				if (abs(solve(Q[id],n,t,ID[name])(Q[id].T)-Q[id].x0).max()>1e-3) l=n+1; else r=n;
			}
			n=l;
			auto begin=chrono::steady_clock::now();
			auto res=solve(Q[id],n,t,ID[name]);
			auto end=chrono::steady_clock::now();
			out<<scientific<<setprecision(5)<<name<<": step = "<<n<<", run time: "<<((chrono::duration<double>)(end-begin)).count()<<" seconds."<<endl;
		}
	};
	auto test=[&]()
	{
		ofstream out("tmp.in");
		string name="Adams-Bashforth(1)";
		int id=1;
		int t=1;
		int n=5e6;
		auto res=solve(Q[id],n,t,ID[name]);
		db step=Q[id].T/n;
		int k=1;
		out<<n*t/k+1<<endl;
		for (int i=0; i<=n*t; i+=k)
		{
			auto tmp=res(step*i);
			out<<tmp[0]<<' '<<tmp[1]<<'\n';
		}

	};
	plot();
	// error(1);
	// error_latex(1);
	// error(2);
	// error_latex(2);
	// bsearch();
	// bsearch2();
	// test();
}