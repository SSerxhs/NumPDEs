#pragma once
#include "basic.h"
#include "multigrid.h"
template <class func_f,class func_boundary>
class FVMOL_1  // fourth-order FV-MOL method
{
	const static db gamma;
	const static valarray<db> b,c;
	const static vector<valarray<db>> AE,AI;
	const static int ns;
	func_f f;
	func_boundary boundary;
	int N,n,step;
	db nu,h,k;
	point u;
	vector<valarray<db>> phi;
	class equations
	{
		func_f f;
		point u;
		func_boundary boundary;
		int N,n,shift;
		db h,nu;
		sparse lhs;
		valarray<db> coef;
		valarray<db> rhs;

	public:
		equations(const point &_u,const func_f &_f,const func_boundary &_boundary,db _nu,int _N,bool flg=1)
			: u(_u),f(_f),boundary(_boundary),N(_N),n(N *N),h(1.0/N),nu(_nu),lhs(flg?n:0),
			coef(flg?0:n),rhs(n)
		{
			if (flg)
			{
				for (int i=0; i<n; i++)
					add(lhs,i,i,-1);
			}
		};
		void add_coef(int x,int y,db d,db t)
		{
			if (x>=0&&x<N&&y>=0&&y<N)
			{
				add(lhs,shift,x*N+y,d);
				return;
			}
			assert((x>=0&&x<N)||(y>=0&&y<N));
			assert((x>=0-2&&x<N+2)&&(y>=0-2&&y<N+2));
			const static vector<valarray<db>> coef_D={
				{4,-13.0/3,5.0/3,-1.0/3},{16,-70.0/3,32.0/3,-7.0/3}};
			const static vector<valarray<db>> coef_N={
				{6.0/5,5.0/10,9.0/10,-5.0/10,1.0/10},{6,-75.0/10,145.0/10,-75.0/10,15.0/10}};
			int dx=1,dy=1,typ,i;
			if (y>=0&&y<N)
			{
				dy=0;
				if (x<0)
				{
					typ=-x-1;
					dx=-dx;
					x=0;
				}
				else
				{
					typ=x-N;
					x=N-1;
				}
			}
			else
			{
				dx=0;
				if (y<0)
				{
					typ=-y-1;
					dy=-dy;
					y=0;
				}
				else
				{
					typ=y-N;
					y=N-1;
				}
			}
			const auto &current=(dx==0?boundary((x+0.5)*h,!!y,t):boundary(!!x,(y+0.5)*h,t));
			const auto &coef=(current.index()==0?coef_D:coef_N)[typ];
			db value=(current.index()==0?get<0>(current).value():get<1>(current).value()*h);
			for (i=1; i<coef.size(); i++)
				add_coef(x-(i-1)*dx,y-(i-1)*dy,coef[i]*d,t);
			rhs[shift]+=coef[0]*d*value;
		}
		void add_rhs(int x,int y,db d,db t)
		{
			if (x>=0&&x<N&&y>=0&&y<N)
			{
				rhs[shift]+=d*coef[x*N+y];
				return;
			}
			assert((x>=0&&x<N)||(y>=0&&y<N));
			assert((x>=0-2&&x<N+2)&&(y>=0-2&&y<N+2));
			const static vector<valarray<db>> coef_D={
				{4,-13.0/3,5.0/3,-1.0/3},{16,-70.0/3,32.0/3,-7.0/3}};
			const static vector<valarray<db>> coef_N={
				{6.0/5,5.0/10,9.0/10,-5.0/10,1.0/10},{6,-75.0/10,145.0/10,-75.0/10,15.0/10}};
			int dx=1,dy=1,typ,i;
			if (y>=0&&y<N)
			{
				dy=0;
				if (x<0)
				{
					typ=-x-1;
					dx=-dx;
					x=0;
				}
				else
				{
					typ=x-N;
					x=N-1;
				}
			}
			else
			{
				dx=0;
				if (y<0)
				{
					typ=-y-1;
					dy=-dy;
					y=0;
				}
				else
				{
					typ=y-N;
					y=N-1;
				}
			}
			const auto &current=(dx==0?boundary((x+0.5)*h,!!y,t):boundary(!!x,(y+0.5)*h,t));
			const auto &coef=(current.index()==0?coef_D:coef_N)[typ];
			db value=(current.index()==0?get<0>(current).value():get<1>(current).value()*h);
			for (i=1; i<coef.size(); i++)
				add_rhs(x-(i-1)*dx,y-(i-1)*dy,coef[i]*d,t);
			rhs[shift]+=coef[0]*d*value;
		}
		void add_coef_dx(int x,int y,db d,db t)
		{
			d/=12*h;
			add_coef(x+1,y,15*d,t);
			add_coef(x,y,-15*d,t);
			add_coef(x+2,y,-d,t);
			add_coef(x-1,y,d,t);
		}
		void add_rhs_dx(int x,int y,db d,db t)
		{
			d/=12*h;
			add_rhs(x+1,y,15*d,t);
			add_rhs(x,y,-15*d,t);
			add_rhs(x+2,y,-d,t);
			add_rhs(x-1,y,d,t);
		}
		void add_coef_dy(int x,int y,db d,db t)
		{
			d/=12*h;
			add_coef(x,y+1,15*d,t);
			add_coef(x,y,-15*d,t);
			add_coef(x,y+2,-d,t);
			add_coef(x,y-1,d,t);
		}
		void add_rhs_dy(int x,int y,db d,db t)
		{
			d/=12*h;
			add_rhs(x,y+1,15*d,t);
			add_rhs(x,y,-15*d,t);
			add_rhs(x,y+2,-d,t);
			add_rhs(x,y-1,d,t);
		}
		void add_coef_diff(int x,int y,db d,db t)
		{
			d*=nu/h;
			add_coef_dx(x,y,d,t);
			add_coef_dx(x-1,y,-d,t);
			add_coef_dy(x,y,d,t);
			add_coef_dy(x,y-1,-d,t);
		}
		void add_rhs_diff(int x,int y,db d,db t)
		{
			d*=nu/h;
			add_rhs_dx(x,y,d,t);
			add_rhs_dx(x-1,y,-d,t);
			add_rhs_dy(x,y,d,t);
			add_rhs_dy(x,y-1,-d,t);
		}
		void add_rhs_x(int x,int y,db d,db t)
		{
			d/=12;
			add_rhs(x,y,d*7,t);
			add_rhs(x+1,y,d*7,t);
			add_rhs(x-1,y,-d,t);
			add_rhs(x+2,y,-d,t);
		}
		void add_rhs_y(int x,int y,db d,db t)
		{
			d/=12;
			add_rhs(x,y,d*7,t);
			add_rhs(x,y+1,d*7,t);
			add_rhs(x,y-1,-d,t);
			add_rhs(x,y+2,-d,t);
		}
		void add_rhs_adv(int x,int y,db d,db t)
		{
			d/=h;
			add_rhs_x(x,y,-d*u.x,t);
			add_rhs_x(x-1,y,d*u.x,t);
			add_rhs_y(x,y,-d*u.y,t);
			add_rhs_y(x,y-1,d*u.y,t);
		}
		void add_rhs_f(int x,int y,db d,db t)
		{
			rhs[shift]+=f((x+.5)*h,(y+.5)*h,t)*d;
		}
		void set_coef(const valarray<db> &_coef)
		{
			coef=_coef;
		}
		void set_shift(int x)
		{
			shift=x;
		}
		void set_rhs(const valarray<db> &_rhs)
		{
			rhs=_rhs;
		}
		sparse &get_lhs()
		{
			return lhs;
		}
		valarray<db> &get_rhs()
		{
			return rhs;
		}
	};
	void iteration()
	{
		int T=phi.size()-1;
		cerr<<"Iteration "<<T<<endl;
		const auto &phi_n=phi.back();
		int i,j,x,y,s;
		vector phi_tmp(ns,valarray<db>(n));
		for (s=0; s<ns; s++)
		{
			if (s==0)
				phi_tmp[0]=phi_n;
			else
			{
				using LHS=sparse;
				using RHS=valarray<db>;
				vector<LHS> lhs_pool(2);
				vector<RHS> rhs_pool(2);
				auto get_equations=[&]()
					{
						int N=1<<lhs_pool.size(),x,y,j;
						equations t(u,f,boundary,nu,N);
						t.set_rhs(reset_size(phi_n,N));
						for (x=0; x<N; x++)
							for (y=0; y<N; y++)
							{
								t.set_shift(x*N+y);
								t.add_coef_diff(x,y,k*gamma,k*(T+c[s]));
							}
						for (j=0; j<s; j++)
						{
							equations q(u,f,boundary,nu,N,0);
							q.set_coef(reset_size(phi_tmp[j],N));
							for (x=0; x<N; x++)
								for (y=0; y<N; y++)
								{
									q.set_shift(x*N+y);
									q.add_rhs_diff(x,y,k*AI[s][j],k*(T+c[j]));
									q.add_rhs_adv(x,y,k*AE[s][j],k*(T+c[j]));
									q.add_rhs_f(x,y,k*AE[s][j],k*(T+c[j]));
								}
							t.get_rhs()+=q.get_rhs();
						}
						lhs_pool.push_back(move(t.get_lhs()));
						rhs_pool.push_back(-t.get_rhs());
						unique_matrix(lhs_pool.back());
						// cerr << "end N = " << N << endl;
					};
				auto get_lhs=[&](int M) -> const LHS &
					{
						assert(M>=4);
						int g=__lg(M);
						while (lhs_pool.size()<=g)
							get_equations();
						return lhs_pool[g];
					};
				auto get_rhs=[&](int M) -> const RHS &
					{
						assert(M>=4);
						int g=__lg(M);
						while (rhs_pool.size()<=g)
							get_equations();
						return rhs_pool[g];
					};
				phi_tmp[s]=multigrid(get_lhs,get_rhs,phi_n).get_result();
			}
		}
		auto phi_np1=phi_n;
		for (s=0; s<ns; s++)
			if (b[s])
			{
				equations q(u,f,boundary,nu,N,0);
				q.set_coef(phi_tmp[s]);
				for (x=0; x<N; x++)
					for (y=0; y<N; y++)
					{
						q.set_shift(x*N+y);
						q.add_rhs_diff(x,y,k*b[s],k*(T+c[s]));
						q.add_rhs_adv(x,y,k*b[s],k*(T+c[s]));
						q.add_rhs_f(x,y,k*b[s],k*(T+c[s]));
					}
				phi_np1+=q.get_rhs();
			}
		phi.push_back(move(phi_np1));
	}

public:
	FVMOL_1(const point &_u,const func_f &_f,const func_boundary &_boundary,const valarray<db> &initial_value,
		int _step,int _N,db _nu,db _k)
		: u(_u),f(_f),boundary(_boundary),phi{initial_value},step(_step),N(_N),n(_N *_N),nu(_nu),h(1.0/_N),
		k(_k)
	{
		for (int i=0; i<step; i++)
			iteration();
	}
	db cell(int x,int y,int t) const
	{
		return phi[t][x*N+y];
	}
};
template <class func_f,class func_boundary>
const db FVMOL_1<func_f,func_boundary>::gamma=0.25;
template <class func_f,class func_boundary>
const valarray<db> FVMOL_1<func_f,func_boundary>::b={
	0.15791629516167136,0,0.18675894052400077,0.6805652953093346,-0.27524053099500667,gamma};
template <class func_f,class func_boundary>
const valarray<db> FVMOL_1<func_f,func_boundary>::c={0,0.5,0.332,0.62,0.85,1};
template <class func_f,class func_boundary>
const vector<valarray<db>> FVMOL_1<func_f,func_boundary>::AE={{0,0,0,0,0,0},
	{2*gamma,0,0,0,0,0},
	{0.221776,0.110224,0,0,0,0},
	{-0.04884659515311857,-0.17772065232640102,0.8465672474795197,0,0,0},
	{-0.15541685842491548,-0.3567050098221991,1.0587258798684427,0.30339598837867193,0,0},
	{0.2014243506726763,0.008742057842904185,0.15993995707168115,0.4038290605220775,0.22606457389066084,0}};
template <class func_f,class func_boundary>
const vector<valarray<db>> FVMOL_1<func_f,func_boundary>::AI={{0,0,0,0,0,0},
	{gamma,gamma,0,0,0,0},
	{0.137776,-0.055776,gamma,0,0,0},
	{0.14463686602698217,-0.22393190761334475,0.4492950415863626,gamma,0,0},
	{0.09825878328356477,-0.5915442428196704,0.8101210538282996,0.283164405707806,gamma,0},
	b};
template <class func_f,class func_boundary>
const int FVMOL_1<func_f,func_boundary>::ns=6;
template <class func_f,class func_boundary>
class FVMOL_1_3d  // fourth-order FV-MOL method
{
	const static db gamma;
	const static valarray<db> b,c;
	const static vector<valarray<db>> AE,AI;
	const static int ns;
	func_f f;
	func_boundary boundary;
	int N,n,step;
	db nu,h,k;
	point_3d u;
	vector<valarray<db>> phi;
	class equations
	{
		func_f f;
		point_3d u;
		func_boundary boundary;
		int N,n,shift;
		db h,nu;
		sparse lhs;
		valarray<db> coef;
		valarray<db> rhs;

	public:
		equations(const point_3d &_u,const func_f &_f,const func_boundary &_boundary,db _nu,int _N,bool flg=1)
			: u(_u),f(_f),boundary(_boundary),N(_N),n(N *N *N),h(1.0/N),nu(_nu),lhs(flg?n:0),
			coef(flg?0:n),rhs(n)
		{
			if (flg)
			{
				for (int i=0; i<n; i++)
					add(lhs,i,i,-1);
			}
		};
		void add_coef(int x,int y,int z,db d,db t)
		{
// cerr<<N<<" * "<<N<<" = "<<n<<", "<<"("<<x<<", "<<y<<", "<<z<<")\n";
			if (x>=0&&x<N&&y>=0&&y<N&&z>=0&&z<N)
			{
				add(lhs,shift,(x*N+y)*N+z,d);
				return;
			}
			assert((x>=0&&x<N)+(y>=0&&y<N)+(z>=0&&z<N)==2);
			assert((x>=0-2&&x<N+2)&&(y>=0-2&&y<N+2)&&(z>=0-2&&z<N+2));
			const static vector<valarray<db>> coef_D={
				{4,-13.0/3,5.0/3,-1.0/3},{16,-70.0/3,32.0/3,-7.0/3}};
			const static vector<valarray<db>> coef_N={
				{6.0/5,5.0/10,9.0/10,-5.0/10,1.0/10},{6,-75.0/10,145.0/10,-75.0/10,15.0/10}};
			int dx=1,dy=1,dz=1,typ,i;
			if (x<0||x>=N)
			{
				dy=dz=0;
				if (x<0)
				{
					typ=-x-1;
					dx=-dx;
					x=0;
				}
				else
				{
					typ=x-N;
					x=N-1;
				}
			}
			else if (y<0||y>=N)
			{
				dx=dz=0;
				if (y<0)
				{
					typ=-y-1;
					dy=-dy;
					y=0;
				}
				else
				{
					typ=y-N;
					y=N-1;
				}
			}
			else
			{
				dx=dy=0;
				if (z<0)
				{
					typ=-z-1;
					dz=-dz;
					z=0;
				}
				else
				{
					typ=z-N;
					z=N-1;
				}
			}
			const auto &current=(dx?boundary(!!x,(y+0.5)*h,(z+0.5)*h,t)
				:(dy?boundary((x+0.5)*h,!!y,(z+0.5)*h,t)
					:boundary((x+0.5)*h,(y+0.5)*h,!!z,t)));
			const auto &coef=(current.index()==0?coef_D:coef_N)[typ];
			db value=(current.index()==0?get<0>(current).value():get<1>(current).value()*h);
			for (i=1; i<coef.size(); i++)
				add_coef(x-(i-1)*dx,y-(i-1)*dy,z-(i-1)*dz,coef[i]*d,t);
			rhs[shift]+=coef[0]*d*value;
		}
		void add_rhs(int x,int y,int z,db d,db t)
		{
			if (x>=0&&x<N&&y>=0&&y<N&&z>=0&&z<N)
			{
				rhs[shift]+=d*coef[(x*N+y)*N+z];
				return;
			}
			assert((x>=0&&x<N)+(y>=0&&y<N)+(z>=0&&z<N)==2);
			assert((x>=0-2&&x<N+2)&&(y>=0-2&&y<N+2)&&(z>=0-2&&z<N+2));
			const static vector<valarray<db>> coef_D={
				{4,-13.0/3,5.0/3,-1.0/3},{16,-70.0/3,32.0/3,-7.0/3}};
			const static vector<valarray<db>> coef_N={
				{6.0/5,5.0/10,9.0/10,-5.0/10,1.0/10},{6,-75.0/10,145.0/10,-75.0/10,15.0/10}};
			int dx=1,dy=1,dz=1,typ,i;
			if (x<0||x>=N)
			{
				dy=dz=0;
				if (x<0)
				{
					typ=-x-1;
					dx=-dx;
					x=0;
				}
				else
				{
					typ=x-N;
					x=N-1;
				}
			}
			else if (y<0||y>=N)
			{
				dx=dz=0;
				if (y<0)
				{
					typ=-y-1;
					dy=-dy;
					y=0;
				}
				else
				{
					typ=y-N;
					y=N-1;
				}
			}
			else
			{
				dx=dy=0;
				if (z<0)
				{
					typ=-z-1;
					dz=-dz;
					z=0;
				}
				else
				{
					typ=z-N;
					z=N-1;
				}
			}
			const auto &current=(dx?boundary(!!x,(y+0.5)*h,(z+0.5)*h,t)
				:(dy?boundary((x+0.5)*h,!!y,(z+0.5)*h,t)
					:boundary((x+0.5)*h,(y+0.5)*h,!!z,t)));
			const auto &coef=(current.index()==0?coef_D:coef_N)[typ];
			db value=(current.index()==0?get<0>(current).value():get<1>(current).value()*h);
			for (i=1; i<coef.size(); i++)
				add_rhs(x-(i-1)*dx,y-(i-1)*dy,z-(i-1)*dz,coef[i]*d,t);
			rhs[shift]+=coef[0]*d*value;
		}
		void add_coef_dx(int x,int y,int z,db d,db t)
		{
			d/=12*h;
			add_coef(x+1,y,z,15*d,t);
			add_coef(x,y,z,-15*d,t);
			add_coef(x+2,y,z,-d,t);
			add_coef(x-1,y,z,d,t);
		}
		void add_rhs_dx(int x,int y,int z,db d,db t)
		{
			d/=12*h;
			add_rhs(x+1,y,z,15*d,t);
			add_rhs(x,y,z,-15*d,t);
			add_rhs(x+2,y,z,-d,t);
			add_rhs(x-1,y,z,d,t);
		}
		void add_coef_dy(int x,int y,int z,db d,db t)
		{
			d/=12*h;
			add_coef(x,y+1,z,15*d,t);
			add_coef(x,y,z,-15*d,t);
			add_coef(x,y+2,z,-d,t);
			add_coef(x,y-1,z,d,t);
		}
		void add_rhs_dy(int x,int y,int z,db d,db t)
		{
			d/=12*h;
			add_rhs(x,y+1,z,15*d,t);
			add_rhs(x,y,z,-15*d,t);
			add_rhs(x,y+2,z,-d,t);
			add_rhs(x,y-1,z,d,t);
		}
		void add_coef_dz(int x,int y,int z,db d,db t)
		{
			d/=12*h;
			add_coef(x,y,z+1,15*d,t);
			add_coef(x,y,z,-15*d,t);
			add_coef(x,y,z+2,-d,t);
			add_coef(x,y,z-1,d,t);
		}
		void add_rhs_dz(int x,int y,int z,db d,db t)
		{
			d/=12*h;
			add_rhs(x,y,z+1,15*d,t);
			add_rhs(x,y,z,-15*d,t);
			add_rhs(x,y,z+2,-d,t);
			add_rhs(x,y,z-1,d,t);
		}
		void add_coef_diff(int x,int y,int z,db d,db t)
		{
			d*=nu/h;
			add_coef_dx(x,y,z,d,t);
			add_coef_dx(x-1,y,z,-d,t);
			add_coef_dy(x,y,z,d,t);
			add_coef_dy(x,y-1,z,-d,t);
			add_coef_dz(x,y,z,d,t);
			add_coef_dz(x,y,z-1,-d,t);
		}
		void add_rhs_diff(int x,int y,int z,db d,db t)
		{
			d*=nu/h;
			add_rhs_dx(x,y,z,d,t);
			add_rhs_dx(x-1,y,z,-d,t);
			add_rhs_dy(x,y,z,d,t);
			add_rhs_dy(x,y-1,z,-d,t);
			add_rhs_dz(x,y,z,d,t);
			add_rhs_dz(x,y,z-1,-d,t);
		}
		void add_rhs_x(int x,int y,int z,db d,db t)
		{
			d/=12;
			add_rhs(x,y,z,d*7,t);
			add_rhs(x+1,y,z,d*7,t);
			add_rhs(x-1,y,z,-d,t);
			add_rhs(x+2,y,z,-d,t);
		}
		void add_rhs_y(int x,int y,int z,db d,db t)
		{
			d/=12;
			add_rhs(x,y,z,d*7,t);
			add_rhs(x,y+1,z,d*7,t);
			add_rhs(x,y-1,z,-d,t);
			add_rhs(x,y+2,z,-d,t);
		}
		void add_rhs_z(int x,int y,int z,db d,db t)
		{
			d/=12;
			add_rhs(x,y,z,d*7,t);
			add_rhs(x,y,z+1,d*7,t);
			add_rhs(x,y,z-1,-d,t);
			add_rhs(x,y,z+2,-d,t);
		}
		void add_rhs_adv(int x,int y,int z,db d,db t)
		{
			d/=h;
			add_rhs_x(x,y,z,-d*u.x,t);
			add_rhs_x(x-1,y,z,d*u.x,t);
			add_rhs_y(x,y,z,-d*u.y,t);
			add_rhs_y(x,y-1,z,d*u.y,t);
			add_rhs_z(x,y,z,-d*u.z,t);
			add_rhs_z(x,y,z-1,d*u.z,t);
		}
		void add_rhs_f(int x,int y,int z,db d,db t)
		{
			rhs[shift]+=f((x+.5)*h,(y+.5)*h,(z+.5)*h,t)*d;
		}
		void set_coef(const valarray<db> &_coef)
		{
			coef=_coef;
		}
		void set_shift(int x)
		{
			shift=x;
		}
		void set_rhs(const valarray<db> &_rhs)
		{
			rhs=_rhs;
		}
		sparse &get_lhs()
		{
			return lhs;
		}
		valarray<db> &get_rhs()
		{
			return rhs;
		}
	};
	void iteration()
	{
		int T=phi.size()-1;
		cerr<<"Iteration "<<T<<endl;
		const auto &phi_n=phi.back();
		int i,j,x,y,z,s;
		vector phi_tmp(ns,valarray<db>(n));
		for (s=0; s<ns; s++)
		{
			if (s==0)
				phi_tmp[0]=phi_n;
			else
			{
				using LHS=sparse;
				using RHS=valarray<db>;
				vector<LHS> lhs_pool(2);
				vector<RHS> rhs_pool(2);
				auto get_equations=[&]()
					{
						int N=1<<lhs_pool.size(),x,y,z,j;
						equations t(u,f,boundary,nu,N);
						t.set_rhs(reset_size_3d(phi_n,N));
						for (x=0; x<N; x++)
							for (y=0; y<N; y++)
								for (z=0; z<N; z++)
								{
									t.set_shift((x*N+y)*N+z);
									t.add_coef_diff(x,y,z,k*gamma,k*(T+c[s]));
								}
						for (j=0; j<s; j++)
						{
							equations q(u,f,boundary,nu,N,0);
							q.set_coef(reset_size_3d(phi_tmp[j],N));
							for (x=0; x<N; x++)
								for (y=0; y<N; y++)
									for (z=0; z<N; z++)
									{
										q.set_shift((x*N+y)*N+z);
										q.add_rhs_diff(x,y,z,k*AI[s][j],k*(T+c[j]));
										q.add_rhs_adv(x,y,z,k*AE[s][j],k*(T+c[j]));
										q.add_rhs_f(x,y,z,k*AE[s][j],k*(T+c[j]));
									}
							t.get_rhs()+=q.get_rhs();
						}
						lhs_pool.push_back(move(t.get_lhs()));
						rhs_pool.push_back(-t.get_rhs());
						unique_matrix(lhs_pool.back());
					};
				auto get_lhs=[&](int M) -> const LHS &
					{
						assert(M>=4);
						int g=__lg(M);
						while (lhs_pool.size()<=g)
							get_equations();
						return lhs_pool[g];
					};
				auto get_rhs=[&](int M) -> const RHS &
					{
						assert(M>=4);
						int g=__lg(M);
						while (rhs_pool.size()<=g)
							get_equations();
						return rhs_pool[g];
					};
				phi_tmp[s]=multigrid_3d(get_lhs,get_rhs,phi_n).get_result();
			}
		}
		auto phi_np1=phi_n;
		for (s=0; s<ns; s++)
			if (b[s])
			{
				equations q(u,f,boundary,nu,N,0);
				q.set_coef(phi_tmp[s]);
				for (x=0; x<N; x++)
					for (y=0; y<N; y++)
						for (z=0; z<N; z++)
							if (b[s])
							{
								q.set_shift((x*N+y)*N+z);
								q.add_rhs_diff(x,y,z,k*b[s],k*(T+c[s]));
								q.add_rhs_adv(x,y,z,k*b[s],k*(T+c[s]));
								q.add_rhs_f(x,y,z,k*b[s],k*(T+c[s]));
							}
				phi_np1+=q.get_rhs();
			}
		phi.push_back(move(phi_np1));
	}

public:
	FVMOL_1_3d(const point_3d &_u,const func_f &_f,const func_boundary &_boundary,const valarray<db> &initial_value,
		int _step,int _N,db _nu,db _k)
		: u(_u),f(_f),boundary(_boundary),phi{initial_value},step(_step),N(_N),n(_N *_N *_N),nu(_nu),
		h(1.0/_N),k(_k)
	{
		for (int i=0; i<step; i++)
			iteration();
	}
	db cell(int x,int y,int z,int t) const
	{
		return phi[t][(x*N+y)*N+z];
	}
};
template <class func_f,class func_boundary>
const db FVMOL_1_3d<func_f,func_boundary>::gamma=0.25;
template <class func_f,class func_boundary>
const valarray<db> FVMOL_1_3d<func_f,func_boundary>::b={
	0.15791629516167136,0,0.18675894052400077,0.6805652953093346,-0.27524053099500667,gamma};
template <class func_f,class func_boundary>
const valarray<db> FVMOL_1_3d<func_f,func_boundary>::c={0,0.5,0.332,0.62,0.85,1};
template <class func_f,class func_boundary>
const vector<valarray<db>> FVMOL_1_3d<func_f,func_boundary>::AE={{0,0,0,0,0,0},
	{2*gamma,0,0,0,0,0},
	{0.221776,0.110224,0,0,0,0},
	{-0.04884659515311857,-0.17772065232640102,0.8465672474795197,0,0,0},
	{-0.15541685842491548,-0.3567050098221991,1.0587258798684427,0.30339598837867193,0,0},
	{0.2014243506726763,0.008742057842904185,0.15993995707168115,0.4038290605220775,0.22606457389066084,0}};
template <class func_f,class func_boundary>
const vector<valarray<db>> FVMOL_1_3d<func_f,func_boundary>::AI={{0,0,0,0,0,0},
	{gamma,gamma,0,0,0,0},
	{0.137776,-0.055776,gamma,0,0,0},
	{0.14463686602698217,-0.22393190761334475,0.4492950415863626,gamma,0,0},
	{0.09825878328356477,-0.5915442428196704,0.8101210538282996,0.283164405707806,gamma,0},
	b};
template <class func_f,class func_boundary>
const int FVMOL_1_3d<func_f,func_boundary>::ns=6;
class FVMOL_2  // fourth-order FV-MOL method
{
	const static db gamma;
	const static valarray<db> b,c;
	const static vector<valarray<db>> AE,AI;
	const static int ns;
	int N,n,step;
	db nu,h,k;
	vector<valarray<db>> phi;
	vector<vector<valarray<db>>> ux,uy;
	class equations
	{
		const vector<valarray<db>> &ux,&uy;
		int N,n,shift;
		db h,nu;
		sparse lhs;
		valarray<db> coef;
		valarray<db> rhs;

	public:
		equations(const vector<valarray<db>> &_ux,const vector<valarray<db>> &_uy,db _nu,int _N,bool flg=1)
			: ux(_ux),uy(_uy),N(_N),n(N *N),h(1.0/N),nu(_nu),lhs(flg?n:0),coef(flg?0:n),rhs(n)
		{
			if (flg)
			{
				for (int i=0; i<n; i++)
					add(lhs,i,i,-1);
			}
		};
		void add_coef(int x,int y,db d)
		{
			if (x>=N)
				x-=N;
			else if (x<0)
				x+=N;
			if (y>=N)
				y-=N;
			else if (y<0)
				y+=N;
			assert(x>=0&&x<N&&y>=0&&y<N);
			add(lhs,shift,x*N+y,d);
		}
		void add_rhs(int x,int y,db d)
		{
			if (x>=N)
				x-=N;
			else if (x<0)
				x+=N;
			if (y>=N)
				y-=N;
			else if (y<0)
				y+=N;
			assert(x>=0&&x<N&&y>=0&&y<N);
			rhs[shift]+=d*coef[x*N+y];
		}
		void add_coef_dx(int x,int y,db d)
		{
			d/=12*h;
			add_coef(x+1,y,15*d);
			add_coef(x,y,-15*d);
			add_coef(x+2,y,-d);
			add_coef(x-1,y,d);
		}
		void add_rhs_dx(int x,int y,db d)
		{
			d/=12*h;
			add_rhs(x+1,y,15*d);
			add_rhs(x,y,-15*d);
			add_rhs(x+2,y,-d);
			add_rhs(x-1,y,d);
		}
		void add_coef_dy(int x,int y,db d)
		{
			d/=12*h;
			add_coef(x,y+1,15*d);
			add_coef(x,y,-15*d);
			add_coef(x,y+2,-d);
			add_coef(x,y-1,d);
		}
		void add_rhs_dy(int x,int y,db d)
		{
			d/=12*h;
			add_rhs(x,y+1,15*d);
			add_rhs(x,y,-15*d);
			add_rhs(x,y+2,-d);
			add_rhs(x,y-1,d);
		}
		void add_coef_diff(int x,int y,db d)
		{
			d*=nu/h;
			add_coef_dx(x,y,d);
			add_coef_dx(x-1,y,-d);
			add_coef_dy(x,y,d);
			add_coef_dy(x,y-1,-d);
		}
		void add_rhs_diff(int x,int y,db d)
		{
			d*=nu/h;
			add_rhs_dx(x,y,d);
			add_rhs_dx(x-1,y,-d);
			add_rhs_dy(x,y,d);
			add_rhs_dy(x,y-1,-d);
		}
		void add_rhs_x(int x,int y,db d)
		{
			d/=12;
			add_rhs(x,y,d*7);
			add_rhs(x+1,y,d*7);
			add_rhs(x-1,y,-d);
			add_rhs(x+2,y,-d);
		}
		void add_rhs_y(int x,int y,db d)
		{
			d/=12;
			add_rhs(x,y,d*7);
			add_rhs(x,y+1,d*7);
			add_rhs(x,y-1,-d);
			add_rhs(x,y+2,-d);
		}
		void add_rhs_ux(int x,int y,db d)
		{
			if (x<=0)
				x+=N;
			if (y<=0)
				y+=N;
			add_rhs_x(x,y,d*(uy[x][y+1]-uy[x][y]));
			d*=((uy[x][y+2]-uy[x][y+1])-(uy[x][y]-uy[x][y-1]))/48;
			add_rhs_x(x,y+1,d);
			add_rhs_x(x,y-1,-d);
		}
		void add_rhs_uy(int x,int y,db d)
		{
			if (x<=0)
				x+=N;
			if (y<=0)
				y+=N;
			add_rhs_y(x,y,d*(ux[x+1][y]-ux[x][y]));
			d*=((ux[x+2][y]-ux[x+1][y])-(ux[x][y]-ux[x-1][y]))/48;
			add_rhs_y(x+1,y,d);
			add_rhs_y(x-1,y,-d);
		}
		void add_rhs_adv(int x,int y,db d)
		{
			d/=h;
			add_rhs_ux(x,y,-d);
			add_rhs_ux(x-1,y,d);
			add_rhs_uy(x,y,-d);
			add_rhs_uy(x,y-1,d);
		}
		void set_coef(const valarray<db> &_coef)
		{
			coef=_coef;
		}
		void set_shift(int x)
		{
			shift=x;
		}
		void set_rhs(const valarray<db> &_rhs)
		{
			rhs=_rhs;
		}
		sparse &get_lhs()
		{
			return lhs;
		}
		valarray<db> &get_rhs()
		{
			return rhs;
		}
	};
	vector<sparse> lhs_pool;
	void iteration()
	{
		int T=phi.size()-1;
		cerr<<"Iteration "<<T<<endl;
		const auto &phi_n=phi.back();
		int i,j,x,y,s;
		using LHS=sparse;
		auto get_equations_lhs=[&]()
			{
				int N=1<<lhs_pool.size(),x,y,j;
				equations t(ux[lhs_pool.size()],uy[lhs_pool.size()],nu,N);
				for (x=0; x<N; x++)
					for (y=0; y<N; y++)
					{
						t.set_shift(x*N+y);
						t.add_coef_diff(x,y,k*gamma);
					}
				lhs_pool.push_back(move(t.get_lhs()));
				unique_matrix(lhs_pool.back());
			};
		vector phi_tmp(ns,valarray<db>(n));
		for (s=0; s<ns; s++)
		{
			if (s==0)
				phi_tmp[0]=phi_n;
			else
			{
				int M=N;
				using RHS=valarray<db>;
				vector<RHS> rhs_pool(2);
				auto get_equations_rhs=[&]()
					{
						int N=1<<rhs_pool.size(),x,y,j;
						valarray t=-reset_size(phi_n,N);
						for (j=0; j<s; j++)
						{
							equations q(ux[rhs_pool.size()],uy[rhs_pool.size()],nu,N,0);
							q.set_coef(reset_size(phi_tmp[j],N));
							for (x=0; x<N; x++)
								for (y=0; y<N; y++)
								{
									q.set_shift(x*N+y);
									q.add_rhs_diff(x,y,k*AI[s][j]);
									q.add_rhs_adv(x,y,k*AE[s][j]);
								}
							t-=q.get_rhs();
						}
						rhs_pool.push_back(move(t));
					};
				auto get_lhs=[&](int M) -> const LHS &
					{
						assert(M>=4);
						int g=__lg(M);
						while (lhs_pool.size()<=g)
							get_equations_lhs();
						return lhs_pool[g];
					};
				auto get_rhs=[&](int M) -> const RHS &
					{
						assert(M>=4);
						int g=__lg(M);
						while (rhs_pool.size()<=g)
							get_equations_rhs();
						return rhs_pool[g];
					};
				phi_tmp[s]=multigrid(get_lhs,get_rhs,phi_n).get_result();
				// for (i = 0; i < n; i++)
				// assert(isfinite(phi_tmp[s][i]));
				// cerr << "max(phi_tmp[" << s << "]) = " << abs(phi_tmp[s]).max() << endl;
			}
		}
		auto phi_np1=phi_n;
		for (s=0; s<ns; s++)
			if (b[s])
			{
	 // assert(ux.size() == 1 + __lg(N));
	 // assert(uy.size() == 1 + __lg(N));
				equations q(ux.back(),uy.back(),nu,N,0);
				q.set_coef(phi_tmp[s]);
				for (x=0; x<N; x++)
					for (y=0; y<N; y++)
					{
						q.set_shift(x*N+y);
						q.add_rhs_diff(x,y,k*b[s]);
						q.add_rhs_adv(x,y,k*b[s]);
					}
				phi_np1+=q.get_rhs();
			}
		phi.push_back(move(phi_np1));
	}

public:
	template <class func_ux,class func_uy>
	FVMOL_2(const func_ux &_ux,const func_uy &_uy,const valarray<db> &initial_value,int _step,int _N,db _nu,db _k)
		: phi{initial_value},step(_step),N(_N),n(_N *_N),nu(_nu),h(1.0/_N),k(_k),ux(__lg(N)+1),
		uy(__lg(N)+1),lhs_pool(2)
	{
		int i,j;
		for (int M=1; M<=N;)
		{
			auto &cur_x=ux[__lg(M)];
			auto &cur_y=uy[__lg(M)];
			db h=1.0/M;
			M<<=1;
			cur_x.resize(M+1,valarray<db>(M+1));
			cur_y.resize(M+1,valarray<db>(M+1));
			for (i=0; i<=M; i++)
				for (j=0; j<=M; j++)
				{
					cur_x[i][j]=_ux(i*h,j*h)/h;
					cur_y[i][j]=_uy(i*h,j*h)/h;
					// assert(isfinite(_ux(i * h, j * h)));
					// assert(isfinite(_uy(i * h, j * h)));
				}
		}
		for (i=0; i<step; i++)
			iteration();
	}
	db cell(int x,int y,int t) const
	{
		return phi[t][x*N+y];
	}
};
const db FVMOL_2::gamma=0.25;
const valarray<db> FVMOL_2::b={
	0.15791629516167136,0,0.18675894052400077,0.6805652953093346,-0.27524053099500667,gamma};
const valarray<db> FVMOL_2::c={0,0.5,0.332,0.62,0.85,1};
const vector<valarray<db>> FVMOL_2::AE={{0,0,0,0,0,0},
	{2*gamma,0,0,0,0,0},
	{0.221776,0.110224,0,0,0,0},
	{-0.04884659515311857,-0.17772065232640102,0.8465672474795197,0,0,0},
	{-0.15541685842491548,-0.3567050098221991,1.0587258798684427,0.30339598837867193,0,0},
	{0.2014243506726763,0.008742057842904185,0.15993995707168115,0.4038290605220775,0.22606457389066084,0}};
const vector<valarray<db>> FVMOL_2::AI={{0,0,0,0,0,0},
	{gamma,gamma,0,0,0,0},
	{0.137776,-0.055776,gamma,0,0,0},
	{0.14463686602698217,-0.22393190761334475,0.4492950415863626,gamma,0,0},
	{0.09825878328356477,-0.5915442428196704,0.8101210538282996,0.283164405707806,gamma,0},
	b};
const int FVMOL_2::ns=6;
