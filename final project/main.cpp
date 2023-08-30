#pragma GCC target("popcnt", "sse3", "sse2", "sse", "avx", "sse4", "sse4.1", "sse4.2", "ssse3", "f16c", "avx2")
#pragma GCC optimize("inline", "fast-math", "unroll-loops", "no-stack-protector", "Ofast")
#define NDEBUG
#include "bits/stdc++.h"
#include "FV.h"
using namespace std;
db sqr(db x)
{
	return x*x;
}
const db pi=acos(-1);
void task_1()
{
	cout<<"1.1 Traveling sinusoidal waves"<<endl;
	const static point u{1,0.5},kp{2*pi,4*pi};
	const static db nu=0.01;
	const static db te=1,Cr=1;
	for (int N:{64,128,256,512})
	{
		time_t start_time=time(0);
		int n=N*N,i,j,ti;
		db x,y,t;
		db h=1.0/N;
		db k=Cr*h/1.5;
		auto f=[&](db x,db y,db t) -> db
			{
				point s={sin(kp.x*x-u.x*t),sin(kp.y*y-u.y*t)};
				return nu*kp.norm_sqr()*s.x*s.y+u.x*(kp.x-1)*cos(kp.x*x-u.x*t)*s.y+
					u.y*(kp.y-1)*cos(kp.y*y-u.y*t)*s.x;
			};
		auto integral_xy_f=[&](db x,db y,db t) -> db
			{
				point c={cos(kp.x*x-u.x*t),cos(kp.y*y-u.y*t)};
				// point s = {sqrt(1 - c.x * c.x)*(), sqrt(1 - c.y * c.y)};
				point s={sin(kp.x*x-u.x*t),sin(kp.y*y-u.y*t)};
				return (nu*kp.norm_sqr()*c.x*c.y-u.x*(kp.x-1)*s.x*c.y-u.y*(kp.y-1)*s.y*c.x)/
					(kp.x*kp.y);
			};
		auto Dxx_f=[&](db x,db y,db t) -> db
			{
				point s={-sqr(kp.x)*sin(kp.x*x-u.x*t),sin(kp.y*y-u.y*t)};
				return nu*kp.norm_sqr()*s.x*s.y-sqr(kp.y)*u.x*(kp.x-1)*cos(kp.x*x-u.x*t)*s.y+
					u.y*(kp.y-1)*cos(kp.y*y-u.y*t)*s.x;
			};
		auto Dyy_f=[&](db x,db y,db t) -> db
			{
				point s={sin(kp.x*x-u.x*t),-sqr(kp.y)*sin(kp.y*y-u.y*t)};
				return nu*kp.norm_sqr()*s.x*s.y+u.x*(kp.x-1)*cos(kp.x*x-u.x*t)*s.y-
					sqr(kp.x)*u.y*(kp.y-1)*cos(kp.y*y-u.y*t)*s.x;
			};
		auto average_f=[&](db x,db y,db t) -> db
			{
				return ((integral_xy_f(x+0.5*h,y+0.5*h,t)+integral_xy_f(x-0.5*h,y-0.5*h,t))-
					(integral_xy_f(x-0.5*h,y+0.5*h,t)+integral_xy_f(x+0.5*h,y-0.5*h,t)))/
					(h*h);
			 // return f(x, y, t) + h * h * (Dxx_f(x, y, t) + Dyy_f(x, y, t)) / 24;
			};
		auto integral_xy_phi=[&](db x,db y,db t) -> db
			{
				return cos(kp.x*x-u.x*t)*cos(kp.y*y-u.y*t)/(kp.x*kp.y);
			};
		auto integral_x_phi=[&](db x,db y,db t) -> db
			{
				return -cos(kp.x*x-u.x*t)*sin(kp.y*y-u.y*t)/(kp.x);
			};
		auto integral_y_phi=[&](db x,db y,db t) -> db
			{
				return -sin(kp.x*x-u.x*t)*cos(kp.y*y-u.y*t)/(kp.y);
			};
		auto integral_x_der_y_phi=[&](db x,db y,db t) -> db
			{
				return -cos(kp.x*x-u.x*t)*cos(kp.y*y-u.y*t)*(kp.y/kp.x);
			};
		auto integral_y_der_x_phi=[&](db x,db y,db t) -> db
			{
				return -cos(kp.x*x-u.x*t)*cos(kp.y*y-u.y*t)*(kp.x/kp.y);
			};
		auto boundary=[&](db x,db y,db t) -> variant<Dirichlet,Neumann>
			{
				if (y==0)
					return Dirichlet((integral_x_phi(x+h/2,y,t)-integral_x_phi(x-h/2,y,t))/h);
				else if (y==1)
					return Neumann((integral_x_der_y_phi(x+h/2,y,t)-integral_x_der_y_phi(x-h/2,y,t))/h);
				else if (x==0)
					return Dirichlet((integral_y_phi(x,y+h/2,t)-integral_y_phi(x,y-h/2,t))/h);
				else if (x==1)
				{
					return Neumann((integral_y_der_x_phi(x,y+h/2,t)-integral_y_der_x_phi(x,y-h/2,t))/h);
				}
				throw;
			};
		int T=round(te/k);
		valarray<db> initial_value(n),final_value(n);
		for (i=0; i<N; i++)
		{
			x=i*h;
			for (j=0; j<N; j++)
			{
				y=j*h;
				initial_value[i*N+j]=((integral_xy_phi(x,y,0)+integral_xy_phi(x+h,y+h,0))-
					(integral_xy_phi(x+h,y,0)+integral_xy_phi(x,y+h,0)))/
					(h*h);
				final_value[i*N+j]=((integral_xy_phi(x,y,te)+integral_xy_phi(x+h,y+h,te))-
					(integral_xy_phi(x+h,y,te)+integral_xy_phi(x,y+h,te)))/
					(h*h);
			}
		}
		auto res=FVMOL_1(u,average_f,boundary,initial_value,T,N,nu,k);
		db norm_2=0,norm_inf=0;
		for (i=0; i<N; i++)
			for (j=0; j<N; j++)
			{
				norm_2+=sqr(res.cell(i,j,T)-final_value[i*N+j]);
				norm_inf=max(norm_inf,abs(res.cell(i,j,T)-final_value[i*N+j]));
			}
		norm_2=sqrt(norm_2)/N;
		ofstream fout("1.1(N="s+to_string(N)+").txt");
		fout<<fixed<<setprecision(15);
		fout<<"h = 1/"<<N<<endl;
		fout<<"norm_2 = "<<norm_2<<endl;
		fout<<"norm_inf = "<<norm_inf<<endl;
		time_t end_time=time(0);
		fout<<"Takes "<<end_time-start_time<<" seconds"<<endl;
		cout<<"Takes "<<end_time-start_time<<" seconds"<<endl;
		for (i=0; i<N; i++)
			for (j=0; j<N; j++)
				fout<<res.cell(i,j,T)<<" \n"[j+1==N];
	}
}
void task_1_3d()
{
	cout<<"1.1.2 Traveling sinusoidal waves (3d)"<<endl;
	const static point_3d u{1,0.5,0.25},kp{2*pi,4*pi,6*pi};
	const static db nu=0.01;
	const static db te=1,Cr=1;
	for (int N:{16,32,64,128})
	{
		time_t start_time=time(0);
		int n=N*N*N,i,j,l,ti;
		db x,y,z,t;
		db h=1.0/N;
		db k=Cr*h/1.75;
		auto f=[&](db x,db y,db z,db t) -> db
			{
				point_3d c={cos(kp.x*x-u.x*t),cos(kp.y*y-u.y*t),cos(kp.z*z-u.z*t)};
				point_3d s={sin(kp.x*x-u.x*t),sin(kp.y*y-u.y*t),sin(kp.z*z-u.z*t)};
				return nu*kp.norm_sqr()*s.x*s.y*s.z+u.x*(kp.x-1)*c.x*s.y*s.z+
					u.y*(kp.y-1)*c.y*s.x*s.z+u.z*(kp.z-1)*c.z*s.x*s.y;
			};
		auto integral_xyz_f=[&](db x,db y,db z,db t) -> db
			{
				point_3d c={cos(kp.x*x-u.x*t),cos(kp.y*y-u.y*t),cos(kp.z*z-u.z*t)};
				point_3d s={sin(kp.x*x-u.x*t),sin(kp.y*y-u.y*t),sin(kp.z*z-u.z*t)};
				return (-nu*kp.norm_sqr()*c.x*c.y*c.z+u.x*(kp.x-1)*s.x*c.y*c.z+
					u.y*(kp.y-1)*s.y*c.x*c.z+u.z*(kp.z-1)*s.z*c.x*c.y)/
					(kp.x*kp.y*kp.z);
			};
		auto average_f=[&](db x,db y,db z,db t) -> db
			{
				return (integral_xyz_f(x+h,y+h,z+h,t)-integral_xyz_f(x+h,y+h,z,t)-
					integral_xyz_f(x+h,y,z+h,t)-integral_xyz_f(x,y+h,z+h,t)+
					integral_xyz_f(x+h,y,z,t)+integral_xyz_f(x,y+h,z,t)+
					integral_xyz_f(x,y,z+h,t)-integral_xyz_f(x,y,z,t))/
					(h*h*h);
			};
		auto integral_xyz_phi=[&](db x,db y,db z,db t) -> db
			{
				return -cos(kp.x*x-u.x*t)*cos(kp.y*y-u.y*t)*cos(kp.z*z-u.z*t)/(kp.x*kp.y*kp.z);
			};
		auto integral_xy_phi=[&](db x,db y,db z,db t) -> db
			{
				return cos(kp.x*x-u.x*t)*cos(kp.y*y-u.y*t)*sin(kp.z*z-u.z*t)/(kp.x*kp.y);
			};
		auto integral_xz_phi=[&](db x,db y,db z,db t) -> db
			{
				return cos(kp.x*x-u.x*t)*sin(kp.y*y-u.y*t)*cos(kp.z*z-u.z*t)/(kp.x*kp.z);
			};
		auto integral_yz_phi=[&](db x,db y,db z,db t) -> db
			{
				return sin(kp.x*x-u.x*t)*cos(kp.y*y-u.y*t)*cos(kp.z*z-u.z*t)/(kp.y*kp.z);
			};
		auto integral_xy_der_z_phi=[&](db x,db y,db z,db t) -> db
			{
				return cos(kp.x*x-u.x*t)*cos(kp.y*y-u.y*t)*cos(kp.z*z-u.z*t)*kp.z/(kp.x*kp.y);
			};
		auto integral_xz_der_y_phi=[&](db x,db y,db z,db t) -> db
			{
				return cos(kp.x*x-u.x*t)*cos(kp.y*y-u.y*t)*cos(kp.z*z-u.z*t)*kp.y/(kp.x*kp.z);
			};
		auto integral_yz_der_x_phi=[&](db x,db y,db z,db t) -> db
			{
				return cos(kp.x*x-u.x*t)*cos(kp.y*y-u.y*t)*cos(kp.z*z-u.z*t)*kp.x/(kp.y*kp.z);
			};
		auto boundary=[&](db x,db y,db z,db t) -> variant<Dirichlet,Neumann>
			{
				if (x==0)
					return Dirichlet(
						((integral_yz_phi(x,y+h/2,z+h/2,t)+integral_yz_phi(x,y-h/2,z-h/2,t))-
							(integral_yz_phi(x,y-h/2,z+h/2,t)+integral_yz_phi(x,y+h/2,z-h/2,t)))/
						(h*h));
				else if (x==1)
					return Neumann(((integral_yz_der_x_phi(x,y+h/2,z+h/2,t)+
						integral_yz_der_x_phi(x,y-h/2,z-h/2,t))-
						(integral_yz_der_x_phi(x,y-h/2,z+h/2,t)+
							integral_yz_der_x_phi(x,y+h/2,z-h/2,t)))/
						(h*h));
				else if (y==0)
					return Dirichlet(
						((integral_xz_phi(x+h/2,y,z+h/2,t)+integral_xz_phi(x-h/2,y,z-h/2,t))-
							(integral_xz_phi(x-h/2,y,z+h/2,t)+integral_xz_phi(x+h/2,y,z-h/2,t)))/
						(h*h));
				else if (y==1)
					return Neumann(((integral_xz_der_y_phi(x+h/2,y,z+h/2,t)+
						integral_xz_der_y_phi(x-h/2,y,z-h/2,t))-
						(integral_xz_der_y_phi(x-h/2,y,z+h/2,t)+
							integral_xz_der_y_phi(x+h/2,y,z-h/2,t)))/
						(h*h));
				else if (z==0)
					return Dirichlet(
						((integral_xy_phi(x+h/2,y+h/2,z,t)+integral_xy_phi(x-h/2,y-h/2,z,t))-
							(integral_xy_phi(x-h/2,y+h/2,z,t)+integral_xy_phi(x+h/2,y-h/2,z,t)))/
						(h*h));
				else if (z==1)
					return Neumann(((integral_xy_der_z_phi(x+h/2,y+h/2,z,t)+
						integral_xy_der_z_phi(x-h/2,y-h/2,z,t))-
						(integral_xy_der_z_phi(x-h/2,y+h/2,z,t)+
							integral_xy_der_z_phi(x+h/2,y-h/2,z,t)))/
						(h*h));
				else
					throw;
			};
		int T=round(te/k);
		valarray<db> initial_value(n),final_value(n);
		for (i=0; i<N; i++)
		{
			x=i*h;
			for (j=0; j<N; j++)
			{
				y=j*h;
				for (l=0; l<N; l++)
				{
					z=l*h;
					initial_value[(i*N+j)*N+l]=
						(integral_xyz_phi(x+h,y+h,z+h,0)-integral_xyz_phi(x+h,y+h,z,0)-
							integral_xyz_phi(x+h,y,z+h,0)-integral_xyz_phi(x,y+h,z+h,0)+
							integral_xyz_phi(x+h,y,z,0)+integral_xyz_phi(x,y+h,z,0)+
							integral_xyz_phi(x,y,z+h,0)-integral_xyz_phi(x,y,z,0))/
						(h*h*h);
					final_value[(i*N+j)*N+l]=
						(integral_xyz_phi(x+h,y+h,z+h,te)-integral_xyz_phi(x+h,y+h,z,te)-
							integral_xyz_phi(x+h,y,z+h,te)-integral_xyz_phi(x,y+h,z+h,te)+
							integral_xyz_phi(x+h,y,z,te)+integral_xyz_phi(x,y+h,z,te)+
							integral_xyz_phi(x,y,z+h,te)-integral_xyz_phi(x,y,z,te))/
						(h*h*h);
				}
			}
		}
		auto res=FVMOL_1_3d(u,average_f,boundary,initial_value,T,N,nu,k);
		db norm_2=0,norm_inf=0;
		for (i=0; i<N; i++)
			for (j=0; j<N; j++)
				for (l=0; l<N; l++)
				{
					norm_2+=sqr(res.cell(i,j,l,T)-final_value[(i*N+j)*N+l]);
					norm_inf=max(norm_inf,abs(res.cell(i,j,l,T)-final_value[(i*N+j)*N+l]));
				}
		norm_2=sqrt(norm_2/n);
		ofstream fout("1.1.2(N="s+to_string(N)+").txt");
		fout<<fixed<<setprecision(15);
		fout<<"h = 1/"<<N<<endl;
		fout<<"norm_2 = "<<norm_2<<endl;
		fout<<"norm_inf = "<<norm_inf<<endl;
		time_t end_time=time(0);
		fout<<"Takes "<<end_time-start_time<<" seconds"<<endl;
		cout<<"Takes "<<end_time-start_time<<" seconds"<<endl;
		for (i=0; i<N; i++)
		{
			for (j=0; j<N; j++)
				for (l=0; l<N; l++)
					fout<<res.cell(i,j,l,T)<<" \n"[l+1==N];
			fout<<"\n";
		}
	}
}
void task_2()
{
	cout<<"1.2 Gaussian patch in vortex shear"<<endl;
	const static db av=0.1;
	const static db Cr=1.0;
	const static db nu=1e-3;
	const static db te=1/av;
	auto u=[&](db x,db y) -> point
		{
			return av*point{sqr(sin(pi*x))*sin(2*pi*y),-sin(2*pi*x)*sqr(sin(pi*y))};
		};
	auto integral_x_u_y=[&](db x,db y) -> db { return av/(2*pi)*sqr(sin(pi*y))*cos(2*pi*x); };
	auto integral_y_u_x=[&](db x,db y) -> db { return -av/(2*pi)*sqr(sin(pi*x))*cos(2*pi*y); };
	// auto integral_x_u = [&](db x, db y) -> db { return av * (x / 2 - sin(2 * pi * x) / (4 * pi)) * sin(2 * pi * y);
	// }; auto integral_y_u = [&](db x, db y) -> db { return av * -sin(2 * pi * x) * (y / 2 - sin(2 * pi * y) / (4 *
	// pi)); };
	const static db factor=-1600*log(10);
	const static db sqrtfactor=sqrt(-factor);
	vector<vector<valarray<db>>> result;
	vector<int> tm;
	const static point O{0.5,0.75};
	for (int N:{64,128,256,512})
	// for (int N : {64, 128, 256})
	{
		time_t start_time=time(0);
		int i,j;
		int n=N*N;
		db h=1.0/N;
		db k=h*Cr/(2*av),x,y;
		auto phi=[&](db x,db y) -> db { return exp(factor*(sqr(x-O.x)+sqr(y-O.y))); };
		/*auto Dxx_phi = [&](db x, db y) -> db {
			return 2 * factor * (2 * factor * sqr(x -  O.x) + 1) * exp(factor * (sqr(x -  O.x) + sqr(y - O.y)));
		};
		auto Dyy_phi = [&](db x, db y) -> db {
			return 2 * factor * (2 * factor * sqr(y -  O.y) + 1) * exp(factor * (sqr(x -  O.x) + sqr(y -  O.y)));
		};*/
		// auto average_phi = [&](db x, db y) -> db { return phi(x, y) + h * h * (Dxx_phi(x, y) + Dyy_phi(x, y)) / 24;
		// };
		auto integral_phi=[&](db x) -> db { return sqrt(pi)/(2*sqrtfactor)*erfl(x*sqrtfactor); };
		auto average_phi=[&](db x,db y) -> db
			{
				return (integral_phi(x+h/2-O.x)-integral_phi(x-h/2-O.x))*
					(integral_phi(y+h/2-O.y)-integral_phi(y-h/2-O.y))/(h*h);
			};
		valarray<db> initial_value(n);
		for (i=0; i<N; i++)
		{
			x=(i+0.5)*h;
			for (j=0; j<N; j++)
			{
				y=(j+0.5)*h;
				initial_value[i*N+j]=average_phi(x,y);
				// assert(isfinite(average_phi(x, y)));
			}
		}
		int T=round(te/k);
		auto res=FVMOL_2(integral_x_u_y,integral_y_u_x,initial_value,T,N,nu,k);
		result.push_back(vector(N,valarray<db>(N)));
		for (i=0; i<N; i++)
			for (j=0; j<N; j++)
				result.back()[i][j]=res.cell(i,j,T);
		time_t end_time=time(0);
		tm.push_back(end_time-start_time);
		cout<<"Takes "<<end_time-start_time<<" seconds"<<endl;
		ofstream fout("1.2(N="s+to_string(N)+")tmp.txt");
		fout<<fixed<<setprecision(15);
		for (i=0; i<N; i++)
			for (j=0; j<N; j++)
				fout<<result.back()[i][j]<<" \n"[j+1==N];
	}
	vector<db> e2(result.size()-1),einf(result.size()-1);
	for (int id=0; id<result.size(); id++)
	{
		int N=1<<id+6,i,j,ii,jj;
		if (id+1<result.size())
		{
			db norm_2=0,norm_inf=0;
			auto &F=result[id],&G=result.back();
			int M=G.size()/N;
			for (i=0; i<N; i++)
				for (j=0; j<N; j++)
				{
					db E=0;
					int ri=(i+1)*M,rj=(j+1)*M;
					for (ii=i*M; ii<ri; ii++)
						for (jj=j*M; jj<rj; jj++)
							E+=G[ii][jj];
					E=E/(M*M)-F[i][j];
					norm_2+=sqr(E);
					norm_inf=max(norm_inf,abs(E));
				}
			norm_2=sqrt(norm_2)/N;
			e2[id]=norm_2;
			einf[id]=norm_inf;
		}
	}
	for (int id=0; id<result.size(); id++)
	{
		int N=1<<id+6,i,j,ii,jj;
		ofstream fout("1.2(N="s+to_string(N)+").txt");
		fout<<fixed<<setprecision(15);
		fout<<"h = 1/"<<N<<endl;
		if (id+1<result.size())
		{
			fout<<"norm_2 = "<<e2[id]<<endl;
			fout<<"norm_inf = "<<einf[id]<<endl;
			if (id+3<result.size())
			{
				db e=e2[id]/e2[id+1];
				fout<<"p (calculate by norm_2) = "<<log(sqrt(e*e+2*e-3)+e-1)/log(2)-1<<endl;
				e=einf[id]/einf[id+1];
				fout<<"p (calculate by norm_2) = "<<log(sqrt(e*e+2*e-3)+e-1)/log(2)-1<<endl;
			}
			else if (id+2<result.size())
			{
				db e=e2[id]/e2[id+1];
				fout<<"p (calculate by norm_2) = "<<log(e-1)/log(2)<<endl;
				e=einf[id]/einf[id+1];
				fout<<"p (calculate by norm_2) = "<<log(e-1)/log(2)<<endl;
			}
		}
		fout<<"Takes "<<tm[id]<<" seconds"<<endl;
		for (i=0; i<N; i++)
			for (j=0; j<N; j++)
				fout<<result[id][i][j]<<" \n"[j+1==N];
	}
}
void task_3()
{
	cout<<"1.3 Taylor vortex"<<endl;
	const static db te=0.5;
	for (db Cr:{0.75,1.5})
		for (int Re:{30,300,3000,30000})
			for (int N:{64,128,256,512})
			{
				db h=1.0/N;
				/*auto u = [&](db x, db y, db t) -> point {
					return point{1, 1} + 2 * exp(-8 * pi * pi * nu * t) * point {
						-cos(2 * pi * (x - t)) * sin(2 * pi * (y - t)), sin(2 *)
					}
				};
				auto p = [&](db x, db y, db t) -> db {
					return -exp(-16 * pi * pi * nu * t) * (cos(4 * pi * (x - t)) + cos(4 * pi * (y - t)));
				};*/
			}
}
int main()
{
	cout<<fixed<<setprecision(15);
	cerr<<fixed<<setprecision(15);
	task_1();
	task_1_3d();
	task_2();
	task_3();
}
