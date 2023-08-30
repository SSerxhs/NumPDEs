#include "bits/stdc++.h"
#include "testlib/testlib.h"
#include "Gauss.h"
#include "geometry.h"
using namespace std;
using db=double;
using geometry::point,geometry::line,geometry::circle,geometry::convex,geometry::half_plane;
using geometry::sgn,geometry::eps,geometry::ll,geometry::segment;
using geometry::intersect,geometry::dis;
class solution
{
public:
	vector<vector<db>> u;
	vector<vector<char>> ban;
	bool in_domain(int x,int y) const { return !ban[x][y]; }
	void shift_to(db x,db y,db f)
	{
		db L=1./(u.size()-1);
		db d=f-u[(int)round(x/L)][(int)round(y/L)];
		for (auto &v:u) for (db &x:v) x+=d;
	}
	template<typename func1,typename func2> solution(const func1 &f,const func2 &v,int n) :u(n+1,vector<db>(n+1)),ban(n+1,vector<char>(n+1))//Lu=f,u=v or u=v' on boundary
	{
		if (n<1)
		{
			cerr<<"Too few points.\n";
			abort();
		}
		db h=1./n;
		++n;
		int m=n*n,i,j;
		vector A(m,vector<db>(m));
		vector<db> b(m);
		bool flg=0;
		ban[0][0]=1; ban[0][n-1]=1; ban[n-1][0]=1; ban[n-1][n-1]=1;
		A[0][0]=1; A[n-1][n-1]=1; A[m-n][m-n]=1; A[m-1][m-1]=1;
		for (i=0; i<n; i++) for (j=0; j<n; j++) if (!ban[i][j])
		{
			db x=1.*i/(n-1),y=1.*j/(n-1);
			auto &c=A[i*n+j];
			if (i&&j&&i+1<n&&j+1<n)
			{
				c[i*n+j]=4;
				c[(i-1)*n+j]=-1;
				c[(i+1)*n+j]=-1;
				c[i*n+j-1]=-1;
				c[i*n+j+1]=-1;
				b[i*n+j]=f(x,y)*h*h;
			}
			else
			{
				auto [type,value]=v(x,y);
				b[i*n+j]=value;
				if (type==1)//ux=v / uy=v
				{
					if (i&&i!=n-1)
					{
						if (j==0)
						{
							// c[i*n+j]-=1./h;
							// c[i*n+j+1]+=1./h;

							c[i*n+j]-=1.5/h;
							c[i*n+j+1]+=2.0/h;
							c[i*n+j+2]-=0.5/h;
						}
						else
						{
							// c[i*n+j]-=1./h;
							// c[i*n+j-1]+=1./h;

							c[i*n+j]-=1.5/h;
							c[i*n+j-1]+=2.0/h;
							c[i*n+j-2]-=0.5/h;
						}
					}
					else
					{
						if (i==0)
						{
							// c[i*n+j]-=1./h;
							// c[(i+1)*n+j]+=1./h;

							c[i*n+j]-=1.5/h;
							c[(i+1)*n+j]+=2.0/h;
							c[(i+2)*n+j]-=0.5/h;
						}
						else
						{
							// c[i*n+j]-=1./h;
							// c[(i-1)*n+j]+=1./h;

							c[i*n+j]-=1.5/h;
							c[(i-1)*n+j]+=2.0/h;
							c[(i-2)*n+j]-=0.5/h;
						}
					}
				}
				else c[i*n+j]=1,flg=1;//u=v
			}
		}
		if (flg==0)
		{
			int x=m-n-2;
			A.erase(A.begin()+1); b.erase(b.begin()+1);
			A.reserve(m+1); b.reserve(m+1);
			vector<db> c(m);
			c[1]=1; b.push_back(0);
			A.push_back(c);
		}
		/*for (i=0; i<m; i++)
		{
			for (j=0; j<m; j++) cout<<A[i][j]<<' ';
			cout<<"= "<<b[i]<<endl;
		}*/
		auto X=Gauss(A,b);
		for (i=0; i<m; i++)
		{
			db tmp=0;
			for (j=0; j<m; j++) tmp+=A[i][j]*X[j];
			assert(abs(tmp-b[i])<1e-4);
		}
		for (i=0; i<n; i++) for (j=0; j<n; j++) u[i][j]=X[i*n+j];
	};
	template<typename func1,typename func2> solution(const func1 &f,const func2 &v,int n,db Ox,db Oy,db R) :u(n+1,vector<db>(n+1)),ban(n+1,vector<char>(n+1))
	{
		if (n<1)
		{
			cerr<<"Too few points.\n";
			abort();
		}
		point<db> O{Ox,Oy};
		if (Ox<=R&&Ox+R>=1)
		{
			if ((dis(O,point<db>(0,0))>=R||dis(O,point<db>(1,0))>=R)&&(dis(O,point<db>(0,1))>=R||dis(O,point<db>(1,1))>=R))
			{
				cerr<<"Not connected\n";
				abort();
			}
		}
		if (Oy<=R&&Oy+R>=1)
		{
			if ((dis(O,point<db>(0,0))>=R||dis(O,point<db>(0,1))>=R)&&(dis(O,point<db>(1,0))>=R||dis(O,point<db>(1,1))>=R))
			{
				cerr<<"Not connected\n";
				abort();
			}
		}
		db h=1./n;
		int i,j,k;
		int cnt=0;
		ban[0][0]=1; ban[0][n]=1; ban[n][0]=1; ban[n][n]=1;
		for (i=0; i<=n; i++) for (j=0; j<=n; j++)
		{
			db x=i*h,y=j*h;
			if (sgn(R-dis(O,point{x,y}))>=0) ban[i][j]=1;
		}
		for (i=1; i<n; i++) for (j=1; j<n; j++) cnt+=!(ban[i][j]|ban[i][j-1]|ban[i-1][j]|ban[i+1][j]|ban[i][j+1]);
		if (cnt<4)
		{
			cerr<<"Too few equation-discretization points\n";
			abort();
		}
		++n;
		int m=n*n;
		vector id(m,vector<int>(m));
		vector<vector<db>> A;
		vector<db> b;
		int q=0;
		const db dx[4]={h,0,-h,0},dy[4]={0,h,0,-h};
		for (i=0; i<n; i++) for (j=0; j<n; j++) if (!ban[i][j]) id[i][j]=q++;
		for (i=0; i<n; i++) for (j=0; j<n; j++) if (!ban[i][j])
		{
			db x=i/(n-1.),y=j/(n-1.);
			db value=f(x,y);
			int status=0;
			for (k=0; k<4; k++) status|=(sgn(R-dis(segment(point{x,y},point{x+dx[k],y+dy[k]}),O))>0)<<k;
			assert(__builtin_popcount(status)<=2);
			vector<int> sid(4);
			vector<point<db>> pos(4);
			for (k=0; k<4; k++) if (status>>k&1) sid[k]=q++;
			if (i&&i+1<n&&j&&j+1<n)
			{
				A.emplace_back(q);
				auto &c=A.back();
				b.push_back(value);
				if ((status&5)==0)
				{
					c[id[i][j]]+=2/h/h;
					c[id[i+1][j]]-=1/h/h;
					c[id[i-1][j]]-=1/h/h;
				}
				else if (status&1)
				{
					assert((status&4)==0);
					pos[0]={Ox-sqrt(R*R-(Oy-y)*(Oy-y)),y};
					db theta=(pos[0].x-x)/h;
					c[id[i][j]]+=(1+theta)*2/(theta*(1+theta)*h*h);
					c[id[i-1][j]]-=2/((1+theta)*h*h);
					c[sid[0]]-=2/(theta*(1+theta)*h*h);
				}
				else
				{
					pos[2]={Ox+sqrt(R*R-(Oy-y)*(Oy-y)),y};
					db theta=(x-pos[2].x)/h;
					c[id[i][j]]+=(1+theta)*2/(theta*(1+theta)*h*h);
					c[id[i+1][j]]-=2/((1+theta)*h*h);
					c[sid[2]]-=2/(theta*(1+theta)*h*h);
				}
				if ((status&10)==0)
				{
					c[id[i][j]]+=2/h/h;
					c[id[i][j+1]]-=1/h/h;
					c[id[i][j-1]]-=1/h/h;
				}
				else if (status&2)
				{
					assert((status&8)==0);
					pos[1]={x,Oy-sqrt(R*R-(Ox-x)*(x-x))};
					db theta=(pos[1].y-y)/h;
					c[id[i][j]]+=(1+theta)*2/(theta*(1+theta)*h*h);
					c[id[i][j-1]]-=2/((1+theta)*h*h);
					c[sid[1]]-=2/(theta*(1+theta)*h*h);
				}
				else
				{
					pos[3]={x,Oy+sqrt(R*R-(Ox-x)*(x-x))};
					db theta=(y-pos[3].y)/h;
					c[id[i][j]]+=(1+theta)*2/(theta*(1+theta)*h*h);
					c[id[i][j+1]]-=2/((1+theta)*h*h);
					c[sid[3]]-=2/(theta*(1+theta)*h*h);
				}
			}
			else
			{
				auto [type,value]=v(x,y);
				b.push_back(value);
				A.emplace_back(q);
				auto &c=A.back();
				if (type==1)//ux=v / uy=v
				{
					if (i&&i!=n-1)
					{
						if (j==0)
						{
							// c[i*n+j]-=1./h;
							// c[i*n+j+1]+=1./h;

							c[id[i][j]]-=1.5/h;
							c[id[i][j+1]]+=2.0/h;
							c[id[i][j+2]]-=0.5/h;
						}
						else
						{
							// c[i*n+j]-=1./h;
							// c[i*n+j-1]+=1./h;

							c[id[i][j]]-=1.5/h;
							c[id[i][j-1]]+=2.0/h;
							c[id[i][j-2]]-=0.5/h;
						}
					}
					else
					{
						if (i==0)
						{
							// c[i*n+j]-=1./h;
							// c[(i+1)*n+j]+=1./h;

							c[id[i][j]]-=1.5/h;
							c[id[i+1][j]]+=2.0/h;
							c[id[i+2][j]]-=0.5/h;
						}
						else
						{
							// c[i*n+j]-=1./h;
							// c[(i-1)*n+j]+=1./h;

							c[id[i][j]]-=1.5/h;
							c[id[i-1][j]]+=2.0/h;
							c[id[i-2][j]]-=0.5/h;
						}
					}
				}
				else c[id[i][j]]=1;
			}
			for (k=0; k<4; k++) if (status>>k&1)
			{
				auto [type,value]=v(pos[k].x,pos[k].y);
				b.push_back(value);
				A.emplace_back(q);
				auto &c=A.back();
				c[sid[k]]=1;
			}
		}
		for (auto &v:A) v.resize(q);
		auto X=Gauss(A,b);
		for (i=0; i<n; i++) for (j=0; j<n; j++) if (!ban[i][j]) u[i][j]=X[id[i][j]];
	}
	db operator()(db x,db y) const
	{
		if (u.size()<2)
		{
			cerr<<"Wrong construction.\n";
			abort();
		}
		if (x<0||x>1||y<0||y>1)
		{
			cerr<<"Wrong coordinates.\n";
			abort();
		}
		db L=1./(u.size()-1);
		return u[(int)round(x/L)][(int)round(y/L)];
	}
};
template<typename func> db norm_1(const solution &f,const func &g,int n)
{
	int i,j;
	db res=0;
	for (i=1; i<n; i++) for (j=1; j<n; j++) if (f.in_domain(i,j)) res+=abs(f((db)i/n,(db)j/n)-g((db)i/n,(db)j/n));
	return res;
}
template<typename func> db norm_2(const solution &f,const func &g,int n)
{
	int i,j;
	db res=0,tmp;
	for (i=1; i<n; i++) for (j=1; j<n; j++) if (f.in_domain(i,j))
	{
		tmp=f((db)i/n,(db)j/n)-g((db)i/n,(db)j/n);
		res+=tmp*tmp;
	}
	return sqrt(res);
}
template<typename func> db norm_infty(const solution &f,const func &g,int n)
{
	int i,j;
	db res=0;
	for (i=1; i<n; i++) for (j=1; j<n; j++) if (f.in_domain(i,j)) res=max(res,abs(f((db)i/n,(db)j/n)-g((db)i/n,(db)j/n)));
	return res;
}
int main(int argc,char *argv[])
{
	auto F=[&](db x,db y) { return exp(y+sin(x)); };
	auto DxF=[&](db x,db y) { assert(x==0||x==1); return exp(y+sin(x))*cos(x); };
	auto DyF=[&](db x,db y) { assert(y==0||y==1); return exp(y+sin(x)); };
	auto n_D2F=[&](db x,db y) { x=sin(x); return exp(y+x)*(x+2)*(x-1); };
	// auto F=[&](db x,db y) { return x*exp(y); };
	// auto DxF=[&](db x,db y) { assert(x==0||x==1); return exp(y); };
	// auto DyF=[&](db x,db y) { assert(y==0||y==1); return x*exp(y); };
	// auto n_D2F=[&](db x,db y) { return x*exp(y); };
	cout<<fixed<<setprecision(5);
	auto Dirichlet=[&](int n)
	{
		solution res(n_D2F,[&](db x,db y) { return pair{0,F(x,y)}; },n);
		return tuple{norm_1(res,F,n),norm_2(res,F,n),norm_infty(res,F,n)};
	};
	auto Neumann=[&](int n)
	{
		solution res(n_D2F,[&](db x,db y)
			{
				assert((x==0||x==1)^(y==0||y==1));
				db r=x==0||x==1?DxF(x,y):DyF(x,y);
				if (x==1||y==1) return pair{1,-r};
				return pair{1,r};
			},n);
		res.shift_to(0.5,0.5,F(0.5,0.5));
		return tuple{norm_1(res,F,n),norm_2(res,F,n),norm_infty(res,F,n)};
	};
	auto Mixed=[&](int n)
	{
		solution res(n_D2F,[&](db x,db y)
			{
				if (x==0||y==1)
				{
					assert((x==0||x==1)^(y==0||y==1));
					db r=x==0||x==1?DxF(x,y):DyF(x,y);
					if (x==1||y==1) return pair{1,-r};
					return pair{1,r};
				}
				return pair{0,F(x,y)};
			},n);
		// res.shift_to(0.5,0.5,F(0.5,0.5));
		return tuple{norm_1(res,F,n),norm_2(res,F,n),norm_infty(res,F,n)};
	};
	auto Dirichlet_Special=[&](int n)
	{
		solution res(n_D2F,[&](db x,db y) { return pair{0,F(x,y)}; },n,0.5,0.5,0.2);
		return tuple{norm_1(res,F,n),norm_2(res,F,n),norm_infty(res,F,n)};
	};
	auto Mixed_Special=[&](int n)
	{
		solution res(n_D2F,[&](db x,db y)
			{
				if (x==0||x==1||y==0||y==1)
				{
					assert((x==0||x==1)^(y==0||y==1));
					db r=x==0||x==1?DxF(x,y):DyF(x,y);
					if (x==1||y==1) return pair{1,-r};
					return pair{1,r};
				}
				return pair{0,F(x,y)};
			},n,0.5,0.5,0.2);
		// res.shift_to(0.5,0.5,F(0.5,0.5));
		return tuple{norm_1(res,F,n),norm_2(res,F,n),norm_infty(res,F,n)};
	};
	registerGen(argc,argv,1);
	string type=opt<string>("type");
	vector<tuple<db,db,db>> res;
	if (type=="All")
	{
		for (string type:{"Dirichlet"s,"Neumann"s,"Mixed"s,"Dirichlet_Special"s,"Mixed_Special"s})
		{
			if (type=="Dirichlet"s) for (int n:{8,16,32,64}) res.push_back(Dirichlet(n));
			else if (type=="Neumann"s) for (int n:{8,16,32,64}) res.push_back(Neumann(n));
			else if (type=="Mixed"s) for (int n:{8,16,32,64}) res.push_back(Mixed(n));
			else if (type=="Dirichlet_Special"s) for (int n:{8,16,32,64}) res.push_back(Dirichlet_Special(n));
			else if (type=="Mixed_Special"s) for (int n:{8,16,32,64}) res.push_back(Mixed_Special(n));
			else cerr<<"Wrong test\n",abort();
			ofstream out(type+".out");
			// ofstream out(type+"2.out");
			for (int i=0; i<res.size(); i++) out<<get<0>(res[i])<<" \n"[i+1==res.size()];
			for (int i=0; i<res.size(); i++) out<<get<1>(res[i])<<" \n"[i+1==res.size()];
			for (int i=0; i<res.size(); i++) out<<get<2>(res[i])<<" \n"[i+1==res.size()];
			for (int i=1; i<res.size(); i++) out<<log(get<0>(res[i])/get<0>(res[i-1]))/log(0.5)<<" \n"[i+1==res.size()];
			for (int i=1; i<res.size(); i++) out<<log(get<1>(res[i])/get<1>(res[i-1]))/log(0.5)<<" \n"[i+1==res.size()];
			for (int i=1; i<res.size(); i++) out<<log(get<2>(res[i])/get<2>(res[i-1]))/log(0.5)<<" \n"[i+1==res.size()];
		}
		return 0;
	}
	if (type=="Dirichlet"s) for (int n:{8,16,32,64}) res.push_back(Dirichlet(n));
	else if (type=="Neumann"s) for (int n:{8,16,32,64}) res.push_back(Neumann(n));
	else if (type=="Mixed"s) for (int n:{8,16,32,64}) res.push_back(Mixed(n));
	else if (type=="Dirichlet_Special"s) for (int n:{8,16,32,64}) res.push_back(Dirichlet_Special(n));
	else if (type=="Mixed_Special"s) for (int n:{8,16,32,64}) res.push_back(Mixed_Special(n));
	else cerr<<"Wrong test\n",abort();
	ofstream out(type+".out");
	// ofstream out(type+"2.out");
	for (int i=0; i<res.size(); i++) out<<get<0>(res[i])<<" \n"[i+1==res.size()];
	for (int i=0; i<res.size(); i++) out<<get<1>(res[i])<<" \n"[i+1==res.size()];
	for (int i=0; i<res.size(); i++) out<<get<2>(res[i])<<" \n"[i+1==res.size()];
	for (int i=1; i<res.size(); i++) out<<log(get<0>(res[i])/get<0>(res[i-1]))/log(0.5)<<" \n"[i+1==res.size()];
	for (int i=1; i<res.size(); i++) out<<log(get<1>(res[i])/get<1>(res[i-1]))/log(0.5)<<" \n"[i+1==res.size()];
	for (int i=1; i<res.size(); i++) out<<log(get<2>(res[i])/get<2>(res[i-1]))/log(0.5)<<" \n"[i+1==res.size()];
}