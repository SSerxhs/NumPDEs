#include "bits/stdc++.h"
using db=double;
#include "heat.h"
#include "advection.h"
using namespace std;
template<typename method,typename func> void solve1(const string &filename,const func &boundary,db h,db r,db nu,int n,const vector<int> &vt)
{
    ofstream out(filename+".txt");
    db k=r*h*h/nu;
    int m=round(1/h);
    auto res=method().solve(boundary,m,n,nu,r);
    for (int i=0;i<=m;i++) out<<i*h<<" \n"[i==m];
    for (int t:vt)
    {
       for (int i=0;i<=m;i++) out<<res(i*h,k*t)<<" \n"[i==m];
    }
}
template<typename method,typename func> void solve2(const string &filename,const func &boundary,db h,db k,db a,db T)
{
    ofstream out(filename+".txt");
    db xl=15,xr=25,tl=0,tr=T;
    auto res=method().solve(boundary,a,xl,xr,tl,tr,h,k);
    int m=round((xr-xl)/h);
    for (int i=0;i<=m;i++) out<<xl+i*h<<" \n"[i==m];
    for (int i=0;i<=m;i++) out<<res(xl+i*h,T)<<" \n"[i==m];
}
int main()
{
    auto task_1=[&]()
    {
        auto boundary=[&](db x,db t) -> db
        {
            if (x==0||x==1)
            {
                return 0;
            }
            return 1-min<db>(1,20*abs(x-0.5));
        };
        db h=1.0/20;
        db nu=1;
        int n=10;
        vector<int> id={1,2,10};
        solve1<Crank_Nicolson>("CN(r=1)",boundary,h,1,nu,n,id);
        solve1<Crank_Nicolson>("CN(r=2)",boundary,h,2,nu,n,id);
        solve1<BTCS>("BTCS(r=1)",boundary,h,1,nu,n,id);
        solve1<BTCS>("BTCS(r=10)",boundary,h,10,nu,n,id);
        solve1<FTCS>("FTCS(r=1)",boundary,h,1,nu,n,id);
        solve1<FTCS>("FTCS(r=0.5)",boundary,h,0.5,nu,n,id);
        solve1<RKM2>("RKM2(r=1)",boundary,h,1,nu,n,id);
        solve1<RKM2>("RKM2(r=10)",boundary,h,10,nu,n,id);
        solve1<GLRKM1>("GLRKM1(r=1)",boundary,h,1,nu,n,id);
        solve1<GLRKM1>("GLRKM1(r=10)",boundary,h,10,nu,n,id);
    };
    auto task_2=[&]()
    {
        auto boundary=[&](db x) -> db
        {
            return exp(-20*(x-2)*(x-2))+exp(-(x-5)*(x-5));
        };
        db h=1.0/20;
        db a=1;
        db T=17;
        solve2<leapfrog>("leapfrog(k=0.8h)",boundary,h,0.8*h,a,T);
        solve2<LF>("LF(k=0.8h)",boundary,h,0.8*h,a,T);
        solve2<LW>("LW(k=0.8h)",boundary,h,0.8*h,a,T);
        solve2<upwind>("upwind(k=0.8h)",boundary,h,0.8*h,a,T);
        solve2<BW>("BW(k=0.8h)",boundary,h,0.8*h,a,T);
        solve2<leapfrog>("leapfrog(k=h)",boundary,h,h,a,T);
        solve2<LW>("LW(k=h)",boundary,h,h,a,T);
    };
    task_1();task_2();
}