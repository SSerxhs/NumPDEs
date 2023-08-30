#pragma once
#include "solution.h"
#include "bits/stdc++.h"
using namespace std;
class advection
{
    vector<pair<int,db>> A0,A1,A2;
    int get_min() const
    {
        int mn=0;
        for (auto [x,y]:A0) mn=min(mn,x);
        for (auto [x,y]:A1) mn=min(mn,x);
        for (auto [x,y]:A2) mn=min(mn,x);
        return mn;
    }
    int get_max() const
    {
        int mx=0;
        for (auto [x,y]:A0) mx=max(mx,x);
        for (auto [x,y]:A1) mx=max(mx,x);
        for (auto [x,y]:A2) mx=max(mx,x);
        return mx;
    }
protected:
    advection(const vector<pair<int,db>> &A0,const vector<pair<int,db>> &A1,const vector<pair<int,db>> &A2):
        A0(A0),A1(A1),A2(A2){}
public:
    template<typename func> solution solve(const func &boundary,db a,db xl,db xr,db tl,db tr,db h,db k) const
    {
        xl+=ceil((tr-tl)/k)*h*(get_min()-1);
        xr+=ceil((tr-tl)/k)*h*(get_max()+1);
        int m=round((xr-xl)/h);
        int n=round((tr-tl)/k);
        db mu=a*k/h;
        vector y(n+1,valarray<db>(m+1)); 
        for (int i=0;i<=m;i++)
        {
            y[0][i]=boundary(xl+i*h);
        }
        ++m;
        for (int i=1;i<=n;i++)
        {
            const auto &y0=y[i-1];
            for (int j=0;j<m;j++)
            {
                db &tmp=y[i][j];
                db tmp1=0,tmp2=0;
                for (auto [x,k]:A0) tmp+=y0[(j+m+x)%m]*k;
                for (auto [x,k]:A1) tmp1+=y0[(j+m+x)%m]*k;
                for (auto [x,k]:A2) tmp2+=y0[(j+m+x)%m]*k;
                tmp+=(tmp2*mu+tmp1)*mu;
            }
        }
        return solution(y,xl,xr,tl,tr);
    }
};
class leapfrog
{
public:
template<typename func> solution solve(const func &boundary,db a,db xl,db xr,db tl,db tr,db h,db k) const
    {
        xl+=ceil((tr-tl)/k)*h*(-2);
        xr+=ceil((tr-tl)/k)*h*2;
        int m=round((xr-xl)/h);
        int n=round((tr-tl)/k);
        db mu=a*k/h;
        vector y(n+1,valarray<db>(m+1)); 
        for (int i=0;i<=m;i++)
        {
            y[0][i]=boundary(xl+i*h);
        }
        ++m;
        for (int i=0;i<m;i++)
        {
            y[1][i]=0.5*(y[0][(i+1)%m]+y[0][(i+m-1)%m])-(mu/2)*(y[0][(i+1)%m]-y[0][(i+m-1)%m]);
        }
        for (int i=1;i<n;i++)
        {
            for (int j=0;j<m;j++)
            {
                y[i+1][j]=y[i-1][j]-mu*(y[i][(j+1)%m]-y[i][(j+m-1)%m]);
            }
        }
        return solution(y,xl,xr,tl,tr);
    }
};
class LF:public advection
{
public:
    LF():advection({{1,1/2.},{-1,1/2.}},{{1,-1/2.},{-1,1/2.}},{}){}
};
class LW:public advection
{
public:
    LW():advection({{0,1.}},{{1,-1/2.},{-1,1/2.}},{{1,1/2.},{0,-1.},{-1,1/2.}}){}
};
class upwind:public advection
{
public:
    upwind():advection({{0,1.}},{{0,-1.},{-1,1.}},{}){}
};
class BW:public advection
{
public:
    BW():advection({{0,1.}},{{0,-3/2.},{-1,2.},{-2,-1/2.}},{{0,1/2.},{-1,-1.},{-2,1/2.}}){}
};