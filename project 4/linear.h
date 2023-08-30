#pragma once
#include "bits/stdc++.h"
using namespace std;
template<typename T> valarray<T> Tridiagonal_Matrix(valarray<T> l,valarray<T> d,valarray<T> u,valarray<T> b)
{
	int n=b.size();
	assert(d.size()==n&&l.size()==n-1&&u.size()==n-1);
	for (int i=0; i<n-1; i++)
	{
		u[i]/=d[i];
		d[i+1]-=l[i]*u[i];
		b[i+1]-=l[i]/d[i]*b[i];
	}
	valarray<T> x=b/d;
	for (int i=n-2; i>=0; i--) x[i]-=u[i]*x[i+1];
	return x;
}
template<typename T> valarray<T> Gauss(vector<valarray<T>> a,valarray<T> b)//ax=b
{
	auto A=a,B=b;
	int n=a[0].size(),m=a.size(),i,j,k;
	assert(n&&n<=m&&m==b.size());
	for (i=0; i<n; i++)
	{
		k=i;
		for (j=i+1; j<m; j++) if (abs(a[k][i])<abs(a[j][i])) k=j;
		swap(a[i],a[k]); swap(b[i],b[k]);
		for (j=0; j<m; j++) if (j!=i)
		{
			T xs=a[j][i]/a[i][i];
			for (k=i; k<n; k++) a[j][k]-=a[i][k]*xs;
			b[j]-=b[i]*xs;
		}
	}
	valarray<T> x(n);
	for (i=0; i<n; i++) x[i]=b[i]/a[i][i];
	return x;
}