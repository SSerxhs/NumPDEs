#ifndef __Header_Gauss
#define __Header_Gauss
#include "bits/stdc++.h"
using namespace std;
template<typename T> vector<T> Gauss(vector<vector<T>> a,vector<T> b)//ax=b
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
	// if (n+1==m) assert(abs(*max_element(a.back().begin(),a.back().end()))<1e-9);
	vector<T> x(n);
	for (i=0; i<n; i++) x[i]=b[i]/a[i][i];
	return x;
}
#endif