#pragma once
#include <type_traits>
template<typename func,typename T> T secant(const func &f,T x0,T eps=1e-13)
{
	static_assert(is_floating_point<T>::value);
	T x1=x0+0.1;
	T f0=f(x0);
	T f1=f(x1);
	while (abs(x1-x0)>eps)
	{
		x0=x1-f1/(f1-f0)*(x1-x0);
		if (abs(x1-x0)<eps) break;
		f0=f(x0);
		x1=x0-f0/(f0-f1)*(x0-x1);
		f1=f(x1);
	}
	return x1;
}