#include "bits/stdc++.h"
using namespace std;
struct A
{
	int y;
};
using t=unique_ptr<A>(*)();
void f(t x)
{

}
int main()
{
	f(unique_ptr<A>());
}