#pragma once
#include "bits/stdc++.h"
using namespace std;
template<typename T> class factory
{
	map<string,shared_ptr<T>> mp;
	factory()=default;
	factory(const factory &)=default;
	factory &operator=(const factory &)=default;
	~factory()=default;
public:
	static factory &create_factory()
	{
		static factory fac;
		return fac;
	};
	bool insert(const string &name,shared_ptr<T> p)
	{
		return mp.insert({name,p}).second;
	}
	int count(const string &name) const
	{
		return mp.count(name);
	}
	bool erase(const string &name)
	{
		if (count(name))
		{
			mp.erase(name);
			return 1;
		}
		return 0;
	}
	shared_ptr<T> at(const string &name)
	{
		if (count(name))
		{
			return mp[name];
		}
		throw invalid_argument("Unknown name \""+name+"\".");
	}
	shared_ptr<T> operator[](const string &name)
	{
		if (count(name))
		{
			return mp[name];
		}
		throw invalid_argument("Unknown name \""+name+"\".");
	}
};