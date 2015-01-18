// Copyright (C) 2013  Anand Bhaskar
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// <http://www.gnu.org/licenses/>
// 
// email: bhaskar@berkeley.edu


// Contains some common functionality needed by the other files

#ifndef __COMMON__
#define __COMMON__

#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <iostream>
using namespace std;

#define FOR(i,a,b)	for(int i=(a);i<(b);++i)
#define REP(iter,v) for(typeof((v).begin()) iter = (v).begin(); iter != (v).end(); ++iter)
#define ALL(v)	((v).begin()), ((v).end())
#define MP make_pair
#define PB push_back
#define SZ size()

typedef long long ll, int64;
typedef long double lld;
typedef vector<int> VI;
typedef vector<vector<int> > VVI;
typedef lld DOUBLE;
typedef vector<DOUBLE> VD;
typedef vector<vector<DOUBLE> > VVD;
typedef vector<vector<vector<DOUBLE> > > VVVD;

template<class T>
void splitstr(const string &s, vector<T> &out)
{
    istringstream in(s);
    out.clear();
    copy(istream_iterator<T>(in), istream_iterator<T>(), back_inserter(out));
}

DOUBLE binom(int n, int k)	{
	DOUBLE ret = 1.;
	if (k > n - k)	k = n - k;
	for (int i = n; i >= n - k + 1; i--)	{
		ret = ret * i / (n - i + 1);
	}
	return ret;
}

#endif
