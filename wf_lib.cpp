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


// Contains some common functionality related to computation under the Wright-Fisher model that is needed by the other files

#ifndef __WF_LIB__
#define __WF_LIB__

#include "common.cpp"
#include <cassert>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
using namespace std;


int getTIdx(int t, const VI& times)	{
	unsigned int idx = (lower_bound(ALL(times), t) - times.begin());
	if (idx >= times.SZ)	{
		assert (false);
	}
	if (times[idx] == t)	{		//on the boundary (t_{idx - 1}, t_idx], so use the population size/Qmat in the next discretization interval
		idx++;
	}
	return idx;
}


// reads the Wright-Fisher model in fileName and stores the changepoints in times and population sizes in Ns
void readWFModelFromFile(string& fileName, VI& times, vector<ll>& Ns)	{
	ifstream fin(fileName.c_str());
	string buff;
	vector<string> arr;
	bool readFirstLine = true;
	cout << "# Reading model from file " << fileName << endl;

	ll Nconst;
	while (!fin.eof())	{
		getline(fin, buff);
		if (buff[0] == '#')	continue;
		if (buff.SZ == 0)	break;
		if (readFirstLine)	{
			Nconst = atol(buff.c_str());	//ancestral population size
			readFirstLine = false;
		}
		else	{
			splitstr(buff, arr);	//time	population_size
			times.PB(atoi(arr[0].c_str()));
			Ns.PB(atol(arr[1].c_str()));
		}
	}
	times.PB(numeric_limits<int>::max());	Ns.PB(Nconst);
	
	fin.close();
}



//Q is the transition matrix for the WF model. Q[i][j] = Probability that i individuals at the current 
//generation (with population curPopSize) have j parents in the previous generation (with population popSize)
//the returned Q matrix has dimensions min(n, curPopSize) rows and min(n, popsize, curPopSize) columns
VVD calcQ(int n, ll popSize, ll curPopSize)	{
	// assert (popSize >= n);
	int imax = (int)min((ll)n, curPopSize);
	int jmax = (int)min((ll)n, popSize);
	jmax = min(imax, jmax);
	VVD Q(imax + 1, VD(jmax + 1, 0.0));

	Q[0][0] = 1.;

	for (int i = 1; i <= imax; i++)	{
		int jmax = (int)min((ll)i, popSize);
		#pragma omp parallel for
		for (int j = 1; j <= jmax; j++)   {
			Q[i][j] = Q[i-1][j-1] * (popSize - j + 1.) / popSize + (j <= (i - 1) ? Q[i-1][j] * j / popSize : 0);
		}
	}

	return Q;
}


//same as calcQ except that the entries where Q would be < eps are removed. For each row i, the left and right 
//index bounds for j where P(i -> j) > eps are returned in the vector v.
VVD calcQWithTruncation(int n, ll popSize, DOUBLE eps, vector<pair<int, int> >& v, ll curPopSize)	{
	//assert (popSize >= n);
	int imax = (int)min((ll)n, curPopSize);
	int jmax = (int)min((ll)n, popSize);
	jmax = min(imax, jmax);
	VVD _Q(imax + 1, VD(jmax + 1, 0.0));

	_Q[0][0] = 1.;

	for (int i = 1; i <= imax; i++)	{
		int jmax = (int)min((ll)i, popSize);
		#pragma omp parallel for
		for (int j = 1; j <= jmax; j++)   {
			_Q[i][j] = _Q[i-1][j-1] * (popSize - j + 1.) / popSize + _Q[i-1][j] * j / popSize;
		}
	}

	VVD Q(imax + 1);
	v.resize(imax + 1);
	
	FOR (i, 0, imax + 1)	{
		int idx = max_element(ALL(_Q[i])) - _Q[i].begin();
		int j;			
		for (j = idx - 1; j >= 1; j--)	{
			assert(_Q[i][j] <= _Q[i][j+1]);
			if (_Q[i][j] < eps * _Q[i][idx])	break;
		}
		int lbd = j + 1;

		int jmax = (int)min((ll)i, popSize);
		for (j = idx + 1; j <= jmax; j++)	{
			assert(_Q[i][j] <= _Q[i][j-1]);
			if (_Q[i][j] < eps * _Q[i][idx])	break;
		}
		int rbd = j - 1;

		v[i] = MP(lbd, rbd);
		
		//set Q[i] to be _Q[i][lbd..rbd]
		Q[i] = vector<DOUBLE>(_Q[i].begin() + lbd, _Q[i].begin() + rbd + 1);
	}

	return Q;
}




//probability of k-drops due to only simultaneous mergers (i.e. no multiple mergers)
VVD calcR(int n, ll popsize, ll curPopSize)	{
	//assert (popsize >= n);
	int imax = (int)min((ll)n, curPopSize);
	int jmax = (int)min((ll)n, popsize);
	jmax = min(imax, jmax);

	VVD R(imax + 1);

	R[1] = VD(2, 0.);
	R[1][1] = 1.;

	for (int j = 2; j <= jmax; j++)	{
		R[j] = VD(jmax + 1, 0.0);
		R[j][j] = R[j-1][j-1] * ((popsize - (j-1)) / (DOUBLE)popsize);
	}
	
	for (int j = 1; j <= jmax; j++)   {
		for (int i = j+1; i <= imax; i++)	{
			R[i][j] = R[i-1][j] / (DOUBLE)popsize * i / (i - j) * (2*j - i + 1) / 2.;
		}
	}
	return R;
}


#endif
