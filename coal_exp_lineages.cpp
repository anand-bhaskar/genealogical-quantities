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

// Computes the expected value and standard deviation of the number of surviving lineages in the coalescent
// as a function of time
//
// Example usage: 
// ./coal_exp_lineages coalModel3.txt 20 10
// Outputs the expected number and standard deviation under the coalescent for the number of ancestors of a 
// sample of size 20 for each of the most recent 10 generations.


#ifndef __COAL_EXP_LINEAGES__
#define __COAL_EXP_LINEAGES__

#include "coal_freq_spectrum_lib.cpp"
using namespace std;

#include <string.h>

vector<EpochInfo> epochs;

VD factor, factor2;

void precomputeFactors(int n)	{
	factor.resize(n + 1, 0.);
	factor2.resize(n + 1, 0.);
	DOUBLE p = 1.;
	FOR (i, 1, n + 1)	{
		p *= 1. * (n - i + 1) / (n + i - 1);
		factor[i] = (2. * i - 1.) * p;
        factor2[i] = factor[i] * (((DOUBLE)i)*((DOUBLE)i) - i + 1.);
	}
}

void goExpLineages(int n, int nGen)	{
	FOR (i, 0, nGen + 1)	{
		DOUBLE t = i;
		DOUBLE exponent = 0., expectedValue = 0., expectedSqValue = 0.;
		unsigned int eidx = 0;
		while (epochs[eidx].t <= t)	{
			if (epochs[eidx].expGrowth == false)	{
				exponent += (epochs[eidx].t - (eidx == 0 ? 0. : epochs[eidx - 1].t)) / epochs[eidx].N;
			}
			else	{
				exponent += (exp(epochs[eidx].beta * (epochs[eidx].t - (eidx == 0 ? 0. : epochs[eidx - 1].t))) -  1.) / epochs[eidx].beta / epochs[eidx].N;
			}
			eidx++;
		}
		assert(eidx < epochs.SZ);
		if (epochs[eidx].expGrowth == false)	{
			exponent += (t - (eidx == 0 ? 0. : epochs[eidx - 1].t)) / epochs[eidx].N;
		}
		else	{
			exponent += (exp(epochs[eidx].beta * (t - (eidx == 0 ? 0. : epochs[eidx - 1].t))) -  1.) / epochs[eidx].beta / epochs[eidx].N;
		}
		FOR (j, 1, n + 1)	{
			expectedValue += exp(-1. * j * (j - 1) / 2. * exponent) * factor[j];
			expectedSqValue += exp(-1. * j * (j - 1) / 2. * exponent) * factor2[j];
		}
		DOUBLE stdDev = expectedSqValue < expectedValue * expectedValue ? 0. : sqrt(expectedSqValue - expectedValue * expectedValue);
		printf("%d\t%.10Lg\t%.10Lg\n", i, expectedValue, stdDev);
	}
}


int main(int argc, char **argv)	{
	assert (argc >= 4);
    
	string modelFile = argv[1];
    int n = atoi(argv[2]);
	int nGens = atoi(argv[3]);

	precomputeFactors(n);

	epochs = initModelFromFile(modelFile);

	goExpLineages(n, nGens);
}

#endif