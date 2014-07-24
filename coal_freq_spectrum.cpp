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

// Computes the expected frequency spectrum under the coalescent
// Example usage: 
// ./coal_freq_spectrum coalModel3.txt 10 5
// Prints the first 5 entries of the expected frequency spectrum under the coalescent for model coalModel3.txt 
// for a sample of size 10.


#ifndef __COAL_FREQ_SPECTRUM__
#define __COAL_FREQ_SPECTRUM__

#include <boost/math/special_functions/expint.hpp>
using namespace boost::math;

#include "coal_freq_spectrum_lib.cpp"
using namespace std;


int main(int argc, char **argv)	{
	assert (argc >= 3);
    
	string modelFile = string(argv[1]);
    int n = atoi(argv[2]);

	//number of leading entries of the frequency spectrum desired
	int maxA = n - 1;
	if (argc >= 4)	{
		maxA = atoi(argv[3]);
		if (maxA >= n)	maxA = n - 1;
		if (maxA < 1)	maxA = 1;
		assert (1 <= maxA && maxA < n);
	}

	cout << "sample size n = " << n << endl;
	cout << endl;

	vector<EpochInfo> epochs = initModelFromFile(modelFile);
	
	DOUBLE Nconst = epochs[epochs.SZ - 1].N;

	DOUBLE normalizationConstant;
    VD coalFreqS = coalFreqSpectrumHelper(epochs, n, maxA, normalizationConstant);
    cout << "Expected # of segregating sites with i derived alleles, with theta = 1, first " << maxA << " entries" << endl;
    FOR (i, 0, coalFreqS.SZ)	{
    	printf("%0.20Lg\t", coalFreqS[i] / 2. / Nconst);
    }
    printf("\n\n");

    cout << "Probability distribution of # derived alleles at a segregating site, first " << maxA << " entries" << endl;
    FOR (i, 0, coalFreqS.SZ)	{
    	printf("%0.20Lg\t", coalFreqS[i]/normalizationConstant);
    }
 	printf("\n");
}

#endif
