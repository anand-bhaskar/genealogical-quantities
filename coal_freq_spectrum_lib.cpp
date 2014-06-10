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


// Contains some common functionality related to computation under the coalescent that is needed by the other files

#ifndef __COAL_FREQ_SPECTRUM_LIB__
#define __COAL_FREQ_SPECTRUM_LIB__

#include "omp.h"
#include "common.cpp"
#include <fstream>
using namespace std;

#include <boost/math/special_functions/expint.hpp>
using namespace boost::math;



//Represents an epoch of exponential growth/decline, or constant population size.
//Several objects of this structure can be stored together in a vector to represent a demographic model.
struct EpochInfo {
	//time (in units of present haploid population size N0) at which this epoch ends (looking from the present going back). 
	//Don't confuse this with the duration/length of the epoch, which can be got by looking at the time the previous epoch ends.
	//Eg. t = 0.5 means this epoch ends at generation 0.5*N0 (and begins at a time dependent on when the previous epoch ended)
    DOUBLE t;
	//population size (in units of present haploid population size N0) at the beginning of the epoch
    DOUBLE N;
	//true if the epoch is an exponential growth/decline epoch, and false if the population size is constant in the epoch
    bool expGrowth;
	//forward-in-time growth rate (in units of present haploid population size N0) for the epoch. This is positive if the population 
	//is undergoing an expansion going forward in time.
    DOUBLE beta;
    
    EpochInfo(DOUBLE _t, DOUBLE _N, bool _expGrowth, DOUBLE _beta) : t(_t), N(_N), expGrowth(_expGrowth), beta(_beta) {
        ;
    }	

	EpochInfo()	{
		;
	}
	
	EpochInfo(const EpochInfo& o)	{
		t = o.t;
		N = o.N;
		expGrowth = o.expGrowth;
		beta = o.beta;
	}
};


void printEpochInfos(const vector<EpochInfo>& epochs, FILE* out = stderr)	{
	FOR (i, 0, epochs.SZ)	{
		EpochInfo epoch = epochs[i];
		if (epoch.expGrowth)	{
			fprintf(out, "# Epoch %d, exponential growth, N = %0.10Lg, time = %0.8Lg, rate = %0.6Lg\n", i, epoch.N, epoch.t, exp(epoch.beta) - 1.);
		}
		else {
			fprintf(out, "# Epoch %d, constant, N = %0.10Lg, time = %0.8Lg\n", i, epoch.N, epoch.t);
		}
	}
}

//Compute exp(x) Ei(-x) by the divergent asymptotic expansion \sum_{n=0}^{order} n!/(-x)^{n+1}
DOUBLE asymptoticExpTimesEi(DOUBLE x)	{
	int order = 5;
	DOUBLE sgn = -1., p = 1. / x, ret = 0.;
	FOR (n, 0, order)	{
		ret += sgn * p;
		p *= (n+1.) / x;
		sgn = -sgn;
	}
	return ret;
}

//Equation 12 of Polanski and Kimmel (Genetics, 2003)
VD computeV(int n)    {
    VD V(n + 1, 0.);
    
    //V[j] = (2j - 1) \fall{n - 1}{j-1} \over \rise{n + 1}{j - 1}   if j is even, 0 if j is odd
    V[2] = (n - 1.) / (n + 1.);
    FOR (j, 4, n + 1)   if (j % 2 == 0) {
        V[j] = V[j - 2] * (n - j + 1) / (n + j - 1) * (n - j + 2) / (n + j - 2);
    }
    FOR (j, 2, n + 1)   {
        V[j] *= 2. * (2*j - 1);
    }
	return V;
}

//Equations 13-15 of Polanski and Kimmel (Genetics, 2003)
VD computeW(int n, int b)    {
    VD W(n + 1, 0.);
    W[2] = 6. / (n + 1.);
	if (n >= 3)	{
    	W[3] = 30. * (n - 2.*b) / (n + 1.) / (n + 2.);
	}
    FOR (j, 2, n - 2 + 1)   {
        W[j+2] = - (1. + j) / j * (3. + 2*j) / (2*j - 1.) * (n - j) / (n + j + 1.) * W[j] + (3. + 2.*j) / j * (n - 2.*b) / (n + j + 1.) * W[j+1];
    }
    return W;
}


//Equation 3 of Polanski and Kimmel (Genetics, 2003)
VD computeEjUnscaled(const vector<EpochInfo>& epochs, int n)   {
    int M = epochs.SZ;	//# of time intervals
    
    VD Omega(M, 0.);
	FOR (i, 0, M)	{
	    if (i > 0)  {
            Omega[i] = Omega[i - 1];
	    }
        if (epochs[i].expGrowth)    {
            Omega[i] += 1. / epochs[i].N / epochs[i].beta * (exp(epochs[i].beta * (epochs[i].t - (i == 0 ? 0. : epochs[i - 1].t))) - 1.);
        }   
        else    {
            Omega[i] += 1. / epochs[i].N * (epochs[i].t - (i == 0 ? 0. : epochs[i - 1].t));
        }
	}
	
	DOUBLE useAsymptoticExpansionConstant = 45.;	//use the asymptotic approximation for exp(x) Ei(-x)
	VD Ej(n + 1, 0.);
	FOR (j, 2, n + 1)   {
		DOUBLE jChoose2 = binom(j, 2);
	    FOR (i, 0, M)   { 
	        if (epochs[i].expGrowth)    {                
                DOUBLE f1 = 1. / epochs[i].N / epochs[i].beta;
                DOUBLE f2 = 1. / epochs[i].N / epochs[i].beta * exp(epochs[i].beta * (epochs[i].t - (i == 0 ? 0. : epochs[i - 1].t)));
                
				assert(f2 >= f1);
				assert(Omega[i] >= (i == 0 ? 0. : Omega[i-1]));
				
				DOUBLE t = 0., t2 = 0.;
				if (jChoose2 * f1 >= useAsymptoticExpansionConstant)	{
					t = 1. / epochs[i].beta * exp(- jChoose2 * (i == 0 ? 0. : Omega[i - 1])) * asymptoticExpTimesEi(jChoose2 * f1);
				}
				else	{
					t = 1. / epochs[i].beta * exp(- jChoose2 * ((i == 0 ? 0. : Omega[i - 1]) - f1)) * expint(- jChoose2 * f1);
				}
								
				if (jChoose2 * f2 >= useAsymptoticExpansionConstant)	{					
					t2 = 1. / epochs[i].beta * exp(- jChoose2 * Omega[i]) * asymptoticExpTimesEi(jChoose2 * f2);
				}
				else	{
					t2 = 1. / epochs[i].beta * exp(- jChoose2 * (Omega[i] - f2)) * expint(- jChoose2 * f2);
				}			

                Ej[j] += t2 - t;
	        }
	        else    {              
		        DOUBLE t2 = 1. * epochs[i].N / jChoose2 * exp(-jChoose2 * (i == 0 ? 0. : Omega[i-1]));
                DOUBLE t = 1. * epochs[i].N / jChoose2 * exp(-jChoose2 * Omega[i]);
                Ej[j] += t2 - t;
            }
	    }
	}
	
	return Ej;
}

//Returns the first maxA entries of the branch lengths (in units of generations) in the coalescent subtending 
//i lineages (1 <= i <= maxA) for a sample of size n under the infinite-sites coalescent model for the demographic 
//model represented by epochs.
//This can be normalized to get the expected frequency spectrum as a probability distribution, and this normalization
//constant is returned in the variable normalizationConstant
//Computes it using equation 8 of Polanski and Kimmel (Genetics, 2003)
VD coalFreqSpectrumHelper(const vector<EpochInfo>& epochs, int n, int maxA, DOUBLE& normalizationConstant)	{
    VD V = computeV(n);
    VD Ej = computeEjUnscaled(epochs, n);
    
    VD ret(min(maxA, n - 1), 0.);

	DOUBLE den = 0.;
	FOR (j, 2, n + 1)   {
    	den += Ej[j] * V[j];
	}
	
    #pragma omp parallel for
    FOR (b, 1, min(maxA, n - 1) + 1)   {
        VD W = computeW(n, b);
        DOUBLE num = 0.;
        FOR (j, 2, n + 1)   {
            num += Ej[j] * W[j];
        }

		//den is the normalizing constant to get a probability distribution from the expected frequency spectrum.
		//To get the unnormalized expected frequency spectrum (in units of generations), do not divide by den.
        // ret[b - 1] = num / den;
		ret[b - 1] = num;
    }
    
 	normalizationConstant = den;
    return ret;
}


//Returns the first maxA entries of the branch lengths (in units of generations) in the coalescent subtending 
//i lineages (1 <= i <= maxA) for a sample of size n under the infinite-sites coalescent model for the demographic 
//model represented by epochs.
VD coalFreqSpectrum(const vector<EpochInfo>& epochs, int n, int maxA)	{
	DOUBLE tmp;
	return coalFreqSpectrumHelper(epochs, n, maxA, tmp);
}

//Returns the first maxA entries of the expected frequency spectrum normalized to be a probability distribution.
VD coalFreqSpectrumNormalized(const vector<EpochInfo>& epochs, int n, int maxA)	{
	DOUBLE normalizationConstant;
	VD ret = coalFreqSpectrumHelper(epochs, n, maxA, normalizationConstant);

	FOR (i, 0, ret.SZ)	{
		ret[i] /= normalizationConstant;
	}

	return ret;
}


//Returns the first maxA entries of the branch lengths (in units of generations) in the coalescent subtending 
//either i or n - i lineages (1 <= i <= maxA) for a sample of size n under the infinite-sites coalescent model 
//for the demographic model represented by epochs.
//This can be normalized to get the expected folded frequency spectrum as a probability distribution.
//Computes it using equation 8 of Polanski and Kimmel (Genetics, 2003)
VD coalFoldedFreqSpectrumHelper(const vector<EpochInfo>& epochs, int n, int maxA, DOUBLE& normalizationConstant)	{
	VD V = computeV(n);
	VD Ej = computeEjUnscaled(epochs, n);
    
    VD ret(min(maxA, n / 2), 0.);

	DOUBLE den = 0.;
	FOR (j, 2, n + 1)   {
    	den += Ej[j] * V[j];
	}
	
    #pragma omp parallel for
    FOR (b, 1, min(maxA, n / 2) + 1)   {
		VI x;
		x.PB(b);
		if (b != n - b)	x.PB(n - b);
		FOR (l, 0, x.SZ)	{
	        VD W = computeW(n, x[l]);
	        DOUBLE num = 0.;
	        FOR (j, 2, n + 1)   {
	            num += Ej[j] * W[j];
	        }
			ret[b - 1] += num;
		}
    }
    
    normalizationConstant = den;

    return ret;
}


//Returns the first maxA entries of the branch lengths (in units of generations) in the coalescent subtending 
//i or n - i lineages (1 <= i <= maxA) for a sample of size n under the infinite-sites coalescent model for the demographic 
//model represented by epochs.
VD coalFoldedFreqSpectrum(const vector<EpochInfo>& epochs, int n, int maxA)	{
	DOUBLE tmp;
	return coalFoldedFreqSpectrumHelper(epochs, n, maxA, tmp);
}

//Returns the first maxA entries of the expected frequency spectrum normalized to be a probability distribution.
VD coalFoldedFreqSpectrumNormalized(const vector<EpochInfo>& epochs, int n, int maxA)	{
	DOUBLE normalizationConstant;
	VD ret = coalFoldedFreqSpectrumHelper(epochs, n, maxA, normalizationConstant);

	FOR (i, 0, ret.SZ)	{
		ret[i] /= normalizationConstant;
	}

	return ret;
}


//read epoch info from the configuration file modelFile
//samplingTime is the time for which the SFS computation is being done
vector<EpochInfo> initModelFromFile(const string& modelFile, int samplingTime = 0)	{
	ifstream fin(modelFile.c_str());
	vector<string> arr;
	string buff;

	getline(fin, buff);
	splitstr(buff, arr);
	DOUBLE NAnc = atof(arr[0].c_str());

	vector<EpochInfo> epochs;
	
	cerr << "# Reading demographic model from file " << modelFile << endl;
	cerr << "# Initializing coalescent model for sampling time " << samplingTime << endl;
	
	vector<bool> epochNIsKnown;
	
	while (! fin.eof())	{
		getline(fin, buff);
		splitstr(buff, arr);
	
		if (arr.SZ == 0)	continue;
	
		EpochInfo epoch;

		assert (arr[0] == string("e") || arr[0] == string("c"));

		if (arr[0] == string("e"))	{	//exponential epoch
			assert (arr.SZ == 4);
	
			if (arr[1] == string("*"))	{	//need to fill in the population size at the start of the epoch later
				epochNIsKnown.PB(false);
			}
			else	{
				epochNIsKnown.PB(true);
				epoch.N = atof(arr[1].c_str());
			}
			epoch.expGrowth = true;
			epoch.t = atof(arr[2].c_str());
			epoch.beta = log(1. + atof(arr[3].c_str()));
		}
		
		if (arr[0] == string("c"))	{	//constant pop size epoch
			assert (arr.SZ == 3);

			epochNIsKnown.PB(true);
			epoch.expGrowth = false;
			epoch.N = atof(arr[1].c_str());
			epoch.t = atof(arr[2].c_str());
		}
		
		epochs.PB(epoch);
	}
	
	//last epoch of ancestral population size
	epochNIsKnown.PB(true);
	epochs.PB(EpochInfo(numeric_limits<DOUBLE>::max(), NAnc, false, -1.0));
	
	int numEpochs = epochs.SZ;
	
	//fill in the population size for the epochs where it is unknown
	for (int i = numEpochs - 1; i >= 0; i--)	if (epochs[i].expGrowth && ! epochNIsKnown[i])	{
		DOUBLE tGrowthInterval = epochs[i].t - (i == 0 ? 0. : epochs[i-1].t);
		epochs[i].N = epochs[i+1].N * exp(epochs[i].beta * tGrowthInterval);
	}

	cerr << endl;
	printEpochInfos(epochs);
	cerr << endl;
	
	//remove the epochs before the sampling time
	int i = 0;
	while (i < numEpochs && epochs[i].t <= samplingTime)	i++;
	//Now epochs[i].t > samplingTime
	if (i < numEpochs && epochs[i].expGrowth)	{	//in the middle of an exponential growth phase, so set the appropriate population size at samplingTime
		epochs[i].N = epochs[i].N * exp(- epochs[i].beta * (samplingTime - (i == 0 ? 0. : epochs[i-1].t)));
	}
	epochs.erase(epochs.begin(), epochs.begin() + i);	//erase older epochs
	FOR (i, 0, epochs.SZ)	{
		epochs[i].t -= samplingTime;
	}
	cerr << "# After changing the sampling time, demography looks like:" << endl;
	printEpochInfos(epochs);
	cerr << endl;

	fin.close();
	
	return epochs;
}


#endif