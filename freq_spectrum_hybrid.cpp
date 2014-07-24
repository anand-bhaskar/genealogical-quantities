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

// Computes the frequency spectrum using the discrete WF model for the first few generations, and then
// switches to the coalescent for the computation further back in time.
//
// Example usage: 
// ./freq_spectrum_hybrid model3.txt coalModel3.txt 2 1e-100 10 20 30 40
// Computes the first 10 entries of the frequency spectrum for model 3 for sample sizes 20, 30, and 40 using the discrete-time Wright-Fisher 
// model for the most recent 2 generations, and then using the coalescent further back in time. The demographic model is given in 
// model3.txt (Wright-Fisher) and coalModel3.txt (coalescent). 
// Uses a parameter of epsilon = 1e-100 for truncating the transition matrix of the Wright-Fisher (this is helpful for very large sample sizes).
//
// ./freq_spectrum_hybrid model3.txt coalModel3.txt 0 0 10 20 30 40
// Computes the first 10 entries of the frequency spectrum under the coalescent for model 3 for sample sizes 20, 30, and 40. By specifying 0 for
// the third argument, the computation is done under the coalescent.

#ifndef __FREQ_SPECTRUM_HYBRID__
#define __FREQ_SPECTRUM_HYBRID__

#include "common.cpp"
#include "coal_freq_spectrum_lib.cpp"
#include "wf_lib.cpp"
#include <numeric>
using namespace std;


DOUBLE adaptiveEpsConstPop = 0., adaptiveEpsVarPop = 0.;
VVD W;
vector<VVD> Qmap;
vector<vector<pair<int, int> > > truncIdxesMap;
vector<pair<int, int> > truncIdxesConstPop;

vector<ll> Ns;
VI times;
string WFModelFile, coalModelFile;
int tUseCoal;

ll Nconst;
VVD fAdaptVarMemo;

DOUBLE uMut;


void goFVarPopSize(int n, int maxA)  {
	int cur = 0, prev = 1;
	VVVD fAdaptVarMemoInternal(2, VVD(maxA + 1, VD(n + 1, 0.0)));
    ll prevNt = -1, prevCurNt = -1;
	VVD Qmat;
	vector<pair<int, int> > truncIdxes;
	
	//create the demographic scenario for the cutoff time tUseCoal
	vector<EpochInfo> epochs = initModelFromFile(coalModelFile, tUseCoal);
	for (int sampleSize = 2; sampleSize <= n; sampleSize++)	{
		//compute first min(maxA, n-1, sampleSize) entries of the frequency spectrum under the coalescent
		int numEntries = min(min(maxA, n - 1), sampleSize - 1);
		VD coalFreqS = coalFreqSpectrum(epochs, sampleSize, numEntries);
		for (int a = 1; a <= maxA && a <= sampleSize - 1 && a <= n - 1; a++)	{
			fAdaptVarMemoInternal[cur][a][sampleSize - a] = coalFreqS[a - 1] / binom(sampleSize, a);
		}
	}
	
	for (int t = tUseCoal - 1; t >= 0; t--)	{
		cur ^= 1;
		prev ^= 1;
		int tIdx = getTIdx(t, times);
	    ll Nt = Ns[tIdx];
		ll curNt = Ns[getTIdx(t-1, times)];
	    if (Nt != prevNt || curNt != prevCurNt)   {
			Qmat = calcQWithTruncation(n, Nt, adaptiveEpsVarPop, truncIdxes, curNt);
		}

		#pragma omp parallel for
		FOR (a, 1, min(maxA, n) + 1)	{
			int aStart = truncIdxes[a].first, aEnd = truncIdxes[a].second;
			FOR (b, 1, n + 1)	if (a + b <= n && a + b <= curNt)	{
				int bStart = truncIdxes[b].first, bEnd = truncIdxes[b].second;
				if (a == 1) {
					if (b == 1)	{
						fAdaptVarMemoInternal[cur][a][b] = 1. + Qmat[2][2-truncIdxes[2].first] * fAdaptVarMemoInternal[prev][a][b];
					}
					else	{
						DOUBLE ret = 1.;
						FOR (m, bStart, bEnd + 1)	if (1 + m <= Nt)	{
					        ret += 1. * (Nt - m) / Nt * Qmat[b][m-bStart] * fAdaptVarMemoInternal[prev][a][m];
					    }
						fAdaptVarMemoInternal[cur][a][b] = ret;
					}
			    }
				else	{
			    	DOUBLE ret = 0.;
			    	FOR (j, aStart, aEnd + 1)   {
						DOUBLE factor = 1.;
						FOR (i, 0, j)	factor *= 1. * (Nt - bStart - i) / (Nt - i);
			        	FOR (k, bStart, bEnd + 1)    if (j + k <= Nt)	{
			            	ret += Qmat[a][j-aStart] * Qmat[b][k-bStart] * factor * fAdaptVarMemoInternal[prev][j][k];
							factor *= 1. * (Nt - j - k) / (Nt - k);
			        	}
			    	}
					fAdaptVarMemoInternal[cur][a][b] = ret;
				}
			}
		}
		prevNt = Nt;
		prevCurNt = curNt;
	}
	fAdaptVarMemo = fAdaptVarMemoInternal[cur];
}


void freqSpectrumVarPopSize(VI& ns, int maxA)	{
	int nmax = *max_element(ALL(ns));
	
    VD ret;
	goFVarPopSize(nmax, maxA);
	
	cout << endl;
    cout << "Expected # of segregating sites with i derived alleles, with theta = " << uMut * 2 * Nconst << endl;
	REP (it, ns)	{
		int n = *it;
		cout << "n = " << n << endl;
	    FOR (k, 1, min(n - 1, maxA) + 1)   {
			DOUBLE val = fAdaptVarMemo[k][n-k] * binom(n, k) * uMut;
			printf("%0.20Lg\t", val);
		}
		cout << endl;
	}
}



int main(int argc, char **argv)	{
	assert (argc >= 7);
	WFModelFile = argv[1];
	coalModelFile = argv[2];
	tUseCoal = atoi(argv[3]);
	adaptiveEpsVarPop = atof(argv[4]);
	int maxA = atoi(argv[5]);

	VI ns;
	FOR (i, 6, argc)	{
		ns.PB(atoi(argv[i]));
	}

	cout << "adaptive truncation eps for variable pop size model = " << adaptiveEpsVarPop << endl;
	cout << "sample size(s) n = ";
	REP (it, ns)	cout << *it << ", ";
	cout << endl;
	int nmax = *max_element(ALL(ns));
	cout << "largest sample size = " << nmax << endl;
	
	readWFModelFromFile(WFModelFile, times, Ns);
	Nconst = Ns[Ns.SZ - 1];
	uMut = 1. / (2. * Nconst);

    cout << "constant population size N = " << Nconst << endl;
    cout << endl;

	assert(Ns.SZ == times.SZ);
	
    cout << "Variable population size" << endl;
    cout << "------------------------" << endl;    
    freqSpectrumVarPopSize(ns, maxA);
	cout << endl;
	
	return 0;
}


#endif
