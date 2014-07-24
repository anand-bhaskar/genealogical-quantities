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

// Computes the expected frequency spectrum using the discrete-time Wright-Fisher model.
//
// Example usage: 
// ./wf_freq_spectrum model3.txt 0 0 9 10 20 30 40
// Computes the expected frequency spectrum for model 3 for sample sizes 10, 20, 30, and 40 using the discrete-time Wright Fisher 
// model with the demography given in model3.txt.
// Uses no truncation (the two 0 values) for truncating the transition matrix of the Wright-Fisher model 
// The first truncation parameter is for the last epoch of constant population size, while the second truncation parameter is 
// for each epoch of the varying population size model.


#ifndef __WF_FREQ_SPECTRUM__
#define __WF_FREQ_SPECTRUM__

#include "wf_lib.cpp"
#include "common.cpp"

#include <fstream>
#include <iostream>
#include <limits>
using namespace std;


DOUBLE adaptiveEpsConstPop = 0., adaptiveEpsVarPop = 0.;
VVD W;
vector<VVD> Qmap;
vector<vector<pair<int, int> > > truncIdxesMap;
vector<pair<int, int> > truncIdxesConstPop;

vector<ll> Ns;
VI times;
string modelFile;

ll Nconst;
VVD Qconst, QconstFullMatrix;

VVD fAdaptMemo, fAdaptVarMemo;
VD hMemo, hVarMemo;

DOUBLE uMut;


// precomputation for the last epoch where the population size is constant further back in time
void goF(int n, int maxA)  {
	fAdaptMemo = VVD(maxA + 1, vector<DOUBLE>(n + 1, 0.0));
	FOR (a, 1, min(maxA, n) + 1)	{
		int aStart = truncIdxesConstPop[a].first, aEnd = truncIdxesConstPop[a].second;
		FOR (b, 1, n + 1)	if (a + b <= n && a + b <= Nconst)	{
			int bStart = truncIdxesConstPop[b].first, bEnd = truncIdxesConstPop[b].second;
			if (a == 1)	{
				DOUBLE factor = 1. - QconstFullMatrix[b + 1][b + 1];
				DOUBLE ans = 1. / factor;
				FOR (m, bStart, min(bEnd, b - 1) + 1)	{
					ans += 1. * (Nconst - m) / Nconst * Qconst[b][m - bStart] / factor * fAdaptMemo[a][m];
				}				
				fAdaptMemo[a][b] = ans;
			}
			else	{
				DOUBLE ans = 0.;
				#pragma omp parallel for reduction (+:ans) num_threads(4)
				FOR (j, aStart, aEnd + 1)	{
					DOUBLE factor = 1.;
					FOR (i, 0, j)	factor *= 1. * (Nconst - bStart - i) / (Nconst - i);
					FOR (k, bStart, bEnd + 1)	if (j != a || k != b)    {
			            ans += Qconst[a][j - aStart] * Qconst[b][k - bStart] * factor / (1. - QconstFullMatrix[a + b][a + b]) * fAdaptMemo[j][k];
						factor *= 1. * (Nconst - j - k) / (Nconst - k);
			        }
			    }
				fAdaptMemo[a][b] = ans;
			}
		}
	}
}

// Computes the total branch length subtending i descendants (1 <= i <= maxA) for a sample of size n for the given 
// demographic model. This is the numerator of the expected SFS probability distribution. 
void goFVarPopSize(int n, int maxA)  {
	int cur = 0, prev = 1;
	VVVD fAdaptVarMemoInternal(2, VVD(maxA + 1, VD(n + 1, 0.0)));
	fAdaptVarMemoInternal[cur] = fAdaptMemo;
	int tLast = -1;
	if (times.SZ >= 2)	{
		tLast = times[times.SZ - 2];
	}
    ll prevNt = -1, prevCurNt = -1;
	VVD Qmat;
	vector<pair<int, int> > truncIdxes;
	for (int t = tLast; t >= 0; t--)	{
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
						FOR (i, 0, j)	factor *= 1. * (Nt - bStart - i) / (Nt - i);	//should be faster to do this
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


// precomputation for the last epoch where the population size is constant further back in time
void goH(int n)	{
	hMemo.resize(n + 1, 0.0);
	hMemo[1] = 0.;
	
	FOR (a, 2, n + 1)	if(a <= Nconst)	{
		DOUBLE ret = 1. * a / (1. - QconstFullMatrix[a][a]);
		FOR (m, 1, a)	{
			ret += QconstFullMatrix[a][m] / (1. - QconstFullMatrix[a][a]) * hMemo[m];
		}
		hMemo[a] = ret;
	}
}


// Computes the total branch length for a sample of size n for the given demographic model.
// This is the denominator of the expected SFS probability distribution. 
void goHVarPopSize(int n)	{
	int cur = 0, prev = 1;
	VVD hVarMemoInternal(2, VD(n + 1, 0.0));
	vector<pair<int, int> > truncIdxes;
	hVarMemoInternal[cur] = hMemo;
	int tLast = -1;
	if (times.SZ >= 2)	{
		tLast = times[times.SZ - 2];
	}
    ll prevNt = -1, prevCurNt = -1;
	VVD Qmat;
	for (int t = tLast; t >= 0; t--)	{
		cur ^= 1;
		prev ^= 1;
		int tIdx = getTIdx(t, times);
		ll Nt = Ns[tIdx];
		ll curNt = Ns[getTIdx(t-1, times)];
	    if (Nt != prevNt || curNt != prevCurNt)   {
			Qmat = calcQWithTruncation(n, Nt, 0., truncIdxes, curNt);	//not using any truncation
		}
		hVarMemoInternal[cur][1] = 0.;
		int aMax = (int)min((ll)n, curNt);
		#pragma omp parallel for
		FOR (a, 2, aMax + 1)	{
			int aStart = truncIdxes[a].first;
			
			DOUBLE ret = a;
			int mMax = (int)min((ll)a, Nt);
			FOR (m, 1, mMax + 1)	{
				ret += Qmat[a][m-aStart] * hVarMemoInternal[prev][m];
			}
			hVarMemoInternal[cur][a] = ret;
		}
        prevNt = Nt;
		prevCurNt = curNt;
	}
	hVarMemo = hVarMemoInternal[cur];
}


void initQAndTruncationMatrices(int n)	{
	Qconst = calcQWithTruncation(n, Nconst, adaptiveEpsConstPop, truncIdxesConstPop, Nconst);
	QconstFullMatrix = calcQ(n, Nconst, Nconst);
}


void freqSpectrumVarPopSize(VI& ns, int maxA)	{
	int nmax = *max_element(ALL(ns));
	
    VD ret;
	goFVarPopSize(nmax, maxA);
	
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

	goHVarPopSize(nmax);
	cout << endl;
	cout << "Probability distribution of # derived alleles at a segregating site, first " << maxA << " entries" << endl;
	REP (it, ns)	{
		int n = *it;
		cout << "n = " << n << endl;
	    FOR (k, 1, min(n - 1, maxA) + 1)   {
			DOUBLE val = fAdaptVarMemo[k][n-k] * binom(n, k) / hVarMemo[n];
			printf("%0.20Lg\t", val);
		}
		cout << endl;
	}
}


// precomputation for the last epoch of constant population size
void freqSpectrum(VI& ns, int maxA)	{
	int nmax = *max_element(ALL(ns));

	goF(nmax, maxA);
	goH(nmax);
}



int main(int argc, char **argv)	{
	assert (argc >= 6);
	modelFile = argv[1];
	adaptiveEpsConstPop = atof(argv[2]);
	adaptiveEpsVarPop = atof(argv[3]);
	int maxA = atoi(argv[4]);

	VI ns;
	FOR (i, 5, argc)	{
		ns.PB(atoi(argv[i]));
	}

	cout << "adaptive truncation eps for constant pop size model = " << adaptiveEpsConstPop << endl;
	cout << "adaptive truncation eps for variable pop size model = " << adaptiveEpsVarPop << endl;
	cout << "sample size(s) n = ";
	REP (it, ns)	cout << *it << ", ";
	cout << endl;
	int nmax = *max_element(ALL(ns));
	cout << "largest sample size = " << nmax << endl;
	
	readWFModelFromFile(modelFile, times, Ns);
	Nconst = Ns[Ns.SZ - 1];
	uMut = 1. / (2. * Nconst);
	
	initQAndTruncationMatrices(nmax);
		
	assert(Ns.SZ == times.SZ);
	
	freqSpectrum(ns, maxA);
		
	cout << endl;
    cout << "Variable population size" << endl;
    cout << "------------------------" << endl;    
    freqSpectrumVarPopSize(ns, maxA);
	
	return 0;
}

#endif
