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

// Computes the expected number of k-mergers in each generation in the discrete-time Wright-Fisher model.
//
// Example usage: 
// ./exp_multiple_mergers model3.txt 20 5
// Outputs the expected number of k-mergers (2 <= k <= 5) in each generation for model 3 for sample size 
// 20 using the discrete-time Wright-Fisher model with the demography given in model3.txt


#ifndef __EXP_MULTIPLE_MERGERS__
#define __EXP_MULTIPLE_MERGERS__

#include <omp.h>
#include "wf_lib.cpp"

vector<ll> Ns;
VI times;

ll Nconst;

VVD kMergersConstPop;

void goExpKMergersConstPop(int K, int n)	{
	VVD Qconst = calcQ(n, Nconst, Nconst);
	kMergersConstPop.resize(K + 1, VD(n + 1, 0.));
	FOR (k, 2, K + 1)	{
		FOR (a, k, n + 1)	if (a <= Nconst)	{
			DOUBLE rak = 0.;
			if (a == k)	rak = Qconst[k][1] / (1. - Qconst[a][a]);
			FOR (m, k + 1, a + 1)	if (m - k + 1 <= Nconst)	{
				DOUBLE p = 1. * (Nconst - m + k) / Nconst;
				rak += Qconst[k][1] * Qconst[a - k][m - k] / (1. - Qconst[a][a]) * p;
			}
			DOUBLE ans = rak * binom(a, k);
			FOR (m, 1, a)	if (m <= Nconst)	{
				ans += Qconst[a][m] / (1. - Qconst[a][a]) * kMergersConstPop[k][m];
			}
			kMergersConstPop[k][a] = ans;
		}
	}
}


void goExpKMergersVarPopFunctionOfTime(int K, int n)	{
	VD expCumKMergers(K + 1, 0.);
	VVD probNAncestors(2, VD(n + 1, 0.0));
	int tLast = 1000;	//default number of generations of output if only 1 epoch of constant population size
	if (times.SZ >= 2)	{
		tLast = times[times.SZ - 2];
	}

	int cur = 0, prev = 1;
	probNAncestors[prev][n] = 1.0;
	VVD Qmat;
	ll prevNt = -1, prevCurNt = -1;;
	for (int t = 0; t <= tLast; t++)	{
		cur ^= 1;
		prev ^= 1;
		int tIdx = getTIdx(t, times);
		ll Nt = Ns[tIdx];
		ll curNt = Ns[getTIdx(t-1, times)];
		if (Nt != prevNt || curNt != prevCurNt)	{
			Qmat = calcQ(n, Nt, curNt);
		}
		
		if (t == tLast)	break;

		#pragma omp parallel for
		FOR (k, 2, K + 1)	{
			FOR (a, k, n + 1)	if (a <= curNt)	{
				DOUBLE rak = 0.;
				if (a == k)	rak = Qmat[k][1];
				FOR (m, k + 1, a + 1)	if (m - k + 1 <= Nt)	{
					DOUBLE p = 1. * (Nt - m + k) / Nt;
					rak += Qmat[k][1] * Qmat[a - k][m - k] * p;
				}
				expCumKMergers[k] += probNAncestors[cur][a] * rak * binom(a, k);				
			}
		}
		printf("%d", t);
		FOR (k, 2, K + 1)	{
			printf("\t%0.12Lg", expCumKMergers[k]);
		}
		printf("\n");
				
		#pragma omp parallel for
		FOR (j, 1, n + 1)	if (j <= Nt)	{
			DOUBLE ans = 0.;
			FOR (m, j, n + 1)	if (m <= curNt)	{
				ans += probNAncestors[cur][m] * Qmat[m][j];
			}
			probNAncestors[prev][j] = ans;
		}
		
		prevNt = Nt;
		prevCurNt = curNt;
	}
	
	#pragma omp parallel for
	FOR (k, 2, K + 1)	{
		FOR (a, k, n + 1)	if (a <= Nconst)	{
			expCumKMergers[k] += probNAncestors[cur][a] * kMergersConstPop[k][a];
		}
	}
	printf("%d", tLast);
	FOR (k, 2, K + 1)	{
		printf("\t%0.12Lg", expCumKMergers[k]);
	}
	printf("\n");
}


int main(int argc, char **argv)	{
	assert (argc >= 4);
	
	string modelFile = argv[1];
	int n = atoi(argv[2]);
	int K = atoi(argv[3]);
	
	cout << "sample size n = " << n << endl;
	
	readWFModelFromFile(modelFile, times, Ns);
	Nconst = Ns[Ns.SZ - 1];

	goExpKMergersConstPop(K, n);
	
	goExpKMergersVarPopFunctionOfTime(K, n);
}


#endif
