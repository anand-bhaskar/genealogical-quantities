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

// Computes the expected number of k-drops and k-pairwise-simultaneous mergers (1 <= k <= 9) per epoch of the discrete-time
// Wright-Fisher model
//
// Example usage: 
// ./exp_mergers model3.txt 20 9
// Computes the expected number of k-drops and k-pairwise-simultaneous mergers (1 <= k <= 9) per epoch of the Wright-Fisher 
// model model3.txt for a sample of size 100.


#ifndef __EXP_MERGERS__
#define __EXP_MERGERS__

#include <omp.h>
#include "wf_lib.cpp"
#include "common.cpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <limits>
using namespace std;


bool precomputeQ = false;
bool printEveryGeneration = true;
VVD W, X, Y;
vector<VVD> Qmap, Rmap, Smap;

vector<ll> Ns;
VI times;

ll Nconst;

VVD kDropsConstPop, kDropsSim2MergersConstPop;

//calculates expected number of mergers which reduce the sample by exactly K, assuming a constant pop size
void goExpKDropsConstPop(int K, int n)	{
	VVD Qconst = calcQ(n, Nconst, Nconst);
	VVD Rconst = calcR(n, Nconst, Nconst);
	
	kDropsConstPop.resize(K + 1, VD(n + 1, 0.0));
	kDropsSim2MergersConstPop.resize(K + 1, VD(n + 1, 0.0));

	#pragma omp parallel for
	FOR (k, 1, K + 1)	{
		FOR (m, k + 1, n + 1)	if (m <= Nconst)	{
			DOUBLE ans = Qconst[m][m - k] / (1. - Qconst[m][m]);
			DOUBLE ans2 = Rconst[m][m - k] / (1. - Qconst[m][m]);
			FOR (j, k + 1, m)	{
				ans += Qconst[m][j] / (1. - Qconst[m][m]) * kDropsConstPop[k][j];
				ans2 += Qconst[m][j] / (1. - Qconst[m][m]) * kDropsSim2MergersConstPop[k][j];
			}
			kDropsConstPop[k][m] = ans;
			kDropsSim2MergersConstPop[k][m] = ans2;
		}
	}
}


void goPerEpochStatisticsVarPop(int printK, int n)	{
	int K = n - 1;
	int tLast = 1000;	//default number of generations of output if only 1 epoch of constant population size
	if (times.SZ >= 2)	{
		tLast = times[times.SZ - 2];
	}
	
	int cur = 0, prev = 1;
	VVD probNAncestors(2, VD(n + 1, 0.0));
	
	VD kDropsEpoch(K + 1, 0.0);
	VD kDropsSim2MergersEpoch(K + 1, 0.0);

	VD kDropsAllEpochs(K + 1, 0.0);
	VD kDropsSim2MergersAllEpochs(K + 1, 0.0);
	
	DOUBLE grandTot = 0., grandTot2 = 0.;
	int lenEpoch = 0, epochIdx = 0;
	
	probNAncestors[cur][n] = 1.0;
    DOUBLE expNInds = n;
	DOUBLE expSqNInds = n*n;
	VVD Qmat, Rmat;
	ll prevNt = -1, prevCurNt = -1;
	for (int t = 0; t <= tLast; t++)	{
		int tIdx = getTIdx(t, times);
		ll Nt = Ns[tIdx];
		ll curNt = Ns[getTIdx(t - 1, times)];
		if (Nt != prevNt || curNt != prevCurNt)	{
			Qmat = calcQ(n, Nt, curNt);
		}

		ll prevEpochNt = Ns[getTIdx(t - 1, times)];
		if (Nt != prevEpochNt || (t != 0 && printEveryGeneration))	{	//epoch has changed
			//print
			DOUBLE tot = 0., tot2 = 0.;
			printf("Epoch index = %d, duration = %d generations, population size = %lld\n", epochIdx, lenEpoch, prevEpochNt);
            printf("Expected sample size at start of epoch = %#18.12Lg\n", expNInds);
			printf("Standard deviation in sample size at start of epoch = %#18.12Lg\n", sqrt(expSqNInds - expNInds * expNInds));
			FOR (k, 1, K + 1)   {
				if (k <= printK)	{
					printf("%6d\t%#20.12Lg\t%#20.12Lg\t%#20.12Lg\n", k, kDropsEpoch[k], kDropsSim2MergersEpoch[k], (kDropsSim2MergersEpoch[k] / kDropsEpoch[k]));
					
				}
				tot += kDropsEpoch[k] * k;
				tot2 += kDropsSim2MergersEpoch[k] * k;
				
				kDropsAllEpochs[k] += kDropsEpoch[k];
				kDropsSim2MergersAllEpochs[k] += kDropsSim2MergersEpoch[k];
			}
			grandTot += tot;
			grandTot2 += tot2;

			printf("sum k*E(k-drops in this epoch) = %#18.12Lg\t%#20.12Lg\t%#20.12Lg\n", tot, tot2, (tot2 / tot));
			printf("sum k*E(k-drops until this epoch) = %#18.12Lg\t%#20.12Lg\t%#20.12Lg\n", grandTot, grandTot2, (grandTot2 / grandTot));
			printf("\n");
			fflush(stdout);
			
			fill(ALL(kDropsEpoch), 0.0);
			fill(ALL(kDropsSim2MergersEpoch), 0.0);

			epochIdx++;
			lenEpoch = 0;
			
            expNInds = 0.;
			expSqNInds = 0.;
			FOR (m, 1, n + 1)	{
                expNInds += probNAncestors[cur][m] * m;
				expSqNInds += probNAncestors[cur][m] * m * m;
            }
		}
		
		if (t == tLast)	break;

		//add to the expected number of k drops in this epoch
		#pragma omp parallel for
		FOR (k, 1, K + 1)	{
			DOUBLE ans = 0., ans2 = 0.;
			DOUBLE Rmat_m_mMinusk = 1.;
			
			FOR (i, 1, 2*k)	{
				Rmat_m_mMinusk = Rmat_m_mMinusk * i / Nt;
				if (i <= k)	{
					Rmat_m_mMinusk = Rmat_m_mMinusk / 2. / i;
				}
				if (i <= k-1)	{
					Rmat_m_mMinusk = Rmat_m_mMinusk * (Nt - i + 1);
				}
			}
			
			//Rmat_m_mMinusk = (2k-1)! / k! / 2^k \fall{N}_{k-1} / N^{2k-1}
			FOR (m, k + 1, n + 1)	if (m <= curNt && m - k <= Nt)	{
				ans += probNAncestors[cur][m] * Qmat[m][m - k];
				//Rmat_m_mMinusk = m! / (m-2k)! / k! / 2^k * \fall{Nt}{m-k} / N^m
				if (m >= 2*k)	{
					Rmat_m_mMinusk = Rmat_m_mMinusk * m * (Nt - (m - k) + 1) / Nt;
					if (m > 2*k)	Rmat_m_mMinusk = Rmat_m_mMinusk / (m - 2*k);

					ans2 += probNAncestors[cur][m] * Rmat_m_mMinusk;
				}
			}
			kDropsEpoch[k] += ans;
			kDropsSim2MergersEpoch[k] += ans2;
		}

		cur ^= 1;
		prev ^= 1;

		//calculate the prob of n ancestors at time t + 1
		#pragma omp parallel for
		FOR (j, 1, n + 1)	if (j <= Nt)	{
			DOUBLE ans = 0.;
			FOR (m, j, n + 1)	if (m <= curNt)	{
				ans += probNAncestors[prev][m] * Qmat[m][j];
			}
			probNAncestors[cur][j] = ans;
		}
		
		prevNt = Nt;
		prevCurNt = curNt;
		
		lenEpoch++;
	}
	
	//constant population size epoch.	probNAncestors[cur] is the current probability of number of ancestors.
	DOUBLE tot = 0., tot2 = 0.;
	printf("Last epoch with constant population size infinitely far back in time, population size = %lld\n", Nconst);
    printf("Expected sample size at start of epoch = %#18.12Lg\n", expNInds);
	printf("Standard deviation in sample size at start of epoch = %#18.12Lg\n", sqrt(expSqNInds - expNInds * expNInds));
	FOR (k, 1, printK + 1)	{
		DOUBLE ans = 0., ans2 = 0.;
		FOR (m, k + 1, n + 1)	{
			ans += probNAncestors[cur][m] * kDropsConstPop[k][m];
			ans2 += probNAncestors[cur][m] * kDropsSim2MergersConstPop[k][m];
		}
		kDropsAllEpochs[k] += ans;
		kDropsSim2MergersAllEpochs[k] += ans2;

		printf("%6d\t%#20.12Lg\t%#20.12Lg\t%#20.12Lg\n", k, ans, ans2, (ans2 / ans));

		tot += ans * k;
		tot2 += ans2 * k;
	}
	grandTot += tot;
	grandTot2 += tot2;

	printf("sum k*E(k-drops in this epoch) = %#18.12Lg\t%#20.12Lg\t%#20.12Lg\n", tot, tot2, (tot2 / tot));
	printf("sum k*E(k-drops until this epoch) = %#18.12Lg\t%#20.12Lg\t%#20.12Lg\n", grandTot, grandTot2, (grandTot2 / grandTot));
	printf("\n\n");
	
	cout << "Variable population scenario" << endl;
	cout << "----------------------------" << endl;
	cout << "# The 6 columns are: k, expected number of k-drops due to all kinds of mergers, \
expected number of k-drops due to only simultaneous mergers (i.e. no multiple mergers), \
proportion of the k-drops due to simultaneous mergers, \
expected number of k-drops due to (k+1)-mergers, \
proportion of the k-drops due to (k+1)-mergers" << endl;
	tot = 0., tot2 = 0.;
	FOR (k, 1, K + 1)	{
		if (k <= printK)	{
			printf("%6d\t%#20.12Lg\t%#20.12Lg\t%#20.12Lg\n", k, kDropsAllEpochs[k], kDropsSim2MergersAllEpochs[k], (kDropsSim2MergersAllEpochs[k] / kDropsAllEpochs[k]));
		}
		tot += kDropsAllEpochs[k] * k;
		tot2 += kDropsSim2MergersAllEpochs[k] * k;
	}
	printf("\nsum k*E(k-drops) = %#18.12Lg\t%#20.12Lg\t%#20.12Lg\n", tot, tot2, (tot2 / tot));
}




void readModelFromFile(string& fileName)	{
	ifstream fin(fileName.c_str());
	string buff;
	vector<string> arr;
	bool readFirstLine = true;
	cout << "# Reading model from file " << fileName << endl;
	while (!fin.eof())	{
		getline(fin, buff);
		if (buff[0] == '#')	continue;
		if (buff.SZ == 0)	break;
		if (readFirstLine)	{
			Nconst = atol(buff.c_str());	//ancestral population size
			// cout << "# Ancestral population size = " << Nconst << endl;
			readFirstLine = false;
		}
		else	{
			splitstr(buff, arr);	//time	population_size
			times.PB(atoi(arr[0].c_str()));
			Ns.PB(atol(arr[1].c_str()));
			// cout << "# " << times[times.SZ - 1] << "\t" << Ns[Ns.SZ - 1] << endl;
		}
	}
	times.PB(numeric_limits<int>::max());	Ns.PB(Nconst);	
}



int main(int argc, char **argv)	{
	assert (argc >= 4);
	
	string modelFile = argv[1];
	int n = atoi(argv[2]);
	int K = atoi(argv[3]);

	cout << "sample size n = " << n << endl;

	readModelFromFile(modelFile);

	//some precomputation needed for the last epoch back in time having a constant population size
	goExpKDropsConstPop(K, n);
	
    cout << endl;
    cout << endl;
    cout << "Statistics for variable population scenario" << endl;
    cout << "-------------------------------------------" << endl;
    cout << "# The 6 columns are: k, expected number of k-drops due to all kinds of mergers in this epoch of constant population size, \
expected number of k-drops due to only simultaneous mergers (i.e. no multiple mergers) in this epoch, \
proportion of the k-drops due to simultaneous mergers in this epoch" << endl;
	goPerEpochStatisticsVarPop(K, n);
    cout << endl;

	return 0;
}

#endif