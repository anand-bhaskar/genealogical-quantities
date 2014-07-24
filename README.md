### About

Contains several programs to compute various genealogical quantities under Kingman's coalescent and the discrete-time Wright-Fisher models of random mating. These programs along with the given text files were used to generate the results in this paper:

[Distortion of genealogical properties when the sample is very large](http://www.pnas.org/content/111/6/2385)  
Bhaskar A., Clark A.G., Song Y.S.  
PNAS, 111(6):2385–2390, 2014


### Files

- freq_spectrum_hybrid.cpp -- computes the expected frequency spectrum using the WF model for the first few generations, and then switches to the coalescent for more ancient times
- wf_freq_spectrum.cpp -- computes the expected frequency spectrum using the WF model
- coal_freq_spectrum.cpp -- computes the expected frequency spectrum using the coalescent
- coal_exp_lineages.cpp -- computes the expected number of ancestors of a sample as a function of time under the coalescent
- exp_multiple_mergers.cpp -- computes the expected number of k-mergers as a function of time under the WF model
- exp_mergers.cpp -- computes the expected number of k-pairwise-simultaneous mergers and the expected number of lineages lost in each generation under the WF model
- wf_lib.cpp, coal_freq_spectrum_lib.cpp, and common.cpp contain common functionality needed by the other programs

Each of these can be compiled separately by doing:
```
make <filename_without_dot_cpp_extension>
```

The demographic models corresponding to Models 1--4 in the paper are in model[1-4].txt (WF) and coalModel[1-4].txt (coalescent). 
The demographic model files model[1-4].txt are meant for computations under the WF model, and contain 1 line each containing the time and the population size at time points where the population size changes from its value at the previous generation. The first line contains the constant population size in the ancient past.

The demographic model files coalModel[1-4].txt are meant for computations in the coalescent. The first line is the constant population size in the ancient past. Each subsequent line has one of these 2 types depending on whether the epoch is an exponential phase or a constant population size: 
"e", pop_size, onset_time, growth_rate - this says that at time onset_time, the population size instantaneously changes to a value x such that it grows at growth_rate (i.e. growth_rate * 100%) per generation until it becomes pop_size at the next demographic event (or time 0, whichever happens first).
"c", pop_size, onset_time - at time onset_time, the population size instantly becomes pop_size.

All population sizes in the model files are in haploids.


### Dependencies

The computation of the frequency spectrum under the coalescent requires the [boost](http://www.boost.org/users/download/) math library. Change this line in the Makefile to point to the directory where the boost header files are installed.
```
INCL = -I/opt/local/include/
```

The various computations can also take advantage of loop-level parallelism through the OpenMP library. If you have a [GNU C++ compiler that supports OpenMP](https://gcc.gnu.org/wiki/openmp) (&geq; 4.2), you can enable parallelism by uncommenting the OPENMP_FLAG line in the Makefile by deleting the ```#``` symbol in the following line:
```
# OPENMP_FLAG = -fopenmp
```
You do not have to change anything in the source code files to turn on parallelism.