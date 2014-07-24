CC = g++

# uncomment the line below to turn on OpenMP support
# OPENMP_FLAG = -fopenmp

CCFLAGS_NO_OFLAG = $(OPENMP_FLAG) -m64 -Wall -Wno-unknown-pragmas
OFLAG = -O3
CCFLAGS = $(CCFLAGS_NO_OFLAG) $(OFLAG)

INCL = -I/opt/local/include/

all: coal_exp_lineages coal_freq_spectrum exp_mergers exp_multiple_mergers freq_spectrum_hybrid wf_freq_spectrum

clean:
	rm -f *.o coal_exp_lineages coal_freq_spectrum exp_mergers exp_multiple_mergers freq_spectrum_hybrid wf_freq_spectrum

coal_freq_spectrum_lib.o: coal_freq_spectrum_lib.cpp common.cpp
	$(CC) $(CCFLAGS) $(INCL) -c $<

wf_lib.o: wf_lib.cpp common.cpp
	$(CC) $(CCFLAGS) $(INCL) -c $<
	
coal_exp_lineages: coal_exp_lineages.cpp coal_freq_spectrum_lib.o
	$(CC) $(CCFLAGS) $(INCL) $< -o $@

coal_freq_spectrum: coal_freq_spectrum.cpp coal_freq_spectrum_lib.o
	$(CC) $(CCFLAGS) $(INCL) $< -o $@

exp_mergers: exp_mergers.cpp wf_lib.o
	$(CC) $(CCFLAGS) $(INCL) $< -o $@

#there seems to be a compiler bug (at least with gcc 4.5.0) when using -O2 or -O3
exp_multiple_mergers: exp_multiple_mergers.cpp wf_lib.o
	$(CC) $(CCFLAGS_NO_OFLAG) $(INCL) $< -o $@

freq_spectrum_hybrid: freq_spectrum_hybrid.cpp coal_freq_spectrum_lib.o wf_lib.o
	$(CC) $(CCFLAGS) $(INCL) $< -o $@

wf_freq_spectrum: wf_freq_spectrum.cpp wf_lib.o
	$(CC) $(CCFLAGS) $(INCL) $< -o $@
