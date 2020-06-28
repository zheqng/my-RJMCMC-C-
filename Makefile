# You may need to edit this file to reflect the type and capabilities of your system.
# The defaults are for a Linux system and may need to be changed for other systems (eg. Mac OS X).


CXX=g++

#CXX=CC
## When using the Sun Studio compiler


# # flags configured by CMake
# ifeq (unix,macos)
#   LIB_FLAGS = -larmadillo -framework Accelerate
# else
#   LIB_FLAGS = -larmadillo
#   ## NOTE: on Ubuntu and Debian based systems you may need to add -lgfortran
#
#   #LIB_FLAGS = -larmadillo -library=sunperf
#   ## When using the Sun Studio compiler
# endif


LIB_FLAGS= -I /home/zheqng/src/armadillo-9.200.7/include -DARMA_DONT_USE_WRAPPER -lopenblas -llapack
# LIB_FLAGS = -larmadillo

OPT = -O2
OPTOP = -O0
## As the Armadillo library uses recursive templates, compilation times depend on the level of optimisation:
##
## -O0: quick compilation, but the resulting program will be slow
## -O1: good trade-off between compilation time and execution speed
## -O2: produces programs which have almost all possible speedups, but compilation takes longer
## -O3: enables auto vectorisation when using gcc

#OPT = -xO4 -xannotate=no
## When using the Sun Studio compiler


#EXTRA_OPT = -fwhole-program
## Uncomment the above line if you're compiling all source files into one program in a single hit


DEBUG = -DARMA_EXTRA_DEBUG
## Uncomment the above line to enable low-level debugging.
## Lots of debugging information will be printed when a compiled program is run.
## Please enable this option when reporting bugs.


FINAL = -DARMA_NO_DEBUG
## Uncomment the above line to disable Armadillo's checks.
## Not recommended unless your code has been first thoroughly tested!
OMPFLAGS = -fopenmp

#export OMP_NUM_THREADS=8

CXXFLAGSNP = -std=c++14  $(DEBUG) $(FINAL) $(OPTNP) $(EXTRA_OPT)
CXXFLAGS = -std=c++14  $(DEBUG) $(FINAL) $(OPT) $(EXTRA_OPT)

#OBJECTS = rj_mix.o RWMH.o HMC.o mix_lib.o alea.o

# .c.o:
	# $(CXX) $(CXXFLAGS) -c -g $@  $^  $(LIB_FLAGS)

all: rj_mix

rj_mix:  nuts.o HMC.o mix_lib.o alea.o
	$(CXX) $(CXXFLAGS)    -o  $@ $(@).cpp $^   -fopenmp

test_kiss: mix_lib.o alea.o
	$(CXX) $(CXXFLAGS)     -o  $@ $(@).cpp $^

test_omp: mix_lib.o alea.o
	$(CXX) $(CXXFLAGS) $(OMPFLAGS)   -o $@ $(@).cpp $^

nuts.o: nuts.cpp nuts.h
HMC.o: HMC.cpp HMC.h
mix_lib.o:mix_lib.cpp mix_lib.h
alea.o: alea.cpp alea.h

nuts.o:
		$(CXX) $(CXXFLAGS)   -c $(@:.o=.cpp) -o $@

HMC.o:
	$(CXX) $(CXXFLAGS)    -c $(@:.o=.cpp) -o $@     -fopenmp

mix_lib.o:
	$(CXX) $(CXXFLAGS)   -c $(@:.o=.cpp) -o $@    -fopenmp


alea.o:
	$(CXX) $(CXXFLAGS)    -c $(@:.o=.cpp) -o $@   


.PHONY: clean

clean:
	rm -f rj_mix *.o
cleanpdf:
	rm *.pdf *.html generator.res .100 .101

