# You may need to edit this file to reflect the type and capabilities of your system.
# The defaults are for a Linux system and may need to be changed for other systems (eg. Mac OS X).


CXX=g++

#CXX=CC
## When using the Sun Studio compiler


# flags configured by CMake
ifeq (unix,macos)
  LIB_FLAGS = -larmadillo -framework Accelerate
else
  LIB_FLAGS = -larmadillo
  ## NOTE: on Ubuntu and Debian based systems you may need to add -lgfortran

  #LIB_FLAGS = -larmadillo -library=sunperf
  ## When using the Sun Studio compiler
endif





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




CXXFLAGSNP = -std=c++14  $(DEBUG) $(FINAL) $(OPTNP) $(EXTRA_OPT)
CXXFLAGS = -std=c++14  $(DEBUG) $(FINAL) $(OPT) $(EXTRA_OPT)

#OBJECTS = rj_mix.o RWMH.o HMC.o mix_lib.o alea.o

# .c.o:
	# $(CXX) $(CXXFLAGS) -c -g $@  $^  $(LIB_FLAGS)

all: rj_mix

rj_mix:  nuts.o HMC.o mix_lib.o alea.o
	$(CXX) $(CXXFLAGSNP)  -g   -o  $@ $(@).cpp $^  $(LIB_FLAGS)

# rj_mix.o:rj_mix.cpp RWMH.h
# RWMH.o: RWMH.cpp RWMH.h HMC.h
nuts.o: nuts.cpp nuts.h
HMC.o: HMC.cpp HMC.h
mix_lib.o:mix_lib.cpp mix_lib.h
alea.o: alea.cpp alea.h

nuts.o:
		$(CXX) $(CXXFLAGSNP) -g   -c $(@:.o=.cpp) -o $@    $(LIB_FLAGS)

 HMC.o mix_lib.o  alea.o:
	$(CXX) $(CXXFLAGSNP)   -g -c $(@:.o=.cpp) -o $@    $(LIB_FLAGS)

# demo : traindata.dat traindata.arg ../rj_mix
	# ../rj_mix traindata
#	../ct_mix galaxy


# clean :
	# rm -f *.sts *.txt 

.PHONY: clean

clean:
	rm -f rj_mix *.o
