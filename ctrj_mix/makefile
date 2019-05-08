# CT/RJ-Mix: Makefile for rj_mix and ct_mix
#
# $Id: makefile,v 3.3 2003/12/02 12:54:05 cappe Exp $

# May need to be modified to fit your local compiler configuration
CC=cc

# This is the default build for rj_mix and ct_mix
all : rj_mix.c mix_lib.c mix_lib.h alea.c alea.h
	${CC}  -O -o rj_mix rj_mix.c mix_lib.c alea.c -lm
#	${CC} -O -o ct_mix ct_mix.c mix_lib.c alea.c -lm

# Targets for specific uses
debug : rj_mix.c ct_mix.c mix_lib.c mix_lib.h alea.c alea.h
	${CC} -g -o rj_mix rj_mix.c mix_lib.c alea.c -lm
	${CC} -g -o ct_mix ct_mix.c mix_lib.c alea.c -lm

mix_lib.o: mix_lib.c mix_lib.h

alea.o: alea.c alea.h

mix_lib.o alea.o :
	${CC} -c -o $@ $(@:.o=.c)

clean :
	rm -f *.o *.out rj_mix ct_mix
