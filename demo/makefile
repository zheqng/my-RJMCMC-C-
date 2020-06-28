# CT/RJ-Mix: Makefile for running the programs on the demo example
#
# $Id: makefile,v 1.1 2003/11/28 16:50:19 cappe Exp $

#demo : galaxy.dat galaxy.arg ../rj_mix ../ct_mix
#	../rj_mix galaxy
#	../ct_mix galaxy
demo : electric_load_2008.dat electric_load_2008.arg ../rj_mix
	../rj_mix electric_load_2008
#	../ct_mix galaxy

demo_omp: traindata.dat traindata.arg ../test_omp
	../test_omp traindata

demo_omp2: traindata.dat traindata.arg ../test_omp2
		../test_omp2 traindata
		
clean :
	rm -f *.sts *.txt
