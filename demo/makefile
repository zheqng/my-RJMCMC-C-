# CT/RJ-Mix: Makefile for running the programs on the demo example
#
# $Id: makefile,v 1.1 2003/11/28 16:50:19 cappe Exp $

#demo : galaxy.dat galaxy.arg ../rj_mix ../ct_mix
#	../rj_mix galaxy
#	../ct_mix galaxy
demo : traindata.dat traindata.arg ../rj_mix
	../rj_mix traindata
#	../ct_mix galaxy

demo_test : traindata.dat traindata.arg ../test_kiss
	../test_kiss traindata

clean :
	rm -f *.sts *.txt 
