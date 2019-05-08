# CT/RJ-Mix: Makefile for running the programs on the demo example
#
# $Id: makefile,v 1.1 2003/11/28 16:50:19 cappe Exp $

demo : galaxy.dat galaxy.arg ../rj_mix 
	../rj_mix galaxy
#	../ct_mix galaxy

check : galaxy.sts galaxy.st2
	md5sum galaxy.sts
	md5sum galaxy.st2

clean : galaxy.sts galaxy.res galaxy.st2 galaxy.rs2
	rm -f galaxy.sts galaxy.res galaxy.st2 galaxy.rs2
