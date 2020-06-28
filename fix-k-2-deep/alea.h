/************************************************************************/
/* CT/RJ-Mix Transdimensional MCMC for Gaussian mixtures		*/
/*									*/
/* As discussed in section 4 of						*/
/* O. Capp? C. Robert, and T. Ryd�n. Reversible jump, birth-and-death	*/
/* and more general continuous time Markov chain Monte Carlo samplers.	*/
/* J. Royal Statist. Soc. Ser. B, 65(3):679-700, 2003.			*/
/* Further information at http://www.tsi.enst.fr/~cappe/ctrj_mix/       */
/*									*/
/* Version: $Id: alea.h,v 1.8 2005/10/06 07:52:36 cappe Exp $		*/
/*									*/
/*									*/
/* Copyright (C) 2003, Olivier Capp? Tobias Ryd�n, Christian P. Robert	*/
/*									*/
/* This program is free software; you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation.					*/
/*									*/
/* This program is distributed in the hope that it will be useful,	*/
/* but WITHOUT ANY WARRANTY; without even the implied warranty of	*/
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	*/
/* GNU General Public License for more details.				*/
/************************************************************************/

/************************************************************************/
/* alea.c: Library of simulation functions by Christian P. Robert,	*/
/*         this is the header file					*/
/************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <stdio.h>
// #include <armadillo>
//#include "lapacke.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>

using namespace boost::random;

using namespace std;
// using namespace arma;
/* How long are unsigned long int? (WE ARE ASSUMING THAT THE ANSWER IS	*/
/* 4 BYTES, if not this would be a problem for the kiss generator, the  */
/* 8 bytes variant has never been tested!)				*/
#if __WORDSIZE == 64
#  define ULONG_NORM  18446744073709551616.0
#else
#  define ULONG_NORM  4294967296.0
#endif

typedef short boolean;
#define FALSE 0
#define TRUE 1

#define PI 3.14159265358979


struct Generator {
        unsigned long int i,j,k;
public:
        Generator& operator=(const Generator &g0){
                if(this == &g0)
                        return *this;

                i = g0.i;
                j = g0.j;
                k = g0.k;

                return *this;
        }
        void print(string s){
                cout<<s<<i<<" "<<j<<" "<<k<<endl;
        }
        void write_file(string s,ofstream & myfile){
                myfile<<s<<i<<" "<<j<<" "<<k<<endl;
        }
};
extern Generator g;


double kiss(Generator & g);
void nor(double* t1, double*t2, Generator & g0);
double Gam_1(double alf, Generator & g0);
double Gam_2(double alf, Generator & g0);
double Gam_3(double alf, Generator & g0);
double Gam(double alf, Generator & g0);
double Beta(double alf, double bet,Generator & g0);
