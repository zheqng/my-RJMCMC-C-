/************************************************************************/
/* CT/RJ-Mix Transdimensional MCMC for Gaussian mixtures		*/
/*									*/
/* As discussed in section 4 of						*/
/* O. Capp�, C. Robert, and T. Ryd�n. Reversible jump, birth-and-death	*/
/* and more general continuous time Markov chain Monte Carlo samplers.	*/
/* J. Royal Statist. Soc. Ser. B, 65(3):679-700, 2003.			*/
/* Further information at http://www.tsi.enst.fr/~cappe/ctrj_mix/      */
/*									*/
/* Version: $Id: mix_lib.h,v 3.9 2003/12/02 12:54:05 cappe Exp $	*/
/*									*/
/*									*/
/* Copyright (C) 2003, Olivier Capp�, Tobias Ryd�n, Christian P. Robert	*/
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
/* mix_lib.h: Header file for mix_lib, rj_mix and ct_mix		*/
/************************************************************************/

/************************************************************************/
/* Below are definitions you may want to CUSTOMIZE			*/
/************************************************************************/
/* Maximal value of n (this is not the actual value of n)		*/
#define MaxN 2000
/* Maximal value of k (this is not the actual value of k)		*/
#define MaxK 20
/* Idem for combinations of pairs					*/
/* (should be larger than MAxK(MAxK-1)/2 + MaxK + 3)			*/
#define MaxKPairs 215
/************************************************************************/
/* End CUSTOMIZABLE section						*/
/************************************************************************/


#include <stdio.h>
#include <math.h>
#include "alea.h"
#include <omp.h>
/************************************************************************/
/* Useful functions aliases and definitions				*/
/************************************************************************/
#define square(x) ((x)*(x))
/* Minimum acceptable negative value (before -inf)			*/
#define NEG_REAL_MIN -1.79e+308

/************************************************************************/
/* Data types								*/
/************************************************************************/
typedef struct {
        int k;
        double w[MaxK];
        double mu[MaxK];
        double var[MaxK];
} mix;

typedef struct {
        int l;
        double dat[MaxN];
} vec;

/************************************************************************/
/* Global variables							*/
/************************************************************************/
/* Seeds for the kiss generator (see alea)				*/
extern unsigned long int SK0, SK1, SK2;
/* Model and sampler parameters (defined in main program)		*/
extern int M;
extern double Kappa, SqrtKappa, AlphaVar, BetaVar;
extern double SqrtEta, Rho, SqrtRho, SqrtNu;
extern double Gamma_S, Rho_S, SqrtRho_S, Nu_S, SqrtNu_S, SplitLogConstFact;

/************************************************************************/
/* Routines								*/
/************************************************************************/
double log_likelihood(vec* Y, mix* m);
double log_likelihood_micro(vec* Y, double* w, double *mu, double *var, int k);
void likelihood_individual_components(vec* Y, mix* m, double gpdf[MaxN][MaxK]);
short MNRWMH_weights(vec* Y, mix* m, double* logl);
short RWMH_means(vec* Y, mix* m, double* logl);
short MRWMH_variances(vec* Y, mix* m, double* logl);
void mixcpy(mix* m_src, mix* m_dest);
void mixprint(mix* m);
double prop_birth(vec* Y, mix* m);
double prop_death(vec* Y, mix* m_old, mix* m_new, int k);
int invert_disc_df(double* p, int l, double ran);
mix draw_from_prior(void);
int draw_initial_model(vec* Y, mix* m, double* logl);
double gammaln(double x);
double prop_split(vec* Y, mix* m, int k);
double prop_merge(vec* Y, mix* m_old, mix* m_new, int k1, int k2);
double compute_log_split_ratio(mix* m, mix* m_split, int k, int k_split1, int k_split2);
