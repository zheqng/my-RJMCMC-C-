/************************************************************************/
/* CT/RJ-Mix Transdimensional MCMC for Gaussian mixtures		*/
/*									*/
/* As discussed in section 4 of						*/
/* O. Cappé, C. Robert, and T. Rydén. Reversible jump, birth-and-death	*/
/* and more general continuous time Markov chain Monte Carlo samplers.	*/
/* J. Royal Statist. Soc. Ser. B, 65(3):679-700, 2003.			*/
/* Further information at http://www.tsi.enst.fr/~cappe/ctrj_mix/      */
/*									*/
/* Version: $Id: ct_mix.c,v 2.6 2003/12/02 12:54:05 cappe Exp $		*/
/*									*/
/*									*/
/* Copyright (C) 2003, Olivier Cappé, Tobias Rydén, Christian P. Robert	*/
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
/* ct_mix: Implementation of continuous time MCMC for mixtures		*/
/*									*/
/* Ouput formats:							*/
/*									*/
/* Stat. file (".st2") is a text file with, on each line:		*/
/* | n iter | k | log-l | move type |					*/
/*				    |					*/
/*       for type 1 (fixed k):      | base2(accpt w,mu,var) | weight |  */
/*       for type 2 (birth/death):  | birth                 | weight |  */
/*       for type 3 (split/merge):  | split                 | weight |  */
/*									*/
/* Res. file (".rs2") is a binary file which contains doubles:		*/
/* weight(0), w(0,1), ..., w(0,k(0)), mu(0,1), ..., mu(0,k(0)),		*/
/* var(0,1), ..., var(0,k(0)), weight(SubSamp), w(SubSamp,1),...,	*/
/* w(SubSamp,k(SubSamp)), mu(SubSamp,1), ...				*/
/*									*/
/* Written on stdout is a copy of the input argument (".arg" file) with */
/* updated values of the seeds of the random generator (SK1, SK2, SK3)  */
/*									*/
/*									*/
/* NB: move type, accpt w, accpt mu, accpt var, birth ans plit refer to	*/
/* iteration "n iter" (it is arbitrary for n iter = 0 which is the	*/
/* inial value). All other values (k, log-l, weight) refer to the value */
/* of the parameters after "n iter" iterations.	 			*/
/************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mix_lib.h"

/* Definitions								*/
#define NArgs 20
#define NStat 2
#define StrLen 160
#define ArgExt ".arg"
#define ResExt ".rs2"
#define StatExt ".st2"

/* Global variables							*/
/* Seeds for the kiss generator (see alea)				*/
unsigned long int SK0, SK1, SK2;
/* Files								*/
char DataFile[StrLen], ResFile[StrLen], StatFile[StrLen];
FILE *ResFP, *StatFP;
/* Prior hyperparameters						*/
int M;
double Kappa, SqrtKappa, AlphaVar, BetaVar;
/* Sampler settings							*/
int NOut, SubSamp, NIt;
/* Fixed k move								*/
double  Eta, SqrtEta, Rho, SqrtRho, Nu, SqrtNu, PFixed;
/* Birth and death (PDeath is not used and is here only for		*/
/* compatibility with rj_mix                                            */
double PBirth, PDeath;
/* Split and merge							*/
double PSplit, Gamma_S, Rho_S, SqrtRho_S, Nu_S, SqrtNu_S, SplitLogConstFact;

/* Function prototypes							*/
int main(int argc, char** argv);
void read_parameters(int argc, char** argv, vec* Y);
void write_data(mix* theta, double wght, double logl, short* stats, int in);
double compute_rates(vec* Y, mix* m, double logl, double* move_prob,
  mix m_dead[MaxK], double* logl_dead, mix m_merge[MaxKPairs],
  double* logl_merge);
void compute_death_rates(vec* Y, mix* m, double logl, double* rates,
  mix m_dead[MaxK], double* logl_dead, double gpdf[MaxN][MaxK]);
void compute_merge_rates(vec* Y, mix* m, double logl, double* rates, mix m_merge[MaxKPairs], double* logl_merge, double gpdf[MaxN][MaxK]);

/************************************************************************/
/* Main									*/
/************************************************************************/
int main(int argc, char** argv) {
  vec Y;
  mix theta, theta_dead[MaxK], theta_merge[MaxKPairs];
  double logl, wght, move_prob[MaxKPairs], logl_dead[MaxK], logl_merge[MaxKPairs];
  double stats_tmp1, stats_tmp2, stats_tmp3;
  short stats[NStat];
  int in, k, move_type, nr;

  read_parameters(argc, argv, &Y);
  /* Initial value							*/
  draw_initial_model(&Y, &theta, &logl);
  /* Compute weight and move type					*/
  wght = compute_rates(&Y, &theta, logl, move_prob, theta_dead, logl_dead,
    theta_merge, logl_merge);
  nr = 3 + theta.k + (theta.k*(theta.k-1))/2;
  move_type = invert_disc_df(move_prob, nr, kiss(&SK0,&SK1,&SK2));
  write_data(&theta, wght, logl, stats, 0);
  /* Main loop								*/
  for (in=1; in<=NIt; in++) {
    if (move_type == 0) {
      /* Fixed k move							*/
      stats[0] = 1;
      stats_tmp1 = MNRWMH_weights(&Y, &theta, &logl);
      stats_tmp2 = RWMH_means(&Y, &theta, &logl);
      stats_tmp3 = MRWMH_variances(&Y, &theta, &logl);
      stats[1] = 4*stats_tmp1+2*stats_tmp2+stats_tmp3;
    }
    else if (move_type == 1){
      /* Birth								*/
      stats[0] = 2; stats[1] = 1;
      logl = prop_birth(&Y, &theta);
    }
    else if (move_type <= theta.k+1) {     
      /* Death								*/
      stats[0] = 2; stats[1] = 0;
      k = move_type - 2;
      mixcpy(&(theta_dead[k]), &theta);
      logl = logl_dead[k];
    }
    else if (move_type == theta.k+2) {
      /* Split								*/
      stats[0] = 3; stats[1] = 1;
      k = (int)floor((double)theta.k * kiss(&SK0,&SK1,&SK2));
      logl = prop_split(&Y, &theta, k);
    }
    else {
      /* Merge								*/
      stats[0] = 3; stats[1] = 0;
      k = move_type - (theta.k+3);
      mixcpy(&(theta_merge[k]), &theta);
      logl = logl_merge[k];
    }
    /* We can safely avoid recomputing the rates and weight if		*/
    /* everything was rejected on a fixed k move			*/
    if (!((stats[0] == 1) && (stats[1] == 0))) {
      wght = compute_rates(&Y, &theta, logl, move_prob, theta_dead, logl_dead,
        theta_merge, logl_merge);
    }
    nr = 3 + theta.k + (theta.k*(theta.k-1))/2;
    move_type = invert_disc_df(move_prob, nr, kiss(&SK0,&SK1,&SK2));
    /* Note: The test below is true for in = 0 and in = NIt		*/
    if ((div(in, SubSamp)).rem == 0)
      write_data(&theta, wght, logl, stats, in);
    }
  return(0);
}

/************************************************************************/
/* Read parameters (including seeds for the random generator) and data  */
/************************************************************************/
void read_parameters(int argc, char** argv, vec* Y) {
  char argfile[StrLen];
  int  argflag[NArgs], i, n;
  double y;
  FILE *argfp, *datafp;

  if (argc < 2) {
    fprintf(stderr, "ct_mix base_name (without .arg extension)\n");
    exit(0);
  }
  else {
    strcpy(ResFile, argv[1]);
    strcpy(argfile, argv[1]);
    strcat(argfile, ArgExt);
  }
  argfp = fopen(argfile, "r");
  if (argfp == NULL) {
    fprintf(stderr, "Can't open %s\n", argfile);
    exit(1);
  }
  else {
    i = -1;
    argflag[++i] = fscanf(argfp, "SK0=%lu\n", &SK0);
    argflag[++i] = fscanf(argfp, "SK1=%lu\n", &SK1);
    argflag[++i] = fscanf(argfp, "SK2=%lu\n", &SK2);
    argflag[++i] = fscanf(argfp, "NOut=%d\n", &NOut);
    argflag[++i] = fscanf(argfp, "SubSamp=%d\n", &SubSamp);
    argflag[++i] = fscanf(argfp, "Data=%s\n", DataFile);
    argflag[++i] = fscanf(argfp, "M=%d\n", &M);
    argflag[++i] = fscanf(argfp, "Kappa=%lf\n", &Kappa);
    argflag[++i] = fscanf(argfp, "AlphaVar=%lf\n", &AlphaVar);
    argflag[++i] = fscanf(argfp, "BetaVar=%lf\n", &BetaVar);
    argflag[++i] = fscanf(argfp, "Eta=%lf\n", &Eta);
    argflag[++i] = fscanf(argfp, "Rho=%lf\n", &Rho);
    argflag[++i] = fscanf(argfp, "Nu=%lf\n", &Nu);
    argflag[++i] = fscanf(argfp, "PFixed=%lf\n", &PFixed);
    argflag[++i] = fscanf(argfp, "PBirth=%lf\n", &PBirth);
    /* Note: PDeath will not be used and is here only to maintain	*/
    /* compatibility with the rj_mix arg file				*/
    argflag[++i] = fscanf(argfp, "PDeath=%lf\n", &PDeath);
    argflag[++i] = fscanf(argfp, "PSplit=%lf\n", &PSplit);
    argflag[++i] = fscanf(argfp, "Gamma_S=%lf\n", &Gamma_S);
    argflag[++i] = fscanf(argfp, "Rho_S=%lf\n", &Rho_S);
    argflag[++i] = fscanf(argfp, "Nu_S=%lf\n", &Nu_S);
  }
  fclose(argfp);
  /* Minimal check							*/
  for (i=0; i<NArgs; i++) {
    if (argflag[i] != 1) {
      fprintf(stderr, "Argument number %d could not be read\n", i);
      exit(1);
    }
  }
  if (M > MaxK) {
    fprintf(stderr, "Too many classes, aborting\n");
    exit(1);
  }
  /* Compute frequently used quantities					*/
  NIt = NOut*SubSamp;
  SqrtKappa = sqrt(Kappa);
  SqrtRho = sqrt(Rho);
  SqrtNu = sqrt(Nu);
  SqrtEta = sqrt(Eta);
  SqrtRho_S = sqrt(Rho_S);
  SqrtNu_S = sqrt(Nu_S);
  SplitLogConstFact = 2.0*gammaln(Gamma_S)-gammaln(2.0*Gamma_S)+log(4.0)
    +0.5*log(Rho_S/Kappa*2.0*PI*Nu_S)+AlphaVar*log(BetaVar)-gammaln(AlphaVar);
  /* Read data								*/
  datafp = fopen(DataFile, "r");
  if (datafp == NULL) {
    fprintf(stderr, "Can't open input file %s\n", DataFile);
    exit(1);
  }
  n = -1;
  while (fscanf(datafp,"%lf",&y)!=EOF) {
    Y->dat[++n] = y;
    if (n >= MaxN-1) {
      fprintf(stderr, "Data is too long, aborting\n");
      exit(1);
    }
  }
  Y->l = n + 1;
  fclose(datafp);
  fprintf(stderr, "Data has %d samples. Running %d x %d iterations of the sampler...\n", Y->l, NOut, SubSamp);
}

/************************************************************************/
/* Write iterations on disk						*/
/************************************************************************/
void write_data(mix* theta, double wght, double logl, short* stats, int in) {

  if (in == 0) {
    /* Initial value, we need to open the files				*/
    /* Open statfile							*/
    strcpy(StatFile, ResFile);
    strcat(StatFile, StatExt);
    StatFP = fopen(StatFile, "w");
    if (StatFP == NULL) {
      fprintf(stderr, "Can't open output file %s\n", StatFile);
      exit(1);
    }
    /* Open resfile							*/
    strcat(ResFile, ResExt);	/* ResFile is corrupted there!		*/ 
    ResFP = fopen(ResFile, "w");
    if (ResFP == NULL) {
      fprintf(stderr, "Can't open output file %s\n", ResFile);
      exit(1);
    }
  }
  /* Write statistics and data						*/
  fprintf(StatFP, "%d %d %.2f %d %d %f\n", in, theta->k, logl,
    stats[0], stats[1], wght);
  fwrite(&wght, sizeof(double), (size_t) 1, ResFP);
  fwrite(theta->w, sizeof(double), (size_t) theta->k, ResFP);
  fwrite(theta->mu, sizeof(double), (size_t) theta->k, ResFP);
  fwrite(theta->var, sizeof(double), (size_t) theta->k, ResFP);
  if (in == NIt) {
    /* That's the end, close result file and write paremeters on stdout	*/
    fclose(StatFP);
    fclose(ResFP);  
    printf("SK0=%lu\n", SK0);
    printf("SK1=%lu\n", SK1);
    printf("SK2=%lu\n", SK2);
    printf("NOut=%d\n", NOut);
    printf("SubSamp=%d\n", SubSamp);
    printf("Data=%s\n", DataFile);
    printf("M=%d\n", M);
    printf("Kappa=%f\n", Kappa);
    printf("AlphaVar=%f\n", AlphaVar);
    printf("BetaVar=%f\n", BetaVar);
    printf("Eta=%f\n", Eta);
    printf("Rho=%f\n", Rho);
    printf("Nu=%f\n", Nu);
    printf("PFixed=%f\n", PFixed);
    printf("PBirth=%f\n", PBirth);
    printf("PDeath=%f\n", PDeath);
    printf("PSplit=%f\n", PSplit);
    printf("Gamma_S=%f\n", Gamma_S);
    printf("Rho_S=%f\n", Rho_S);
    printf("Nu_S=%f\n", Nu_S);
  }
}

/************************************************************************/
/* Compute the rate of all possible moves				*/
/************************************************************************/
double compute_rates(vec* Y, mix* m, double logl, double* move_prob,
  mix m_dead[MaxK], double* logl_dead, mix m_merge[MaxKPairs], double* logl_merge) {
  double weight = 0.0;
  double gpdf[MaxN][MaxK];
  int i, nr;

  /* Fixed move rate							*/
  move_prob[0] = PFixed;
  /* Birth and split rate (no possible birth when we have M classes)    */
  if (m->k == M) {
    move_prob[1] = 0.0;
    move_prob[2+m->k] = 0.0;
  }
  else {
    move_prob[1] = PBirth;
    move_prob[2+m->k] = PSplit;
  }
  /* Compute likelihood if needed					*/
  if ((m->k > 1) && ((PBirth > 0.0) || (PSplit > 0.0)))
    likelihood_individual_components(Y, m, gpdf);
  /* Death and merge rates						*/
  /* (they are zero when we only have 1 class)				*/
  if (m->k > 1) {
    compute_death_rates(Y, m, logl, &(move_prob[2]), m_dead, logl_dead, gpdf);
    compute_merge_rates(Y, m, logl, &(move_prob[3+m->k]), m_merge, logl_merge, gpdf);
  }
  else {
    move_prob[2] = 0.0;
    move_prob[3+m->k] = 0.0;
  }
  /* Normalize rates to a probability vector				*/
  nr = 3 + m->k + (m->k * (m->k - 1))/2;
  for (i=0; i<nr; i++)
    weight += move_prob[i];  
  for (i=0; i<nr; i++)
    move_prob[i] /= weight;
  /* Return the expected holding time rather than the rate		*/
  weight = 1.0/weight;
  return weight;
}

/************************************************************************/
/* Compute rates of death moves						*/
/************************************************************************/
void compute_death_rates(vec* Y, mix* m, double logl, double* rates, mix m_dead[MaxK], double* logl_dead, double gpdf[MaxN][MaxK]) {
  int i, k, t;
  double tmp, sum;

  if (PBirth > 0.0) {
    /* Compute all possible dead models					*/
    for (k=0; k<(m->k); k++) {
      m_dead[k].k = (m->k)-1;
      sum = 0.0;
      /* 1 New parameters						*/
      for (i=0; i<k; i++) {
	m_dead[k].w[i] = m->w[i];
	sum += m->w[i];
	m_dead[k].mu[i] = m->mu[i];
	m_dead[k].var[i] = m->var[i];
      }
      for (i=1+k; i<(m->k); i++) {
	m_dead[k].w[i-1] = m->w[i];
	sum += m->w[i];
	m_dead[k].mu[i-1] = m->mu[i];
	m_dead[k].var[i-1] = m->var[i];
      }
      /* We could divide by (1-m->w[k]) but this is perhaps safer if	*/
      /*  one of the weights is very close to 1				*/
      for (i=0; i<m_dead[k].k; i++) {
	m_dead[k].w[i] /= sum;
      }
      /* 2 Compute log-likelihood for the dead model			*/
      logl_dead[k] = 0.0;
      for (t=0; t<(Y->l); t++){
	tmp = 0.0;
	for (i=0; i<k; i++) {
	  tmp += gpdf[t][i];
	}
	for (i=1+k; i<(m->k); i++) {
	  tmp += gpdf[t][i];
	}
	logl_dead[k] += log(tmp);
      }
      /* Perhaps faster than dividing "tmp" by "sum" above		*/
      logl_dead[k] += -(Y->l)*log(sum);
      rates[k] = exp(logl_dead[k]-logl) * PBirth/(double)m->k;
    }
  }
  else {
    /* If the birth rate is zero then the death rate also is zero	*/
    for (k=0; k<(m->k); k++)
      rates[k] = 0.0;
  }
}

/************************************************************************/
/* Compute rates of merge moves						*/
/************************************************************************/
void compute_merge_rates(vec* Y, mix* m, double logl, double* rates, mix m_merge[MaxKPairs], double* logl_merge, double gpdf[MaxN][MaxK]) {
  int i, k, l1, l2, lm=0, im, t;
  double tmp, std_merge, pdf_merge[MaxN];

  if (PSplit > 0.0) {
    k = (m->k)-1;
    /* Compute all possible merged models				 */
    for (l1=0; l1<(m->k)-1; l1++) {
      for (l2=l1+1; l2<(m->k); l2++) {
        m_merge[lm].k = k;
        /* 1 New parameters						*/
        im = 0;
        for (i=0; i<(m->k); i++) {
	  if ((i != l1) && (i != l2)) {
	    m_merge[lm].w[im] = m->w[i];
	    m_merge[lm].mu[im] = m->mu[i];
	    m_merge[lm].var[im] = m->var[i];
            im++;
	  }
        }
        /* Put merged class in last position				*/
        m_merge[lm].w[k-1] = m->w[l1] + m->w[l2];
        m_merge[lm].mu[k-1] = 0.5*(m->mu[l1] + m->mu[l2]);
        m_merge[lm].var[k-1] = sqrt(m->var[l1] * m->var[l2]);
        /* 2 Compute log-likelihood for the merged model		*/
        /* Compute weighted pdf for the merged model			*/
        std_merge = sqrt(m_merge[lm].var[k-1]);
        for (t=0; t<(Y->l); t++)
          pdf_merge[t] = m_merge[lm].w[k-1] *
            exp(-square(Y->dat[t]-(m_merge[lm].mu[k-1]))/(2.0*m_merge[lm].var[k-1]))
            /std_merge;
        logl_merge[lm] = 0.0;
        for (t=0; t<(Y->l); t++){
	   tmp = pdf_merge[t];
           for (i=0; i<(m->k); i++)
	     if ((i != l1) && (i != l2))
               tmp += gpdf[t][i];
	   logl_merge[lm] += log(tmp);
        }
        rates[lm] = exp(logl_merge[lm]-logl-compute_log_split_ratio(&(m_merge[lm]), m, k-1, l1, l2)) * 2.0*PSplit/(double)(k*(k+1));
        lm++;
      }
    }
  }
  else {
    /* If the split rate is zero then the merge rate also is zero	*/
    for (k=0; k<(m->k*((m->k)-1))/2; k++)
      rates[k] = 0.0;
  }
}
