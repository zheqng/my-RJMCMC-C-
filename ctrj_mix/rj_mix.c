/************************************************************************/
/* CT/RJ-Mix Transdimensional MCMC for Gaussian mixtures		*/
/*									*/
/* As discussed in section 4 of						*/
/* O. Capp�, C. Robert, and T. Ryd�n. Reversible jump, birth-and-death	*/
/* and more general continuous time Markov chain Monte Carlo samplers.	*/
/* J. Royal Statist. Soc. Ser. B, 65(3):679-700, 2003.			*/
/* Further information at http://www.tsi.enst.fr/~cappe/ctrj_mix/      */
/*									*/
/* Version: $Id: rj_mix.c,v 3.8 2003/12/02 12:54:05 cappe Exp $		*/
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
/* rj_mix: Implementation of reversible jump MCMC for mixtures		*/
/*									*/
/* Ouput formats:							*/
/*									*/
/* Stat. file (".sts") is a text file with, on each line:		*/
/* | n iter | k | log-l | move type |					*/
/*				    |					*/
/*       for type 1 (fixed k):      | accpt w | accpt mu | accpt var  | */
/*       for type 2 (birth/death):  | birth   | accpt    | unused (-1)| */
/*       for type 3 (split/merge):  | split   | accpt    | unused (-1)| */
/*									*/
/* Res. file (".res") is a binary file which contains doubles:		*/
/* w(0,1), ..., w(0,k(0)), mu(0,1), ..., mu(0,k(0)), var(0,1), ...,	*/
/* var(0,k(0)), w(SubSamp,1),..., w(SubSamp,k(SubSamp)), mu(SubSamp,1),	*/
/* ...									*/
/*									*/
/* Written on stdout is a copy of the input argument (".arg" file) with */
/* updated values of the seeds of the random generator (SK1, SK2, SK3)  */
/*									*/
/* NB: move type, accpt w, accpt mu, accpt var, birth... refer to	*/
/* iteration "n iter" (it is arbitrary for n iter = 0 which is the	*/
/* inial value). All other values (k, log-l) refer to the value of the	*/
/* parameters after "n iter" iterations.				*/
/************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mix_lib.h"
#include <time.h>
/* Definitions								*/
#define NArgs 20
#define NStat 4
#define StrLen 160
#define ArgExt ".arg"
#define ResExt ".res"
#define StatExt ".sts"

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
double Eta, SqrtEta, Rho, SqrtRho, Nu, SqrtNu, PFixed;
/* Birth and death							*/
double PBirth, PDeath, PFixed_or_BD;
/* Split and merge							*/
double PSplit, PMerge, Gamma_S, Rho_S, SqrtRho_S, Nu_S, SqrtNu_S, SplitLogConstFact;

/* Function prototypes							*/
int main(int argc, char** argv);
void read_parameters(int argc, char** argv, vec* Y);
void write_data(mix* theta, double logl, short* stats, int in);
void RJMH_birth_or_death(vec* Y, mix* m, double* logl, short* stats);
void RJMH_split_or_merge(vec* Y, mix* m, double* logl, short* stats);


/************************************************************************/
/* Main									*/
/************************************************************************/
int main(int argc, char** argv) {
        vec Y;
        mix theta;
        double logl, ran;
        short stats[NStat];
        int in;
        time_t t_start, t_end;
        // ofstream TimeFP("time.txt");
        double DiffTime;
        t_start = time(NULL);




        read_parameters(argc, argv, &Y);
        /* Initial value							*/
        draw_initial_model(&Y, &theta, &logl);
        write_data(&theta, logl, stats, 0);
        /* Main loop								*/
        for (in=1; in<=NIt; in++) {
                ran = kiss(&SK0,&SK1,&SK2);
                if (ran < PFixed) {
                        /* Fixed k move							*/
                        stats[0] = 1;
                        stats[1] = MNRWMH_weights(&Y, &theta, &logl);
                        stats[2] = RWMH_means(&Y, &theta, &logl);
                        stats[3] = MRWMH_variances(&Y, &theta, &logl);
                }
                else if (ran < PFixed_or_BD) {
                        /* Birth or death							*/
                        stats[0] = 2;
                        RJMH_birth_or_death(&Y, &theta, &logl, &(stats[1]));
                        stats[3] = -1;
                }
                else {
                        /* Split or merge							*/
                        stats[0] = 3;
                        RJMH_split_or_merge(&Y, &theta, &logl, &(stats[1]));
                        stats[3] = -1;
                }
                /* Note: The test below is true for in = 0 and in = NIt		*/
                if ((div(in, SubSamp)).rem == 0)
                        write_data(&theta, logl, stats, in);
                // printf("%d\n",in);

        }
        t_end = time(NULL);
        DiffTime = difftime(t_end, t_start);
        printf("%f\n",DiffTime);
        return(0);
}


/************************************************************************/
/* Read parameters (including seeds for the random generator) and data  */
/************************************************************************/
void read_parameters(int argc, char** argv, vec* Y) {
        char argfile[StrLen];
        int argflag[NArgs], i, n;
        double y;
        FILE *argfp, *datafp;

        if (argc < 2) {
                fprintf(stderr, "rj_mix base_name (without .arg extension)\n");
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
        PMerge = 1.0 - (PFixed + PBirth + PDeath + PSplit);
        /* Beware of incorrect probabilities					*/
        if ((PFixed < 0.0) || (PBirth < 0.0) || (PDeath < 0.0) || (PSplit < 0.0)
            || (PMerge > 1.0)) {
                fprintf(stderr, "Something is wrong with the probabilities\n");
                exit(1);
        }
        PFixed_or_BD = PFixed + PBirth + PDeath;
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
void write_data(mix* theta, double logl, short* stats, int in) {

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
                strcat(ResFile, ResExt); /* ResFile is corrupted there!		*/
                ResFP = fopen(ResFile, "w");
                if (ResFP == NULL) {
                        fprintf(stderr, "Can't open output file %s\n", ResFile);
                        exit(1);
                }
        }
        /* Write statistics and data						*/
        fprintf(StatFP, "%d %d %.2f %d %d %d %d\n", in, theta->k, logl,
                stats[0], stats[1], stats[2], stats[3]);
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
/* Implements the RJ MH move based on birth/death proposal		*/
/************************************************************************/
void RJMH_birth_or_death(vec* Y, mix* m, double* logl, short* stats) {
        double prop_ratio, logl_new;
        mix m_new;
        int k;
        short birth, accept;

        /* Take care of case 1, 2, M-1 and M components (this could be done	*/
        /* more simply)							*/
        if (m->k == 1) {
                birth = 1;
                prop_ratio = PDeath/(PBirth+PDeath);
        }
        else if (m->k == M) {
                birth = 0;
                prop_ratio = PBirth/(PBirth+PDeath);
        }
        else {
                if (kiss(&SK0,&SK1,&SK2) < PBirth/(PBirth+PDeath)) {
                        /* Birth								*/
                        birth = 1;
                        if (m->k == (M-1))
                                prop_ratio = (PBirth+PDeath)/PBirth;
                        else
                                prop_ratio = PDeath/PBirth;
                }
                else {
                        /* Death								*/
                        birth = 0;
                        if (m->k == 2)
                                prop_ratio = (PBirth+PDeath)/PDeath;
                        else
                                prop_ratio = PBirth/PDeath;
                }
        }
        if (birth) {
                mixcpy(m, &m_new);
                logl_new = prop_birth(Y, &m_new);
        }
        else {
                k = (int)floor((double)m->k * kiss(&SK0,&SK1,&SK2));
                logl_new = prop_death(Y, m, &m_new, k);
        }
        /* Accept/reject							*/
        if (kiss(&SK0,&SK1,&SK2) < (exp(logl_new-(*logl))*prop_ratio)) {
                accept = 1;
                /* Modify the parameters and log likelihood				*/
                mixcpy(&m_new, m);
                *logl = logl_new;
        }
        else
                accept = 0;
        stats[0] = birth;
        stats[1] = accept;
}


/************************************************************************/
/* Implements the RJ MH move based on split/merge proposal		*/
/************************************************************************/
void RJMH_split_or_merge(vec* Y, mix* m, double* logl, short* stats) {
        double prop_ratio, logl_new, add_logratio;
        mix m_new;
        int k, k1, k2;
        short split, accept;

        /* Take care of case 1, 2, M-1 and M components (this could be done	*/
        /* more simply)							*/
        if (m->k == 1) {
                split = 1;
                prop_ratio = PMerge/(1.0-PFixed_or_BD);
        }
        else if (m->k == M) {
                split = 0;
                prop_ratio = PSplit/(1.0-PFixed_or_BD);
        }
        else {
                if (kiss(&SK0,&SK1,&SK2) < PSplit/(1.0-PFixed_or_BD)) {
                        /* Split								*/
                        split = 1;
                        if (m->k == (M-1))
                                prop_ratio = (1.0-PFixed_or_BD)/PSplit;
                        else
                                prop_ratio = PMerge/PSplit;
                }
                else {
                        /* Merge								*/
                        split = 0;
                        if (m->k == 2)
                                prop_ratio = (1.0-PFixed_or_BD)/PMerge;
                        else
                                prop_ratio = PSplit/PMerge;
                }
        }
        if (split) {
                mixcpy(m, &m_new);
                k = (int)floor((double)m->k * kiss(&SK0,&SK1,&SK2));
                logl_new = prop_split(Y, &m_new, k);
                add_logratio = compute_log_split_ratio(m, &m_new, k, k, m->k);
        }
        else {
                /* Draw two distinct indices using modulo m->k addition		*/
                k1 = (int)floor((double)m->k * kiss(&SK0,&SK1,&SK2));
                k2 = 1+(int)floor((double)(m->k - 1) * kiss(&SK0,&SK1,&SK2));
                k2 = div(k1+k2, m->k).rem;
                logl_new = prop_merge(Y, m, &m_new, k1, k2);
                add_logratio = -compute_log_split_ratio(&m_new, m, m_new.k-1, k1, k2);
        }
        /* Accept/reject							*/
        if (kiss(&SK0,&SK1,&SK2) < (exp(logl_new-(*logl)+add_logratio)*prop_ratio)) {
                accept = 1;
                /* Modify the parameters and log likelihood				*/
                mixcpy(&m_new, m);
                *logl = logl_new;
        }
        else
                accept = 0;
        stats[0] = split;
        stats[1] = accept;
}
