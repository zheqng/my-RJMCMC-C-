/************************************************************************/
/* CT/RJ-Mix Transdimensional MCMC for Gaussian mixtures		*/
/*									*/
/* As discussed in section 4 of						*/
/* O. Capp�, C. Robert, and T. Ryd�n. Reversible jump, birth-and-death	*/
/* and more general continuous time Markov chain Monte Carlo samplers.	*/
/* J. Royal Statist. Soc. Ser. B, 65(3):679-700, 2003.			*/
/* Further information at http://www.tsi.enst.fr/~cappe/ctrj_mix/      */
/*									*/
/* Version: $Id: mix_lib.c,v 3.7 2003/12/02 12:54:05 cappe Exp $	*/
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
/* mix_lib: Library of functions used by both rj_mix and ct_mix.	*/
/************************************************************************/

#include "mix_lib.h"

/************************************************************************/
/* Computes the loglikelihood up to constant factors			*/
/************************************************************************/
double log_likelihood(vec* Y, mix* m) {
        double std[MaxK], tmp[MaxK], norm;
        double logl=0.0;
        int t, i;
        // #pragma omp parallel for
        for (i=0; i<(m->k); i++) {
                std[i] = sqrt(m->var[i]);
        }
        // #pragma omp parallel for
        for (t=0; t<(Y->l); t++) {
                norm = 0;
                for (i=0; i<(m->k); i++) {
                        tmp[i] = m->w[i]*exp(-square(Y->dat[t]-(m->mu[i]))/(2.0*m->var[i]))/std[i];
                        norm += tmp[i];
                }
                logl += log(norm);
        }
        return logl;
}

/************************************************************************/
/* Computes the loglikelihood up to constant factors:			*/
/* This is the same function as above but with more detailed args.,	*/
/* (usefull for the fixed k moves)					*/
/************************************************************************/
double log_likelihood_micro(vec* Y, double* w, double *mu, double *var, int k) {
        double std[MaxK], tmp[MaxK], norm;
        double logl=0.0;
        int t, i;

        for (i=0; i<k; i++) {
                std[i] = sqrt(var[i]);
        }
        for (t=0; t<(Y->l); t++) {
                norm = 0;
                for (i=0; i<k; i++) {
                        tmp[i] = w[i]*exp(-square(Y->dat[t]-(mu[i]))/(2.0*var[i]))/std[i];
                        norm += tmp[i];
                }
                logl += log(norm);
        }
        return logl;
}

/************************************************************************/
/* Computes the likelihood for each component of the mixture		*/
/************************************************************************/
void likelihood_individual_components(vec* Y, mix* m, double gpdf[MaxN][MaxK]) {
        double std[MaxK];
        int t, i;

        /* 0. Compute standard deviations once				*/
        for (i=0; i<(m->k); i++) {
                std[i] = sqrt(m->var[i]);
        }
        /* 1. Compute all Gaussian pdf					*/
        for (t=0; t<(Y->l); t++) {
                for (i=0; i<(m->k); i++) {
                        gpdf[t][i] = m->w[i]*exp(-square(Y->dat[t]-(m->mu[i]))/(2.0*m->var[i]))/std[i];
                }
        }
}

/************************************************************************/
/* MNRW MH on the weights with log-normal proposal			*/
/************************************************************************/
short MNRWMH_weights(vec* Y, mix* m, double* logl) {
        double ran1, ran2, prop=0.0, sum=0.0, ratio, logl_new;
        double w_new[MaxK];
        int i;
        short accept;

        /* Propose new weights with log-normal perturbation			*/
        for (i=0; i<(m->k); i++) {
                /* New veights and accumulate for normalizing			*/
                nor(&ran1, &ran2, &SK0, &SK1, &SK2);
                w_new[i] = m->w[i]*exp(ran1*SqrtEta);
                sum += w_new[i];
        }
        for (i=0; i<(m->k); i++) {
                /* Normalize the weigths						*/
                w_new[i] /= sum;
                /* Accumulate for acceptance ratio					*/
                prop += log(w_new[i]/m->w[i]);
        }
        /* Compute acceptance ratio */
        logl_new = log_likelihood_micro(Y, w_new, m->mu, m->var, m->k);
        ratio = exp(logl_new-(*logl)+prop);
        if (kiss(&SK0,&SK1,&SK2) < ratio) {
                accept = 1;
                /* Update model and log-likelihood					*/
                for (i=0; i<(m->k); i++)
                        m->w[i] = w_new[i];
                *logl = logl_new;

        }
        else
                accept = 0;
        return accept;
}

/************************************************************************/
/* RW MH on the means with normal proposal                              */
/************************************************************************/
short RWMH_means(vec* Y, mix* m, double* logl) {
        double ran1, ran2, prop=0.0, ratio, std, logl_new;
        double mu_new[MaxK];
        int i;
        short accept;

        /* Propose new means with normal perturbation				*/
        /* The variance of the proposal is proportional to 1/k		*/
        std = SqrtRho/sqrt((double) m->k);
        for (i=0; i<(m->k); i++) {
                /* New mean								*/
                nor(&ran1, &ran2, &SK0, &SK1, &SK2);
                mu_new[i] = m->mu[i]+ran1*std;
                /* Accumulate for acceptance ratio					*/
                prop += square(mu_new[i])-square(m->mu[i]);
        }
        /* Compute acceptance ratio */
        logl_new = log_likelihood_micro(Y, m->w, mu_new, m->var, m->k);
        ratio = exp(logl_new-(*logl)-prop/(2*Kappa));
        if (kiss(&SK0,&SK1,&SK2) < ratio) {
                accept = 1;
                /* Update model and log-likelihood					*/
                for (i=0; i<(m->k); i++)
                        m->mu[i] = mu_new[i];
                *logl = logl_new;
        }
        else
                accept = 0;
        return accept;
}

/************************************************************************/
/* MRW MH on the variances with log-normal proposal                     */
/************************************************************************/
short MRWMH_variances(vec* Y, mix* m, double* logl) {
        double ran1, ran2, prop=0.0, ratio, logl_new;
        double var_new[MaxK];
        int i;
        short accept;

        /* Propose new variances with log-normal perturbation			*/
        for (i=0; i<(m->k); i++) {
                /* New variances							*/
                nor(&ran1, &ran2, &SK0, &SK1, &SK2);
                var_new[i] = m->var[i]*exp(ran1*SqrtNu);
                /* Accumulate for acceptance ratio					*/
                prop += -AlphaVar*log(var_new[i]/m->var[i])
                        -BetaVar*(1/var_new[i]-1/m->var[i]);
        }
        /* Compute acceptance ratio */
        logl_new = log_likelihood_micro(Y, m->w, m->mu, var_new, m->k);
        ratio = exp(logl_new-(*logl)+prop);
        if (kiss(&SK0,&SK1,&SK2) < ratio) {
                accept = 1;
                /* Update model and log-likelihood					*/
                for (i=0; i<(m->k); i++)
                        m->var[i] = var_new[i];
                *logl = logl_new;
        }
        else
                accept = 0;
        return accept;
}

/************************************************************************/
/* Duplicates a mixture model						*/
/************************************************************************/
void mixcpy(mix* m_src, mix* m_dest) {
        int i;

        m_dest->k = m_src->k;
        for (i=0; i<(m_src->k); i++) {
                m_dest->w[i] = m_src->w[i];
                m_dest->mu[i] = m_src->mu[i];
                m_dest->var[i] = m_src->var[i];
        }
}

/************************************************************************/
/* Prints the details of a mixture model on STDOUT (for debugging	*/
/* purposes								*/
/************************************************************************/
void mixprint(mix* m){
        int i;

        printf("k=%d\n", m->k);
        /* Weights								*/
        printf("w=[");
        for (i=0; i<(m->k)-1; i++)
                printf("%f ", m->w[i]);
        printf("%f]\n", m->w[(m->k)-1]);
        /* Means								*/
        printf("mu=[");
        for (i=0; i<(m->k)-1; i++)
                printf("%f ", m->mu[i]);
        printf("%f]\n", m->mu[(m->k)-1]);
        /* Variances								*/
        printf("var=[");
        for (i=0; i<(m->k)-1; i++)
                printf("%f ", m->var[i]);
        printf("%f]\n", m->var[(m->k)-1]);
}

/************************************************************************/
/* Proposes a birth move, put new class in last position and returns	*/
/* log-likelihood							*/
/************************************************************************/
double prop_birth(vec* Y, mix* m) {
        int i;
        double ran1, ran2, logl;

        if (m->k >= M) {
                fprintf(stderr, "prop_birth called with %d classes while maximum allowed is %d, aborting\n", m->k, M);
                exit(1);
        }
        /* Propose new class, which is stored in position m->k		*/
        /* (since the first one is stored at position 0 !!!)			*/
        /* Beta for the weight						*/
        m->w[m->k] = Beta(1.0, (double)m->k, &SK0, &SK1, &SK2);
        /* Normalize the weights						*/
        for (i=0; i<(m->k); i++)
                m->w[i] = m->w[i]*(1-(m->w[m->k]));
        /* Normal for the mean						*/
        nor(&ran1, &ran2, &SK0, &SK1, &SK2);
        m->mu[m->k] = ran1*SqrtKappa;
        /* Inverse gamma for the variance					*/
        m->var[m->k] = BetaVar/Gam(AlphaVar, &SK0, &SK1, &SK2);
        /* Add 1 to the number of classes					*/
        m->k++;
        logl = log_likelihood(Y, m);
        return logl;
}

/************************************************************************/
/* Proposes a death move and returns log-likelihood for new parameters  */
/************************************************************************/
double prop_death(vec* Y, mix* m_old, mix* m_new, int k) {
        int i;
        double sum=0.0, logl;

        if (m_old->k <= 1) {
                fprintf(stderr, "prop_death called with %d classe(s), aborting\n", m_old->k);
                exit(1);
        }
        if ((k < 0) || (k >= m_old->k)) {
                fprintf(stderr, "prop_death cannot delete class %d when k=%d, aborting\n", k, m_old->k);
                exit(1);
        }
        m_new->k = (m_old->k)-1;
        for (i=0; i<k; i++) {
                m_new->w[i] = m_old->w[i];
                sum += m_old->w[i];
                m_new->mu[i] = m_old->mu[i];
                m_new->var[i] = m_old->var[i];
        }
        for (i=1+k; i<(m_old->k); i++) {
                m_new->w[i-1] = m_old->w[i];
                sum += m_old->w[i];
                m_new->mu[i-1] = m_old->mu[i];
                m_new->var[i-1] = m_old->var[i];
        }
        /* We could divide by (1-m_old->w[k]) but this is perhaps safer if	*/
        /*  one of the weights is very close to 1				*/
        for (i=0; i<m_new->k; i++) {
                m_new->w[i] /= sum;
        }
        logl = log_likelihood(Y, m_new);
        return logl;
}

/************************************************************************/
/* Applies the inverse of the discrete df corresponding to probability  */
/* vector p to a number 0 <= ran < 1					*/
/* (p sums to 1 and ran >= 0 are not checked!)				*/
/************************************************************************/
int invert_disc_df(double* p, int l, double ran) {
        int i;
        double F[MaxKPairs];

        /* Avoid infinite loops						*/
        if (ran >= 1.0) {
                fprintf(stderr, "Uniform random deviate %f is not stricly less than 1\n", ran);
                exit(1);
        }
        /* Construct distribution function					*/
        F[0] = p[0];
        /* We can safely stop if the df is already equal to 1			*/
        for (i=1; i<(l-1) && (F[i-1] < 1.0); i++)
                F[i] = F[i-1] + p[i];
        /* Last value set to 1 in case the prob. do not sum exactly to 1	*/
        F[l-1] = 1.0;
        i = 0;
        while (ran >= F[i])
                i++;
        return i;
}

/************************************************************************/
/* Draws a model from the prior						*/
/************************************************************************/
mix draw_from_prior(void) {
        mix m;
        double ran1, ran2, sum=0.0;
        int i;

        m.k = 1 + (int)floor((double)M * kiss(&SK0,&SK1,&SK2));
        for (i=0; i<(m.k); i++) {
                /* Expon. for the weights (making sure we dont take the log of 0)   */
                m.w[i] = -log(1.0-kiss(&SK0, &SK1, &SK2));
                sum += m.w[i];
                /* Normal for the mean						*/
                nor(&ran1, &ran2, &SK0, &SK1, &SK2);
                m.mu[i] = ran1*SqrtKappa;
                /* Inverse gamma for the variance					*/
                m.var[i] = BetaVar/Gam(AlphaVar, &SK0, &SK1, &SK2);
        }
        /* Make the weights Dirichlet(1, ..., 1)				*/
        for (i=0; i<(m.k); i++) {
                m.w[i] /= sum;
        }
        return m;
}

/************************************************************************/
/* Initializes the parameters of the model (and returns the scaling	*/
/* needed  to obtain a non -Inf log-likelihood, this should normaly be	*/
/* 1)									*/
/************************************************************************/
int draw_initial_model(vec* Y, mix* m, double* logl) {
        int j, f=1;

        /* Start by drawing from prior					*/
        *m = draw_from_prior();
        *logl = log_likelihood(Y, m);
        /* This can only ever happen when all variances are too small,	*/
        /* scale them all by a factor 2 untill this is OK			*/
        while (*logl < NEG_REAL_MIN) {
                for (j=0; j<(m->k); j++)
                        m->var[j] *= 2.0;
                *logl = log_likelihood(Y, m);
                f <<= 1;
        }
        /* Print out a warning message if this happens			*/
        if (f > 1)
                fprintf(stderr, "Warning: Scaled initial variances by %d\n", f);
        return f;
}

/************************************************************************/
/* Computes the log of the gamma function				*/
/************************************************************************/
/* Adapted from GAUSS code by Paul L. Fackler available at		*/
/* http://www.american.edu/academic.depts/cas/econ/gaussres/pdf/	*/
/*     loggamma.src							*/
/* Original comments follow:						*/
/* Source: Pike, M.C., and I.D. Hill, Algorithm 291, Communications	*/
/* of the ACM, 9,9:p.684 (Sept, 1966). Accepts X>0, accuracy to 10	*/
/* decimal places.  Uses Sterling's formula.				*/
double gammaln(double x) {
        double z;

        x = x+6;
        z = 1.0/(square(x));
        z = (((-0.000595238095238*z+0.000793650793651)
              *z-0.002777777777778)*z+0.083333333333333)/x;
        z = (x-0.5)*log(x)-x+0.918938533204673+z
            -log(x-1)-log(x-2)-log(x-3)-log(x-4)-log(x-5)-log(x-6);
        return z;
}

/************************************************************************/
/* Proposes a split move and returns log-likelihood			*/
/* (put new classes in current and last	position)			*/
/************************************************************************/
double prop_split(vec* Y, mix* m, int k) {
        double logl, eps_w, eps_mu, eps_var;

        if (m->k >= M) {
                fprintf(stderr, "prop_split called with %d classes while maximum allowed is %d, aborting\n", m->k, M);
                exit(1);
        }
        if ((k < 0) || (k >= m->k)) {
                fprintf(stderr, "prop_split called with argument %d while the number of classes is %d, aborting\n", k, m->k);
                exit(1);
        }
        /* Splits class k and put resulting classes in position k and m->k	*/
        eps_w = Beta(Gamma_S, Gamma_S, &SK0, &SK1, &SK2);
        nor(&eps_mu, &eps_var, &SK0, &SK1, &SK2);
        eps_mu = SqrtRho_S*eps_mu;
        eps_var = exp(SqrtNu_S*eps_var);
        m->w[m->k] = eps_w*m->w[k];
        m->w[k] = m->w[k]-(m->w[m->k]); /* Makes sure the sum stays 1		*/
        m->mu[m->k] = m->mu[k]-eps_mu;
        m->mu[k] = m->mu[k]+eps_mu;
        m->var[m->k] = m->var[k]/eps_var;
        m->var[k] = m->var[k]*eps_var;
        /* Add 1 to the number of classes					*/
        m->k++;
        logl = log_likelihood(Y, m);
        return logl;
}

/************************************************************************/
/* Proposes a merge move and returns log-likelihood for new parameters  */
/* (puts the merged class in last position)				*/
/************************************************************************/
double prop_merge(vec* Y, mix* m_old, mix* m_new, int k1, int k2) {
        int i_old, i_new;
        double logl;

        if (m_old->k <= 1) {
                fprintf(stderr, "prop_merge called with %d classe(s), aborting\n", m_old->k);
                exit(1);
        }
        if ((k1 < 0) || (k2 < 0) || (k1 >= m_old->k) || (k2 >= m_old->k) || (k1 == k2)) {
                fprintf(stderr, "prop_merge cannot merge classes %d and %d when k=%d, aborting\n", k1, k2, m_old->k);
                exit(1);
        }
        m_new->k = (m_old->k)-1;
        /* Copy parameters for all classes but k1 and k2			*/
        i_new = 0;
        for (i_old=0; i_old<(m_old->k); i_old++) {
                if (!((i_old == k1) || (i_old == k2))) {
                        m_new->w[i_new] = m_old->w[i_old];
                        m_new->mu[i_new] = m_old->mu[i_old];
                        m_new->var[i_new] = m_old->var[i_old];
                        i_new++;
                }
        }
        /* Merge classes k1 and k2 in the last one of the new model		*/
        m_new->w[(m_new->k)-1] = m_old->w[k1] + m_old->w[k2];
        m_new->mu[(m_new->k)-1] = 0.5 * (m_old->mu[k1] + m_old->mu[k2]);
        m_new->var[(m_new->k)-1] = sqrt(m_old->var[k1]*m_old->var[k2]);
        /* Compute new log-likelihood						*/
        logl = log_likelihood(Y, m_new);
        return logl;
}

/************************************************************************/
/* Computes the part of the split ratio that does not depend upon the	*/
/* data									*/
/************************************************************************/
double compute_log_split_ratio(mix* m, mix* m_split, int k, int k_split1, int k_split2) {
        double lograt;

        /* Fixed part								*/
        lograt = SplitLogConstFact;
        /* Weights								*/
        lograt += log((double)m->k)
                  - (Gamma_S-1)*log(m_split->w[k_split1]*m_split->w[k_split2]/square(m->w[k]))
                  + log(m->w[k]);
        /* Means								*/
        lograt += (square(m->mu[k])-(square(m_split->mu[k_split1])+square(m_split->mu[k_split2])))/(2.0*Kappa)
                  + square(m_split->mu[k_split1]-m_split->mu[k_split2])/(8.0*Rho_S);
        /* Variances								*/
        lograt += -AlphaVar*log(m->var[k])
                  + BetaVar*(1.0/m->var[k]-(1.0/m_split->var[k_split1]+1.0/m_split->var[k_split2]))
                  + square(log(m_split->var[k_split2]/m_split->var[k_split1]))/(8.0*Nu_S);
        return lograt;
}
