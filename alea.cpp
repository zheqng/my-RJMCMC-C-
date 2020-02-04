/************************************************************************/
/* CT/RJ-Mix Transdimensional MCMC for Gaussian mixtures		*/
/*									*/
/* As discussed in section 4 of						*/
/* O. Capp? C. Robert, and T. Ryd�n. Reversible jump, birth-and-death	*/
/* and more general continuous time Markov chain Monte Carlo samplers.	*/
/* J. Royal Statist. Soc. Ser. B, 65(3):679-700, 2003.			*/
/* Further information at http://www.tsi.enst.fr/~cappe/ctrj_mix/       */
/*									*/
/* Version: $Id: alea.c,v 1.10 2005/10/06 07:53:27 cappe Exp $		*/
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
/* alea.c: Library of simulation functions by Christian P. Robert	*/
/************************************************************************/

#include "alea.h"

/*__________________Kiss, The generator____________________________*/
/* (This is basically a version of the most basic KISS generator,  */
/*  with a 2^32 period, as proposed by George Marsaglia. See       */
/*  http://mathforum.org/discuss/sci.math/a/m/473973/474544 for a  */
/*  recent update on this issue. The floating point value returned */
/*  by this version is in [0,(2^32-1)/2^32], that is x == 0.0 may  */
/*  happen but not x == 1.0.)                                      */
/*_________________________________________________________________*/


double kiss(Generator & g) {

								double x;

								(g.j) = (g.j) ^ ((g.j)<<17);
								(g.k) = ( (g.k) ^ ( (g.k) <<18)) & 0x7FFFFFFF;
								(g.i) = 69069UL*((g.i))+23606797UL;
								(g.j) ^= ((g.j)>>15);
								(g.k) ^= ((g.k)>>13);
								x = (double)((g.i)+(g.j)+(g.k)) / ULONG_NORM;
								// cout<<"kiss(g):"<<x<<" ";
								return(x);
}

/*__________________Normal (Box-Muller)____________________________*/
void nor(double* t1, double*t2, Generator & g0) {
								double z,u2;
								Generator g;

								g = g0;
								u2=kiss(g);
								if (!((u2>0.)&&(u2<1.)))
																u2=kiss(g);
								z = sqrt( (-2.0*( log( u2 ) ) ) );
								u2 = 2.0*PI*kiss(g);
								g0 = g;
								*t1 = (double) ( z*cos(u2) );
								*t2 = (double) ( z*sin(u2) );
}

/*__________________Gamma(a,1), a>1____________________________*/
double Gam_1(double alf, Generator & g0) {
								double b,d,u1,v,w;
								Generator g;
								boolean stoop;

								g = g0;
								b = alf -1;
								d = 3*alf - 0.75;
								stoop = FALSE;

								while (!stoop)
								{  u1 = kiss(g);
											w = u1*(1-u1);
											v = b + sqrt(d/w)*(u1-0.5);
											stoop = (v>0) && (kiss(g)<(exp(b-v+b*log(v/b))/sqrt(64*w*w*w)));}
								g0 = g;
								return(v);
}

/*__________________Gamma(a,1) a<1____________________________*/
double Gam_2(double alf, Generator & g0) {
								double x;
								Generator g;

								g = g0;
								x =exp(log(kiss(g))/alf)*Gam_1(alf+1,g);
								g0 = g;
								return (x);
}

/*__________________Gamma(a,1), a>>1____________________________*/
double Gam_3(double alf, Generator & g0) {
								double x,y=-1;

								while (y<0) {
																nor(&x,&y,g0);
																y=x*sqrt(alf)+alf;
								}
								return(y);
}

/*__________________Gamma(a,1), summary____________________________*/
double Gam(double alf, Generator & g0) {
								double x;
								Generator g;

								g = g0;
								if (alf>1.)
								{  if (alf<500.)
																			x = Gam_1(alf,g);
											else
																			x = Gam_3(alf,g); }
								else
								{ if (alf==1.)
																		x = -log(kiss(g));
										else
																		x = Gam_2(alf,g); }

								g0 = g;
								return (x);
}

/*__________________Beta_________________________________________________*/
double Beta(double alf, double bet,Generator & g0) {
								double x,y;
								Generator g;

								g = g0;
								x = Gam (alf,g);
								y = Gam (bet,g);
								g0 = g;
								return ( x/(x+y) );
}
