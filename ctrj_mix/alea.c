/************************************************************************/
/* CT/RJ-Mix Transdimensional MCMC for Gaussian mixtures		*/
/*									*/
/* As discussed in section 4 of						*/
/* O. Cappé, C. Robert, and T. Rydén. Reversible jump, birth-and-death	*/
/* and more general continuous time Markov chain Monte Carlo samplers.	*/
/* J. Royal Statist. Soc. Ser. B, 65(3):679-700, 2003.			*/
/* Further information at http://www.tsi.enst.fr/~cappe/ctrj_mix/       */
/*									*/
/* Version: $Id: alea.c,v 1.10 2005/10/06 07:53:27 cappe Exp $		*/
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
double kiss(unsigned long int* i, unsigned long int*j, unsigned long int*k) {
  double x;

  *j = *j ^ (*j<<17);
  *k = (*k ^ (*k<<18)) & 0x7FFFFFFF;
  *i = 69069UL*(*i)+23606797UL;
  *j ^= (*j>>15);
  *k ^= (*k>>13);
  x = (double)(*i+*j+*k) / ULONG_NORM;
  return(x);
}

/*__________________Normal (Box-Muller)____________________________*/
void nor(double* t1, double*t2, unsigned long int* i0, unsigned long int* j0, unsigned long int* k0) {
  double z,u2;
  unsigned long int i,j,k;

  i=*i0; j=*j0; k=*k0;
  u2=kiss(&i,&j,&k);
  if (!((u2>0.)&&(u2<1.)))
	u2=kiss(&i,&j,&k);
  z = sqrt( (-2.0*( log( u2 ) ) ) );
  u2 = 2.0*PI*kiss(&i,&j,&k);
  *i0=i; *j0=j; *k0=k;
  *t1 = (double) ( z*cos(u2) );
  *t2 = (double) ( z*sin(u2) );
}

/*__________________Gamma(a,1), a>1____________________________*/
double Gam_1(double alf, unsigned long int* i0, unsigned long int* j0, unsigned long int* k0) {
  double b,d,u1,v,w;
  unsigned long int i,j,k;
  boolean stoop;
  
  i=*i0; j=*j0; k=*k0;
  b = alf -1;
  d = 3*alf - 0.75;
  stoop = FALSE;

  while (! stoop)
  {  u1 = kiss(&i,&j,&k);
     w = u1*(1-u1);
     v = b + sqrt(d/w)*(u1-0.5);
     stoop = (v>0) && (kiss(&i,&j,&k)<(exp(b-v+b*log(v/b))/sqrt(64*w*w*w)));
  }
  *i0=i; *j0=j; *k0=k;
  return(v);
}

/*__________________Gamma(a,1) a<1____________________________*/
double Gam_2(double alf, unsigned long int* i0, unsigned long int* j0, unsigned long int* k0) {
  double x;
  unsigned long int i,j,k;

  i=*i0; j=*j0; k=*k0;
  x =exp(log(kiss(&i,&j,&k))/alf)*Gam_1(alf+1,&i,&j,&k);
  *i0=i; *j0=j; *k0=k;
  return (x);
}

/*__________________Gamma(a,1), a>>1____________________________*/
double Gam_3(double alf, unsigned long int* i0, unsigned long int* j0, unsigned long int* k0) {
  double x,y=-1;

  while (y<0){
    nor(&x,&y,i0,j0,k0);
    y=x*sqrt(alf)+alf;
  }
  return(y);
}

/*__________________Gamma(a,1), summary____________________________*/
double Gam(double alf, unsigned long int* i0, unsigned long int* j0, unsigned long int* k0) {
  double x;
  unsigned long int i,j,k;

  i=*i0; j=*j0; k=*k0;
  if (alf>1.)
  {  if (alf<500.)
	x = Gam_1(alf,&i,&j,&k);
     else
	x = Gam_3(alf,&i,&j,&k);
  }
  else
  { if (alf==1.)
	x = -log(kiss(&i,&j,&k));
    else
	x = Gam_2(alf,&i,&j,&k);
  }

  *i0 =i; *j0=j; *k0=k;
  return (x);
}

/*__________________Beta_________________________________________________*/
double Beta(double alf, double bet, unsigned long int* i0, unsigned long int* j0, unsigned long int* k0) {
  double x,y;
  unsigned long int i,j,k;

  i=*i0; j=*j0; k=*k0; 
  x = Gam (alf,&i,&j,&k);
  y = Gam (bet,&i,&j,&k);
  *i0=i; *j0=j; *k0=k;
  return ( x/(x+y) );
}
