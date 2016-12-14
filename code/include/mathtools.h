#ifndef __MATHSTOOLS_H
#define __MATHSTOOLS_H

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <iostream> 

#include <typescalaire.h>

inline int factoriel(const int n)
// return n!
{ int fact=n;
  for(int i=2; i<n;i++)
    fact=fact*i;
  return(fact);
}

inline int binom(const int n, const int k)
// return n!/((n-k)!k!)
{
  return( (factoriel(n)/(factoriel(n-k)*factoriel(k))) );
}
inline double pione(const double cosphi, const double sinphi, 
		     const double A10, const double A01)
{
  return(A10*cosphi+A01*sinphi);
}
inline double pitwo(const double cosphi, const double sinphi, 
		     const double A20, const double A11, const double A02)
{
  return(A20*cosphi*cosphi + A11*cosphi*sinphi + A02*sinphi*sinphi);
}

inline double modulo(const double x,const double y)// return x modulo(y)
{
  return(x - floor(x/y)*y);
}

void Moments(double *M,int n, double *f, double *x, int p, int q );
// Compute the Moments of order "n" of the functions f(i,.)
// in the form M(i)= sum_j f_i(j)x_i(j)^n (but use a trapezoidal rule) 
// f_i(j) is supposed to be stored at f(i*q+j) j=0,...,q-1 i=0,...,p-1
// x_i(j) is supposed to be stored at x(i*q+j) j=0,...,q-1 i=0,...,p-1
// (for the Data Consitency Conditions of order n)

Scalaire DistancePoly(Scalaire *t,Scalaire *f,int n, int degree, int verbose=0);
// measure the norme 2 of the difference of the fonction "f" to the polynomial
// set of degree "degree"
// f is described by it sampled graph (t[i], f[i]), i=0,...,n-1
// the function is given only for degree =0 or 1 or 2 or 3


#endif
