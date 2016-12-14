// Paquetage Pseudofct
// -------------------
//
// Fonctions utiles a l'implementation de la fonction
// de tomographie pseudo-locale 
//

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <pseudofct.h>
#include <integ.h>

using namespace std;

int fact(int n)
//
// n!
//
{
 int i, res ;

 res = 1 ;
 for (i=1;i<=n;i++)
     res = res * i ;
 return res ;
}

double sigma(double p, double rho)
// 
// sigma_rho(p)
// 
{
 double u,res1,res2,eps1;

 eps1 = 1e-6 ;
 if (fabs(p)>rho + eps1)
    {
     printf("******** pseufofct/sigma : sigma a son support dans [-%lf,%lf]\n",rho,rho) ;
     exit(0) ;
    }
 u = p/rho ;
 res1 = 1-u*u*u*u ;
 res2 = res1*res1*res1*res1 ;
 return res2 ;
}

double sigma_sur_t(double p, double rho)
// 
// sigma_rho(p)/p
// 
{
 double res;

 if (p==0.)
    {
     printf("******** pseufofct/sigma_sur_t : t=0 \n") ;
     exit(0) ;
    }
 else
    res = sigma(p,rho)/p ;
 return res ;
}

double omega_sur_t(double r,double m)
// 
// w1(r)/r
// 
{
 double res;

 if (r==0.)
    {
     printf("******** pseufofct/omega_sur_t : r=0 \n") ;
     exit(0) ;
    }
 else
    res = pow(1-r*r,m+1.)/r ;
 return res ;
}


//  
// CALCUL DU FACTEUR A_{sigma,q}
//
double base_a_sq(double r, double q, double m, int n1, double epsint)
// 
// Fonction de r dont l'integrale entre -1 et 1 (a la somme d'une constante pres
// et au produit de la fonction w1(r)/r pres) definit le facteur A{sigma,q}
// n1 est le nombre de points grace auxquels le resultat de cette fonction
// est calcule.
// Cette fonction est definie pour tout r, meme en 0.
// epsint : Parametre servant a definir l'integrale autour du point singulier t = 0
// 
{
 double res,res1,res2,ti,tinter,tsuiv ;
 double *f ;
 int i ;

 f = (double *) calloc(2000,sizeof(double)) ; 

 for (i=0;i<=n1-1;i++)
    {
     ti = position_abscisse_f(r-q,-epsint,n1,i) ;
     tsuiv = position_abscisse_f(r-q,-epsint,n1,i+1) ;
     tinter = (ti + tsuiv)/2. ;
     f[2*i] = sigma(r-ti,q)/ti ;
     f[2*i+1] = sigma(r-tinter,q)/tinter ;
    }
 f[2*n1] = sigma(r-tsuiv,q)/tsuiv ;
 res1 = simpson_f(f,r-q,-epsint,n1) ;
 for (i=0;i<=n1-1;i++)
    {
     ti = position_abscisse_i(epsint,r+q,n1,i) ;
     tsuiv = position_abscisse_i(epsint,r+q,n1,i+1) ;
     tinter = (ti + tsuiv)/2. ;
     f[2*i] = sigma(r-ti,q)/ti ;
     f[2*i+1] = sigma(r-tinter,q)/tinter ;
    }
 f[2*n1] = sigma(r-tsuiv,q)/tsuiv ;
 res2 = simpson_i(f,epsint,r+q,n1) ;
 res = res1 + res2 ;
 free(f) ;
 return res ;
}

double a_sq(double q, double m, double eps, int n1, int n2, double epsint, double epsext)
//
// Facteur A_{sigma,q}
// n represente le nombre de points utilises pour
// l'evaluation de l'integrale qui definit a_sq
// n2,epsext : Definissent l'integrale / r -- n2 doit etre pair
// n1,epsint : Definissent l'integrale / t 
// 
{
 double res,res1,res2,ti,tsuiv,tinter ;
 double *f,*g ;
 int i ;

 f = (double *) calloc(2000,sizeof(double)) ; 

 for (i=0;i<=n2;i++)
    {
     ti = -1+i*(1.-epsext)/n2 ;
     f[i] = omega_sur_t(ti,m) * base_a_sq(ti,q,m,n1,epsint) ;
    }
 res1 = simpson(f,-1.,-epsext,n2) ;

 for (i=0;i<=n2;i++)
    {
     ti = epsext+i*(1.-epsext)/n2 ;
     f[i] = omega_sur_t(ti,m) * base_a_sq(ti,q,m,n1,epsint) ;
    }
 res2 = simpson(f,epsext,1.,n2) ;

 free(f) ;
 res = res1 + res2 ;
 return res ;
}


//
// CALCUL DU FILTRE DE CONVOLUTION
//
double base_omegaprime(double u, double m, double s, double rho, double eps)
// 
// Fonction de u dont l'integrale entre -rho et
// rho definit (a une constante multiplicative pres)
// le filtre omegaprime_{sigma_rho,eps}(s)
// 
{
 double res,p,rapport ;

 if (fabs(s-u) > eps)
    res = 0. ;
 else
    {
     rapport = (s-u)/eps ;
     p = 1. -(2*m+1.)*rapport*rapport ;
     res = sigma_sur_t(u,rho)*pow(1.-rapport*rapport,m-1.)*p ;
    }
 return res ;
}


double filtre_omegaprime(double m, double s, double rho, double eps, int n)
// 
// Fonction omegaprime_{sigma_rho,eps}(s)
// n+1 = nb de points servant a evaluer l'integrale
// 
{
 double Cprime, res, res1, res2, ti, tsuiv, tinter ;
 double borne_sup, borne_inf ; // les bornes de l'integrale a evaluer
 double *f ;
 double eps3, borne ;       // eps3 : precision de calcul des integrales
 int i,m1 ;

 eps3 = 1e-5 ;
 m1 = (int)(m+0.5) ;
 f = (double *) calloc(2000,sizeof(double)) ; 
 Cprime = -2*(m+1.) / (eps*eps*eps) ;
// Cprime = -2*(m+1.) ;

 if ((rho<=s-eps)||(s+eps<=-rho))
    res = 0. ;
 else
    {
     if (rho<s+eps)
        borne_sup = rho ;
     else
        borne_sup = s+eps ;
     if (-rho<s-eps)
        borne_inf = s-eps ;
     else 
        borne_inf = -rho ;
     //
     // Calcul de la partie "negative" de l'integrale
     //
     if (borne_inf>=-eps3)
         res1 = 0.;
     else if (borne_sup >= -eps3)
        {
         borne = -eps3 ;
         for (i=0;i<=n-1;i++)
            {
             ti = position_abscisse_f(borne_inf,borne,n,i) ;
             f[2*i] = base_omegaprime(ti,m,s,rho,eps) ;
//             printf("ti=%lf, fi=%lf\n",ti,f[2*i]);
             tsuiv = position_abscisse_f(borne_inf,borne,n,i+1) ;
             tinter = (ti+tsuiv)/2. ;
             f[2*i+1] = base_omegaprime(tinter,m,s,rho,eps) ;
//             printf("ti=%lf, fi=%lf\n",tinter,f[2*i+1]);
            }
         f[2*n] = base_omegaprime(tsuiv,m,s,rho,eps) ;
//         printf("ti=%lf, fi=%lf\n",tsuiv,f[2*n]);
         res1 = simpson_f(f,borne_inf,borne,n) ;
        }
     else 
        {
         borne = borne_sup ;
         for (i=0;i<=n;i++)
            {
             ti = borne_inf+i*(borne-borne_inf)/n ;
             f[i] = base_omegaprime(ti,m,s,rho,eps) ;
            }
         res1 = simpson(f,borne_inf,borne,n) ;
        }
     //
     // Calcul de la partie "positive" de l'integrale
     //
     if (borne_sup<=eps3)
        res2 = 0.;
     else if (borne_inf <= eps3)
        {
         borne = eps3 ;
         for (i=0;i<=n-1;i++)
            {
             ti = position_abscisse_i(borne,borne_sup,n,i) ;
             f[2*i] = base_omegaprime(ti,m,s,rho,eps) ;
//             printf("ti=%lf, fi=%lf\n",ti,f[2*i]);
             tsuiv = position_abscisse_i(borne,borne_sup,n,i+1) ;
             tinter = (ti+tsuiv)/2. ;
             f[2*i+1] = base_omegaprime(tinter,m,s,rho,eps) ;
//             printf("ti=%lf, fi=%lf\n",tinter,f[2*i+1]);
            }
         f[2*n] = base_omegaprime(tsuiv,m,s,rho,eps) ;
//       printf("ti=%lf, fi=%lf\n",tsuiv,f[2*n]);
         res2 = simpson_i(f,borne,borne_sup,n) ;
        }
     else
        {
         borne = borne_inf ;
         for (i=0;i<=n;i++)
            {
             ti = borne+i*(borne_sup-borne)/n ;
             f[i] = base_omegaprime(ti,m,s,rho,eps) ;
            }
         res2 = simpson(f,borne,borne_sup,n) ;
        }
     res = Cprime*(res1+res2) ;
    }
// printf("res=%lf\n",res);
 free(f) ;
 return res ;
}

