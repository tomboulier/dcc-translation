// Paquetage Integ
// ---------------
//
// Implementation de la methode des trapezes et de la methode de Simpson
// pour le calcul des integrales des fonctions reelles
//

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;

double position_abscisse_i(double init, double fin, int n, int i) 
//
// Calcule la position du ieme point pour un maillage "exponentiel"
// avec des points accumules sur init
// Precondition : init et fin sont du meme signe et init < fin
//
{
 double res, exposant ;

 if (init>0.)
   {
    exposant = ((double)(i)) / ((double)(n)) * log(fin/init) ;
    res = init * exp(exposant) ;
//    res = init * pow(10.,exposant) ;
   }
 else
   {
    printf("********* integ/position_abscisse_i : Cas non considere \n") ;
    exit(0) ;
//    res = fin * pow(10.,exposant) ;
   }
 return res ;
}


double position_abscisse_f(double init, double fin, int n, int i) 
//
// Calcule la position du ieme point pour un maillage "exponentiel"
// avec des points accumules sur fin
// Precondition : init et fin sont du meme signe et init < fin
//
{
 double res,exposant ;

 if (init<0.)
   {
    exposant = ((double)(n-i)) / ((double)(n)) * log(init/fin) ;
    res = fin * exp(exposant) ;
//    res = init * pow(10.,exposant) ;
   }
 else
   {
    printf("********* integ/position_abscisse_f : Cas non considere \n") ;
    exit(0) ;
   }
 return res ;
}



double trapeze(double *f, double a, double b, int n)
// 
// Le calcul d'integrale utilise n+1 points, dont a et b
// a = f[0],b = f[n]
// 
{
 double pas, res ;
 int i ;

 res = 0. ;
 pas = (b-a)/n ;
 for (i=0;i<=n-1;i++)
    res = res + pas * (f[i]+f[i+1])/2. ;
 return res ;
}

double simpson(double *f, double a, double b, int n)
// 
// n doit etre pair, le calcul utilise n+1 points, dont a et b
// 
{
 double res1, res2, res3, res, pas ; 
 int n1, i ; // n/2

 if ((2*(n/2)) != n)  // n est impair
   {
    printf("******** integ/simpson : Nb de points doit etre pair \n") ;
    exit(0) ;
   }
 else
   n1 = n/2 ;

 if (b <= a)   
   {
    printf("******** integ/simpson : On doit avoir a < b \n") ;
    exit(0) ;
   }
 res1 = f[0]+f[n] ; 
 res2 = 0. ;
 for (i=0;i<=n1-1;i++)
    res2 = res2 + f[2*i+1] ;
 res3 = 0. ;
 for (i=0;i<=n1-2;i++)
    res3 = res3 + f[2*i+2] ;

 pas = (b-a)/n ;
 res = pas/3. * (res1+4*res2+2*res3) ;
 return res ;
}

double simpson_i(double *f, double a, double b, int n) 
//
// Calcul de l'integrale par la methode de Simpson
// avec un maillage comportant des points s'accumulant sur a
// Precondition : a et b sont du meme signe et a < b
// Le calcul s'appuie sur 2n+1 points f[0],f[1],...,f[2n]
//
{
 double res, ti, tsuiv, pas_i ; 
 int i ;

 if ((2*(n/2)) != n)  // n est impair
   {
    printf("******** integ/simpson_i : Nb de points doit etre pair \n") ;
    exit(0) ;
   }
 if (a*b <= 0.)   // a et b ne sont pas du meme signe
   {
    printf("******** integ/simpson_i : a et b ne sont pas du meme signe \n") ;
    exit(0) ;
   }
 if (b <= a)   
   {
    printf("******** integ/simpson_i : On doit avoir a < b \n") ;
    exit(0) ;
   }
 res = 0. ;
 for (i=0;i<=n-1;i++)
   {
    ti = position_abscisse_i(a,b,n,i) ;
    tsuiv = position_abscisse_i(a,b,n,i+1) ;
    pas_i = tsuiv - ti ;
    res = res + pas_i/6. * (f[2*i]+4*f[2*i+1]+f[2*i+2]) ;
   }
 return res ;
}

double simpson_f(double *f, double a, double b, int n) 
//
// Calcul de l'integrale par la methode de Simpson
// avec un maillage comportant des points s'accumulant sur b
// Precondition : a et b sont du meme signe et a < b
// Le calcul s'appuie sur 2n+1 points f[0],f[1],...,f[2n]
//
{
 double res, ti, tsuiv, pas_f ; 
 int i;

 if ((2*(n/2)) != n)  // n est impair
   {
    printf("******** integ/simpson_f : Nb de points doit etre pair \n") ;
    exit(0) ;
   }
 if (a*b <= 0.)   // a et b ne sont pas du meme signe
   {
    printf("******** integ/simpson_f : a et b ne sont pas du meme signe \n") ;
    exit(0) ;
   }
 if (b <= a)
   {
    printf("******** integ/simpson_f : On doit avoir a < b \n") ;
    exit(0) ;
   }
 res = 0. ;
 for (i=0;i<=n-1;i++)
   {
    ti = position_abscisse_f(a,b,n,i) ;
    tsuiv = position_abscisse_f(a,b,n,i+1) ;
    pas_f = tsuiv - ti ;
    res = res + pas_f/6. * (f[2*i]+4*f[2*i+1]+f[2*i+2]) ;
   }
 return res ;
}

