// Paquetage Integ
// ---------------
//
// Implementation de la methode des trapezes
// pour le calcul des integrales des fonctions reelles
//

#ifndef __INTEG_H
#define __INTEG_H

#include <math.h>

double position_abscisse_i(double init, double fin, int n, int i) ;
// 
// Calcule la position du ieme point pour un maillage "exponentiel"
// avec des points accumules sur init
// Precondition : init et fin sont du meme signe et init < fin
// 

double position_abscisse_f(double init, double fin, int n, int i) ;
// 
// Calcule la position du ieme point pour un maillage "exponentiel"
// avec des points accumules sur fin
// Precondition : init et fin sont du meme signe et init < fin
// 

double trapeze(double *f, double a, double b, int n) ;
// 
// Le calcul d'integrale utilise n+1 points, dont a et b
// et suppose que le maillage est regulier
// 

double simpson(double *f, double a, double b, int n) ;
// 
// Le calcul d'integrale utilise 2n+1 points, dont a et b
// et suppose que le maillage est regulier
// 

double simpson_i(double *f, double a, double b, int n) ;
// 
// Calcul de l'integrale par la methode de Simpson 
// avec un maillage comportant des points s'accumulant sur a
// Precondition : a et b sont du meme signe et a < b
// Le calcul s'appuie sur 2n+1 points f[0],f[1],...,f[2n]
// 

double simpson_f(double *f, double a, double b, int n) ;
// 
// Calcul de l'integrale par la methode de Simpson 
// avec un maillage comportant des points s'accumulant sur b
// Precondition : a et b sont du meme signe et a < b
// Le calcul s'appuie sur 2n+1 points f[0],f[1],...,f[2n]
// 

#endif
