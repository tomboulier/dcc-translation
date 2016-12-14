#ifndef __PSEUDOFCT_H
#define __PSEUDOFCT_H

double a_sq(double q, double m, double eps, int n1, int n2, double epsint, double epsext);
// 
// Facteur A_{sigma,q} 
// n represente le nombre de points utilises pour
// l'evaluation de l'integrale qui definit a_sq
// epsint et epsext, n1 et n2 : parametres definissant les calculs des integrales
// 

double filtre_omegaprime(double m, double s, double ro, double eps, int n);
// 
// Fonction omegaprime_{sigma_ro,eps}(s)
// n+1 = nb de points servant a evaluer l'integrale
// 

#endif

