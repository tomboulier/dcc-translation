//
// Parallel sinogramm in case of local reconstruction
//

#ifndef __PSCONVOL_H
#define __PSCONVOL_H

#include <sinogramme.h>
#include <math.h>

void PsLitSino(Sinogramme *a,Sinogramme *b,double xima,double yima,double rayon,int nbtransloc);
//
// Construit le sinogramme b a partir des donnees tronquees de a :
// On ne considere que des donnees correspondant a des droites qui 
// coupent la boule de centre (xima,yima) et de rayon 'rayon'
// 

void CopyPar(Sinogramme *a,Sinogramme *b);
// copie le sinogramme a dans b

#endif
