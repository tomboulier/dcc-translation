// FanBeam sinogramme
//

#ifndef __FBCONVOL_H
#define __FBCONVOL_H

#include <sinogramme.h>
#include <math.h>

//class FBSinogramme : public Sinogramme //
//
// Operators
// ==========
//
void LambdaFBConvolution(Projection *u,FBSinogramme *a,FBSinogramme *b);
// b(i,.) <- u convol a(i,.) (for all i)
// The dimension of u must be 2*a->NQ()


void LambdaLitSino(FBSinogramme *a,FBSinogramme *b,Scalaire xima,Scalaire yima,Scalaire h,int ntrans);
// ne prend que les donnees de a a partir de la translation n1 jusqu'a la 
// translation n2=b->NQ()-n1

void copy(FBSinogramme *a,FBSinogramme *b);
//copy le sinogramme a dans b

#endif
