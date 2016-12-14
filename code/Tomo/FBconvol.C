// Operation sur les FBsinogrammes 
// en particulier retroprojection
//
// LD Dec 93

// (c) Copyright TIMC 1993


// #include <sinogramme.h>
#include <FBsinogramme.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;

// Operateurs
// =============================
//
void LambdaFBConvolution(Projection *u,FBSinogramme *a,FBSinogramme *b)
// b(i,.) <- u convol a(i,.) (for all i)
// The dimension of u must be 2*a->NQ()
{
  int i,j,l,n;
  Scalaire h,somme;
  Scalaire fangle;
#if CHECK
  if ( (BadDim(a,b)) || (BadDim2QmoinsUn(a,u)) )
    fputs("Err : Bad dimmenssions in Convolution\n",stderr);
#endif
  n=a->NQ();
  fangle=a->FA();
  if(n==1)
    {
	 puts("Error : Cannot filter with only 1 fan!");
      return;
    }
  else
    {
      h=fangle/(n-1);
    }

  for (i=0;i<a->NP();i++)
    {
      for (j=5;j<=n-5;j++)
	{
	  somme=0.;
	  for (l=j-4;l<=j+4;l++)
	    somme=somme+(*u)(j-l+n)*(*a)(i,l)*(Scalaire)cos(-fangle/2+
							    (l-1)*(fangle/(n-1)));
	  (*b)(i,j)=somme*h;
	}
    }	
}



void LambdaLitSino(FBSinogramme *a,FBSinogramme *b,Scalaire xima,Scalaire yima,Scalaire h,int ntrans)
// ne prend que les donnees de a a partir de la translation n1 jusqu'a la 
// translation n2=b->NQ()-n1

{
  int ip,j,p,q;
  int n1,numtrans;
  Scalaire sx,sy; // abs and ord of the source
  Scalaire smpx,smpy; // abs and ord of (the source - the current pixel)
  Scalaire nsmp,scal; // Norm**2 of source-pixel, dummy
  Scalaire gamma;
  Scalaire radius,fanangle;

  q=a->NQ();
  p=a->NP();
  radius=a->R();
  fanangle=a->FA();
  for (ip=0;ip<p;ip++)
    {
      sx=radius*cos(a->A(ip));
      sy=radius*sin(a->A(ip));
      smpx=sx-xima;
      smpy=sy-yima;
      scal=smpx*sx+smpy*sy;
      nsmp=smpx*smpx+smpy*smpy;
      gamma = acos(scal/sqrt(nsmp*(sx*sx+sy*sy)));
      sx=radius*cos(a->A(ip)+M_PI/2);
      sy=radius*sin(a->A(ip)+M_PI/2);
      if((sx*xima+sy*yima)>0)
	gamma=-gamma;
      numtrans=(int) ((gamma+(fanangle/2))/h);
      n1=numtrans-ntrans;


      for (j=1;j<=b->NQ();j++)
	(*b)(ip,j)=(*a)(ip,j+n1);
    }
}

void copy(FBSinogramme *a,FBSinogramme *b)
//copy le sinogramme a dans b
{
  int i,j,p,q;

  q=a->NQ();
  p=a->NP();
  for (i=0;i<p;i++)
    {
      for (j=1;j<=q;j++)
	(*b)(i,j)=(*a)(i,j);
    }
}
