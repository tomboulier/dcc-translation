//
// Parallel sinogramm in case of local reconstruction
//

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sinogramme.h>

using namespace std;

void PsLitSino(Sinogramme *a,Sinogramme *b,double xima,double yima,double rayon,int nbtransloc)
//
// Construit le sinogramme b a partir des donnees tronquees de a :
// On ne considere que des donnees correspondant a des droites qui 
// coupent la boule de centre (xima,yima) et de rayon 'rayon'
// 
{
 int p,q,i,ja,jb ;
 double theta,pja,ctheta,stheta ;
 double hx,hy,C ;

 q=a->NQ();
 p=a->NP();
 for (i=0;i<p;i++)
   {
    theta = a->A(i) ;
    ctheta = cos(theta) ;
    stheta = sin(theta) ;
    jb = 0 ;
    for (ja=1;ja<=q;ja++)
      {
       //
       // Calcul de hx et hy
       //
       pja = a->X(i,ja) ;
       C = yima*ctheta-xima*stheta ;
       hx = pja*ctheta-C*stheta ;
       hy = C*ctheta+pja*stheta ;
       //
       // Fin du calcul de hx et hy
       //
       if ((hx-xima)*(hx-xima)+(hy-yima)*(hy-yima) <= (rayon*rayon))
         {
          jb = jb+1 ;
          (*b)(i,jb)=(*a)(i,ja);
         }
      } // end for ja
   //
   // Verification 
   //
    printf("*");
   } // end for i
 printf(" \n");
}

void CopyPar(Sinogramme *a,Sinogramme *b)
// copie le sinogramme a dans b
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
