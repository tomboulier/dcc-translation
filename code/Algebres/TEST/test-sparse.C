 /* to test the  SparseColMat class:
	LD Nov 94.
 */
 
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream> 

#include <alp.h>
#include <sparsecolmat.h>

using namespace std;

main()
{
  int i,j,dim1,dim2,n;
  Scalaire x;
  SparseColMat  b(7,7,3);
  SparseColMat *a;

  dim1=4;
  dim2=4;
  n=3;
  a= new SparseColMat(dim1,dim2,n);

  (*a)(1,1)=1.;

//  a->Nul(); 
  b.Nul();

  for(i=1;i<=b.DIM1();i++)
    for(j=i;(j<=b.DIM2()) && (j<(i+b.N()));j++) {
	b.Set((Scalaire) 2*j,i,j);
      }
//     b(i,j)=(Scalaire) j;

  std::cout << " Impression du stockage de la matrice : " << std::endl;
  for(i=1;i<=b.DIM1();i++)
    {
      std::cout << " ligne " << i<< std::endl;
      for(j=1;j<=b.N();j++)
	  std::cout << b.JA(i,j)<<" " ;
      std::cout << std::endl;
      for(j=1;j<=b.N();j++)
	  std::cout << b(i,j)<<" " ;
      std::cout << std::endl;
      std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
    }

  std::cout << " Impression de la matrice : " << std::endl;
  for(i=1;i<=b.DIM1();i++) 
    {
      std::cout << " ligne " << i << " : ";
      for(j=1; j<=b.DIM2() ; j++)
	  std::cout << j<<":" << b.Get(i,j) << " " ;
      std::cout << std::endl; 
    }

}
 
