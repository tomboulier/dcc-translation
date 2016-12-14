 /* to test the  SparseColMat class:
	LD Nov 94.
 */
 
#include <cstdlib>
#include <cstdio>
#include <iostream> 
#include <alp.h>
#include <sparsecolmat.h>

using namespace std;


main()
{
  int i,j,dim1,dim2,n;
//  Scalaire s;
  SparseColMat *a;
  Vector *x,*y,*y1,*y2;
// pour le gradient conjugue
  Scalaire arret;
  int itermin,itermax,iter;
  Vector *r;
  int verbose=1;
//
  dim1=8; // doit être pair pour le test de ProdSubSparseColMatVector
  dim2=8;
  n=2;
  a= new SparseColMat(dim1,dim2,n);
  a->Nul();

  SparseColMat *a1, *a2; // test de ProdSubSparseColMatVector et de ProdtransSubSparseColMatVector 
  a1= new SparseColMat(dim1/2,dim2,n); // dim1 est pair !!
  a2= new SparseColMat(dim1/2,dim2,n); // dim1 est pair !!
  a1->Nul();a2->Nul();


  for(i=1;i<=a->DIM1();i++)
    for(j=i;(j<=a->DIM2()) && (j<(i+a->N()));j++)
      a->Set((Scalaire) j,i,j);
  for(i=1;i<=a1->DIM1();i++)
    for(j=i;(j<=a1->DIM2()) && (j<(i+a->N()));j++)
      a1->Set((Scalaire) j,i,j);
  int lineshift=dim1/2;
  for(i=1;i<=a2->DIM1();i++)
    for(j=i+lineshift;(j<=a2->DIM2()) && (j<(i+lineshift+a->N()));j++)
      a2->Set((Scalaire) j,i,j);

  std::cout << " Impression du stockage de la matrice a: " << std::endl;
  for(i=1;i<=a->DIM1();i++)    {
    std::cout << " ligne " << i<< std::endl;
    for(j=1;j<=a->N();j++)
      std::cout << a->JA(i,j)<<" " ;
    std::cout << std::endl;
    for(j=1;j<=a->N();j++)
      std::cout << (*a)(i,j)<<" " ;
    std::cout << std::endl;
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
  }

  std::cout << " Impression de la matrice a : " << std::endl;
  for(i=1;i<=a->DIM1();i++)  {
    std::cout << " ligne " << i << " : ";
    for(j=1; j<=a->DIM2() ; j++)
      std::cout << a->Get(i,j) << " " ;
      //      std::cout << j<<":" << a->Get(i,j) << " " ;
    std::cout << std::endl; 
  }

  Vector *sumoflines;
  Vector *sumofcols;
  sumofcols=new Vector(dim2);
  sumoflines=new Vector(dim1);
  LineSumOfSparseColMat(a,sumoflines);
  ColSumOfSparseColMat(a,sumofcols);

  std::cout << " Somme des colonnes de la matrice a : " << std::endl;
  for(j=1; j<=a->DIM2() ; j++)
    std::cout << (*sumofcols)(j) << " " ;
  std::cout << std::endl; 
  std::cout << " Somme des lignes de la matrice a : " << std::endl;
  for(i=1;i<=a->DIM1();i++)  
    std::cout << (*sumoflines)(i) << std::endl;
 

  std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
  std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;

  std::cout << " Impression du stockage de la matrice a1: " << std::endl;
  for(i=1;i<=a1->DIM1();i++)    {
    std::cout << " ligne " << i<< std::endl;
    for(j=1;j<=a1->N();j++)
      std::cout << a1->JA(i,j)<<" " ;
    std::cout << std::endl;
    for(j=1;j<=a1->N();j++)
      std::cout << (*a1)(i,j)<<" " ;
    std::cout << std::endl;
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
  }

  std::cout << " Impression de la matrice a1 : " << std::endl;
  for(i=1;i<=a1->DIM1();i++)  {
    std::cout << " ligne " << i << " : ";
    for(j=1; j<=a1->DIM2() ; j++)
      std::cout << a1->Get(i,j) << " " ;
      //      std::cout << j<<":" << a1->Get(i,j) << " " ;
    std::cout << std::endl; 
  }
  std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
  std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;

  std::cout << " Impression du stockage de la matrice a2: " << std::endl;
  for(i=1;i<=a2->DIM1();i++)    {
    std::cout << " ligne " << i<< std::endl;
    for(j=1;j<=a2->N();j++)
      std::cout << a2->JA(i,j)<<" " ;
    std::cout << std::endl;
    for(j=1;j<=a2->N();j++)
      std::cout << (*a2)(i,j)<<" " ;
    std::cout << std::endl;
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
  }

  std::cout << " Impression de la matrice a2 : " << std::endl;
  for(i=1;i<=a2->DIM1();i++)  {
    std::cout << " ligne " << i << " : ";
    for(j=1; j<=a2->DIM2() ; j++)
      std::cout << a2->Get(i,j) << " " ;
      //      std::cout << j<<":" << a2->Get(i,j) << " " ;
    std::cout << std::endl; 
  }

  std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
  std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
  std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;

  x=new Vector(dim2);
  y=new Vector(dim1);
  y1=new Vector( dim1/2);
  y2=new Vector( dim1/2);

  for(j=1;j<=a->DIM2();j++) 
    (*x)(j)=1.;
  std::cout << " initialisation de *x a 1" << std::endl;
  ProdSparseColMatVector(a,x,y);
  std::cout << " Somme sur les colonnes de a (utilisant ProdSparseColMatVector) : ";
  std::cout << std::endl;
  for(i=1;i<=a->DIM1();i++) 
    std::cout << " y =SUM_j a_{i,j} (" << i << ") : "<< (*y)(i) << std::endl;

  ProdSparseColMatVector(a1,x,y1);
  ProdSparseColMatVector(a2,x,y2);
  std::cout << " Somme sur les colonnes de a1 et a2 donc a (utilisant ProdSparseColMatVector) : ";
  std::cout << std::endl;
  for(i=1;i<=a1->DIM1();i++) 
    std::cout << " y1 =SUM_j a_{i,j} (" << i << ") : "<< (*y1)(i) << std::endl;
  for(i=1;i<=a2->DIM1();i++) 
    std::cout << " y2 =SUM_j a_{i,j} (" << i+lineshift << ") : "<< (*y2)(i) << std::endl;
  
  for(i=1;i<=a->DIM1();i++) 
    (*y)(i)=1;
  std::cout << " initialisation de *y a 1" << std::endl;
  ProdtransSparseColMatVector(a,y,x);
  std::cout << " Somme sur les colonnes (utilisant ProdtransSparseColMatVector) : ";
  std::cout << std::endl;
  for(j=1;j<=a->DIM2();j++) 
      std::cout << " x =SUM_i a_{i,j} (" << j << ") : "<< (*x)(j) << std::endl;

  for(i=1;i<=a1->DIM1();i++) 
    (*y1)(i)=1;
  for(i=1;i<=a2->DIM1();i++) 
    (*y2)(i)=1;
  ProdtransSparseColMatVector(a1,y1,x);
  AddProdtransSparseColMatVector(a2,y2,x);
  std::cout << " Somme sur les colonnes (utilisant y1 et y2 et AddProdtransSparseColMatVector) : ";
  std::cout << std::endl;
  for(j=1;j<=a->DIM2();j++) 
      std::cout << " x =SUM_i a_{i,j} (" << j << ") : "<< (*x)(j) << std::endl;
  
  exit(0);

  std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
  std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
  std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;

  std::cout << " Test du Gradient Conjugue " << std::endl;
  
  delete a;
  n=dim1;
  a= new SparseColMat(dim1,dim2,n);
  a->Nul(); 
  
  for(i=1;i<=a->DIM1();i++)
    for(j=1;(j<=a->DIM2()) && (j<(i+a->N()));j++)  {
      if(i==j)
	a->Set((Scalaire) (Scalaire) a->DIM1(),i,j);
      else
	a->Set((Scalaire) -1.,i,j);
    }
  
  for(i=1;i<=a->DIM1() ;i++)
    (*y)(i)=(Scalaire) 1.;
  for(j=1;j<=a->DIM2();j++)
    (*x)(j)=(Scalaire) 0.;
  
  itermin=0;
  itermax=a->DIM2();
  r= new Vector(dim2);
  r->Zero(); 


  std::cout << " BEFORE SLSCG : " << std::endl;

  SLSCG(a,x,y,arret,itermin,itermax,iter,r,verbose);

  std::cout << " Impression de la matrice : " << std::endl;
  for(i=1;i<=a->DIM1();i++) 
    {
      std::cout << " ligne " << i << " : ";
      for(j=1; j<=a->DIM2() ; j++)
	  std::cout << j<<":" << a->Get(i,j) << " " ;
      std::cout << std::endl; 
    }
  std::cout << " Impression du second menbre: " << std::endl;
  for(i=1;i<=a->DIM1();i++) 
    std::cout << " ligne " << i << " : " << (*y)(i)<<std::endl;
  std::cout << " Impression de la solution : " << std::endl;
  for(i=1;i<=a->DIM2();i++) 
    std::cout << " ligne " << i << " : " << (*x)(i)<<std::endl;
  

  std::cout << " ********************************** " << std::endl;
  std::cout << " ********************************** " << std::endl;
  std::cout << " ********************************** " << std::endl;
  std::cout << " ********************************** " << std::endl;

  a->Nul(); 

  for(i=1;i<=a->DIM1();i++)
    for(j=1;(j<=a->DIM2()) && (j<(i+a->N()));j++)  {
      if(i==j){
	a->Set((Scalaire) a->DIM1()-i+1,i,j);
      }
      else{
	if(i<j)
	  a->Set((Scalaire) -1.,i,j);
      }
    }
  
  for(i=1;i<=a->DIM1() ;i++)
    (*y)(i)=(Scalaire) 1.;
  for(j=1;j<=a->DIM2();j++)
    (*x)(j)=(Scalaire) 0.;
  
  itermin=0;
  itermax=a->DIM2();
  r= new Vector(dim2);
  r->Zero(); 
  

  std::cout << " BEFORE SLSCG : " << std::endl;

  SLSCG(a,x,y,arret,itermin,itermax,iter,r,verbose);

  std::cout << " Impression de la matrice : " << std::endl;
  for(i=1;i<=a->DIM1();i++)   {
    std::cout << " ligne " << i << " : ";
    for(j=1; j<=a->DIM2() ; j++)
      std::cout << j<<":" << a->Get(i,j) << " " ;
    std::cout << std::endl; 
  }
  std::cout << " Impression du second menbre: " << std::endl;
  for(i=1;i<=a->DIM1();i++) 
    std::cout << " ligne " << i << " : " << (*y)(i)<<std::endl;
  std::cout << " Impression de la solution : " << std::endl;
  for(i=1;i<=a->DIM2();i++) 
    std::cout << " ligne " << i << " : " << (*x)(i)<<std::endl;
  std::cout << " ********************************** " << std::endl;
  std::cout << " ********************************** " << std::endl;
  std::cout << " ********************************** " << std::endl;
  std::cout << " ********************************** " << std::endl;
  std::cout << "         AVEC SLSCGreg3d            " << std::endl; 
  std::cout << " ********************************** " << std::endl;
  std::cout << " ********************************** " << std::endl;
  a->Nul(); 

  for(i=1;i<=a->DIM1();i++)
    for(j=1;(j<=a->DIM2()) && (j<(i+a->N()));j++)  {
      if(i==j)
	a->Set((Scalaire) a->DIM1()-i+1,i,j);
      else
	if(i<j)
	  a->Set((Scalaire) -1.,i,j);
    }

  for(i=1;i<=a->DIM1() ;i++)
    (*y)(i)=(Scalaire) 1.;
  for(j=1;j<=a->DIM2();j++)
    (*x)(j)=(Scalaire) 0.;
  
  itermin=0;
  itermax=a->DIM2();
  r= new Vector(dim2);
  r->Zero(); 


  std::cout << " BEFORE SLSCGreg3 : " << std::endl;

  int *isout;
  isout = new int[8];
  SLSCGreg3d(a,x,y,2,2,2,0.,isout,0,arret,itermin,itermax,iter,r,verbose);

  std::cout << " Impression de la matrice : " << std::endl;
  for(i=1;i<=a->DIM1();i++) {
    std::cout << " ligne " << i << " : ";
    for(j=1; j<=a->DIM2() ; j++)
      std::cout << j<<":" << a->Get(i,j) << " " ;
    std::cout << std::endl; 
  }
  std::cout << " Impression du second menbre: " << std::endl;
  for(i=1;i<=a->DIM1();i++) 
    std::cout << " ligne " << i << " : " << (*y)(i)<<std::endl;
  std::cout << " Impression de la solution : " << std::endl;
  for(i=1;i<=a->DIM2();i++) 
    std::cout << " ligne " << i << " : " << (*x)(i)<<std::endl;
      
}
  


