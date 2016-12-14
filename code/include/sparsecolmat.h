// Algebre en dimension n
// Partie II : 
//
// LD Nov 94 : creation

// (c) Copyright TIMC 1993

#ifndef __ALGNSPARSECOLMAT_H
#define __ALGNSPARSECOLMAT_H

#include <cmath>
#include <typescalaire.h>
#include <alp.h> 
// #include <alp.h> modif 26/02/98

//
// ========
// SparseColMat
// ========
// ELPACK ITPACK sparse matrix storage format.
//
// we denote in the following A the "true" complete matrix
// and *a the array for its storage on the computer
// We use the following simple storage: all nonnul entries are stored
// if A(i,j)<>0 then
//    there exists a single index k such that
//     *ja(i,k)=j and *a(i,k)=A(i,j)
//
// we use the following rules:
// *ja(i,0) is the number of non zero entries of A in line number i. 
// let k be the first index such that *ja(,k)=0, then A has k-1 non zero 
//     entries on the ligne i.
//
class SparseColMat
{
   protected:
     int dim1,dim2,n;	  // true dimension A(dim1,dim2)
     Scalaire *a;	// matrix entries
     int *ja;	// indirection
                // *a(i,k)= A(i,*ja(i,k))
     Scalaire **lina;	// pointer inside a to accelerate ()
     int **linja;	// pointer inside ja to accelerate ()
		// 
     void Allocate();	// Allocate *a, *ja with dimension dim1*n
 public:
	//
	// Constructors et destructors
	// =============================
	//
    SparseColMat(const int n1=1,const int n2=1,const int width=1);
    // allocate  a sparse matrix of true dimensions n1 x n2   
    // and of real dimension n1*width
    //
//    virtual ~SparseColMat(); // free a SparseColMat
    ~SparseColMat(); // free a SparseColMat
       //
       // structure acces
       // ===============
       //
     inline Scalaire & operator () (const int i,const int j)
     // Acces to *a(i,j) ( i must be >=1 and <=dim1,
     //                    j must be >=1 and <=n,
     // SHOULD NEVER BE USED BY USER (bug to correct TO BE PROTECTED...)
    {
     return lina[i][j];
    }
    inline const Scalaire  operator () (const int i,const int j) const
     // Acces to *a(i,j) ( i must be >=1 and <=dim1,
     //                    j must be >=1 and <=n,
     {
      return lina[i][j];
     }
	void Set(const Scalaire x,const int i,const int j); 
                     //  i must be >=1 and <=dim1,  must be >=1 and <=dim2
                     // A(i,j)=x
                     // i.e. let k be the first integer
                     // such that *ja(i,k)=0 then
		     //         *ja(i,k)=l and *a(i,k)=x
		     //  (if *ja(i,k) <>0 until k=n then ERROR and STOP.
       Scalaire Get(const int i,const int j); 
		     // return A(i,j)
		     // i.e., *a(i,k) if \exists k such that *ja(i,k)=j
                     //       0  otherwise
	int JA(const int i,const int k);
                     // return *ja(i,k) 
	int N();// Dimension 
	int DIM1 ();// first dimension
	int DIM2 ();// second dimension

//
// Initialisations
// ===============
//
    void Nul();  // nul. matrix [both *a (0,0,...,0) and *ja]
    void Indentity();//  *a (1,1,...,1,1,1,...,1) 
                     //  *ja(i,1)=i and *ja(i,k)=0 if k<>1
};

//
// Operators
// ==========
//

void ProdSparseColMatVector(SparseColMat *a,Vector *x,Vector *y);
// y <- A x
void LineSumOfSparseColMat(SparseColMat *a,Vector *sumoflines);
// sumoflines<- A (1,1,...,1)^t 
// A is dim1xdim2 : the result is the vector sumoflines of dim1 such that 
//    $sumoflines(i)=\sum_j a_{i,j}$ 
void ProdtransSparseColMatVector(SparseColMat *a,Vector *x,Vector *y);
// y <- A^t x
void ColSumOfSparseColMat(SparseColMat *a,Vector *sumofcols);
// sumofcols<- A^t (1,1,...,1)^t 
// A is dim1xdim2 : the result is the vector sumofcols of dim2 such that 
//    $sumofcols(j)=\sum_i a_{i,j}$ 
void AddProdtransSparseColMatVector(SparseColMat *a,Vector *x,Vector *y);
// y <- y+A^t x
void ProdSubSparseColMatVector(SparseColMat *a,Vector *x,Vector *y, int lineshift);
// A is a submatrix of large matrix M of nbline starting at line 1+lineshift finishing at  lineshift+a->DIM1()
//  y[lineshift : lineshift+a->DIM1()] <- A x   or 
//      <=>    y[lineshift : lineshift+a->DIM1()] = M[lineshift : lineshift+a->DIM1(), :] x 
void ProdtransSubSparseColMatVector(SparseColMat *a,Vector *x,Vector *y, int lineshift);
// A is a submatrix of large matrix M of nbline starting at line 1+lineshift finishing at  lineshift+a->DIM1()
//  y <- A^t x   or 
//      <=>    y = M[lineshift : lineshift+a->DIM1(), :]^t x[lineshift : lineshift+a->DIM1()]
// y <- A^t x

void SLSCG(SparseColMat *a,Vector *x,Vector *b,
	   Scalaire arret,int itermin,int itermax,int & iter,
	   Vector *r,int verbose);
     // Sparse Least Squares Cojugate Gradient.
     // allows us to solve a sparse linear least squares 
     // Min_{x} || A x - b ||
     //
     // *a is a SparseColMat of dimension DIM1 x DIM2
     // *x is a Vector of dimension DIM2
     // *b is a Vector  of dimension DIM1
     // *r,is a  vector of dimension DIM2 (residual A^t A x - A^t b )
     // arret : STOP criterion
     //         if ||*r(k)|| <= arret or if ||*r(k)||/||*r(0)||
     //           where *r(k) is the residual at iteration k.
     // verbose ==1 then the program is verbose...
     // 

void SLSCGreg3d(SparseColMat *a,Vector *x,Vector *b, 
		int dim1, int dim2, int dim3, Scalaire tau,
		int *out, int proj,
                Scalaire arret,int itermin,int itermax,int & iter,
		Vector *r,int verbose);
     // Sparse Least Squares Conjugate Gradient with Regularisation
     // allows us to solve a sparse linear least squares 
     // Min_{x} || A P x - b || +tau x^t PDP x
     // where D is the 3d Laplacian discretisation
     // x is suppose to be a 3D discretisation stored in the following way:
     // x(i*dim2*dim3+j*dim3+k)=x(i,j,k)
     // Dx=6x(i,j,k)-(x(i+1,j,k)+x(i-1,j,k)+ .....)
     // P is a projector onto support information 
     //       (equal to Identity if noproj=1)
     //
     // *a is a SparseColMat of dimension DIM1 x DIM2
     // *x is a Vector of dimension DIM2
     // *out is an integer array of dimension DIM2 for support information 
     //     introduction: if proj is equal to 0 then *out is suppose 
     //     to be zero. Else out[i-1]=1 means that (*x)(i-1) is nul
     //            (*out) is the complementary of the projector P
     // *b is a Vector of dimension DIM1
     // *r,is a  vector of dimension DIM2 [residual (A^t A +tau D)x - A^t b ]
     // arret : STOP criterion
     //         if ||*r(k)|| <= arret or if ||*r(k)||/||*r(0)||
     //           where *r(k) is the residual at iteration k.
     // verbose ==1 then the program is verbose...

void AddLapla3D(Vector *x, int dim1,int dim2,int dim3, Scalaire tau,Vector *y);
//  *y=*y+tau D *x


void SLSCGreg3dpg(SparseColMat *a,Vector **x,Vector **b, 
		int dim1, int dim2, int dim3, Scalaire tau,
                int *out, int proj,
		Scalaire arret,int itermin,int itermax,int & iter,
		  Vector **r,int verbose, Scalaire ParallelRatio=1);
     // FOR PARALLEL GEOMETRIES !!!
     // Sparse Least Squares Conjugate Gradient with Regularisation
     // allows us to solve a sparse linear least squares 
     // Min_{x} || A P x - b || +tau x^t PDP x
     // where D is the 3d Laplacian discretisation
     // x is suppose to be a 3D discretisation stored in the following way:
     // (*(x[j]))(k*dim1+i+1)=x(i,j,k)
     // Dx=6x(i,j,k)-(x(i+1,j,k)+x(i-1,j,k)+ .....)
     // P is a projector onto support information 
     //       (equal to Identity if noproj=1)
     //
     // *a is a SparseColMat of dimension DIM1A x DIM2A
     // **x are dim1  Vectors of dimension DIM2A=dim3*dim2
     // *out is an integer array of dimension dim1*DIM2A (dim1*dim2*dim3)
     //     for support  information 
     //     introduction: if proj is equal to 0 then *out is suppose 
     //     to be zero. Else out[j*dim2+i-1]=1 means that (*(x[j]))(i-1) is nul
     //            (*out) is the complementary of the projector P
     //
     // Thu Oct  9 19:04:38 CEST 2008
     // *b are NbSclice  Vector of dimension DIM1A
     // *r are NbSclice  Vector of dimension DIM1A [residual (A^t A +tau D)x - A^t b ]
     //        where NbSlice=dim1*ParallelRatio (see below)
     // arret : STOP criterion
     //         if ||*r(k)|| <= arret or if ||*r(k)||/||*r(0)||
     //           where *r(k) is the residual at iteration k.
     // verbose ==1 then the program is verbose...
     //
     // Thu Oct  9 19:04:38 CEST 2008
     // ParallelRatio is the sampling ratio between NbSlice/dim1


void AddLapla3Dgp(Vector **x, int dim1,int dim2,int dim3, 
		Scalaire tau,Vector **y);
     // FOR PARALLEL GEOMETRIES !!!
     //  *y=*y+tau D *x

void SLSCGreg2d(SparseColMat *a,Vector *x,Vector *b, 
		int dim1, int dim2, Scalaire tau,
		int *out, int proj,
                Scalaire arret,int itermin,int itermax,int & iter,
		Vector *r,int verbose);
     // Sparse Least Squares Conjugate Gradient with Regularisation
     // allows us to solve a sparse linear least squares 
     // Min_{x} || A P x - b || +tau x^t PDP x
     // where D is the 2d Laplacian discretisation
     // x is suppose to be a 2D image: discretisation stored in the following way:
     // x(i*dim2+j)=x(i,j)
     // Dx=4x(i,j)-(x(i+1,j)+x(i-1,j)+x(i,j+1)+x(i,j-1))
     // P is a projector onto support information 
     //       (equal to Identity if noproj=1)
     //
     // *a is a SparseColMat of dimension DIM1 x DIM2
     // *x is a Vector of dimension DIM2
     // *out is an integer array of dimension DIM2 for support information 
     //     introduction: if proj is equal to 0 then *out is suppose 
     //     to be zero. Else out[i-1]=1 means that (*x)(i-1) is nul
     //            (*out) is the complementary of the projector P
     // *b is a Vector of dimension DIM1
     // *r,is a  vector of dimension DIM2 [residual (A^t A +tau D)x - A^t b ]
     // arret : STOP criterion
     //         if ||*r(k)|| <= arret or if ||*r(k)||/||*r(0)||
     //           where *r(k) is the residual at iteration k.
     // verbose ==1 then the program is verbose...

void AddLapla2D(Vector *x, int dim1,int dim2, Scalaire tau,Vector *y);
//  *y=*y+tau D *x

void InnerOSEM(SparseColMat *a,Vector *x, Vector *y,
	       Vector *ax, Vector *ws, Vector *atransws, 
	       Vector *sumofcols, int verbose);
// compute an inner loop OSEM update operation i.e. 
//   a projection ax= A x
//   a back projection x=x.* A^t(y./ax)  ./ A^t (1,1,...,1)
// a is a sparse rectangular (sub) matrix A
// x is a positive vector to be updated
// y contain the (sub)data
// ax is a working space of dimension of y: ax= A x
// sumofcol is a working space of dimension of x for the sum of column of A
//       sumofcol =  A^t (1,1,...,1)
// 
//  

#endif

