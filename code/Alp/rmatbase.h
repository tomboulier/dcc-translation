// Linear algebra with C++
// Real matrices
//
// $Header: /home/bainvil/Modules/alp/RCS/rmatbase.h,v 2.2 1995/06/10 13:29:42 bainvil Exp bainvil $

#ifndef __RMATBASE_H
#define __RMATBASE_H

#include <alp/rvecbase.h>
#include <alp/rsvecbase.h>
#include <iostream>

#define min(a,b) ((a)<(b))?(a):(b)

class Matrix
    {
  protected:
    double *x;			// Base address
    int n1,n2;			// n1 rows by n2 columns
    int s1,s2;			// Increments (min(s1,s2)=1)
    //
    // Memory allocation
    // =================
    void Allocate();            // Allocates x from n1 and n2
  public:
    //
    // Constructors
    // ============
    Matrix(const int rows,const int columns);
    Matrix(const int n=3);
    Matrix(const Matrix& a);
    Matrix(const Vector& u);

    //
    // Destructeur
    // ===========
    virtual ~Matrix();

    //
    // Structure access
    // ================
    inline double & operator () (int i1,int i2)
	{ return x[(i1-1)*s1+(i2-1)*s2]; }
    inline double operator () (int i1,int i2) const
	{ return x[(i1-1)*s1+(i2-1)*s2]; }
    inline int R() const { return n1; }
    inline int C() const { return n2; }

    //
    // "conversion" to Vector (to see the matrix as a vector)
    // ======================
    // non const Matrix
    inline SharedVector ShareAsVector()
	{ return SharedVector(x,n1*n2,1); }
    inline friend SharedVector ShareAsVector(Matrix& a)
	{ return SharedVector(a.x,a.n1*a.n2,1); }
    // const matrix
    inline SharedVector ShareAsVector() const
	{ return SharedVector(x,n1*n2,1); }
    inline friend SharedVector ShareAsVector(const Matrix& a)
	{ return SharedVector(a.x,a.n1*a.n2,1); }
    // NOTE : in that case, the vector is NOT allocated, and SHARES the matrix coefficients

    //
    // Diagonal, rows and columns operators
    // ====================================
    // Sharing diagonal, rows and columns of the matrix through vectors
    // const Matrix
    inline friend SharedVector ShareDiag(const Matrix& a)
	{ return SharedVector(a.x,min(a.n1,a.n2),a.s1+a.s2); }
    inline friend SharedVector ShareRow(const Matrix& a,int i)
	{ return SharedVector(a.x+(i-1)*a.s1,a.n2,a.s2); }
    inline friend SharedVector ShareColumn(const Matrix& a,int j)
	{ return SharedVector(a.x+(j-1)*a.s2,a.n1,a.s1); }
    // non const Matrix
    inline friend Vector ShareDiag(Matrix& a)
  { return SharedVector(a.x,min(a.n1,a.n2),a.s1+a.s2); }
    inline friend SharedVector ShareRow(Matrix& a,int i)
	{ return SharedVector(a.x+(i-1)*a.s1,a.n2,a.s2); }
    inline friend SharedVector ShareColumn(Matrix& a,int j)
	{ return SharedVector(a.x+(j-1)*a.s2,a.n1,a.s1); }
    // Copying diagonal, rows and columns to vectors
    friend void GetDiag(const Matrix& a,Vector &u);
    friend void GetRow(const Matrix& a,int i,Vector &u);
    friend void GetColumn(const Matrix& a,int j,Vector &u);
    // Copying diagonal, rows and columns from vectors
    friend void SetDiag(Matrix& a,const Vector &u);
    friend void SetRow(Matrix& a,int i,const Vector &u);
    friend void SetColumn(Matrix& a,int j,const Vector &u);
    
    //
    // Diagonal, rows and columns member functions
    // ===========================================
    // Sharing diagonal, rows and columns of the matrix through vectors
    // non const Matrix
    inline SharedVector ShareDiag()
	{ return SharedVector(x,min(n1,n2),s1+s2); }
    inline SharedVector ShareRow(int i)
	{ return SharedVector(x+(i-1)*s1,n2,s2); }
    inline SharedVector ShareColumn(int j)
	{ return SharedVector(x+(j-1)*s2,n1,s1); }
    // const Matrix
    inline SharedVector ShareDiag() const
	{ return SharedVector(x,min(n1,n2),s1+s2); }
    inline SharedVector ShareRow(int i) const
	{ return SharedVector(x+(i-1)*s1,n2,s2); }
    inline SharedVector ShareColumn(int j) const
	{ return SharedVector(x+(j-1)*s2,n1,s1); }
    // Copying diagonal, rows and columns to vectors
    void GetDiag(Vector &u) const;
    void GetRow(int i,Vector &u) const;
    void GetColumn(int j,Vector &u) const;
    // Copying diagonal, rows and columns from vectors
    void SetDiag(const Vector &u);
    void SetRow(int i,const Vector &u);
    void SetColumn(int j,const Vector &u);
    
    //
    // Initializations
    // ===============
    void Zero();
				// this(i,j) <- 0
    void Id();
				// this(i,j) <- 1 where i==j, and 0 elsewhere
    void Scaling(const double k);
				// this(i,j) <- k where i==j, and 0 elsewhere
    void Rotation2D(const double a);
				// 2D Rotation of angle a
    void Rotation3D(const Vector &u,const double angle);
				// 3D Rotation of UNIT axis u in R3
				// and given angle
/*    void Rotation3D(const Vector &r);
				// 3D Rotation of rotation vector r in R3
*/
    void Rotation3DFromUnitQuaternion(const Vector &q);
				// 3D Rotation associated with
				// the UNIT quaternion q in R4

    //
    // Memory operators
    // ================
    Matrix & operator = (const Matrix &a);
				// this <- a
    friend void Copy(const Matrix &a,Matrix &b);
				// b <- a
    friend void Swap(Matrix &a,Matrix &b);
				// b <-> a

    //
    // Bounds operators
    // ================
    friend void MinAbsIndex(const Matrix &a,int &i,int &j);
    friend void MaxAbsIndex(const Matrix &a,int &i,int &j);
    friend void MinIndex(const Matrix &a,int &i,int &j);
    friend void MaxIndex(const Matrix &a,int &i,int &j);
    friend double MinAbs(const Matrix &a);
    friend double MaxAbs(const Matrix &a);
    friend double Min(const Matrix &a);
    friend double Max(const Matrix &a);

    //
    // Bounds member functions
    // =======================
    void MinAbsIndex(int &i,int &j) const;
    void MaxAbsIndex(int &i,int &j) const;
    void MinIndex(int &i,int &j) const;
    void MaxIndex(int &i,int &j) const;
    double MinAbs() const;
    double MaxAbs() const;
    double Min() const;
    double Max() const;

    //
    // Norms and dot product operators
    // ===============================
    friend double AbsSum(const Matrix &a);
				// <- sum(i=1..r,j=1..c) of |a(i,j)|
    friend double Sum(const Matrix &a);
				// <- sum(i=1..r,j=1..c) of a(i,j)
    friend double Norm2(const Matrix &a);
				// <- sum(i=1..r,j=1..c) of a(i,j)^2
    friend double Norm(const Matrix &a);
				// <- sqrt(Norm2(a))
    friend double Dot(const Matrix &a,const Matrix &b);
				// <- sum(i=1..r,j=1..c) of a(i,j)*b(i,j)
				//    = Trace of a^T * b
    friend double Trace(const Matrix& a);
				// <- sum(i=1..min(r,c)) of a(i,i)

    //
    // Norms and dot product member functions
    // ======================================
    double AbsSum() const;
				// <- sum(i=1..r,j=1..c) of |this(i,j)|
    double Sum() const;
				// <- sum(i=1..r,j=1..c) of this(i,j)
    double Norm2() const;
				// <- sum(i=1..r,j=1..c) of this(i,j)^2
    double Norm() const;
				// <- sqrt(this.Norm2())
    double Trace() const;
				// <- sum(i=1..min(r,c)) of this(i,i)

    //
    // Arithmetic matrix-vector operators
    // ==================================
    friend void Image(const Vector &u,const Matrix &a,Vector &v);
				// v <- a * u
    friend void Image(const Vector &u,
		      double alpha,const Matrix &a,
		      double beta,const Vector &b,Vector &v);
				// v <- alpha * a * u + beta * b
    friend void TransImage(const Vector &u,const Matrix &a,Vector &v);
				// v <- a^T * u
    friend void TransImage(const Vector &u,
			   double alpha,const Matrix &a,
			   double beta,const Vector &b,Vector &v);
				// v <- alpha * a^T * u + beta * b
    friend void Update(const Matrix &a,
		       double alpha,const Vector &x,const Vector &y,Matrix &b);
				// b <- a + alpha * x * y^T
    friend double Quad(const Matrix &a,const Vector& x,const Vector& y);
				// <- x^T * a * y

    //
    // Arithmetic matrix-vector member functions
    // =========================================
    void Image(const Vector& u,Vector& v);
				// v <- this * u
    void Image(const Vector &u,double alpha,
	       double beta,const Vector &b,Vector &v);
				// v <- alpha * this * u + beta * b
    void TransImage(const Vector &u,Vector &v);
				// v <- this^T * u
    void TransImage(const Vector &u,double alpha,
		    double beta,const Vector &b,Vector &v);
				// v <- alpha * this^T * u + beta * b
    void Update(double alpha,const Vector &x,const Vector &y);
				// this <- this + alpha * x * y^T
    double Quad(const Vector& x,const Vector& y);
				// <- x^T * this * y
    // Functions to avoid a warning (see the end of rsvcbase.h for comments)
    inline void Image(const Vector& u,SharedVector v)
	{ Image(u,(Vector&)v); }
    inline void Image(const Vector &u,double alpha,
		      double beta,const Vector &b,SharedVector v)
	{ Image(u,alpha,beta,b,(Vector&)v); }
    inline void TransImage(const Vector &u,SharedVector v)
	{ TransImage(u,(Vector&)v); }
    inline void TransImage(const Vector &u,double alpha,
			   double beta,const Vector &b,SharedVector v)
	{ TransImage(u,alpha,beta,b,(Vector&)v); }

    //
    // Transposition (takes the time of swapping 2 integers : fast)
    // =============
    friend void Transpose(const Matrix &a,Matrix &b);
				// b <- a^T
    void Transpose();
				// this <- this^T

    //
    // Matrix additions operators
    // ==========================
    friend void Add(const Matrix &a,const Matrix &b,Matrix &c);
				// c <- a + b
    friend void Sub(const Matrix &a,const Matrix &b,Matrix &c);
				// c <- a - b
    friend void Add(double alpha,const Matrix &a,
		    double beta,const Matrix &b,Matrix &c);
				// c <- alpha * a + beta * b
    friend void TransAdd(const Matrix &a,const Matrix &b,Matrix &c);
				// c <- a + b^T
    friend void TransAdd(double alpha,const Matrix &a,
			 double beta,const Matrix &b,Matrix &c);
				// c <- alpha * a + beta * b^T

    //
    // Matrix additions member functions
    // =================================
    Matrix& operator += (const Matrix &a);
				// this <- this + a
    Matrix& operator -= (const Matrix &a);
				// this <- this - a
    void Add(const Matrix &a);
				// this <- this + a
    void Sub(const Matrix &a);
				// this <- this - a
    void Add(double k,const Matrix &a);
				// this <- this + k * a
    void Combine(double k,double l,const Matrix &a);
				// this <- k * this + l * a

    //
    // Matrix-scalar products operators
    // ================================
    friend void Mul(double k,const Matrix &a,Matrix &b);
				// b <- k * a
    friend void Div(double k,const Matrix &a,Matrix &b);
				// b <- k * a

    //
    // Matrix-scalar products member functions
    // =======================================
    Matrix& operator *= (double k);
				// this <- k * this
    Matrix& operator /= (double k);
				// this <- (1/k) * this
    void Mul(double k);
				// this <- k * this
    void Div(double k);
				// this <- (1/k) * this

    //
    // Matrix-matrix multiplication operators
    // ======================================
    friend void Mul(const Matrix &a,const Matrix &b,Matrix &c);
				// c <- a * b
    friend void Mul(double alpha,const Matrix &a,const Matrix &b,
		    double beta,const Matrix &c,Matrix &m);
				// m <- alpha * a * b + beta * c
    friend void RightTransMul(const Matrix &a,const Matrix &b,Matrix &c);
				// c <- a * b^T
    friend void LeftTransMul(const Matrix &a,const Matrix &b,Matrix &c);
				// c <- a^T * b
    friend void RightTransMul(double alpha,const Matrix &a,const Matrix &b,
			      double beta,const Matrix &c,Matrix &m);
				// m <- alpha * a * b^T + beta * c
    friend void LeftTransMul(double alpha,const Matrix &a,const Matrix &b,
			     double beta,const Matrix &c,Matrix &m);
				// m <- alpha * a * b^T + beta * c

    //
    // Matrix-matrix multiplication member functions
    // =============================================
    void RightMul(const Matrix &a);
				// this <- this * a
    void LeftMul(const Matrix &a);
				// this <- a * this
    void RightMul(double alpha,const Matrix &a,double beta,const Matrix &b);
				// this <- alpha * this * a + beta * b
    void LeftMul(double alpha,const Matrix &a,double beta,const Matrix &b);
				// this <- alpha * a * this + beta * b

    //
    // IO operators
    // ============

    friend std::ostream & operator << (std::ostream &o,const Matrix &a);
    // ascii
    friend std::ostream & TeX(std::ostream &o,const Matrix &a);
    // ascii to be included in a latex documents
    };

// See the comments in rsvecbase.h for the following functions :

inline void GetDiag(const Matrix& a,SharedVector u)
    { GetDiag(a,(Vector&)u); }
inline void GetRow(const Matrix& a,int i,SharedVector u)
    { GetRow(a,i,(Vector&)u); }
inline void GetColumn(const Matrix& a,int j,SharedVector u)
    { GetColumn(a,j,(Vector&)u); }
inline void Image(const Vector &u,const Matrix &a,SharedVector v)
    { Image(u,a,(Vector&)v); }
inline void Image(const Vector &u,
		  double alpha,const Matrix &a,
		  double beta,const Vector &b,SharedVector v)
    { Image(u,alpha,a,beta,b,(Vector&)v); }
inline void TransImage(const Vector &u,const Matrix &a,SharedVector v)
    { TransImage(u,a,(Vector&)v); }
inline void TransImage(const Vector &u,
		       double alpha,const Matrix &a,
		       double beta,const Vector &b,SharedVector v)
    { TransImage(u,alpha,a,beta,b,(Vector&)v); }

#undef min
#endif
