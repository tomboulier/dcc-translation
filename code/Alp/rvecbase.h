// Linear algebra with C++
// Real vectors
//
// EB 93 94 95

#ifndef __RVECBASE_H
#define __RVECBASE_H

#include <cmath>
#include <iostream>
#include <alp/basic.h>

class SharedVector;

class Vector
    {
  protected:
    //
    // Memory representation
    // =====================
    double *x;                  // Base address
    int n;			// Size of the vector
    int s;			// Increment
    int shared;			// Flag true if x is not allocated
    // 1st component is *x
    // 2nd component is *(x+s)
    // ...
    // nth component is *(x+(n-1)*s)
    //
    // Memory allocation
    // =================
    void Allocate();            // Allocates x knowing n
  public:
    //
    // Constructors
    // ============
    Vector(int n=3);
				// new unitialized vector of dimension n
    Vector(double x1,double x2);
				// new vector of dim 2 <- (x1,x2)
    Vector(double x1,double x2,double x3);
				// new vector of dim 3 <- (x1,x2,x3)
    Vector(double x1,double x2,double x3,double x4);
				// new vector of dim 4 <- (x1,x2,x3,x4)
    Vector(const Vector &u);
				// new vector of dim dim(u) <- u
    //
    // Destructor
    // ==========
    virtual ~Vector();
				// delete this

    //
    // Initializations and structure access
    // ====================================
    void Zero();
				// this <- (0,0,...,0)
    void Set(double x1,double x2);
    void Set(double x1,double x2,double x3);
    void Set(double x1,double x2,double x3,double x4);
				// Set the 2,3,4 first coordinates
    void Get(double& x1,double& x2) const;
    void Get(double& x1,double& x2,double& x3) const;
    void Get(double& x1,double& x2,double& x3,double& x4) const;
				// Get the 2,3,4 first coordinates
    inline double& operator () (int i)
      { return (s==1)?(x[i-1]):(x[(i-1)*s]); }
				// <- reference to this(i)
    inline double operator () (int i) const
      { return (s==1)?(x[i-1]):(x[(i-1)*s]); }
				// <- this(i) for const vectors
    inline int N() const { return n; }
				// <- dim(this)
    SharedVector SubVector(int first,int size);
				// returns the subvector of given indices
    SharedVector SubVector(int first,int size) const;
				// returns the subvector of given indices
    //
    // Bounds operators
    // ================
    // These functions find the index and/or value
    // maximizing the given criterion among the values a(i)
    friend int MinAbsIndex(const Vector& a);
    friend double MinAbs(const Vector& a);
    friend void MinAbsIndexAndValue(const Vector& a,int& i,double& v);
				// -|a(i)|
    friend int MaxAbsIndex(const Vector& a);
    friend double MaxAbs(const Vector& a);
    friend void MaxAbsIndexAndValue(const Vector& a,int& i,double& v);
				// +|a(i)|
    friend int MinIndex(const Vector& a);
    friend double Min(const Vector& a);
    friend void MinIndexAndValue(const Vector& a,int& i,double& v);
				// -a(i)
    friend int MaxIndex(const Vector& a);
    friend double Max(const Vector& a);
    friend void MaxIndexAndValue(const Vector& a,int& i,double& v);
				// +a(i)

    //
    // Bounds member functions
    // =======================
    int MinAbsIndex() const;
    double MinAbs() const;
    void MinAbsIndexAndValue(int& i,double& v) const;
    int MaxAbsIndex() const;
    double MaxAbs() const;
    void MaxAbsIndexAndValue(int& i,double& v) const;
    int MinIndex() const;
    double Min() const;
    void MinIndexAndValue(int& i,double& v) const;
    int MaxIndex() const;
    double Max() const;
    void MaxIndexAndValue(int& i,double& v) const;

    //
    // Memory operators
    // ================
    Vector& operator = (const Vector& u);
				// this <- u
    friend void Copy(const Vector& a,Vector& b);
				// b <- a
    friend void Swap(Vector& a,Vector& b);
				// b <-> a

    //
    // Norms, distances and dot product operators
    // ==========================================
    friend double AbsSum(const Vector& a);
				// <- |a(1)|+...+|a(n)|
    friend double Sum(const Vector& a);
				// <- a(1)+...+a(n)
    friend double Norm2(const Vector& a);
				// <- a(1)^2+...+a(n)^2
    friend double Norm(const Vector& a);
				// <- Sqrt(Norm2(a))
    friend double Dist2(const Vector& a,const Vector& b);
				// <- Norm2(b-a)
    friend double Dist(const Vector& a,const Vector& b);
				// <- Norm(b-a)
    friend double Dot(const Vector &u,const Vector &v);
				// <- <u,v>
    friend void Normalize(const Vector &u,Vector &v);
				// v <- (1/Norm(u)) * u

    //
    // Norms member functions
    // ======================
    double AbsSum() const;
				// <- |this(1)|+...+|this(n)|
    double Sum() const;
				// <- this(1)+...+this(n)
    double Norm2() const;
				// <- this(1)^2+...+this(n)^2
    double Norm() const;
				// <- sqrt(this.Norm2())
    void Normalize();
				// this <- (1/this.Norm()) * this
    
    //
    // Vector space operators
    // ======================
    friend void Mul(const double k,const Vector &u,Vector &v);
				// v <- k * u
    friend void Div(const double k,const Vector &u,Vector &v);
				// v <- (1/k) * u
    friend void Add(const Vector &u,const Vector &v,Vector &w);
				// w <- u + v
    friend void Sub(const Vector &u,const Vector &v,Vector &w);
				// w <- u - v
    friend void Add(const Vector &u,double k,const Vector &v,Vector &w);
				// w <- u + k * v
    friend void Combine(double k1,const Vector &u1,
			double k2,const Vector &u2,Vector &v);
				// v <- k1 * u1 + k2 * u2
    friend void Combine(double k1,const Vector &u1,
			double k2,const Vector &u2,
			double k3,const Vector &u3,Vector &v);
				// v <- k1 * u1 + k2 * u2 + k3 * u3
    friend void Middle(const Vector& a,const Vector& b,Vector& m);
				// m <- Middle([ab])
    friend void Barycenter(double k1,const Vector& a1,
			   double k2,const Vector& a2,Vector& a);
				// a <- 1/(k1+k2) * (k1*a1+k2*a2)
    friend void Barycenter(double k1,const Vector& a1,
			   double k2,const Vector& a2,
			   double k3,const Vector& a3,Vector& a);
				// a <- 1/(k1+k2+k3) * (k1*a1+k2*a2+k3*a3)

    //
    // Vector space member functions
    // =============================
    Vector& operator += (const Vector& u);
				// this <- this + u
    Vector& operator -= (const Vector& u);
				// this <- this - u
    Vector& operator *= (double k);
				// this <- k * this
    Vector& operator /= (double k);
				// this <- (1/k) * this
    void Add(const Vector& u);
				// this <- this + u
    void Sub(const Vector& u);
				// this <- this - u
    void Mul(double k);
				// this <- k * this
    void Div(double k);
				// this <- (1/k) * this
    void Add(double k,const Vector& u);
				// this <- this + k * u
    void Combine(double k,double l,const Vector& u);
				// this <- k * this + l * u

    //
    // Vector order (added Aug 95)
    // ============
    friend int
    IsSup(const Vector& p,const Vector& ref);
				// <-
				// +1 if forall i, p(i)>ref(i)
				// -1 if there exists i, p(i)<ref(i)
				//  0 in other cases, i.e.
				//    if forall i, p(i)>=ref(i) and
				//    there exists i, p(i)=ref(i)

    friend int
    IsInf(const Vector& p,const Vector& ref);
				// <-
				// +1 if forall i, p(i)>ref(i)
				// -1 if there exists i, p(i)<ref(i)
				//  0 in other cases, i.e.
				//    if forall i, p(i)>=ref(i) and
				//    there exists i, p(i)=ref(i)
    
    //
    // Dimension specific operators
    // ============================
    friend double Det(const Vector &u,const Vector &v);
				// <- Det(u,v), vectors in R^2
    friend double Det(const Vector &u,const Vector &v,const Vector &w);
				// <- Det(u,v,w), vectors in R^3
    friend void Cross(const Vector &u,const Vector &v,Vector &w);
				// w <- u x v, vectors in R^3

    //
    // IO operators
    // ============
    /*
    friend std::ostream& operator << (std::ostream& o,const Vector& u);
				// ascii <x1,...,xN>
    friend std::ostream& TeX(std::ostream& o,const Vector& u);
				// ascii to be used in latex documents
    */
    };

#endif
