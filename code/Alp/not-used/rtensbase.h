// Linear algebre with C++
// Real tensors
// $Header: /home/bainvil/Modules/alp/RCS/rtensbase.h,v 1.1 1995/06/10 13:29:42 bainvil Exp bainvil $

#ifndef __RTENSBASE_H
#define __RTENSBASE_H

#include <alp/raggregate.h>  // Real aggregation functions
#include <iostream.h>
#include <math.h>

class SharedTensor;  // Needs to be declared here

// Multi-dimensional array of reals (orders 1 to 6, is 6 enough ?)
class Tensor
    {
  protected:
    int d;	     // Order of the tensor
    int size;	     // Total size = n[0]*n[1]*...*n[d-1]
    double *x;       // Base address
    int *n;	     // Sizes in the d dimensions (n[0]..n[d-1])
    int *s;	     // Increments in the d dimensions (s[0]..s[d-1])
    int shared;	     // Flag : true if shared tensor (x is not allocated by the class)
    //
    // Coordinates <-> array index
    // ===========================
    int CtoI(const int *c) const;
				// c[0..d-1] -> index in x
    void ItoC(int i,int *c) const;
				// index in x -> c[0..d-1]
    //
    // Memory allocation
    // =================
    void AllocD();              // Given d, allocates n and s
    void AllocX();		// Given n[i], allocates x and sets s[i] and size
    //
    // Protected constructor (for shared tensors)
    // ==========================================
    Tensor(int order,char a);   // The char avoids confusion with other constructors
  public:
    //
    // Constructors (each order has its constructor)
    // ============
    Tensor(int n1);             // Vector
    Tensor(int n1,int n2);	// Matrix
    Tensor(int n1,int n2,int n3);
    Tensor(int n1,int n2,int n3,int n4);
    Tensor(int n1,int n2,int n3,int n4,int n5);
    Tensor(int n1,int n2,int n3,int n4,int n5,int n6);
    Tensor(const Tensor& t);	// copy constructor
    //
    // Destructor
    // ==========
    ~Tensor();
    //
    // Memory operations
    // =================
    friend void Copy(const Tensor& a,Tensor& b);
				// b <- a (same order & dimensions)
    friend void Swap(Tensor& a,Tensor& b);
				// a <-> b (same order & dimensions, calls Copy())
    Tensor& operator = (const Tensor& a);
				// this <- a (same order & dimensions, calls Copy())
    //
    // Initialisations
    // ===============
    void Zero();
				// this <- 0
    //
    // Structure access
    // ================
    inline double & Value(int *c) { return x[CtoI(c)]; }
    inline double Value(int *c) const { return x[CtoI(c)]; }
    inline double & operator () (int i1)
	{ return x[(i1-1)*s[0]]; }
    inline double & operator () (int i1,int i2)
	{ return x[(i1-1)*s[0]+(i2-1)*s[1]]; }
    inline double & operator () (int i1,int i2,int i3)
	{ return x[(i1-1)*s[0]+(i2-1)*s[1]+(i3-1)*s[2]]; }
    inline double & operator () (int i1,int i2,int i3,int i4)
	{ return x[(i1-1)*s[0]+(i2-1)*s[1]+(i3-1)*s[2]+(i4-1)*s[3]]; }
    inline double & operator () (int i1,int i2,int i3,int i4,int i5)
	{ return x[(i1-1)*s[0]+(i2-1)*s[1]+(i3-1)*s[2]+(i4-1)*s[3]+(i5-1)*s[4]]; }
    inline double & operator () (int i1,int i2,int i3,int i4,int i5,int i6)
	{ return x[(i1-1)*s[0]+(i2-1)*s[1]+(i3-1)*s[2]+(i4-1)*s[3]+(i5-1)*s[4]+(i6-1)*s[5]]; }
    inline double operator () (int i1) const
	{ return x[(i1-1)*s[0]]; }
    inline double operator () (int i1,int i2) const
	{ return x[(i1-1)*s[0]+(i2-1)*s[1]]; }
    inline double operator () (int i1,int i2,int i3) const
	{ return x[(i1-1)*s[0]+(i2-1)*s[1]+(i3-1)*s[2]]; }
    inline double operator () (int i1,int i2,int i3,int i4) const
	{ return x[(i1-1)*s[0]+(i2-1)*s[1]+(i3-1)*s[2]+(i4-1)*s[3]]; }
    inline double operator () (int i1,int i2,int i3,int i4,int i5) const
	{ return x[(i1-1)*s[0]+(i2-1)*s[1]+(i3-1)*s[2]+(i4-1)*s[3]+(i5-1)*s[4]]; }
    inline double operator () (int i1,int i2,int i3,int i4,int i5,int i6) const
	{ return x[(i1-1)*s[0]+(i2-1)*s[1]+(i3-1)*s[2]+(i4-1)*s[3]+(i5-1)*s[4]+(i6-1)*s[5]]; }
    inline int Order() const { return d; }
    inline int Dimension(int i) const { return n[i-1]; }
    inline int Size() const { return size; }
    //
    // Sharing parts of the tensor
    // ===========================
    SharedTensor Share(int k,int *si,int *c);
				// Shares a k-order sub-tensor of (*this)
				// index j (0<=j<k) of the result is index
				// si[j] of (*this) c gives the other coordinates
				// of the sub-tensor in (*this)
    const SharedTensor Share(int k,int *si,int *c) const;
				// const Tensor version
    SharedTensor Share1(int i);
				// Shares the (d-1)-order of (*this) having the
				// first index i
    const SharedTensor Share1(int i) const;
				// const Tensor version
    //
    // Total aggregation
    // =================
    double Aggregate(AggregationFunction &F) const;
    //
    // Partial aggregation
    // ===================
    void Aggregate(int k,int *si,int *di,AggregationFunction &F,Tensor& a) const;
				// a  is a (d-k) order tensor (the result)
				// si is the array of the k indexes over
				// which the aggregation is performed
				// di gives the destination index di[i] in a
				// of index i in (*this)
    //
    // Real-tensor products
    // ====================
    Tensor& operator *= (double k);                      // this <- k * this
    Tensor& operator /= (double k);			 // this <- 1/k * this
    void Mul(double k);					 // this <- k * this
    void Div(double k);					 // this <- 1/k * this
    friend void Mul(double k,const Tensor& a,Tensor& b); // b <- k * a
    friend void Div(double k,const Tensor& a,Tensor& b); // b <- 1/k * a
    //
    // Tensor-tensor sums
    // ==================
    Tensor& operator += (const Tensor& a);		 // this <- this + a
    Tensor& operator -= (const Tensor& a);		 // this <- this - a
    void Add(const Tensor& a);				 // this <- this + a
    void Sub(const Tensor& a);				 // this <- this - a
    void Add(double k,const Tensor& a);			 // this <- this + k * a
    void Sub(double k,const Tensor& a);			 // this <- this - k * a
    void Combine(double k1,double k2,const Tensor& a);	 // this <- k1 * this + k2 * a
    friend void Add(const Tensor& a,
		    const Tensor& b,Tensor &c);          // c <- a + b
    friend void Sub(const Tensor& a,
		    const Tensor& b,Tensor &c);          // c <- a - b
    friend void Add(const Tensor& a,double k,
		    const Tensor& b,Tensor &c);          // c <- a + k * b
    friend void Sub(const Tensor& a,double k,
		    const Tensor& b,Tensor &c);          // c <- a - k * b
    friend void Combine(double k1,const Tensor& a1,
			double k2,const Tensor& a2,
			Tensor& c);                      // c <- k1 * a1 + k2 * a2
    //
    // Total product
    // =============
    friend double Mul(const Tensor& a,const Tensor& b);
				// <- sum of a(i) * b(i) for all i
				// a,b : same order and dimensions
    //
    // Partial product
    // ===============
    friend void Mul(int k,
		    int *asi,int *adi,const Tensor& a,
		    int *bsi,int *bdi,const Tensor& b,
		    Tensor& c);
				// c is a (a.d-k+b.d-k)-order tensor
				// adi[j] is the index in c of the index j in a
				// bdi[j] is the index in c of the index j in b
				// k is the number of indexes to sum over
				// asi are the k indexes to sum in a
				// bsi are the k corresponding in b
    //
    // Text output
    // ===========
    friend ostream& operator << (ostream& o,const Tensor& a);
    };

#endif
