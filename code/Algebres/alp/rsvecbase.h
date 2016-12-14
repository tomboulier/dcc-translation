// Linear algebra with C++
// Real shared vectors
//
// $Header: /home/bainvil/Modules/alp/RCS/rsvecbase.h,v 1.1 1995/06/10 13:29:42 bainvil Exp bainvil $

#ifndef __RSVECBASE_H
#define __RSVECBASE_H

#include <alp/rvecbase.h>

// A SharedVector doesn't allocate its array of values.
// The object is thus "attached" to another object (matrix, vector, ...).

class SharedVector : public Vector
    {
  public:
    //
    // Constructors
    // ============
    SharedVector();
    SharedVector(double *ptr,int size,int increment);
    SharedVector(const SharedVector &u);

    //
    // Assignment
    // ==========
    inline SharedVector& operator = (const SharedVector& u)
	{ x=u.x;n=u.n;s=u.s;return *this; }
    //
    // Attachment to an array (or a portion of an array)
    // ======================
    inline void Attach(double *ptr,int size,int increment)
	{ x=ptr;n=size;s=increment; }
    };

// Redefinition of some functions to avoid a warning such as
//
// rmatbase.C: In method `void Matrix::SetColumn(int, const class Vector &)':
// rmatbase.C:73: warning: initialization of non-const `Vector &'
//                from rvalue `SharedVector'
//
// when using functions such as
//
// void Copy(const Vector& a,Vector& b)
//
// where b is a SharedVector

inline void Copy(const Vector& a,SharedVector b)
    { Copy(a,(Vector&)b); }
inline void Swap(SharedVector a,SharedVector b)
    { Swap((Vector&)a,(Vector&)b); }
inline void Normalize(const Vector &u,SharedVector v)
    { Normalize(u,(Vector&)v); }
inline void Mul(const double k,const Vector &u,SharedVector v)
    { Mul(k,u,(Vector&)v); }
inline void Div(const double k,const Vector &u,SharedVector v)
    { Div(k,u,(Vector&)v); }
inline void Add(const Vector &u,const Vector &v,SharedVector w)
    { Add(u,v,(Vector&)w); }
inline void Sub(const Vector &u,const Vector &v,SharedVector w)
    { Sub(u,v,(Vector&)w); }
inline void Add(const Vector &u,double k,const Vector &v,SharedVector w)
    { Add(u,k,v,(Vector&)w); }
inline void Combine(double k1,const Vector &u1,
		    double k2,const Vector &u2,SharedVector v)
    { Combine(k1,u1,k2,u2,(Vector&)v); }
inline void Combine(double k1,const Vector &u1,
		    double k2,const Vector &u2,
		    double k3,const Vector &u3,SharedVector v)
    { Combine(k1,u1,k2,u2,k3,u3,(Vector&)v); }
inline void Middle(const Vector& a,const Vector& b,SharedVector m)
    { Middle(a,b,(Vector&)m); }
inline void Barycenter(const double k1,const Vector& a1,
		       const double k2,const Vector& a2,SharedVector a)
    { Barycenter(k1,a1,k2,a2,(Vector&)a); }
inline void Barycenter(const double k1,const Vector& a1,
		       const double k2,const Vector& a2,
		       const double k3,const Vector& a3,SharedVector a)
    { Barycenter(k1,a1,k2,a2,k3,a3,(Vector&)a); }

#endif
