// Linear algebra with C++
// Affine transformations
//
// $Header: /home/bainvil/Modules/alp/RCS/raffbase.h,v 1.2 1995/02/13 21:40:01 bainvil Exp bainvil $

#ifndef __RAFFBASE_H
#define __RAFFBASE_H

#include <alp/rvecbase.h>
#include <alp/rmatbase.h>

class Transform
    {
  protected:
    // The use of two "origins" makes things easier
    int ndep,narr; // Dimensions of spaces
    Vector *deporg,*arrorg;
    Matrix *m;	   // narr rows and ndep columns
    void Allocate();	   // allocates deporg,arrorg,m knowing ndep and narr
  public:
    // Constructors and destructor
    // ===========================
    Transform(int dim=3);               // Depart dim. = arrival dim.
    Transform(int dep_dim,int arr_dim);	// General case
    Transform(const Vector &image_of_zero,const Matrix &m);
    Transform(const Vector &a,const Vector &image_of_a,const Matrix &m);
    ~Transform();

    // Initializations and structure access
    // ====================================
    void Identity();
    // this <- Id
    void Scaling(const Vector &center,double k);
    // this <- Homothety with given center and coefficient
    void Translation(const Vector &u);
    // this <- Translation of vector u (given in Arr space)
    void Rotation2d(const Vector &center,double angle);
    // this (2d->2d) <- rotation with given characteristics
    void Rotation3d(const Vector &axis_point,const Vector &axis_direction,double angle);
    // this (3d->3d) <- rotation with given characteristics

    // Acces aux membres
    inline Vector &DepOrigin() { return *deporg; }
    inline Vector &ArrOrigin() { return *arrorg; }
    inline Matrix &LinearPart() { return *m; }
    // Version const
    inline const Vector &DepOrigin() const { return *deporg; }
    inline const Vector &ArrOrigin() const { return *arrorg; }
    inline const Matrix &LinearPart() const { return *m; }

    inline int DepDimension() const { return ndep; }
    inline int ArrDimension() const { return narr; }

    // Modification operators
    // ======================
    void PreConcatenate(const Transform &t);
    // this <- this o t, this'(x)=this(t(x))
    void PostConcatenate(const Transform &t);
    // this <- t o this, this'(x)=t(this(x))
    };

// Memory operators
// ================

void Copy(const Transform &f,Transform &g);
// g <- f
void Swap(Transform &f,Transform &g);
// g <-> f

// Operators
// =========

void Compose(const Transform &f,const Transform &g,Transform &r);
// r <- f o g, r(x) = f(g(x))
void Inverse(const Transform &f,Transform &r);
// r <- f^-1 (DepDim == ArrDim required)
void Image(const Vector &u,const Transform &f,Vector& v);
// v <- f(u)
void LinearImage(const Vector &u,const Transform &f,Vector &v);
// v <- linear_part(f)(u)

#endif







