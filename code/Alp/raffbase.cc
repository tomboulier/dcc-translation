// Linear algebra with C++
// Affine transformations
//
// $Header: /home/bainvil/Modules/alp/RCS/raffbase.C,v 1.2 1995/02/13 21:40:01 bainvil Exp bainvil $

#include <alp.h>

// Protected functions
// ===================

void Transform::Allocate()
    {
    deporg=new Vector(ndep);
    arrorg=new Vector(narr);
    m=new Matrix(narr,ndep);
    }

// Constructors and destructor
// ===========================
Transform::Transform(int dim) : ndep(dim),narr(dim)
    { Allocate(); }

Transform::Transform(int dep_dim,int arr_dim) : ndep(dep_dim),narr(arr_dim)
    { Allocate(); }

Transform::Transform(const Vector &image_of_zero,const Matrix &mm)
    {
    ndep=mm.C();
    narr=mm.R();
    Allocate();
    Copy(mm,*m);
    deporg->Zero();
    Copy(image_of_zero,*arrorg);
    }

Transform::Transform(const Vector &a,const Vector &image_of_a,const Matrix &mm)
    {
    ndep=mm.C();
    narr=mm.R();
    Allocate();
    Copy(mm,*m);
    Copy(a,*deporg);
    Copy(image_of_a,*arrorg);
    }

Transform::~Transform()
    {
    delete deporg;delete arrorg;delete m;
    }

// Initializations and structure access
// ====================================
void Transform::Identity()
    {
    deporg->Zero();arrorg->Zero();m->Id();
    }

void Transform::Scaling(const Vector &center,double k)
    {
    Copy(center,*deporg);
    Copy(center,*arrorg);
    m->Scaling(k);
    }

#if 0
void Transform::Translation(const Vector &u)
    {
    }

void Transform::Rotation2d(const Vector &center,double angle)
    {
    }

void Transform::Rotation3d(const Vector &axis_point,const Vector &axis_direction,double angle)
    {
    }

// Modification operators
// ======================
void Transform::PreConcatenate(const Transform &t)
    {
    }
void Transform::PostConcatenate(const Transform &t)
    {
    }

// Memory operators
// ================

void Copy(const Transform &f,Transform &g)
    {
    }
void Swap(Transform &f,Transform &g)
    {
    }

// Operators
// =========

void Compose(const Transform &f,const Transform &g,Transform &r)
    {
    }
void Inverse(const Transform &f,Transform &r)
    {
    }
void Image(const Vector &u,const Transform &f,Vector& v)
    {
    }
void LinearImage(const Vector &u,const Transform &f,Vector &v)
    {
    }
#endif
