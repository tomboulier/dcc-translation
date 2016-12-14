// Linear algebra with C++
// Real shared vectors
//
// EB Juin 95

#include <alp.h>

SharedVector::SharedVector() : Vector(0)
   { x=0;n=s=0; }
SharedVector::SharedVector(double *ptr,int size,int increment) : Vector(0)
   { x=ptr;n=size;s=increment; }
SharedVector::SharedVector(const SharedVector &u) : Vector(0)
   { x=u.x;n=u.n;s=u.s; }
