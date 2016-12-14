// Linear algebre with C++
// Real shared tensors
// $Header: /home/bainvil/Modules/alp/RCS/rstensbase.C,v 1.1 1995/06/10 13:29:42 bainvil Exp bainvil $

#include <alp.h>

//
// Constructors
// ============
SharedTensor::SharedTensor(int order,double *data,int *dimensions,int *increments) : Tensor(order,'0')
    {
    int j;
    for (j=0;j<d;j++) { n[j]=dimensions[j];s[j]=increments[j]; }
    x=data;
    }
SharedTensor::SharedTensor(const SharedTensor& t) : Tensor(t.d,'0')
    {
    int j;
    for (j=0;j<d;j++) { n[j]=t.n[j];s[j]=t.s[j]; }
    x=t.x;
    }
//
// operator =
// ==========
SharedTensor& SharedTensor::operator = (const SharedTensor& t)
    {
    int j;
    for (j=0;j<d;j++) { n[j]=t.n[j];s[j]=t.s[j]; }
    x=t.x;
    return *this;
    }
