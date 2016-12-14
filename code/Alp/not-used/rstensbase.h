// Linear algebre with C++
// Real shared tensors
// $Header: /home/bainvil/Modules/alp/RCS/rstensbase.h,v 1.1 1995/06/10 13:29:42 bainvil Exp bainvil $

#ifndef __RSTENSBASE_H
#define __RSTENSBASE_H

#include <alp/raggregate.h>  // Real aggregation functions
#include <alp/rtensbase.h>   // Real tensors

class SharedTensor : public Tensor
    {
  public:
    //
    // Constructors
    // ============
    SharedTensor(int order,double *data,int *dimensions,int *increments);
    SharedTensor(const SharedTensor& t);
    //
    // operator =
    // ==========
    SharedTensor& operator = (const SharedTensor& t);
    };

#endif
