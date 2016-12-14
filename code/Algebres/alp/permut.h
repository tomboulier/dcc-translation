// Linear algebra with C++
// Permutations
//
// $Header: /home/bainvil/Modules/alp/RCS/permut.h,v 1.1 1995/03/16 18:15:21 bainvil Exp bainvil $

#ifndef __PERMUT_H
#define __PERMUT_H

#include <iostream>

class Permutation
    {
  protected:
    int *x;			// Values
    int n;			// Size
  public:
    // Constructor - destructor
    // ========================
    Permutation(int n);
    ~Permutation();

    // Initialization
    // ==============
    void Id();                  // Identity

    // Modification
    // ============
    void Exch(int i,int j);     // Exchange F(i) and F(j)

    // Use
    // ===
    inline int operator () (int i) const
	{ if (i<1 || i>n) return i; else return x[i-1]; }

    // I/O
    // ===
    friend std::ostream& operator << (std::ostream& o,Permutation& p);
    };

#endif
