// Linear algebra with C++
//
// EB 93 94 95

#ifndef __ALP_H
#define __ALP_H

#include <alp/coptions.h>	// Display of compilation options
#include <alp/basic.h>		// Basic definitions

// Auxiliary objects
// =================
#if 0
#include <alp/raggregate.h>	// Real aggregation functions
#endif
#include <alp/permut.h>		// Permutations

// Main objects
// ============
#if 0
#include <alp/rtensbase.h>	// Real tensors
#include <alp/rstensbase.h>	// Real shared tensors
#endif

#include <alp/rvecbase.h>	// Real vectors
#include <alp/rsvecbase.h>      // Real shared vectors
#include <alp/rmatbase.h>	// Real matrices
#include <alp/raffbase.h>	// Real affine applications

// Non-trivial matrix operators
// ============================
#include <alp/rmatops.h>	// Real matrices algorithms

// Geometry
// ========
#include <alp/ibox.h>           // isothetic boxes

#endif
