// Basic definitions
//
// EB Aug 95

#ifndef __basic_h
#define __basic_h

//Begin LD Mar  6 déc 2016 10:12:30 CET
//#include <values.h>  // LD Mar  6 déc 2016 10:12:30 CET
//#include <floats.h>
#include <limits.h> 
//end LD Mar  6 déc 2016 10:12:30 CET

const double AlpEpsilon=1.0e-8;	    // A value for 0
const double AlpInfty=1.0e100;      // A value for +oo
				    // -AlpInfty can always be used as -oo
				    // since MINDOUBLE < -MAXDOUBLE
int
AlpCompare(double x,double y);
				// "epsilon comparison" ; returns :
				// -1 if x < y-AlpEpsilon
				// +1 if x > y+AlpEpsilon
				//  0 if |x-y| <= AlpEpsilon

#endif
