// Basic definitions
//
// EB Aug 95

#include <alp.h>

int
AlpCompare(double x,double y)
    {
    if (x>y+AlpEpsilon) return 1;
    if (x<y-AlpEpsilon) return -1;
    return 0;
    }
