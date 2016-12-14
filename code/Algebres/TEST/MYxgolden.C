#include <iostream>
#include <iomanip>
#include <cmath>
#include "nr.h"

extern "C" {double MYcfunc(double x);}

using namespace std;

// Driver for routine golden

DP func(const DP x)
{
  //        return NR::bessj0(x);
  //  return (x-3)*(x-3);
  return DP(MYcfunc(x));
}

int main(void)
{
        const DP TOL=1.0e-6,EQL=1.0e-3;
        bool newroot;
        int i,j,nmin=0;
        DP ax,bx,cx,fa,fb,fc,xmin;

        cout << "Minima of the function x^2" << endl;
        cout << setw(10) << "min. #" << setw(9) << "x";
        cout << setw(18) << "(x-3)*(x-3)" <<  endl;
        cout << fixed << setprecision(6);
	ax=DP(-5);
	bx=DP(5);
	NR::mnbrak(ax,bx,cx,fa,fb,fc,func);
	NR::golden(ax,bx,cx,func,TOL,xmin);
	cout << setw(7) << nmin << setw(16) << xmin;
	cout << setw(13) << func(xmin);
	cout << endl;

        return 0;
}
