// Linear algebra with C++
// Permutations
//
// $Header: /home/bainvil/Modules/alp/RCS/permut.C,v 1.1 1995/03/16 18:15:21 bainvil Exp bainvil $

#include "permut.h"

//using namespace std;


Permutation::Permutation(int size) : n(size)
    { x=new int[n];Id(); }
Permutation::~Permutation()
    { delete [] x; }
void Permutation::Id()
    { for (int i=1;i<=n;i++) x[i-1]=i; }
void Permutation::Exch(int i,int j)
    { int m=x[i-1];x[i-1]=x[j-1];x[j-1]=m; }
std::ostream& operator << (std::ostream& o,Permutation& p)
    {
    int i;
    o<<'(';
    for (i=1;i<=p.n;i++)
	{
	o<<p.x[i-1];
	o<<((i==p.n)?')':' ');
	}
    return o;
    }
