 /* to test the  Interpolation methods
	LD Nov 2012.
 */
 
#include <cstdlib>
#include <cstdio>
#include <iostream> 
#include <fstream>

#include <sparsecolmat.h>
#include <mathtools.h>
main()
{
  // function x^2 on x \in [-3 3]
  double tdeb=-3.0,tfin=3.0;
  int n=50, i;
  double h=(tfin-tdeb)/(n-1);
  double t [n];
  t[0]=tdeb;
  for(i=1;i<n;i++)
    t[i]=t[i-1]+h;
  std::cout<<" Sampling t "<< std::endl;   
  for(i=0;i<n;i++)
    std::cout<<t[i]<< " "<<std::flush;
  std::cout<< std::endl;   
  double f [n];
  for(i=0;i<n;i++)
    f[i]=3*t[i]*t[i]*t[i]+2*t[i]*t[i]+t[i]+1;
    //        f[i]=t[i]*t[i]+t[i]+1;
    //    f[i]=t[i]+1;
  std::cout<<" Sampling f(t)=3*t*t*t + 2*t*t + t + 1 "<< std::endl;   
  //  std::cout<<" Sampling f(t)=t*t + t + 1 "<< std::endl;   
  //  std::cout<<" Sampling f(t)= t + 1 "<< std::endl;   

  for(i=0;i<n;i++)
    std::cout<<f[i]<< " "<<std::flush;
  std::cout<< std::endl;   
  int degree;
  std::cout<<" Entrez le degré de l'espace des Polynome : "<<std::flush;   
  std::cin>>degree;
  Scalaire distance = DistancePoly(t,f,n, degree,1);
  std::cout<<" distance ="<< distance<<std::endl;   

}
 
