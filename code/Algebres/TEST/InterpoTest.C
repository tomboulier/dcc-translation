 /* to test the  Interpolation methods
	LD Nov 2012.
 */
 
#include <cstdlib>
#include <cstdio>
#include <iostream> 
#include <fstream>

#include <interpolation.h>
main()
{
  // sampling of the function x+y+z 
  Scalaire vx0y0z0=0, vx1y0z0=1,vx0y1z0=1, vx1y1z0=2,vx0y0z1=1, vx1y0z1=2,vx0y1z1=2, vx1y1z1=3;
  Scalaire   x, y, z  ; // must be between 0 and 1

  x=0.1;y=0.1;z=0.1;
  std::cout<<" x= "<<x<<" ; y= "<<y<<" ; z= "<<z<<std::flush;
  std::cout<<" : TrilinearInterpo = "<< 
    TrilinearInterpo(vx0y0z0,vx1y0z0,vx0y1z0,vx1y1z0,vx0y0z1,vx1y0z1,vx0y1z1,vx1y1z1,x,y,z)<< " et x+y+z= "<<x+y+z<< std::endl;
  //  std::cout<<" : 1D NearInterpo = "<<  NearInterpo(vx0y0z0,vx1y0z0,x)<< std::endl;
  std::cout<<"                          : NearInterpo = "<< 
    NearInterpo(vx0y0z0,vx1y0z0,vx0y1z0,vx1y1z0,vx0y0z1,vx1y0z1,vx0y1z1,vx1y1z1,x,y,z)<< " et rint(x)+rint(y)+rint(z)= "<<rint(x)+rint(y)+rint(z)<<std::endl;
  std::cout<<" ######################################### "<<std::endl;
  x=0.9;y=0.1;z=0.1;
  std::cout<<" x= "<<x<<" ; y= "<<y<<" ; z= "<<z<<std::flush;
  std::cout<<" : TrilinearInterpo = "<< 
    TrilinearInterpo(vx0y0z0,vx1y0z0,vx0y1z0,vx1y1z0,vx0y0z1,vx1y0z1,vx0y1z1,vx1y1z1,x,y,z)<< " et x+y+z= "<<x+y+z<< std::endl;
  std::cout<<"                         : NearInterpo = "<<  
    NearInterpo(vx0y0z0,vx1y0z0,vx0y1z0,vx1y1z0,vx0y0z1,vx1y0z1,vx0y1z1,vx1y1z1,x,y,z)<< " et rint(x)+rint(y)+rint(z)= "<<rint(x)+rint(y)+rint(z)<<std::endl;
}
 
