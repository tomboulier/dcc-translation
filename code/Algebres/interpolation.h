#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <typescalaire.h>

// ========================================
//  interpolation
// ========================================

inline Scalaire LinearInterpo(Scalaire vx0,Scalaire vx1,
				Scalaire x)
// linear interpolation on [0,1] at x from 
// v(x0),v(x1)
// x must be in [0,1]
{  return( (x*vx1+(1-x)*vx0));}

inline Scalaire BilinearInterpo(Scalaire vx0y0,Scalaire vx1y0,
				Scalaire vx0y1,Scalaire vx1y1,
				Scalaire x,Scalaire y)
// bilinear interpolation on [0,1]x[0,1] at (x,y) from 
// v(x0,y0),v(x1,y0), v(x0,y1),v(x1,y1)
// (x,y) must be in [0,1]x[0,1]
{  return( (1-y)*(x*vx1y0+(1-x)*vx0y0) + y*(x*vx1y1+(1-x)*vx0y1) );}

inline Scalaire TrilinearInterpo(Scalaire vx0y0z0,Scalaire vx1y0z0,
				 Scalaire vx0y1z0,Scalaire vx1y1z0,
				 Scalaire vx0y0z1,Scalaire vx1y0z1,
				 Scalaire vx0y1z1,Scalaire vx1y1z1,
				 Scalaire x,Scalaire y,Scalaire z)
// trilinear interpolation on [0,1]x[0,1]x[0,1] at (x,y,z) from 
// v(x0,y0,z0),v(x1,y0,z0), v(x0,y1,z0),v(x1,y1,z0),v(x0,y0,z1),v(x1,y0,z1), v(x0,y1,z1),v(x1,y1,z1)
// (x,y,z) must be in [0,1]x[0,1]x[0,1]
{  return( (1-z)*
	   ( (1-y)*(x*vx1y0z0+(1-x)*vx0y0z0) + y*(x*vx1y1z0+(1-x)*vx0y1z0) )
	   +z * ( (1-y)*(x*vx1y0z1+(1-x)*vx0y0z1) + y*(x*vx1y1z1+(1-x)*vx0y1z1) )
	   );
}


inline Scalaire NearInterpo(Scalaire vx0,Scalaire vx1,
			    Scalaire x)
// nearest interpolation on [0,1] at x from 
// v(x0),v(x1)
// x must be in [0,1]
{ Scalaire ix=rint(x);return( ((1-ix)*vx0+ix*vx1));}

inline Scalaire NearInterpo(Scalaire vx0y0,Scalaire vx1y0,
				Scalaire vx0y1,Scalaire vx1y1,
				Scalaire x,Scalaire y)
// nearest interpolation on [0,1]x[0,1] at (x,y) from 
// v(x0,y0),v(x1,y0), v(x0,y1),v(x1,y1)
// (x,y) must be in [0,1]x[0,1]
{   Scalaire ix=rint(x), iy=rint(y);
  return( (1-iy)*((1-ix)*vx0y0+ix*vx1y0) + iy*((1-ix)*vx0y1+ix*vx1y1) );
}


inline Scalaire NearInterpo(Scalaire vx0y0z0,Scalaire vx1y0z0,
				 Scalaire vx0y1z0,Scalaire vx1y1z0,
				 Scalaire vx0y0z1,Scalaire vx1y0z1,
				 Scalaire vx0y1z1,Scalaire vx1y1z1,
				 Scalaire x,Scalaire y,Scalaire z)
// Nearest interpolation on [0,1]x[0,1]x[0,1] at (x,y,z) from 
// v(x0,y0,z0),v(x1,y0,z0), v(x0,y1,z0),v(x1,y1,z0),v(x0,y0,z1),v(x1,y0,z1), v(x0,y1,z1),v(x1,y1,z1)
// (x,y,z) must be in [0,1]x[0,1]x[0,1]
{   Scalaire ix=rint(x), iy=rint(y),iz=rint(z);
  return((1-iz)*
	 ((1-iy)*((1-ix)*vx0y0z0+ix*vx1y0z0) + iy*((1-ix)*vx0y1z0+ix*vx1y1z0))
	 +iz*
	 ((1-iy)*((1-ix)*vx0y0z1+ix*vx1y0z1) + iy*((1-ix)*vx0y1z1+ix*vx1y1z1))
	 );
}


