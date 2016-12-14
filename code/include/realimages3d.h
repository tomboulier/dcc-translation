// manipulation des donnees 3D
// vues comme une succession d'images 2D
//
// LD Nov 94 : creation

// (c) Copyright TIMC 1994

#ifndef __REALIMAGES3D_H
#define __REALIMAGES3D_H


#include <typescalaire.h>
#include <realimage.h>
// #include <geom.h>


class RealImages3D
{
private:
  RealImage **cube;
  int nx,ny,nz;
  Scalaire xbegin,xend,ybegin,yend,zbegin,zend;
  void Allocate();
  // allocate nz RealImage*
  // allocate nz RealImage(s) of dimension nx x ny
  // (*cube)(n) is the first element of the nth RealImage
  // (*cube)(n,i,j) is the value of pixel (i,j) of the nth RealImage
  //                i.e. voxel n,i,j
public:
  //
  // Constructors et destructors
  // =============================
  //
  RealImages3D(int nnx=32,int nny=32,int nnz=32);
  // create a collection of nz real images(nx,ny)
  // by default 32 (32x32)
  RealImages3D(int nnx,int nny,int nnz,
	       Scalaire xb,Scalaire xe,Scalaire yb,Scalaire ye,
	       Scalaire zb,Scalaire ze);
  // create a collection of nz real images(nx,ny)
  // with corners....
  RealImages3D(const RealImages3D & image3d);
  // constructor by copy
  ~RealImages3D();
  // free the nz images and the 3D structure
  
  //
  // Initializations
  // ===============
  //
  void Zero();
  // set the image to 0.
  void AddEllipsoide(Scalaire xaxe,Scalaire yaxe,Scalaire zaxe,
		     Scalaire xcenter,Scalaire ycenter,Scalaire zcenter,
		     Scalaire density,
		     Scalaire eulerx=0, Scalaire eulery=0, Scalaire eulerz=0);
  // add "density" times the indicator of the ellipsis 
  // ((Xr-xcenter)/xaxe)^2+((Yr-ycenter)/yaxe)^2+((Zr-zcenter)/zaxe)^2<1
  // where (Xr,Yr,Zr) Rot[eulerx,eulery,eulerz] (X,Y,Z) 
  // Euler angles are given in degrees

  //  void AddNoise(Scalaire mean,Scalaire sigma);
  // add a normal noise to this.
 
  //
  // Structure Acces 
  // ===============
  //
  inline RealImage & operator () (const int k) 
  // Access to the image number k
  {
#if CHECK
    if ( k<0 || k>=nz)
      fputs("Err : Adress out of RealImage\n",stderr);
#endif
    return (*(cube[k]));
  }
  inline const RealImage  operator () (const int k) const
  // Access to the image number k
  {
#if CHECK
    if ( k<0 || k>=nz)
      fputs("Err : Adress out of RealImage\n",stderr);
#endif
    return (*(cube[k]));
  }

  inline  Scalaire & operator () (const int i,const int j,const int k)
  // Access to the value of pixel i,j
  {
#if CHECK
    if (i<0 || j<0 || k<0 ||i>=nx || j>=ny || k>=nz)
      fputs("Err : Adress out of RealImages3D\n",stderr);
#endif
    return (*(cube[k]))(i,j);
  }

  inline  const Scalaire  operator () (const int i,const int j,const int k) const
  // Acces to pixel i,j.
    {
#if CHECK
      if (i<0 || j<0 || k<0 ||i>=nx || j>=ny || k>=nz)
	fputs("Err : Adress out of RealImages3D\n",stderr);
#endif
      return (*(cube[k]))(i,j);
    }

  int NX() const;
  // Pixel number in X.
  int NY() const;
  // Pixel number in Y.
  int NZ() const;
  // Pixel number in Z.
  void GetCorners(Scalaire *x1,Scalaire *x2,
		  Scalaire *y1,Scalaire *y2,
		  Scalaire *z1,Scalaire *z2) const;
  // Get the 6 corners of the icube (image3D): 
  // blc (x1,y1), brc (x2,y1), tlc (x1,y2) ,trc (x2,y2) ...
  // x1 get the value of xbegin, x2 get the valu of xend, etc....
  void SetCorners(Scalaire x1,Scalaire x2,
		  Scalaire y1,Scalaire y2,
		  Scalaire z1,Scalaire z2);
  // Set the 6 corners of the cube (image3D): 
    Scalaire Interpo0(Scalaire x,Scalaire y, Scalaire z);
  // return the nearest value (images3D(x,y,z))
    void Write(char *filename);
  // Write the image in a file.
    void Examine() const;
  // display the contents of the images
Scalaire Min();
// computes the min
Scalaire Max();
// computes the max
void Plus(const Scalaire val);
// add the constant val to each voxel
void Mul(const Scalaire m);
// multiply the constant m to each voxel
void Log();
// Perform a Log Transform
void JWeight(double *weight);
// multiply each column (J direction) term by term by 
// the weighting vector weight
 void JConvol(double *filter);
// Convolution of each column (J direction) term by term by 
// the filter filter
// filter is vector of length 2 ny - 1 such that
// filter(ny-1) is the filter function at 0
 void JConvol(double *filter,const RealImages3D & image3d);
 // idem as previous but *this <- Jconvol of *filter by *image;
void JConvol(double *filter, const int halfwidth, const RealImages3D & image3d);
// idem as previous but *this <- Jconvol of filter by image;
// but the convolution width is on [-halfwidth halfwidth]
// filter is of size 2*halfwith+1
// filter(halfwidth) is the filter function at 0

/* WARNING praxim code is TRASH now

 void GetSlice(RealImage& slice, CVector3D u, CVector3D v, CPoint3D o);
  // Computes an oblique slice in the real image 3D
  // 

 void GetResample(RealImages3D& resample, 
		  CVector3D u, CVector3D v, CVector3D w, 
		  CPoint3D o);
// Get a resampled cube centeres on the point o 
// and according to the directions u, v and w

*/

//
// IO methods
  void PGMWrite(char *filename);
  // Write the image in pgm format in the file filename.
  void PGMWrite(char *filename, float imagemin, float imagemax);
  // Write the image in pgm format in the file filename.
  // In RAWBITE
  // imagemin and imagemax are used for scaling the data 
  // for producing the image
  void PGMsWrite(char *filename);
  // Write the image in pgm RawBit format in the file filename.
  // computes the maximum and the minimum of the 3D images.
  // these values are used in each 2D slice for normalization 
  // so that the same scale is used for all slices
  void CharWrite(char *filename);
  // Write the image in char format in the file filename.
  void CharWrite(char *filename, int begz, int endz,
	     int sx, int sy, int sz);
  // Write the image in char format in the file filename.
  // keep only the slices k, begz<= k < endz
  // sx,sy,sz are srides
  void MHDWrite(char *filename);
  // Write the image in a bin  file filename : k times (j,i)
  void AVSWrite(char *filename);
  // Write the image in a "easy for AVS" format in the bin file filename 
  // k times (i,j)
  void SPECTWrite(char *filename );
  // Write the image to a file having nuclear medicine (SPECT) format..
  void SPECTRead(char *filename );
  // Read  the image from a file from the gamma camera (SPECT)
  void TRIXELLRead(char *filename );
  // Read  the image from a file from the TRIXELL detector
  void CharRead(char *filename, int trans=0 );
     // Read  the image from a file written by CharWrite.
     // trans = 1-> (*(cube[k]))(j,i)=proji[k*nx*ny+i*ny+j];
     //         else (DEFAULT) -> (*(cube[k]))(i,j)=proji[k*nx*ny+i*ny+j];
  void AVSRead(char *filename );
  // Read  the image from a file written by AVSWrite.
  void INTRead(char *filename );
  // Read  the image from a file written by Joe Aoun (int)
  void ShortRead(char *filename );
     // Read  the image from a file from the file image stored in Short (CT)..
  void UShortWrite(char *filename);
  // Write the image in unsigned short format in the file filename.
  // computes the maximum and the minimum of the 3D images.
  // these values are used in each 2D slice for normalization 
  // so that the same scale is used for all slices
  void PGMRead(char *filename);
  // Read the image in pgm format in the file filename.
  // OK for RAWBIT or ASCII....
};
//
// Operators
// ==========
//
  //
  // Arithmetic RealImages3D operators
  // ===========================
  void Mul(const Scalaire r,const RealImages3D &u,RealImages3D &v);
  // v <- r.u
  void Add(const RealImages3D &u,const RealImages3D &v, 
		  RealImages3D &w);
  // w <- u+v
  void Sub(const RealImages3D &u,const RealImages3D &v,
		  RealImages3D &w);
  // w <- u-v
  void Mul(const RealImages3D &u,const RealImages3D &v,
		  RealImages3D &w);
  // w <- u*v
  void Div(const RealImages3D &u,const RealImages3D &v,
		  RealImages3D &w);
  // w <- u/v
//  void RAdd(const RealImages3D &u,const double r,
//		     const RealImages3D &v,RealImages3D &w);
  // w <- u+r.v
  Scalaire Norm2(const RealImages3D &a);
  // <- a(1)^2+...+a(n)^2
  Scalaire  Dist2(const RealImages3D &a,const RealImages3D &b);
    // <- Norm2(b-a)





    
#endif
