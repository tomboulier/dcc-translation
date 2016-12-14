// Operation sur les images reelles
// en particulier pour la reconstruction en tomographie
// 
// LD Dec 93

// (c) Copyright TIMC 1993

#ifndef __REALIMAGE_H
#define __REALIMAGE_H


#include <math.h>
#include <tomotypes.h>

class RealImage
{
 private:
  Scalaire *ima;	// Coefficients or values
  Scalaire **line; // pointer inside ima to accelerate
  int nx,ny;	// dimensions of an 2D image
  Scalaire xbegin,xend,ybegin,yend;
  // Bottom Left Corner=(xbegin,ybegin).
  void Allocate();	// Allocate ima(nx,ny)
// We suppose that images are regularly sampled in both x and y directions.
// Thus, only the value of nx,ny, xbegin,xend,ybegin,yend are necessary 
// to know where the pixel ima[i,j] is.
//
 public:
//
// Constructeurs et destructeurs
// =============================
//
		
  RealImage(int nnx=32,int nny=32);
// create an real image(nx,ny)
// by default 32x32
// by default [-1, 1]x[-1, 1];
//
   RealImage(int nnx,int nny,
		      Scalaire xb,Scalaire xe,Scalaire yb,Scalaire ye);
// create an real image(nx,ny)
// with corners (xb,yb) (xb,ye) (xe,yb) (xe,ye)

  RealImage(const RealImage & image);
// constructor by copy
  ~RealImage();
// free the image
//
// Initialisations
// ===============
//
  void Zero();
  // set the image to 0.
  void AddGaussienne(Scalaire sigma=1,Scalaire xmean=0,Scalaire ymean=0);
  // Add a Gaussian to an image
  // Gaussian is e^{((x-xmean)**2+(y-ymean)**2)/sigma**2}
  void AddEllipse(Scalaire xaxe,Scalaire yaxe,
		  Scalaire xcenter,Scalaire ycenter,
		  Scalaire density,Scalaire angle=0);
  //add an ellipse of grey level "density"
  // centered on xcenter, ycenter 
  // of (long, at least should be) semi axis xaxe 
  //  with angle direction "angle"
  // and orthogonal short semi axis yaxe
  void AddPhantomDef(char *filename);
  //add a Phantom defined in the definition file filaname
  // (ellipse, disk, traingle (not yet implemented)


  void TVdenoising(RealImage &denoisedima, RealImage &Zx,RealImage &Zy,
		   int nbiter, Scalaire alpha, Scalaire gamma=0.2);
  // denoisedima is the denoised image result
  // Z is a variable homogeneous to grad(*this) must be zero 
  // gamma = 0.8*1/2^n : for two dim 0.8*.25=0.2 (gamma is a progression step)
  // alpha is an hyper parameter (we solve argmin(alpha TV(u)+1/2|||x-u|^2)
  // nbiter is the iteration number 

  //  void AddNoise(Scalaire mean,Scalaire sigma);
  // add a normal noise to this.
  void UnderSampling(RealImage *image,
		     int undrow,int undcol,int iofset, int jofset);
  // set the RealImage to the UnderSampling of *image
  // (*this)(i,j) <- Local Mean Close to (*image)(iofset+undrow*i,jofset=undcol*j)
  void EvenInterlacedUnderSampling(RealImage *image,
					      int undrow,int undcol,
					      int iofset, int jofset);
  // set the RealImage to the Interlaced UnderSampling of *image
  // (*this)(i,j) <- Local Mean Close to 
  // (*image)(iofset+undrow*i,jofset + i%2 + undcol*2*j)

void OddInterlacedUnderSampling(RealImage *image,
					   int undrow,int undcol,
				int iofset, int jofset);
  // set the RealImage to the Interlaced UnderSampling of *image
  // (*this)(i,j) <- Local Mean Close to 
  // (*image)(iofset+undrow*i,jofset + (i+1)%2 + undcol*2*j)


//
// Acces structure
// ===============
//
  inline  Scalaire & operator () (const int i,const int j)
    // Access to the value of pixel i,j
      {
#if CHECK
	if (i<0 || j<0 || i>=nx || j>=ny)
	  fputs("Err : Adress out of RealImage\n",stderr);
#endif
	return line[i][j];
      }
  inline const Scalaire  operator () (const int i,const int j) const
    // Acces to pixel i,j.
    // Acces to *a(i,j) : i must be >=1 and <=dim1,
    //                    j must be >=1 and <=n,
      {
	return line[i][j];
      }
  void GetCorners(Scalaire *x1,Scalaire *x2,
		  Scalaire *y1,Scalaire *y2) const;
// Get the 4 corners of the image: 
// blc (x1,y1), brc (x2,y1), tlc (x1,y2) ,trc (x2,y2)
  void SetCorners(Scalaire x1,Scalaire x2,
		  Scalaire y1,Scalaire y2);
// Set the 4 corners of the image: 
  int NX() const;
// Pixel number in X.
  int NY() const;
// Pixel number in Y.
     Scalaire Min();
// Minimum of an image
     Scalaire Max();
// Maximum of an image

 void Plus(const Scalaire val);
// add the constant val to each pixel
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
 void JConvol(double *filter, const RealImage & image);
 // idem as previous but *this <- Jconvol of filter by image;
 void JConvol(double *filter, const int halfwidth, const RealImage & image);
// idem as previous but *this <- Jconvol of filter by image;
// but the convolution width is on [-halfwidth halfwidth]
// filter is of size 2*halfwith+1
// filter(halfwidth) is the filter function at 0

  void Examine() const;
// display the contents of the image.

  void Write(char *filename);
// Write the image in the file filename (formated)
//  "fprintf(pfile,"%e",ima[i*ny+j]);"

  void MHDWrite(char *filename);
// Write the image in a bin file. (j,i)
  void AVSWrite(char *filename);
// Write the image in a bin file.(i,j)

// Read the image in an AVS file.
  void AVSRead(char *filename);

  void INTRead(char *filename );
  // Read  the image from a file written by Joe Aoun

// Read  the image from a file from the gamma camera (SPECT)..
     void SPECTRead(char *filename);

// Read  the image from a file from the TRIXELL detector..
     void TRIXELLRead(char *filename, int headsize=0);

// Write the image in a file.
     void PGMaWrite(char *filename);
// Write the image in pgm format in the file filename.
// format P2 : ASCII 
     void PGMWrite(char *filename, float imagemin, float imagemax);
// Write the image in pgm format in the file filename.
// In RAWBITE
// imagemin and imagemax are used for scaling the data 
// for producing the image
  void PGMWrite(char *filename);
// Write the image in pgm RawBit format in the file filename.

  void RAWrbWrite(char *filename);
// Write the image in rawbit (pgm without header) format in the file filename.
  void UShortWrite(char *filename, float imagemin=0, float imagemax=0);
// Write the image in Unsigned Short format in the file filename.
// imagemin and imagemax are used for scaling the data 
// for producing the image
  // if imagemin==imagemax==0 then the max and min are computed

  void PGMRead(char *filename);
// Read the image in pgm format in the file filename.
// OK for RAWBIT or ASCII....

  void VMIWrite(char *filename);
// Write the image in vmi format in the file filename.

};

//
// Operateurs
// ==========
//

  void Mul(const Scalaire r,const RealImage &u,RealImage &v);
  // v <- k.u
  void Add(const RealImage &u,const RealImage &v, 
		  RealImage &w);
  // w <- u+v
  void Sub(const RealImage &u,const RealImage &v,
		  RealImage &w);
// w <- u-v
   void Mul(const RealImage &u,const RealImage &v,
		  RealImage &w);
// w <- u*v
   void Div(const RealImage &u,const RealImage &v,
		  RealImage &w);
// w <- u/v
//  void RAdd(const RealImage &u,const double r,
//		     const RealImage &v,RealImage &w);
  // w <- u+r.v
  Scalaire Norm2(const RealImage &a);
  // <- a(1)^2+...+a(n)^2
  Scalaire  Dist2(const RealImage &a,const RealImage &b);
    // <- Norm2(b-a)

// to be written....     extern "C" void Imio_WriteImage(char *chemin,char *data,short lig,short col,int size_erg);
/* 
   just for calling the vmi write image proc.
   the file opetrator.h in $(TIMC)/include does not work...????
*/


int LineInterImage(const RealImage &image,
		   const Scalaire cosphi,const Scalaire sinphi,const Scalaire s,
		   Scalaire *intersection, int *iline, int *jcol );
// Compute the intersection lengths between the image "image" and the line  
// $\left\{x\in\bR^2, (cosphi,sinphi)\cdot x= s\right\}
// i.e the Radon Line (\phi,s) [ rq that this line is equal to (\phi+\pi, -s) ]
//     
// output
// LineInterImage is the number of intersected pixels (if 0, no intersection)
// else for l=0,l<LineInterImage intersection[l] is the length of the intersection of the line with the pixel iline[l] iline[l]
// WARNING : the dimensions of the 3 vectors intesection, iline, jcol should be
// nx+ny-1 where (nx,ny) are the dimensions of the image
// WARNING : cosphi must be >0 and sinphi!=0
//     => if cosphi<0 consider the line (\phi+\pi,-s)) 
//         if cosphi < 0 then compute  -(cosphi,sinphi)\cdot x= - s
//           i.e. replace \phi by \phi+\pi
//     => if sinphi=0 (ortho orizontal lines) or cosphi=0 (orthp vertical lines)
//           cf infra
//////////

int HorizonLineInterImage(const RealImage &image,
			  const Scalaire s,
			  Scalaire *intersection, int *iline, int *jcol);
// case where cosphi=0 \vzeta is (1,0) or (-1,0) sinphi= +or- 1
// WARNING we assume \phi=-\pi/2 thus -sinphi=1
int VertLineInterImage(const RealImage &image,
		       const Scalaire s,
		       Scalaire *intersection, int *iline, int *jcol);
// case where sinphi=0 \vzeta is (0,1) or (0,-1) cosphi= +or- 1
// WARNING we assume \phi=0 thus cosphi=1

#endif
