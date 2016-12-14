// manipulation des donnees 3D
// vues comme une succession d'images 2D
//
// LD Nov 94 : creation
// (c) Copyright TIMC 1994
//
// A FAIRE : la gestion des coins!!!
//

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <realimages3d.h>
// #include <fstream.h>

using namespace std;

void RealImages3D::Allocate()
     // allocate nz RealImage*
     // allocate nz RealImage(s) of dimension nx x ny
     // (*cube)(n) is the first element of the nth RealImage
     // (*cube(n,i,j) is the value of pixel (i,j) of the nth RealImage
     //                i.e. voxel n,i,j
{
  cube= new  RealImage* [nz];
}

RealImages3D::RealImages3D(int nnx,int nny,int nnz)
     // create a collection of nz real images(nx,ny)
     // by default 32 (32x32)
{
  int k;
  nx=nnx;ny=nny;nz=nnz;
  xbegin=-1;
  xend=1;
  ybegin=-1;
  yend=1;
  zbegin=-1;
  zend=1;
  Allocate();
  for(k=0;k<nz;k++)
      cube[k]= new RealImage(nx,ny,xbegin,xend,ybegin,yend);
//      RealImage (cube[k]) (nx,ny,xbegin,xend,ybegin,yend);
}
RealImages3D::RealImages3D(int nnx,int nny,int nnz,
	       Scalaire xb,Scalaire xe,Scalaire yb,Scalaire ye,
	       Scalaire zb,Scalaire ze)
{
  int k;
  nx=nnx;ny=nny;nz=nnz;
  xbegin=xb;
  xend=xe;
  ybegin=yb;
  yend=ye;
  zbegin=zb;
  zend=ze;
  Allocate();
  for(k=0;k<nz;k++)
      cube[k]= new RealImage(nx,ny,xbegin,xend,ybegin,yend);
//      RealImage (cube[k]) (nx,ny,xbegin,xend,ybegin,yend);
}

RealImages3D::RealImages3D(const RealImages3D & image3d)
     // constructor by copy
{
  int i,j,k;
  Scalaire xb,xe,yb,ye,zb,ze;
  nx=image3d.NX();
  ny=image3d.NY();
  nz=image3d.NZ();
  image3d.GetCorners(&xb,&xe,&yb,&ye,&zb,&ze);
  xbegin=xb;
  xend=xe;
  ybegin=yb;
  yend=ye;
  zbegin=zb;
  zend=ze;
  
  Allocate();
  
  for(k=0; k<nz; k++)
       cube[k] = new RealImage(*image3d.cube[k]);
//      RealImage (cube[k]) (nx,ny,xbegin,xend,ybegin,yend);
}

RealImages3D::~RealImages3D()
     // free the images and the 3D structure
{
  int k;
  // free the RealImage(s)
  for(k=0;k<nz;k++)
      delete  cube[k];
  delete  cube;
}

//
// Initializations
// ===============
//
void RealImages3D::Zero()
     // set all the images to 0.
{
 int k;
 for (k=0;k<nz;k++)
   {
     cube[k]->Zero();
     // i.e.		ima[i*ny+j]=0.0;
//     cout << "dans RealImages3D::Zero examen du bloc " << k << endl;
//    cube[k].Examine();
    }
}

void RealImages3D::AddEllipsoide(Scalaire xaxe,Scalaire yaxe,Scalaire zaxe,
				 Scalaire xcenter,Scalaire ycenter,Scalaire zcenter,
				 Scalaire density,
				 Scalaire eulerx,Scalaire eulery,Scalaire eulerz)
  // add "density" times the indicator of the ellipsis 
  // ((Xr-xcenter)/xaxe)^2+((Yr-ycenter)/yaxe)^2+((Zr-zcenter)/zaxe)^2<1
  // where (Xr,Yr,Zr) Rot[eulerx,eulery,eulerz] (X,Y,Z) 
  // Euler angles are given in degrees

{
 int i,j,k;
 Scalaire xcrs,ycrs,zcrs,hx,hy,hz;

 if((eulerx!=0)||(eulerx!=0)||(eulerx!=0)){
   std::cout<< " Euler angle must be zero NOW " << std::endl;
   exit(0);
 }
 /*
 CPoint3D X(0,0,0);
 CPoint3D Xr(0,0,0);
 CVector3D translation(0 ,0 , 0);
 // CVector3D translationrot(0 ,0 ,0 );
 CRotation3D rotation(eulerx*M_PI/180, eulery*M_PI/180, eulerz*M_PI/180);
 CTransform3D transfo(rotation , translation);
 */
 Scalaire xr,yr,zr,x,y,z;

 hx=(xend-xbegin)/nx;
 hy=(yend-ybegin)/ny;
 hz=(zend-zbegin)/nz;
 for (k=0;k<nz;k++){
   //   X.Z()=zbegin+hz*(k+.5)-zcenter;
   z=zbegin+hz*(k+.5)-zcenter;
   for (i=0;i<nx;i++){
     //     X.X()=xbegin+hx*(i+.5)-xcenter;
     x=xbegin+hx*(i+.5)-xcenter;
     for (j=0; j<ny;j++){
       // X.Y()=ybegin+hy*(j+.5)-ycenter;
       y=ybegin+hy*(j+.5)-ycenter;
       //Xr=transfo*X;
       // on aimerait une rotation mais ICI identité
       xr=x;yr=y;zr=z;
       zcrs=(zr/zaxe);
       zcrs=zcrs*zcrs;
       xcrs=(xr/xaxe);
       xcrs=xcrs*xcrs;
       ycrs=(yr/yaxe);
       ycrs=ycrs*ycrs;
       if((xcrs+ycrs+zcrs)<1)
	 (*(cube[k]))(i,j)=(*(cube[k]))(i,j)+density;
     }
   }
 }
}


//
// Structure Acces 
// ===============
//
//Scalaire & RealImages3D::operator () (const int i,const int j,const int k)
//     // Access to the value of pixel i,j
//const Scalaire  RealImages3D::operator () (const int i,const int j, const int k) const
// Acces to pixel i,j.

int RealImages3D::NX() const
// Pixel number in X.
{
  return(nx);
}
int RealImages3D::NY() const
// Pixel number in Y.
{
  return(ny);
}
int RealImages3D::NZ() const
// Pixel number in Z.
{
  return(nz);
}


void RealImages3D::GetCorners(Scalaire *x1,Scalaire *x2,
			      Scalaire *y1,Scalaire *y2,
			      Scalaire *z1,Scalaire *z2)  const
			      // Get the 8 corners of the icube (image3D): 
			      // blc (x1,y1,z1), brc (x2,y1,z1),...
{
  *x1=xbegin;
  *x2=xend;
  *y1=ybegin;
  *y2=yend;
  *z1=zbegin;
  *z2=zend;
}

void RealImages3D::SetCorners(Scalaire x1,Scalaire x2,
		  Scalaire y1,Scalaire y2,
		  Scalaire z1,Scalaire z2)
// Set the 8 corners of the cube (image3D):
{
  xbegin=x1;
  xend=x2;
  ybegin=y1;
  yend=y2;
  zbegin=z1;
  zend=z2;
} 


Scalaire RealImages3D::Interpo0(Scalaire x,Scalaire y, Scalaire z)
// return the nearest value (images3D(x,y,z))
{
  int i,j,k;
  i=(int)((x-xbegin)/(xend-xbegin)*nx);
  j=(int)((y-ybegin)/(yend-ybegin)*ny);
  k=(int)((z-zbegin)/(zend-zbegin)*nz);
  if (i<0 || j<0 || k<0 ||i>=nx || j>=ny || k>=nz)
    return (0);
  else
    return (*(cube[k]))(i,j);
} 

void RealImages3D::Write(char *filename)
     // Write the image in a file.
{
  FILE *pfile;
  int i,j,k;
  float val;

/* open the file for writing in */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }

/* write the number of lines and columns */
  fprintf(pfile,"%d %d\n",nx,ny);
/* lecture des lignes suivantes */
  for (k=0;k<nz;k++)
    for (i=0;i<nx;i++)
      for (j=0; j<ny;j++)
	{
	  fprintf(pfile,"%e",(*(cube[k]))(i,j));
	}
  fclose(pfile);
}


void RealImages3D::MHDWrite(char *filename)
// Write the image in a bin  file filename : k times (j,i)
{
  FILE *pfile;
  int i,j,k;
  float val,min,max;
  float *uc;

  uc = new float[nx*ny];

/* open the file for writing in */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
  min=max=(*(cube[0]))(0,0);
  for (k=0;k<nz;k++)
    for (i=0;i<nx;i++)
      for (j=0; j<ny;j++)
	{
	  val=(*(cube[k]))(i,j);
	  if(min>val)min=val;
	  if(max<val)max=val;
	}
  cout << " file " << filename << " min, max = " << min <<" ; "<< max << endl;

  for (k=0;k<nz;k++){
    for (j=0; j<ny;j++)
      for (i=0;i<nx;i++)
	uc[j*nx+i]=(*(cube[k]))(i,j);
    fwrite(uc,sizeof(float),nx*ny,pfile);
  }
  fclose(pfile);
  delete [] uc;
}

void RealImages3D::AVSWrite(char *filename )
  // Write the image in a "easy for AVS" format in the bin file filename 
  // k times (i,j)
{
  FILE *pfile;
  int i,j,k;
  float val,min,max;
  float *uc;

  uc = new float[nx*ny];

/* open the file for writing in */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
  min=max=(*(cube[0]))(0,0);
  for (k=0;k<nz;k++)
    for (i=0;i<nx;i++)
      for (j=0; j<ny;j++)
	{
	  val=(*(cube[k]))(i,j);
	  if(min>val)min=val;
	  if(max<val)max=val;
	}
  cout << " file " << filename << " min, max = " << min <<" ; "<< max << endl;

  for (k=0;k<nz;k++)
    {
      for (i=0;i<nx;i++)
	for (j=0; j<ny;j++)
	  uc[i*ny+j]=(*(cube[k]))(i,j);
      fwrite(uc,sizeof(float),nx*ny,pfile);
    }
  fclose(pfile);
  delete [] uc;
}

void RealImages3D::SPECTWrite(char *filename )
     // Write the image to a file having nuclear medicine (SPECT) format..
{
  FILE *pfile; 
  int i,j,k;
  float val,min,max;
  unsigned short *proji;

  proji= new unsigned short [nx*ny*nz];

// open the file for writing in 
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
  min=max=(*(cube[0]))(0,0);
  for (k=0;k<nz;k++)
    for (i=0;i<nx;i++)
      for (j=0; j<ny;j++)
	{
	  val=(*(cube[k]))(i,j);
	  if(min>val)min=val;
	  if(max<val)max=val;
	}
  cout << " file " << filename << " min, max = " << min <<" ; "<< max << endl;

/* write the image in binary  */

  for (k=0;k<nz;k++)
    {
      for (j=0; j<ny;j++)
	for (i=0;i<nx;i++)
	  proji[k*nx*ny+i*ny+j]=(unsigned short) (*(cube[k]))(j,i);
    }
  fwrite (proji, sizeof(short), nz*ny*nx, pfile);

  fclose(pfile);
  delete [] proji;
}

void RealImages3D::SPECTRead(char *filename )
     // Read  the image from a file from the gamma camera (SPECT)..
{
  FILE *pfile; 
  int i,j,k;
  float val,min,max;
  unsigned short *proji;

  proji= new unsigned short [nx*ny*nz];

// open the file for reading in 
  if ((pfile=fopen(filename,"r"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
// lecture des lignes suivantes 
  fread (proji, sizeof(short), nz*ny*nx, pfile);

  for (k=0;k<nz;k++)
    {
      for (j=0; j<ny;j++)
	for (i=0;i<nx;i++)
	  (*(cube[k]))(j,i)=proji[k*nx*ny+i*ny+j];
    }
  fclose(pfile);
  delete [] proji;
}

void RealImages3D::ShortRead(char *filename )
     // Read  the image from a file from the file image stored in Short (CT)..
{
  FILE *pfile; 
  int i,j,k;
  float val,min,max;
  short *proji;

  proji= new short [nx*ny*nz];

// open the file for reading in 
  if ((pfile=fopen(filename,"r"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
// lecture des lignes suivantes 
  fread (proji, sizeof(short), nz*ny*nx, pfile);

  for (k=0;k<nz;k++)
    {
      for (j=0; j<ny;j++)
	for (i=0;i<nx;i++)
	  (*(cube[k]))(j,i)=proji[k*nx*ny+i*ny+j];
    }
  fclose(pfile);
  delete [] proji;
}

void RealImages3D::TRIXELLRead(char *filename )
     // Read  the image from a file written by imageJ.
{
  FILE *pfile;
  int i,j,k;
  //  float val,min,max;
  unsigned short *proji;

  proji = new unsigned short [nx*ny];

/* open the file for reading in */
  if ((pfile=fopen(filename,"r"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
/* lecture des lignes suivantes */
  for (k=0;k<nz;k++)
    {
      fread (proji, sizeof(short), ny*nx, pfile);
      for (i=0;i<nx;i++)
	for (j=0; j<ny;j++)
	  (*(cube[k]))(i,j)=proji[i*ny+j];
    }
  fclose(pfile);
  delete [] proji;
}

void RealImages3D::CharRead(char *filename, int trans)
     // Read  the image from a file written by CharWrite.
     // trans = 1-> (*(cube[k]))(j,i)=proji[k*nx*ny+i*ny+j];
     //         else (DEFAULT) -> (*(cube[k]))(i,j)=proji[k*nx*ny+i*ny+j];
{
  FILE *pfile; 
  int i,j,k;
  float val,min,max;
  unsigned char *proji;

  proji= new unsigned char [nx*ny*nz];

/* open the file for reading in */
  if ((pfile=fopen(filename,"r"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
/* lecture des lignes suivantes */
  fread (proji,  sizeof(char), nz*ny*nx, pfile);

  if(trans!=1)
    for (k=0;k<nz;k++)
      for (i=0;i<nx;i++)
	for (j=0; j<ny;j++)
	  (*(cube[k]))(i,j)=proji[k*nx*ny+i*ny+j];
  else
    for (k=0;k<nz;k++)
      for (i=0;i<nx;i++)
	for (j=0; j<ny;j++)
	  (*(cube[k]))(j,i)=proji[k*nx*ny+i*ny+j];
    
  fclose(pfile);
  delete [] proji;
}

void RealImages3D::AVSRead(char *filename )
     // Read  the image from a file written by AVSWrite.
{
  FILE *pfile;
  int i,j,k;
  //  float val,min,max;
  float *uc;

  uc = new float[nx*ny];

/* open the file for reading in */
  if ((pfile=fopen(filename,"r"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
/* lecture des lignes suivantes */
  for (k=0;k<nz;k++)
    {
      fread(uc,sizeof(float),nx*ny,pfile);
      for (i=0;i<nx;i++)
	for (j=0; j<ny;j++)
	  (*(cube[k]))(i,j)=uc[i*ny+j];
    }
  fclose(pfile);
  delete [] uc;
}

void RealImages3D::INTRead(char *filename )
     // Read  the image from a file written by AVSWrite.
{
  FILE *pfile;
  int i,j,k;
  int *uc;

  uc = new int[nx*ny];

/* open the file for reading in */
  if ((pfile=fopen(filename,"r"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
/* lecture des lignes suivantes */
  for (k=0;k<nz;k++)
    {
      fread(uc,sizeof(int),nx*ny,pfile);
      for (i=0;i<nx;i++)
	for (j=0; j<ny;j++)
	  (*(cube[k]))(i,j)=uc[i*ny+j];
    }
  fclose(pfile);
  delete [] uc;
}

void RealImages3D::Examine() const
// display the contents of the images
{
  int k;
  for (k=0;k<nz;k++)
    {
      cout << " Examining slice " << k << endl;
      cube[k]->Examine();
      cout << " ++++++++++++++++++++++++ " << endl;
    }
      cout << " ========================================" << endl;
      cout << " ========================================" << endl;
}

void RealImages3D::Plus(const Scalaire val)
// add the constant val to each voxel
{
  int k;
  for (k=0;k<nz;k++)
    cube[k]->Plus(val);
}

void  RealImages3D::Mul(const Scalaire m)
// multiply the constant m to each voxel
{
  for (int k=0;k<nz;k++)
    cube[k]->Mul(m);
}

void  RealImages3D::Log()
// Perform a Log Transform
{
  for (int k=0;k<nz;k++)
    cube[k]->Log();
}

void RealImages3D::JWeight(double *weight)
// multiply each column (J direction) term by term by 
// the weighting vector weight
{
  int k;
  for (k=0;k<nz;k++)
    cube[k]->JWeight(weight);
}
void RealImages3D::JConvol(double *filter)
// Convolution of each column (J direction) term by term by 
// the filter filter
// filter is vector of length 2 ny - 1 such that
// filter(ny-1) is the filter function at 0
{
  int k;
  for (k=0;k<nz;k++)
    cube[k]->JConvol(filter);
}
void RealImages3D::JConvol(double *filter,const RealImages3D & image3d)
// idem as previous but *this <- Jconvol of *filter by *image;
{
  int k;
  for (k=0;k<nz;k++)
    cube[k]->JConvol(filter,image3d(k));
}
void RealImages3D::JConvol(double *filter,const int halfwidth, const RealImages3D & image3d)
// idem as previous but *this <- Jconvol of *filter by *image;
// but the convolution width is on [-halfwidth halfwidth]
// filter is of size 2*halfwith+1
// filter(halfwidth) is the filter function at 0
{
  int k;
  for (k=0;k<nz;k++)
    cube[k]->JConvol(filter,halfwidth,image3d(k));
}

Scalaire RealImages3D::Min()
// computes the min
{
  Scalaire mini,estcemini;
  int k;
  mini=cube[0]->Min();
  for (k=1;k<nz;k++){
    estcemini=cube[k]->Min();
    if(estcemini<mini)mini=estcemini;
  }
  return(mini);
}

Scalaire RealImages3D::Max()
// computes the max
{
  Scalaire maxi,estcemaxi;
  int k;
  maxi=cube[0]->Max();
  for (k=1;k<nz;k++){
    estcemaxi=cube[k]->Max();
    if(estcemaxi>maxi)maxi=estcemaxi;
  }
  return(maxi);
}

/* WARNING praxim code is TRASH now

void RealImages3D::GetSlice(RealImage& slice, CVector3D u, CVector3D v, CPoint3D o)
  // Computes an oblique slice in the real image
  // 
  //
{
  Scalaire x1,x2,y1,y2;
  // LD 2003 : Dot(u,u)=> u*u
   Scalaire invnormu=1/sqrt(u*u);
   Scalaire invnormv=1/sqrt(v*v);
   CVector3D vo(o.X(),o.Y(),o.Z());
   x1=(vo*u)*invnormu;
   y1=(vo*v)*invnormv;
   CVector3D end;
   end=vo+slice.NX()*u+slice.NY()*v;
   x2=(end*u)*invnormu;
   y2=(end*v)*invnormv;
   slice.SetCorners(x1, x2, y1, y2);
  
  Scalaire hx=(xend-xbegin)/nx,hy=(yend-ybegin)/ny,hz=(zend-zbegin)/nz;

  Scalaire decax=(o.X()-xbegin)/hx,decay=(o.Y()-ybegin)/hy,decaz=(o.Z()-zbegin)/hz;
  Scalaire dux=u.X()/hx,duy=u.Y()/hy,duz=u.Z()/hz;
  Scalaire dvx=v.X()/hx,dvy=v.Y()/hy,dvz=v.Z()/hz;

  double imdeb,jmdeb,kmdeb;
  double dim,djm,dkm;
  int im,jm,km;

  double l_Dx,l_Dx1,l_Dy,l_Dy1,l_Dz,l_Dz1;

  for(register int j=0; j<slice.NY(); j++){
    imdeb=decax+j*dvx;
    jmdeb=decay+j*dvy;
    kmdeb=decaz+j*dvz;
    
    for(register int i=0; i<slice.NX(); i++)  {
      dim=imdeb+i*dux;
      djm=jmdeb+i*duy;
      dkm=kmdeb+i*duz;
      im=(int) floor(dim-0.5);
      jm=(int) floor(djm-0.5);
      km=(int) floor(dkm-0.5);

      Scalaire interpo = 0;

      if(im>=0 && im<(nx-1) && jm>=0 && jm<(ny-1) && km>=0 && km<(nz-1))
	{
	  // point inside the cube 
	  // mode d'interpolation tri-linéaire
	  // valeur du pixel = interpolation tri-lineaire autour du voxel le plus proche
	  l_Dx = dim-im-0.5;
	  l_Dy = djm-jm-0.5;
	  l_Dz = dkm-km-0.5;
	  l_Dx1 = 1 - l_Dx;
	  l_Dy1 = 1 - l_Dy;
	  l_Dz1 = 1 - l_Dz;

	  interpo=l_Dz1*l_Dy1*l_Dx1 * (*(cube[km]))(im,jm);
	  interpo+=l_Dz1*l_Dy1*l_Dx * (*(cube[km]))(im+1,jm);
	  interpo+=l_Dz1*l_Dy*l_Dx1 * (*(cube[km]))(im,jm+1);
	  interpo+=l_Dz1*l_Dy*l_Dx * (*(cube[km]))(im+1,jm+1);
	  interpo+=l_Dz*l_Dy1*l_Dx1 * (*(cube[km+1]))(im,jm);
	  interpo+=l_Dz*l_Dy1*l_Dx * (*(cube[km+1]))(im+1,jm);
	  interpo+=l_Dz*l_Dy*l_Dx1 * (*(cube[km+1]))(im,jm+1);
	  interpo+=l_Dz*l_Dy*l_Dx * (*(cube[km+1]))(im+1,jm+1);

	  slice(i,j) = interpo;
	      }
	  }
      }
}


void RealImages3D::GetResample(RealImages3D& resample, 
			       CVector3D u, CVector3D v, CVector3D w, 
			       CPoint3D o)
  // Get a resampled cube according to the direction u, v and u vectorial v
{
  Scalaire x1,x2,y1,y2;
  Scalaire z1,z2;
  // LD 2003 modif 
  // Scalaire invnormu=1/sqrt(Dot(u,u));
  // Scalaire invnormv=1/sqrt(Dot(v,v));
  // Scalaire invnormw=1/sqrt(Dot(w,w));
  Scalaire invnormu=1/sqrt((u*u));
  Scalaire invnormv=1/sqrt((v*v));
  Scalaire invnormw=1/sqrt((w*w));

  CVector3D vo(o.X(),o.Y(),o.Z());
   
  //  x1=Dot(o,u)*invnormu;
  //y1=Dot(o,v)*invnormv;
  //z1=Dot(o,w)*invnormw;
  x1=(vo*u)*invnormu;
  y1=(vo*v)*invnormv;
  z1=(vo*w)*invnormw;

  //  CPoint3D end;
  CVector3D end;
  end=vo+resample.NX()*u+resample.NY()*v+resample.NZ()*w;
  //  x2=Dot(end,u)*invnormu;
  //  y2=Dot(end,v)*invnormv;
  //  z2=Dot(end,w)*invnormw;
  x2=(end*u)*invnormu;
  y2=(end*v)*invnormv;
  z2=(end*w)*invnormw;
  resample.SetCorners(x1, x2, y1, y2, z1,z2);

  CPoint3D oz=o;
  int k;
  for (k=0;k<resample.NZ();k++){
    this->GetSlice(resample(k),u,v,oz);
    oz+=w;
  }
}
*/
//
// IO methods
void RealImages3D::PGMWrite(char *filename)
     // Write the image in pgm format in the file filename.
{
  int k;
  char buffer[400];
  for (k=0;k<nz;k++)
    {
      if(k >= 100)
	sprintf(buffer,"%s.%d",filename,k);
      else
	{
	  if(k>=10)
	    sprintf(buffer,"%s.0%d",filename,k);
	  else
	    sprintf(buffer,"%s.00%d",filename,k);
	}
      //      Cout << " writing slice " << k <<" in file "<< buffer << endl;
      cube[k]->PGMWrite(buffer);
    }
}

void RealImages3D::PGMWrite(char *filename, float imagemin, float imagemax)
// Write the image in pgm format in the file filename.
// In RAWBITE
// imagemin and imagemax are used for scaling the data 
// for producing the image
{
  int k;
  char buffer[400];
  for (k=0;k<nz;k++)
    {
      if(k >= 100)
	sprintf(buffer,"%s.%d",filename,k);
      else
	{
	  if(k>=10)
	    sprintf(buffer,"%s.0%d",filename,k);
	  else
	    sprintf(buffer,"%s.00%d",filename,k);
	}
      //      Cout << " writing slice " << k <<" in file "<< buffer << endl;
      cube[k]->PGMWrite(buffer,imagemin,imagemax);
    }
}


void RealImages3D::PGMsWrite(char *filename)
  // Write the image in pgm RawBit format in the file filename.
  // computes the maximum and the minimum of the 3D images.
  // these values are used in each 2D slice for normalization 
  // so that the same scale is used for all slices
{
  int k;
  float maximum;
  float minimum;
  maximum=(float) this->Max();
  minimum=(float) this->Min();
  char buffer[400];
  for (k=0;k<nz;k++)
    {
      if(k >= 100)
	sprintf(buffer,"%s.%d",filename,k);
      else
	{
	  if(k>=10)
	    sprintf(buffer,"%s.0%d",filename,k);
	  else
	    sprintf(buffer,"%s.00%d",filename,k);
	}
      cout << " writing slice " << k <<" in file "<< buffer << endl;
      cube[k]->PGMWrite(buffer,minimum,maximum);
    }
}


void RealImages3D::UShortWrite(char *filename)
  // Write the image in pgm RawBit format in the file filename.
  // computes the maximum and the minimum of the 3D images.
  // these values are used in each 2D slice for normalization 
  // so that the same scale is used for all slices
{
  int k;
  float maximum;
  float minimum;
  maximum=(float) this->Max();
  minimum=(float) this->Min();
  char buffer[400];
  cout << " MIN et MAX du volume 3D "<< minimum<<" ;  "<<maximum<<endl;
  for (k=0;k<nz;k++){
    if(k >= 100)
      sprintf(buffer,"%s.%d",filename,k);
    else{
      if(k>=10)
	sprintf(buffer,"%s.0%d",filename,k);
      else
	sprintf(buffer,"%s.00%d",filename,k);
    }
    cout << " writing slice " << k <<" in file "<< buffer << endl;
    cube[k]->UShortWrite(buffer,minimum,maximum);
  }
}


void RealImages3D::CharWrite(char *filename)
     // Write the image in char format in the file filename.
{
  this->CharWrite(filename,0,nz,1,1,1);
}

void RealImages3D::CharWrite(char *filename, int begz, int endz,
			     int sx, int sy, int sz)
     // Write the image in char format in the file filename.
     // keep only the slices k, begz<= k < endz
     // sx,sy,sz are srides
{
  FILE *pfile;
  int i,j,k;
  int dimx,dimy,dimz;
  Scalaire min=(*(cube[0]))(0,0),max=(*(cube[0]))(0,0);
  Scalaire scale;

  dimx=nx/sx;
  dimy=ny/sy;
  dimz=(endz-begz)/sz;
  for (k=0;k<dimz;k++)
    for (i=0;i<dimx;i++)
      for (j=0;j<dimy;j++)
	{
	  if((*(cube[begz+k*sz]))(i*sx,j*sy)<min)
	    min=(*(cube[begz+k*sz]))(i*sx,j*sy);
	  else
	    if((*(cube[begz+k*sz]))(i*sx,j*sy)>max)
	      max=(*(cube[begz+k*sz]))(i*sx,j*sy);
	}
  if(max!=min) 
    scale=255./(max-min);
  else
    scale=1.;
  
  /* open the file for writing in */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
  unsigned char *uc;
  uc = new unsigned char[dimx*dimy*(endz-begz)];
  for (k=0;k<dimz;k++)
    for (i=0;i<dimx;i++)
      for (j=0;j<dimy;j++)
	{
	  uc[k*dimx*dimy+i*dimy+j]=
	    (unsigned char)(scale*((*(cube[begz+k*sz]))(i*sx,j*sy)-min)) ;
	}
  fwrite(uc,1,dimx*dimy*(endz-begz),pfile);
  cout << " The result is a cube of dim dimx="<< dimx 
    << "; dimy="<< dimy <<"; dimz=" <<dimz <<endl;
  fclose(pfile);
  delete [] uc;
}



void RealImages3D::PGMRead(char *filename)
     // read the image in pgm  format in the file filename.
     // OK for RAWBIT or ASCII....
{
  int k;
  char buffer[400];
  for (k=0;k<nz;k++)
    {
      sprintf(buffer,"%s%d",filename,k);
      cout << " writing slice " << k <<" in file "<< buffer << endl;
      cube[k]->PGMRead(buffer);
    }
}



void Mul(const Scalaire r,const RealImages3D &u,
			    RealImages3D &v)
// v <- r.u
{
  int k,dimz;
  
  dimz=u.NZ();
  if(dimz!=v.NZ())
    {
      cout << " ERR: Mul bad dimensions NZ " << endl;
      return;
    }
  for (k=0;k<dimz;k++)
      Mul(r,u(k),v(k));
}
void Add(const RealImages3D &u,const RealImages3D &v, 
	 RealImages3D &w)
     // w <- u+v
{
  int k,dimz;
  
  dimz=u.NZ();
  if( (dimz!=v.NZ()) || (dimz!=w.NZ()) )
    {
      cout << " ERR: Add bad dimensions NZ " << endl;
      return;
    }
  for (k=0;k<dimz;k++)
      Add(u(k),v(k),w(k));
}

void Sub(const RealImages3D &u,const RealImages3D &v,
	 RealImages3D &w)
     // w <- u-v
{
  int k,dimz;
  
  dimz=u.NZ();
  if( (dimz!=v.NZ()) || (dimz!=w.NZ()) )
    {
      cout << " ERR: Sub bad dimensions NZ " << endl;
      return;
    }
  for (k=0;k<dimz;k++)
      Sub(u(k),v(k),w(k));
} 

void Mul(const RealImages3D &u,const RealImages3D &v,
	 RealImages3D &w)
     // w <- u*v
{
  int k,dimz;
  
  dimz=u.NZ();
  if( (dimz!=v.NZ()) || (dimz!=w.NZ()) )
    {
      cout << " ERR: Sub bad dimensions NZ " << endl;
      return;
    }
  for (k=0;k<dimz;k++)
      Mul(u(k),v(k),w(k));
} 
 
void Div(const RealImages3D &u,const RealImages3D &v,
	 RealImages3D &w)
     // w <- u*v
{
  int k,dimz;
  
  dimz=u.NZ();
  if( (dimz!=v.NZ()) || (dimz!=w.NZ()) )
    {
      cout << " ERR: Sub bad dimensions NZ " << endl;
      return;
    }
  for (k=0;k<dimz;k++)
      Div(u(k),v(k),w(k));
} 
 
Scalaire Norm2(const RealImages3D &a)
     // <- a(1)^2+...+a(n)^2
{
  int k,dimz;
  Scalaire norme2;
  dimz=a.NZ();
  norme2=0.;
  for (k=0;k<dimz;k++)
      norme2=norme2+Norm2(a(k));
  return(norme2);
}
Scalaire  Dist2(const RealImages3D &a,const RealImages3D &b)
     // <- Norm2(b-a)
{
  int k,dimz;
  Scalaire dist2;
  dimz=a.NZ();
  if( (dimz!=b.NZ())  )
    {
      cout << " ERR: Mul bad dimensions NZ " << endl;
      return(-1);
    }
  dist2=0.;
  for (k=0;k<dimz;k++)
      dist2=dist2+Dist2(a(k),b(k));
  return(dist2);
}

