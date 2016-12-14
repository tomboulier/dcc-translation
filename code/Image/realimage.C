// operation on realimage used for reconstruction algorithmes
// such as thoose appearing in tomography
//
// LD Dec 93.
//
// (c) Copyright TIMC 1993


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include <realimage.h>
#include <random.h>
// #include <vmi.h>
// #include <operator.h>

using namespace std;


void RealImage::Allocate()
// Allocate ima
	{
	ima=new Scalaire [nx*ny];
	line= new Scalaire* [nx];
	}

//
// Constructeurs et destructeurs
// =============================
//

RealImage::RealImage(int nnx,int nny)
     // Alloue une RealImage of dimension nnx X nny
{
  int i;
  nx=nnx;ny=nny;
  xbegin=-1;
  xend=1;
  ybegin=-1;
  yend=1;
  Allocate();
  for (i=0;i<nx;i++)
    line[i]=&(ima[i*ny]);
  Zero();
}

RealImage::RealImage(int nnx,int nny,
		     Scalaire xb,Scalaire xe,Scalaire yb,Scalaire ye)
     // Alloue une RealImage of dimension nnx X nny
{
  int i;
  nx=nnx;ny=nny;
  xbegin=xb;
  xend=xe;
  ybegin=yb;
  yend=ye;
  Allocate();
  for (i=0;i<nx;i++)
    line[i]=&(ima[i*ny]);
  Zero();
}

RealImage::RealImage(const RealImage & image)
     // constructor by copy
{
  int i,j;
  Scalaire xb,xe,yb,ye;
  nx=image.NX();
  ny=image.NY();
  image.GetCorners(&xb,&xe,&yb,&ye);
  xbegin=xb;
  xend=xe;
  ybegin=yb;
  yend=ye;
  Allocate();
  for (i=0;i<nx;i++)
    {
      line[i]=&(ima[i*ny]);
      for (j=0;j<ny;j++)
	line[i][j]=image(i,j);
    }
}
RealImage::~RealImage()
// Libere la RealImage
	{
	delete [] ima;
	delete [] line;
	}

//
// Initialisations
// ===============
//
void RealImage::Zero()
// set the image content to 0.
{
  int i,j;
  for (i=0;i<nx;i++)
    for (j=0;j<ny;j++)
      line[i][j]=(Scalaire)0.0;
  // i.e.		ima[i*ny+j]=0.0;
}

void RealImage::AddGaussienne(Scalaire sigma,Scalaire xmean,Scalaire ymean)
// Add a Gaussian to an image
// Gaussian is e^{((x-xmean)**2+(y-ymean)**2)/sigma**2}
{
  int i,j;
  Scalaire sigma2=sigma*sigma;
  Scalaire x2,y2,hx,hy,x,y;
  hx=(xend-xbegin)/nx;
  hy=(yend-ybegin)/ny;
  x=xbegin+hx/2;
  for (i=0;i<nx;i++){
    y=ybegin+hy/2;
    x2=(x-xmean)*(x-xmean);
    for (j=0;j<ny;j++){
      y2=(y-ymean)*(y-ymean);
      line[i][j]+=exp(- ((x2+y2)/sigma2) );
      y=y+hy;
    }
    x=x+hx;
  }
}

void RealImage::AddEllipse(Scalaire xaxe,Scalaire yaxe,
		  Scalaire xcenter,Scalaire ycenter,
		  Scalaire density,Scalaire angle)
     //add an ellipse of grey level "density"
     // centered on xcenter, ycenter 
     // of (long, at least should be) semi axis xaxe 
     //  with angle direction "angle"
     // and orthogonal short semi axis yaxe
{
  int i,j;
  Scalaire x,y;
  Scalaire xrot,yrot,xcrs,ycrs,hx,hy;
  Scalaire c=cos((double) angle);
  Scalaire s=sin((double) angle);

  hx=(xend-xbegin)/nx;
  hy=(yend-ybegin)/ny;
  x=xbegin+hx*.5-xcenter;
  for (i=0;i<nx;i++){
    y=ybegin+hy*.5-ycenter;
    for (j=0; j<ny;j++){
      //      xrot=c*x-s*y;
      //      yrot=s*x+c*y;
      xrot=c*x+s*y;
      yrot=-s*x+c*y;
      xcrs=xrot/xaxe;
      ycrs=yrot/yaxe;
      xcrs*=xcrs;
      ycrs*=ycrs;
      if((xcrs+ycrs)<1)
	 line[i][j]+=density;
      y+=hy;
    }
    x+=hx;
  }
}

void  RealImage::AddPhantomDef(char *filename)
//add a Phantom defined in the definition file filaname
// (ellipse, disk, traingle (not yet implemented)
/* open the file for reading in */
{
  FILE *pfile;
  if ((pfile=fopen(filename,"r"))==NULL){
    std::cout<< "file "<<filename<<" not found "<< std::endl;
    exit(1);
  }
  int nbell,nbdisk,nbtri;
  fscanf(pfile,"%d",&nbell);
  fscanf(pfile,"%d",&nbdisk);
  fscanf(pfile,"%d",&nbtri);
  if(nbtri!=0)
    std::cout<<" triangles are not considered ..."<<std::endl;

  Scalaire er, etheta, epsi, ea, eb, edensity; // ellipse parameters
  Scalaire xcenter, ycenter;
  for(int i=0;i<nbell;i++){
    fscanf(pfile,"%lf %lf %lf %lf %lf %lf ",
	   &er,&etheta,&ea,&eb,&epsi,&edensity);
    printf("r, theta, psi, a, b, density = :%e %e %e %e %e %e\n ", 
	   er, etheta, epsi, ea, eb,edensity);
    
    xcenter=er*cos(etheta);
    ycenter=er*sin(etheta);
    this->AddEllipse(ea,eb,xcenter,ycenter,edensity,epsi);
    //    imageIN.AddEllipse(LengthSemiAxisHorizon,LengthSemiAxisVert,
    //         Xcenter,Ycenter,density,PsiInDegree/180*M_PI);
  }
  for(int i=0;i<nbdisk;i++){
    fscanf(pfile,"%lf %lf %lf %lf ",
	   &er,&etheta,&ea,&edensity);
    epsi=0.; eb=ea; 
    printf("r, theta, psi, a, b, density = :%e %e %e %e %e %e\n ", 
	   er, etheta, epsi, ea, eb,edensity);
    xcenter=er*cos(etheta);
    ycenter=er*sin(etheta);
    this->AddEllipse(ea,eb,xcenter,ycenter,edensity,epsi);
  }
  fclose(pfile);
}

void RealImage::TVdenoising(RealImage &denoisedima, RealImage &Zx,RealImage &Zy,
		   int nbiter, Scalaire alpha, Scalaire gamma)
  // denoisedima is the denoised image result
  // Z is a variable homogeneous to grad(*this) must be zero 
  // gamma = 0.8*1/2^n : for two dim 0.8*.25=0.2 (gamma is a progression step)
  // alpha is an hyper parameter (we solve argmin alpha TV(u)+1/2|||x-u|^2
  // nbiter is the iteration number 
{
  int i,j,iter;
  Scalaire normeZ;
  for (i=0;i<nx;i++)
    for (j=0;j<ny;j++){
      Zx(i,j)=0;
      Zy(i,j)=0;
    }
  //  std::cout<<"Zx et Zy = 0  "<<std::endl;
  for(iter=0;iter<nbiter;iter++){
    if((iter-((iter/10)*10))==0)std::cout<<"*"<<std::flush;
    // denoisedima = this - div(Z)=dZx/dx + dZy/dy
    for (i=0;i<nx;i++)
      for (j=0;j<ny;j++)
	denoisedima(i,j)=line[i][j] ;
    // div... WARNING div must be the transpose of grad (in the l2 sense)
    for (i=1;i<nx;i++)
      for (j=0;j<ny;j++)
	denoisedima(i,j)-=  (Zx(i,j)-Zx(i-1,j));
    // i=0
    for (j=0;j<ny;j++)
      denoisedima(0,j)-=  (Zx(0,j));
    // // i=nx-1
    // for (j=0;j<ny;j++)
    //   denoisedima(nx-1,j)-=  (-Zx(nx-2,j));
    // // for (i=0;i<(nx-1);i++)
    // //   for (j=0;j<ny;j++)
    // // 	denoisedima(i,j)-=  (Zx(i+1,j)-Zx(i,j));
    // // // for (i=(nx-1))
    // //   for (j=0;j<ny;j++)
    // // 	denoisedima(nx-1,j)-=  (Zx(0,j)-Zx(nx-1,j));
    for (i=0;i<nx;i++)
      for (j=1;j<ny;j++)
	denoisedima(i,j)-= (Zy(i,j)-Zy(i,j-1));
    //j=0
    for (i=0;i<nx;i++)
      denoisedima(i,0)-= (Zy(i,0));
    // //j=ny-1
    // for (i=0;i<nx;i++)
    //   denoisedima(i,ny-1)-= (Zy(i,ny-2));
    // // for (i=0;i<nx;i++)
    // //   for (j=0;j<(ny-1);j++)
    // // 	denoisedima(i,j)-= (Zy(i,j+1)-Zy(i,j));
    // // for (i=0;i<nx;i++)
    // //   //  for (j=(ny-1))
    // //   denoisedima(i,ny-1)-= (Zy(i,0)-Zy(i,ny-1));

    // Z=Z - gamma grad(denoisedima)
    for (i=0;i<(nx-1);i++)
      for (j=0;j<ny;j++)
	Zx(i,j)-= gamma*(denoisedima(i+1,j)-denoisedima(i,j));
    //for (i=(nx-1))
    for (j=0;j<ny;j++)
      Zx(nx-1,j)-= 0;
    //     Zx(nx-1,j)-= gamma*(denoisedima(0,j)-denoisedima(nx-1,j));
    for (i=0;i<nx;i++)
      for (j=0;j<(ny-1);j++)
	Zy(i,j)-= gamma*(denoisedima(i,j+1)-denoisedima(i,j));
    for (i=0;i<nx;i++)
      Zy(i,ny-1)-= 0;
      //      for (j=(ny-1))
      //Zy(i,ny-1)-= gamma*(denoisedima(i,0)-denoisedima(i,ny-1));
 
    // N(i)=||Z(i)|| 
    for (i=0;i<nx;i++)
      for (j=0;j<ny;j++){
	normeZ=sqrt(Zx(i,j)*Zx(i,j)+Zy(i,j)*Zy(i,j));
	if(normeZ>alpha){
	  Zx(i,j)*=(alpha/normeZ);
	  Zy(i,j)*=(alpha/normeZ);
	}
      }
  }
}
  //    for(i=0;i<nx*ny;i++)
  //if(ima[i]<min)min=ima[i];

/*
void RealImage::AddNoise(Scalaire mean,Scalaire sigma)
     // add a normal noise to this.
{
  int i,j;
  if((mean==0) && (sigma==1))
    {
      for (i=0;i<nx;i++)
	for (j=0;j<ny;j++)
	  line[i][j]+=rndnorm();
    }
  else
    {
      for (i=0;i<nx;i++)
	for (j=0;j<ny;j++)
	  line[i][j]+=rndnorm()*sigma+mean;
    }
}
*/

void RealImage::UnderSampling(RealImage *image,
			      int undrow,int undcol,int iofset, int jofset)
  // set the RealImage to the UnderSampling of *image
  // (*this)(i,j) <- Local Mean Close to (*image)(iofset+undrow*i,jofset+undcol*j)
{
  int i,j,l,c;
  Scalaire mean;
  if( ((iofset+(nx-1)*undrow)>(image->NX())) || 
      ((jofset+(ny-1)*undcol)>(image->NY())) ){
    cout << " Error UnderSampling : Bad dimensions " << endl;
    return;
  }
  int iloc,jloc;
  iloc=iofset;
  for (i=0;i<nx;i++){
    //    cout << endl<<" i "<<i<<endl ;
    //    cout << " j "<<flush ;
    //    iloc=iofset+i*undrow;    
    jloc=jofset;
    for (j=0;j<ny;j++){
      //      cout <<j<<","<<flush;
      mean=0.;
      //      jloc=jofset+j*undcol;
      for(l=0;l<undrow;l++)
	for(c=0;c<undcol;c++)
	  mean+=(*image)(iloc+l,jloc+c);
      line[i][j]=mean;
      jloc+=undcol;
    }
    iloc+=undrow;    
  }
}

void RealImage::EvenInterlacedUnderSampling(RealImage *image,
					    int undrow,int undcol,int iofset, int jofset)
  // set the RealImage to the Interlaced UnderSampling of *image
  // (*this)(i,j) <- Local Mean Close to 
  // (*image)(iofset+undrow*i,jofset + i%2 + undcol*2*j)
{
  int i,j,l,c;
  Scalaire mean;
  if( ((iofset+(nx-1)*undrow)>(image->NX())) || 
      ((jofset+(ny-1)*2*undcol)>(image->NY())) ){
    cout << " Error UnderSampling : Bad dimensions " << endl;
    return;
  }
  int iloc,jloc;
  iloc=iofset;
  for (i=0;i<nx;i++){
    //    cout << endl<<" i "<<i<<endl ;
    //    cout << " j "<<flush ;
    //    iloc=iofset+i*undrow;    
    jloc=jofset + i%2;
    for (j=0;j<ny;j++){
      //      cout <<j<<","<<flush;
      mean=0.;
      //      jloc=jofset+j*undcol;
      for(l=0;l<undrow;l++)
	for(c=0;c<undcol;c++)
	  mean+=(*image)(iloc+l,jloc+c);
      line[i][j]=mean;
      jloc+=2*undcol;
    }
    iloc+=undrow;    
  }
}

void RealImage::OddInterlacedUnderSampling(RealImage *image,
					    int undrow,int undcol,int iofset, int jofset)
  // set the RealImage to the Interlaced UnderSampling of *image
  // (*this)(i,j) <- Local Mean Close to 
  // (*image)(iofset+undrow*i,jofset + (i+1)%2 + undcol*2*j)
{
  int i,j,l,c;
  Scalaire mean;
  if( ((iofset+(nx-1)*undrow)>(image->NX())) || 
      ((jofset+(ny-1)*2*undcol)>(image->NY())) ){
    cout << " Error UnderSampling : Bad dimensions " << endl;
    return;
  }
  int iloc,jloc;
  iloc=iofset;
  for (i=0;i<nx;i++){
    //    cout << endl<<" i "<<i<<endl ;
    //    cout << " j "<<flush ;
    //    iloc=iofset+i*undrow;    
    jloc=jofset + (i+1)%2;
    for (j=0;j<ny;j++){
      //      cout <<j<<","<<flush;
      mean=0.;
      //      jloc=jofset+j*undcol;
      for(l=0;l<undrow;l++)
	for(c=0;c<undcol;c++)
	  mean+=(*image)(iloc+l,jloc+c);
      line[i][j]=mean;
      jloc+=2*undcol;
    }
    iloc+=undrow;    
  }
}
//
// Acces structure
// ===============
//

void RealImage::GetCorners(Scalaire *x1,Scalaire *x2,Scalaire *y1,Scalaire *y2) const
// Get the 4 corners of the image: 
// blc (x1,y1), brc (x2,y1), tlc (x1,y2) ,trc (x2,y2)
{
  *x1=xbegin;
  *x2=xend;
  *y1=ybegin;
  *y2=yend;
}

void RealImage::SetCorners(Scalaire x1,Scalaire x2,Scalaire y1,Scalaire y2) 
// Set the 4 corners of the image: 
{
  xbegin=x1;
  xend=x2;
  ybegin=y1;
  yend=y2;
}
int  RealImage::NX() const
// Pixel number in X.
{
  return(nx);
}
int  RealImage::NY() const
// Pixel number in Y.
{
  return(ny);
}

Scalaire RealImage::Min()
// Minimum of an image
{
  int i;
  Scalaire min=ima[0];
  for(i=0;i<nx*ny;i++)
    if(ima[i]<min)min=ima[i];
  return(min);
}
Scalaire RealImage::Max()
// Maximum of an image
{
  int i;
  Scalaire max=ima[0];
  for(i=0;i<nx*ny;i++)
    if(ima[i]>max)max=ima[i];
  return(max);
}
void RealImage::Plus(const Scalaire val)
// add the constant val to each pixel
{
  int i;
  for(i=0;i<nx*ny;i++)
    ima[i]=ima[i]+val;
}
void  RealImage::Mul(const Scalaire m)
// multiply the constant m to each voxel
{
  for(int i=0;i<nx*ny;i++)
    ima[i]=ima[i]*m;
}
void  RealImage::Log()
// Perform a Log Transform (and put to 0 negative or null values)
{
  for(int i=0;i<nx*ny;i++)
    if(ima[i]>0)ima[i]=log(ima[i]);
    else ima[i]=0;
}


void RealImage::JWeight(double *weight)
// multiply each column (J direction) term by term by 
// the weighting vector weight
{
  for (int i=0;i<nx;i++)
    for (int j=0; j<ny;j++)
      ima[i*ny+j]*=weight[j];
}

void RealImage::JConvol(double *filter)
// Convolution of each column (J direction) term by term by 
// the filter filter
// filter is vector of length 2 ny - 1 such that
// filter(ny-1) is the filter function at 0
{
  Scalaire v[ny];
  for (int i=0;i<nx;i++){
    for (int j=0; j<ny;j++){
      v[j]=line[i][j];
      ima[i*ny+j]=0;
    }
    for (int j=0; j<ny;j++)
      for (int l=0; l<ny;l++)
	ima[i*ny+j]+=filter[j-l+ny-1]*v[l];
  }
}

void RealImage::JConvol(double *filter, const RealImage & image)
// idem as previous but *this <- Jconvol of *filter by *image;
{
  Scalaire v[ny];
  for (int i=0;i<nx;i++){
    for (int j=0; j<ny;j++){
      v[j]=image(i,j);
      ima[i*ny+j]=0;
    }
    for (int j=0; j<ny;j++)
      for (int l=0; l<ny;l++)
	ima[i*ny+j]+=filter[j-l+ny-1]*v[l];
  }
}

void RealImage::JConvol(double *filter, const int halfwidth, const RealImage & image)
// idem as previous but *this <- Jconvol of *filter by *image;
// but the convolution width is on [-halfwidth halfwidth]
// filter is of size 2*halfwith+1
// filter(halfwidth) is the filter function at 0
{
  Scalaire v[ny];
  int begl,endl;
  for (int i=0;i<nx;i++){
    for (int j=0; j<ny;j++){
      v[j]=image(i,j);
      ima[i*ny+j]=0;
    }
    for (int j=0; j<ny;j++){
      begl=max(0,j-halfwidth);endl=min(ny,j+halfwidth);
      for (int l=begl; l<endl;l++)
	ima[i*ny+j]+=filter[j-l+halfwidth]*v[l];
    }
  }
}


void RealImage::Examine() const
// display the contents of the image.
	{
        int i,j;
        cout << " ----------------------------------------- " << endl;
        cout << " ----------------------------------------- " << endl;
        cout << " IMAGE (" << nx <<","<<ny<< "):" << endl;
	cout << " x1:" << xbegin << " x2:" << xend
		<< " y1:" << ybegin<< " y2:" << yend << endl;
        cout << " ----------------------------------------- " << endl;
	for (i=0;i<nx;i++)
		{
		for (j=0;j<ny;j++)
		   cout << ima[i*ny+j] << " " ;
		cout  << endl;
		}
	cout << " =========================================== " << endl;
       	cout << " =========================================== " << endl;
	}

void RealImage::Write(char *filename)
// Write the image in the file filename.
{
  FILE *pfile;
  int i,j;
  float val;

/* open the file for writing in */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }

/* write the number of lines and columns */
  fprintf(pfile,"%d %d\n",nx,ny);
/* lecture des lignes suivantes */
  for (i=0;i<nx;i++)
    for (j=0; j<ny;j++)
      {
	fprintf(pfile,"%e",ima[i*ny+j]);
      }
  fclose(pfile);
}

void RealImage::MHDWrite(char *filename )
// Write the image in a bin file. (j,i)
{
  FILE *pfile;
  int i,j;
  float val,min,max;
  float *uc;

  uc = new float[nx*ny];

/* open the file for writing in */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
  min=max=ima[0];
  for (i=0;i<nx;i++)
    for (j=0; j<ny;j++)
      {
	val=ima[i*ny+j];
	if(min>val)min=val;
	if(max<val)max=val;
      }
  //  cout << " file " << filename << " min, max = " << min <<" ; "<< max << endl;

/* ecriture des lignes suivantes */
  for (j=0; j<ny;j++)
    for (i=0;i<nx;i++)
      uc[i*ny+j]=ima[i*ny+j];
  fwrite(uc,sizeof(float),nx*ny,pfile);
  fclose(pfile);
  delete [] uc;
}

void RealImage::AVSWrite(char *filename )
// Write the image in a bin file.(i,j)
{
  FILE *pfile;
  int i,j;
  float val,min,max;
  float *uc;

  uc = new float[nx*ny];

/* open the file for writing in */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
  min=max=ima[0];
  for (i=0;i<nx;i++)
    for (j=0; j<ny;j++)
      {
	val=ima[i*ny+j];
	if(min>val)min=val;
	if(max<val)max=val;
      }
  //  cout << " file " << filename << " min, max = " << min <<" ; "<< max << endl;

/* write the number of lines and columns */
/*  fprintf(pfile,"AVS\n");
  fprintf(pfile,"# done by RealImage::AVSWrite \n");
  fprintf(pfile,"# Corners x1 %f   x2 %f   y1 %f   y2 %f  \n",
	  (float) xbegin,(float) xend,(float) ybegin,(float) yend );
  fprintf(pfile,"# Min = %f\n",min);
  fprintf(pfile,"# Max = %f\n",max);
  fprintf(pfile,"%d %d %d\n",nx,ny);
*/
/* ecriture des lignes suivantes */
  for (i=0;i<nx;i++)
    for (j=0; j<ny;j++)
      uc[i*ny+j]=ima[i*ny+j];
  fwrite(uc,sizeof(float),nx*ny,pfile);
  fclose(pfile);
  delete [] uc;
}

void RealImage::AVSRead(char *filename )
     // Read  the image from a file written by AVSWrite.
{
  FILE *pfile;
  int i,j;
  float *uc;

/* open the file for reading in */
  if ((pfile=fopen(filename,"r"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
  uc = new float[nx*ny];
/* lecture des lignes suivantes */
  fread(uc,sizeof(float),nx*ny,pfile);
  fclose(pfile);
  for (i=0;i<nx;i++)
    for (j=0; j<ny;j++)
      ima[i*ny+j]=uc[i*ny+j];
  delete [] uc;
}

void RealImage::INTRead(char *filename )
     // Read  the image from a file written by Joe Aoun
{
  FILE *pfile;
  int i,j;
  int *uc;

/* open the file for reading in */
  if ((pfile=fopen(filename,"r"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
  uc = new int[nx*ny];
/* lecture des lignes suivantes */
  fread(uc,sizeof(int),nx*ny,pfile);
  fclose(pfile);
  for (i=0;i<nx;i++)
    for (j=0; j<ny;j++)
      ima[i*ny+j]=uc[i*ny+j];
  delete [] uc;
}

void RealImage::SPECTRead(char *filename )
     // Read  the image from a file from the gamma camera (SPECT)..
{
  FILE *pfile; 
  int i,j;
  float val,min,max;
  unsigned short *proji;

  proji= new unsigned short [nx*ny];

// open the file for reading in 
  if ((pfile=fopen(filename,"r"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
// lecture des lignes suivantes 
  fread (proji, sizeof(short), ny*nx, pfile);

  for (j=0; j<ny;j++)
    for (i=0;i<nx;i++)
      ima[i*ny+j]=proji[i*ny+j];
  fclose(pfile);
  delete [] proji;

}

void RealImage::TRIXELLRead(char *filename, int headsize)
     // Read  the image from a file from the TRIXELL detector..
{
  FILE *pfile; 
  int i,j;
  float val,min,max;

  unsigned short *proji;
  unsigned char *header;

  proji= new unsigned short [nx*ny];
  header= new unsigned char [headsize];

// open the file for reading in 
  if ((pfile=fopen(filename,"r"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
  cout << " Lecture du header ..." << endl;
// lecture des lignes suivantes 
  fread (header, sizeof(char), headsize, pfile);
  cout << " Lecture de l'image ..." << endl;
// lecture des lignes suivantes 
  fread (proji, sizeof(short), ny*nx, pfile);

  for (i=0;i<nx;i++)
    for (j=0; j<ny;j++)
      ima[i*ny+j]=proji[i*ny+j];
  fclose(pfile);
  delete [] proji;
  delete [] header;
}


void RealImage::PGMaWrite(char *filename)
  // Write the image in pgm format in the file filename.
  // format P2 : ASCII 
{
  FILE *pfile;
  int i,j,row,col,rows,cols;
  float val;
  int ival;
  float min,max,scale;
  
  rows=nx;
  cols=ny;
  /* open the file for writing in */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
  max=ima[0];
  min=ima[0];
  for (i=0;i<nx;i++) {
    for (j=0; j<ny;j++){
      val=ima[i*ny+j];
      if(max<val) max = val;
      if(min>val) min = val;
    }
  }
  if(max!=min) 
    scale=255./(max-min);
  else
    scale=1.;
  
  fprintf(pfile,"P2\n");
  fprintf(pfile,"# done by RealImage::PGMaWrite \n");
  fprintf(pfile,"# Corners x1 %f   x2 %f   y1 %f   y2 %f \n",
	  (float) xbegin,(float) xend,(float) ybegin,(float) yend);
  fprintf(pfile,"# Min = %f\n",min);
  fprintf(pfile,"# Max = %f\n",max);
  fprintf(pfile,"%d %d \n",ny,nx);
  fprintf(pfile,"%d\n",255);
  for ( row = 0; row < rows; ++row)
    {
      for ( col = 0; col < cols; ++col)
	{
	  ival = (int) ( scale * ( ima[row*cols+col] - min) );
	  fprintf(pfile,"%d ",ival);
	}
      fprintf(pfile," \n");
    }
  fclose(pfile);
}

void RealImage::PGMWrite(char *filename, float imagemin, float imagemax)
// Write the image in pgm format in the file filename.
// In RAWBITE
// imagemin and imagemax are used for scaling the data 
// for producing the image
{
  FILE *pfile;
  int i,j,row,col,rows,cols;
  float val;
  int ival;
  float min,max,scale;
  unsigned char *uc;

  rows=nx;
  cols=ny;
  /* open the file for writing in */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
  max=ima[0];
  min=ima[0];
  for (i=0;i<nx;i++)
    {
      for (j=0; j<ny;j++)
	{
	  val=ima[i*ny+j];
	  if(max<val) max = val;
	  if(min>val) min = val;
	}
    }
  if(imagemax!=imagemin) 
    scale=255./(imagemax-imagemin);
  else
    scale=1.;
  
  uc = new unsigned char[rows*cols];
  
  fprintf(pfile,"P5\n");
//  fprintf(pfile,"P2\n");
  fprintf(pfile,"# done by RealImage::PGMWrite \n");
  fprintf(pfile,"# Corners x1 %f   x2 %f   y1 %f   y2 %f \n",
	  (float) xbegin,(float) xend,(float) ybegin,(float) yend);
  fprintf(pfile,"# true Min = %f\n",min);
  fprintf(pfile,"# true Max = %f\n",max);
  fprintf(pfile,"# treshold : imagemin = %f; imagemax= %f\n",imagemin,imagemax);
  fprintf(pfile,"%d %d \n",ny,nx);
  fprintf(pfile,"%d\n",255);
  for ( row = 0; row < rows; ++row)
    {
      for ( col = 0; col < cols; ++col)
	{
	  if((ima[row*cols+col])>imagemax)
	    ival=255;
	  else
	    {
	      if((ima[row*cols+col])<imagemin)
		ival=0;
	      else
		ival = (int) ( scale * ( ima[row*cols+col] - imagemin) );
	    }
	  uc[row*cols+col] =  (unsigned char) ival ;
	  //		fprintf(pfile,"%d ",ival);
	}
      //	    fprintf(pfile," \n");
    }
  fwrite(uc,sizeof(unsigned char),rows*cols,pfile);
  fclose(pfile);
  delete [] uc;
}

void RealImage::UShortWrite(char *filename, float imagemin, float imagemax)
// Write the image in Unsigned Short format in the file filename.
// imagemin and imagemax are used for scaling the data 
// for producing the image
  // if imagemin==imagemax==0 then the max and min are computed
{
  FILE *pfile;
  int i,j,row,col,rows,cols;
  float val;
  float min,max,scale;
  unsigned short *uc, usval;

  rows=nx;
  cols=ny;
  /* open the file for writing in */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
  if((imagemax==0)&&(imagemin==0)){
    max=ima[0];
    min=ima[0];
    for (i=0;i<nx;i++){
      for (j=0; j<ny;j++){
	val=ima[i*ny+j];
	if(max<val) max = val;
	if(min>val) min = val;
      }
    }
  }
  else{
    min=imagemin;
    max=imagemax;
  }

  if(max!=min) 
    scale=32767./(max-min);
  else
    scale=1.;

  uc = new unsigned short[rows*cols];
  
  for ( row = 0; row < rows; ++row){
    for ( col = 0; col < cols; ++col){
      if((ima[row*cols+col])>max)
	usval= 32767;
      else{
	if((ima[row*cols+col])<min)
	  usval=0;
	else
	  usval = (unsigned short) floor(scale * (ima[row*cols+col] - min) );
      }
      uc[row*cols+col] = usval ;
      //		fprintf(pfile,"%d ",ival);
    }
    //	    fprintf(pfile," \n");
  }
  fwrite(uc,sizeof(unsigned short),rows*cols,pfile);
  fclose(pfile);
  delete [] uc;
}

void RealImage::PGMRead(char *filename)
// Read the image in pgm format in the file filename.
// OK for RAWBIT or ASCII....

	{
	FILE *pfile;
	int i,j,row,col,rows,cols;
	float val;
	float x1,x2,y1,y2;
	int ival,imaxi;
	float min=0,max=0,scale;
	char magic[2];
	int dummy;
	char dummystring[1000];
	int ascii,rawbit; /* logical */

	rows=nx;
	cols=ny;
/* open the file for reading in */
	if ((pfile=fopen(filename,"r"))==NULL)
	     { printf("file %s non found \n",filename);
	       return;
	     }

	ascii=0;rawbit=0;
	fscanf(pfile,"%s \n",magic);
	//	printf(" magic %s \n",magic);
	if( (magic[0]=='P') && (magic[1]=='2') )
	  ascii=1;
	if( (magic[0]=='P') && (magic[1]=='5') )
	  rawbit=1;
	if( (ascii==0) && (rawbit==0))
	  {
	    printf(" Error Magic Number ! not pgm \n");
	    fclose(pfile);
	    return;
	  }
//	fscanf(pfile,"%s",dummystring);
//	fgets(s1,1,pfile);
	/*
	dummy=fgetc(pfile);
	fgets(dummystring,1000,pfile);
	fscanf(pfile,"%*s %*s %*s %f %*s %f %*s %f %*s %f ",
	       &x1,&x2,&y1,&y2);

	xbegin= (Scalaire)x1;
	xend=(Scalaire)x2;
	ybegin=(Scalaire)y1;
	yend=(Scalaire)y2;

	printf(" xbegin=%f  xend=%f ybegin=%f yend=%f \n",x1,x2,y1,y2);

	fscanf(pfile,"%*s %*s %*s %f",&min);
	fscanf(pfile,"%*s %*s %*s %f",&max); 
	printf(" min=%f  max=%f \n",min,max);
	fscanf(pfile,"%d %d",&cols,&rows);
	*/

	
	fgets(dummystring,999,pfile);
	while(dummystring[0]=='#')
	  {
	    //	    printf("%s \n",dummystring);
	    fgets(dummystring,999,pfile);
	  }
	  
	sscanf(dummystring,"%d %d",&cols,&rows);
	

	if((rows!=nx) || (cols!=ny))
	  {
	    printf(" Error Bad dimensions \n");
	    fclose(pfile);
	    return;
	  }
	else
	  {
	    nx=rows;
	    ny=cols;
	  }
	
	//	printf(" j'ai lu nx=%d ny=%d \n",rows,cols);
	fscanf(pfile,"%d \n",&imaxi);
	//printf(" imaxi=%d \n",imaxi);

	
	if(max!=min) 
	  scale=(max-min)/imaxi;
	else
	   scale=1.;

	if(ascii)
	  {
	    //	    printf(" ASCII \n");
	    for ( row = 0; row < rows; ++row)
	      {
		for ( col = 0; col < cols; ++col)
		  {
		    fscanf(pfile,"%d ",&ival);
		    ima[row*cols+col]=  scale * ival + min;
		  }
	      }
	  }
	else
	  {	
	    //	    printf(" RAWBIT \n");
	    // dummy=fgetc(pfile); // read a RC
	    unsigned char *uc;
	    uc = new unsigned char[rows*cols];
	    fread(uc,1,rows*cols,pfile);
	    for ( row = 0; row < rows; ++row)
	      for ( col = 0; col < cols; ++col)
		ima[row*cols+col]=uc[row*cols+col]*scale+min ;
	    delete [] uc;
	  }
	fclose(pfile);
	}

void RealImage::PGMWrite(char *filename)
// Write the image in pgm RawBit format in the file filename.
{
  FILE *pfile;
  int i,j,row,col,rows,cols;
  float val;
  int ival;
  float min,max,scale;
  unsigned char *uc;
  
  rows=nx;
  cols=ny;
  /* open the file for writing in */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("file %s non found \n",filename);
    exit(1);
    }
  max=ima[0];
  min=ima[0];
  for (i=0;i<nx;i++)
    {
      for (j=0; j<ny;j++)
	{
	  val=ima[i*ny+j];
	  if(max<val) max = val;
	  if(min>val) min = val;
	}
    }
  
  /*
  cout << "-----------------------------------"<< endl;
  cout << " the minimum is : " << min << endl;
  cout << " the maximum is : " << max << endl;
  cout << "-----------------------------------"<< endl;
  */  
  if(max!=min) 
    scale=255./(max-min);
  else
    scale=1.;

  uc = new unsigned char[rows*cols];
  
  fprintf(pfile,"P5\n");
  fprintf(pfile,"# done by RealImage::PGMWrite \n");
  fprintf(pfile,"# Corners x1 %f   x2 %f   y1 %f   y2 %f \n",
	  (float) xbegin,(float) xend,(float) ybegin,(float) yend );
  fprintf(pfile,"# Min = %f\n",min);
  fprintf(pfile,"# Max = %f\n",max);
  fprintf(pfile,"%d %d \n",ny,nx);
  fprintf(pfile,"%d\n",255);   /* NO WHITE SPACE !!! */
  for ( row = 0; row < rows; ++row)
    for ( col = 0; col < cols; ++col)
      uc[row*cols+col] = 
	(unsigned char) ( scale * ( ima[row*cols+col] - min) );
  fwrite(uc,sizeof(unsigned char),rows*cols,pfile);
  fclose(pfile);
  delete [] uc;
}

// Catherine : mars 99
void RealImage::RAWrbWrite(char *filename)
// Write the image in RawBit (pgm without header) format in the file filename.
  
{
  FILE *pfile;
  int i,j,row,col,rows,cols;
  float val;
  int ival;
  float min,max,scale;
  unsigned char *uc;
  
  rows=nx;
  cols=ny;
/* open the file for writing in */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("file %s non found \n",filename);
    exit(1);
    }
  max=ima[0];
  min=ima[0];
  for (i=0;i<nx;i++)
    {
      for (j=0; j<ny;j++)
	{
	  val=ima[i*ny+j];
	  if(max<val) max = val;
	  if(min>val) min = val;
	}
    }
  cout << "-----------------------------------"<< endl;
  cout << " the minimum is : " << min << endl;
  cout << " the maximum is : " << max << endl;
  cout << "-----------------------------------"<< endl;
  if(max!=min) 
    scale=255./(max-min);
  else
    scale=1.;
  
  uc = new unsigned char[rows*cols];
  
  for ( row = 0; row < rows; ++row)
    for ( col = 0; col < cols; ++col)
      uc[row*cols+col] = 
	(unsigned char) ( scale * ( ima[row*cols+col] - min) );
  fwrite(uc,sizeof(unsigned char),rows*cols,pfile);
  fclose(pfile);
  delete [] uc;
}

void RealImage::VMIWrite(char *filename)
// Write the image in vmi format in the file filename.
{
  FILE *pfile;
  int i,j,row,col,rows,cols;
  float val,min,max,scale;
  char* data;
  
  min=ima[0];
  max=ima[0];
  for (i=0;i<nx;i++)
    for (j=0; j<ny;j++)
      {
	val=ima[i*ny+j];
	if(max<val) max = val;
	if(min>val) min = val;
      }
  cout << "-----------------------------------"<< endl;
  cout << " the minimum is : " << min << endl;
  cout << " the maximum is : " << max << endl;
  cout << "-----------------------------------"<< endl;
  if(max!=min) 
    scale=127./(max-min);
  else
    scale=1.;
  
  data = new char(nx*ny);
  for (i=0;i<nx;i++)
    for (j=0; j<ny;j++)
      *(data++)=(char) (scale*(ima[i*ny+j]- min));
  printf(" no more possible: ask  F. Leitner");
  //     	Imio_WriteImage(filename,data,(short)nx,(short) ny,sizeof(char));
}

void Mul(const Scalaire r,const RealImage &u,RealImage &v)
// v <- r.u
{
  int i,j,dimx,dimy;
  dimx=u.NX();
  dimy=u.NY();
  if( (dimx!=v.NX())||(dimy!=v.NY()) )
    {
      cout << " ERR: Mul bad dimensions NX or NY " << endl;
      return;
    }  
  for (i=0;i<dimx;i++)
    for (j=0; j<dimy;j++)
      v(i,j)=r*u(i,j);
}

void Add(const RealImage &u,const RealImage &v, 
		  RealImage &w)
     // w <- u+v
{
  int i,j,dimx,dimy;
  dimx=u.NX();
  dimy=u.NY();
  if( (dimx!=v.NX())||(dimy!=v.NY()) || (dimx!=w.NX())||(dimy!=w.NY()) )
    {
      cout << " ERR: Add bad dimensions NX or NY " << endl;
      return;
    }  
  for (i=0;i<dimx;i++)
    for (j=0; j<dimy;j++)
      w(i,j)=v(i,j)+u(i,j);
}

void Sub(const RealImage &u,const RealImage &v,
		  RealImage &w)
     // w <- u-v
{
  int i,j,dimx,dimy;
  dimx=u.NX();
  dimy=u.NY();
  if( (dimx!=v.NX())||(dimy!=v.NY()) || (dimx!=w.NX())||(dimy!=w.NY()) )
    {
      cout << " ERR: Sub bad dimensions NX or NY " << endl;
      return;
    }  
  for (i=0;i<dimx;i++)
    for (j=0; j<dimy;j++)
      w(i,j)=u(i,j)-v(i,j);
}

void Mul(const RealImage &u,const RealImage &v,
		  RealImage &w)
     // w <- u*v
{
  int i,j,dimx,dimy;
  dimx=u.NX();
  dimy=u.NY();
  if( (dimx!=v.NX())||(dimy!=v.NY()) || (dimx!=w.NX())||(dimy!=w.NY()) )
    {
      cout << " ERR: Sub bad dimensions NX or NY " << endl;
      return;
    }  
  for (i=0;i<dimx;i++)
    for (j=0; j<dimy;j++)
      w(i,j)=u(i,j)*v(i,j);
}

void Div(const RealImage &u,const RealImage &v,
		  RealImage &w)
     // w <- u/v
{
  int i,j,dimx,dimy;
  dimx=u.NX();
  dimy=u.NY();
  if( (dimx!=v.NX())||(dimy!=v.NY()) || (dimx!=w.NX())||(dimy!=w.NY()) )
    {
      cout << " ERR: Sub bad dimensions NX or NY " << endl;
      return;
    }  
  for (i=0;i<dimx;i++)
    for (j=0; j<dimy;j++)
      if(v(i,j)!=0) w(i,j)=u(i,j)/v(i,j);
}

Scalaire Norm2(const RealImage &a)
     // <- a(1)^2+...+a(n)^2
{
  int i,j,dimx,dimy;
  Scalaire norme;
  dimx=a.NX();
  dimy=a.NY();
  norme=0.;
  for (i=0;i<dimx;i++)
    for (j=0; j<dimy;j++)
      norme=norme + a(i,j)*a(i,j);
  return(norme);
}
Scalaire  Dist2(const RealImage &a,const RealImage &b)
     // <- Norm2(b-a)
{
  int i,j,dimx,dimy;
  Scalaire dist2;
  dimx=a.NX();
  dimy=a.NY();
  if( (dimx!=b.NX())||(dimy!=b.NY()) )
    {
      cout << " ERR: Dist2 bad dimensions NX or NY " << endl;
      return(-1.);
    }  
  dist2=0;
  for (i=0;i<dimx;i++)
    for (j=0; j<dimy;j++)
      dist2=dist2 + (a(i,j)-b(i,j))*(a(i,j)-b(i,j));
  return(dist2);
}  

int LineInterImage(const RealImage &image,
		   const Scalaire cosphi, const Scalaire sinphi, const Scalaire s,
		   Scalaire *intersection, int *iline, int *jcol )
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
{
  int nx=image.NX(), ny=image.NY();
  Scalaire xb,xe,yb,ye;
  image.GetCorners(&xb,&xe,&yb,&ye);
  Scalaire dx=(xe-xb)/nx,dy=(ye-yb)/ny;
  
  //compute the two extrem points
  /*  yphisxb=(s-cosphi*xb)/sinphi; 
  yphisxe=(s-cosphi*xe)/sinphi; 
  xphisyb=(s-sinphi*yb)/cosphi;
  xphisye=(s-sinphi*ye)/cosphi;*/
  Scalaire txb=(s*cosphi-xb)/sinphi;
  Scalaire txe=(s*cosphi-xe)/sinphi;
  Scalaire tyb=(yb-s*sinphi)/cosphi;
  //  Scalaire tye=(ye-s*sinphi)/cosphi;
  //
  int sensx; // sens of progression with dx (sign of -sinphi) 
  Scalaire dty=dy/cosphi; // dy progression (j++) makes t progress of dty 
                         // (always >0)
  Scalaire dtx=-dx/sinphi; // dx progression(i++) makes t progress of dtx 
                             // dtx SIGN depends on sinphi sign )
  int i,j; // (i,j) coordinate of the current pixel
  Scalaire t,txproc, typroc ;// t intersection with the (i,j) pixel
    // typroc is t intersection with (i,j+1) 
    // txproc is t intersection with (i+1,j)  (or (i-1,j)) pixel 
  if((-sinphi)>0){
    sensx=1;
    if(txb<tyb){
      t=tyb;
      j=0;
      Scalaire xphisyb=(s-sinphi*yb)/cosphi;
      i=floor((xphisyb-xb)/dx);
      typroc=t+dty;
      txproc=txb+(i+1)*dtx;
    }
    else{
      t=txb;
      i=0;
      Scalaire yphisxb=(s-cosphi*xb)/sinphi;
      j=floor((yphisxb-yb)/dy);
      txproc=t+dtx;
      typroc=tyb+(j+1)*dty;
    }
  }
  else {
    sensx=-1;
    // dtx is negative
    //    dtx=dx/sinphi;
    if(txe<tyb){
      t=tyb;
      j=0;
      Scalaire xphisyb=(s-sinphi*yb)/cosphi;
      i=floor((xphisyb-xb)/dx);
      typroc=t+dty;
      txproc=txb+i*dtx; // sensx is -1 and dtx negative here
    }
    else{
      t=txe;
      i=nx-1;
      Scalaire yphisxe=(s-cosphi*xe)/sinphi;
      j=floor((yphisxe-yb)/dy);
      txproc=t-dtx; // sensx is -1 and dtx negative here
      typroc=tyb+(j+1)*dty;
    }
  }
  // //  std::cout<<" t="<<t<<std::flush;
  // invariant : on entre dans le pixel i,j
  // t, txproc and t
  int l=0;
  if(sensx==1)
    while ( (i>=0 && i<nx) && (j>=0 && j<ny)) {
      iline[l]=i;
      jcol[l]=j;
      if(txproc<typroc){
	intersection[l]=txproc-t;
	t=txproc;
	i++;
	txproc+=dtx;
      }
      else{
	intersection[l]=typroc-t;
	t=typroc;
	j++;
	typroc+=dty;
      }
      l++;
      // //      std::cout<<" t="<<t<<std::flush;
    }
  else
    while ( (i>=0 && i<nx) && (j>=0 && j<ny)) {
      iline[l]=i;
      jcol[l]=j;
      if(txproc<typroc){
	intersection[l]=txproc-t;	
	t=txproc;
	i--;     // sensx is -1 
	txproc-=dtx; // and dtx negative here
      }
      else{
	intersection[l]=typroc-t;
	t=typroc;
	j++;
	typroc+=dty;
      }
      l++;
      // //      std::cout<<" t="<<t<<std::flush;
    }
  // //  std::cout<<std::endl;
  return(l);
}

int HorizonLineInterImage(const RealImage &image,
		   const Scalaire s,
		   Scalaire *intersection, int *iline, int *jcol )
// case where cosphi=0 \vzeta is (1,0) or (-1,0) sinphi= +or- 1
// WARNING we assume \phi=-\pi/2 thus -sinphi=1
{
  int nx=image.NX(), ny=image.NY();
  Scalaire xb,xe,yb,ye;
  image.GetCorners(&xb,&xe,&yb,&ye);
  Scalaire dx=(xe-xb)/nx,dy=(ye-yb)/ny;
  int i,j; // (i,j) coordinate of the current pixel
  if ((s>yb)&&(s<ye)){
    j=floor((s-yb)/dy);
    for(i=0;i<nx;i++){
      intersection[i]=dx;
      iline[i]=i;
      jcol[i]=j;
    }  
    return(nx);
  }
  else
    return(0);
}

int VertLineInterImage(const RealImage &image,
		   const Scalaire s,
		   Scalaire *intersection, int *iline, int *jcol )
// case where sinphi=0 \vzeta is (0,1) or (0,-1) cosphi= +or- 1
// WARNING we assume \phi=0 thus cosphi=1
{
  int nx=image.NX(), ny=image.NY();
  Scalaire xb,xe,yb,ye;
  image.GetCorners(&xb,&xe,&yb,&ye);
  Scalaire dx=(xe-xb)/nx,dy=(ye-yb)/ny;
  int i,j; // (i,j) coordinate of the current pixel
  if ((s>xb)&&(s<xe)){
    i=floor((s-xb)/dx);
    for(j=0;j<ny;j++){
      intersection[j]=dy;
      iline[j]=i;
      jcol[j]=j;
    }  
    return(ny);
  }
  else
    return(0);
}
