// Tomographie 
// Partie I : projections
//
// LD  Nov 93
// (c) Copyright TIMC 1993

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <tomo.h>
#include <pseudofct.h> 

using namespace std;

//
// Constructeurs et destructeurs
// =============================
//

void Projection::Alloue()
// Alloue n Scalaires dans x et y (privee)
	{
	x=new Scalaire [n];
	y=new Scalaire [n];
	}

Projection::Projection(int d)
// Alloue un Projection initialise de dimension d
	{
	n=d;
	Alloue();
        cerr << "Construction (" << n<< ")" << endl;
	Un();
	}


Projection::~Projection()
// Libere la Projection
	{
	delete [] x;
	delete [] y;
        cerr << "Destruction (" << n<< ")" << endl;
	}

//
// Initialisations
// ===============
//

void Projection::Zero()
// Recoit (0,0,...,0)
	{
	int i;
	phi=0.0;
	for (i=0;i<n;i++)
           {
	     x[i]=(double)(-1.+(double)i*(double)2./(n-1));
		y[i]=(double)0.0;
           }
	}

void Projection::Un()
// Recoit (0,0,...,0)
{
  int i;
  phi=0.0;
  for (i=0;i<n;i++)
    {
      x[i]=(double)(-1+(double)i*(double)2./(n-1));
      y[i]=(double)1.0;
    }
}

void Projection::Gaussian(Scalaire sigma, Scalaire mean, 
			  Scalaire sbegin, Scalaire send)
// Recoit constant * exp(-(x(i)-mean)^2/sigma^2), 
{
  int i;
  Scalaire constant=1./(sigma*sqrt(2*M_PI)), h=(send-sbegin)/(n-1);
  phi=0.0;
  Scalaire sigmadeux=sigma*sigma;
  Scalaire s=sbegin;
  for (i=0;i<n;i++,s+=h)   {
    //    x[i]=(double)(-1+(double)i*(double)2./(n-1));
    x[i]=s;
    y[i]=constant*exp(- (x[i]-mean)*(x[i]-mean)/(2*sigmadeux));
  }
}
  
void Projection::FourierIdealFilter(Scalaire b)
// becomes the Fourier ideal filter in the fourier domain, i.e.,
// |sigma| if |sigma| < b, 0 else.... See Natterer 86 p. 18
//
{
  int i,middle;
  Scalaire sigma;
  middle=(n+1)/2; // "indice" corresponding to the middle... in the FOURIER space
  phi=0.0;
  for (i=0;i<middle;i++)
    {
      x[i]=(double)(-1+(double)i*(double)2./(n-1));
      sigma=i*M_PI/2;
      if(sigma>b)
	y[i]=0.0 ;
      else
	y[i]=sigma;		   
    }
  for (i=middle;i<n;i++)
    {
      x[i]=(double)(-1+(double)i*(double)2./(n-1));
      sigma=(n-i)*M_PI/2;
      if(sigma>b)
	y[i]=0.0 ;
      else
	y[i]=sigma;		   
    }
}


void Projection::FourierSigneFilter()
  // becomes the Fourier SIGNE filter in the fourier domain, i.e.,
  // Related to the Hilbert Transform
  // HILBERT transform in the direct domain when used with the 
  // method Filtre and if the results is multiplied by (-i)
  //
{
  int i,middle;
  Scalaire sigma;
  middle=(n+1)/2; // "indice" corresponding to the middle...in the FOURIER space
  phi=0.0;
  for (i=0;i<middle;i++)
    {
      x[i]=(double)(-1+(double)i*(double)2./(n-1));
      sigma=i*M_PI/2;
      y[i]=1;		   
    }
  for (i=middle;i<n;i++)
    {
      x[i]=(double)(-1+(double)i*(double)2./(n-1));
      sigma=(n-i)*M_PI/2;
      y[i]=-1;
    }
  y[0]=0;
  if((2*middle)==n)y[middle]=0;
}  

void Projection::FourierExponentialIdealFilter(Scalaire b,Scalaire mu)
// becomes the Fourier ideal filter in the fourier domain, i.e.,
// |sigma| if abs(mu) < |sigma| < b, 0 else.... See Natterer 86 p. 18
//
{
  int i,middle;
  Scalaire sigma;
  middle=(int)((n+1)*.5); // "indice" corresponding to the middle...
  // in the Fourier space (see TFD fourt)
  phi=0.0;
  for (i=0;i<middle;i++)
    {
      x[i]=(double)(-1+(double)i*(double)2./(n-1));
      sigma=i*M_PI/2;
      if((sigma> b) || (sigma< ((Scalaire)fabs((double)mu))))
	y[i]=0.0 ;
      else
	y[i]=sigma;		   
    }
  for (i=middle;i<n;i++)
    {
      x[i]=(double)(-1+(double)i*(double)2./(n-1));
      sigma=(n-i)*M_PI/2;
      if((sigma> b) || (sigma< ((Scalaire)fabs((double)mu))))
	y[i]=0.0 ;
      else
	y[i]=sigma;		   
    }
}


void Projection::SheppLoganFilter(Scalaire b,Scalaire sbegin,Scalaire  send)
// becomes the shepp and Logan filter in the Direct domain
// see Natterer 86, p.111.
// b is the bandwith
{
  int i;
  Scalaire middle, h=(send-sbegin)/(n-1);
  middle=(n-1)*.5; // "indice" corresponding to the middle...
  phi=0.0;
  for (i=0;i<n;i++) {
      x[i]=sbegin+i*h;
      y[i]= b*b/(M_PI*M_PI*M_PI*M_PI*(1-4*(i-middle)*(i-middle))) ;
    }
}

void Projection::FBSheppLoganFilter(Scalaire b,Scalaire h)
// becomes the shepp and Logan filter in the Direct domain for Fan Beam
// see Natterer 86, p.111. and 113
// b is the bandwith, h is the angle between 2 successive measurement  
// lines for one fixed position of the source.
{
  int i;
  Scalaire middle,s;
  middle=(n-1)*.5; // "indice" corresponding to the middle...
  // do not forget that the indices of a
  // projection are 1,2,...n
  phi=0.0;
  for (i=0;i<n;i++)
    {
      x[i]=(double)(-2+(double)i*(double)4./(n-1));
      y[i]= 0.;
    }
  //  for (i=(int)(middle/2);i<(int)(n-middle/2);i++)
  for (i=0;i<n;i++)    {
    s=b*sin((i-middle)*h);
    y[i]= b*b/(2.*M_PI*M_PI*M_PI);
    y[i]=y[i]*(M_PI/2-s*sin(s))/((M_PI/2)*(M_PI/2)-s*s) ;
  }
}

void Projection::LowPassHilbertFilter(int q, Scalaire h)
// compute the LowPassHilbertFilter projection of length 2n-1 for a parallel proj of length q
// of the form $(1-cos(b s)) / (\pi s)$ or equiv $(1-cos(2 \pi c s)) / (\pi s) $ with $b=2\pi c$
// with the cutoff frequency c=1/(2h) or $b=\pi/h$ !!!!
// at ih i=-n+1:n-1, LPHF(i)=(1-\cos(\pi i)) / (\pi i h) and is 0 at i=0.
// and the indices are translated for starting at 0 (+n-1)
{
  int i;
  if((2*q-1)!=n){
    cout << "error in  LowPassHilbertFilterr(int q, Scalaire h) " ;
    cout << " the projection dim must be 2q-1 !! " ;
    cout << endl;
    exit(0);
  }
  Scalaire s=- (q-1)*h;
  for (i=0;i<n;i++)   { 
    x[i]=s; s+=h; y[i]=0;
  }
  Scalaire deuxh=2*h;
  Scalaire ih=h;
  for (i=q;i<n;i+=2) {
    y[i]=2. /(M_PI * ih);
    ih+=deuxh;
  }
  ih=-h; //  for (i=0;i<q-1;i+=2)
  for (i=q-2;i>-1;i-=2){
    y[i]=2. /(M_PI * ih);
    ih-=deuxh;
  }
}
void Projection::LowPassHilbertFilter(Scalaire b, 
				      Scalaire sbegin, Scalaire send)
// build low pass hilbert  filter in the Direct domain
// see Natterer 2001 inverse problems 17,pp. 113-119
// b is the cut off bandwith (with the Natterer FT def... WARNING)
// return        x[i]=sbegin + i*h;      bxi=b*x[i];
//                  y[i]= ((1- cos( b*xi  ))/ x[i])/M_PI;
//      instead of ( (1-cos(2M_PI c x[i]) / x[i] ) / M_PI  for the classical FT (with 2\pi factor)
// in fact, is the same with 2M_PI c=b
// sbegin and send define the interval (support) of the filter
// warning : sbegin and send are sampled 
{
  int i;
  int middle;
  Scalaire h; // sampling step
  double bxi;
  if((n/2)*2==n){
    cout << "error in  LowPassHilbertFilter " ;
    cout << " the projection dim must be odd " ;
    cout << endl;
    exit(0);
  }
  if(n<2){
    cout << "error in  LowPassHilbertFilter " ;
    cout << " the projection dim must be >2 " ;
    cout << endl;
    exit(0);
  }
  middle=(n-1)/2; // "indice" corresponding to the middle...  q-1 with n=2q-1
  phi=0.0;
  h=(send-sbegin)/(n-1);
  for (i=0;i<middle;i++)
    {
      x[i]=sbegin + i*h;
      bxi=b*x[i];
      y[i]= ((1- cos( bxi  ))/ x[i])/M_PI;
    }
  for (i=middle+1; i<n;i++)
    {
      x[i]=sbegin + i*h;
      bxi=b*x[i];
      y[i]= ((1- cos( bxi  ))/ x[i])/M_PI;
    }
  x[middle]=0.;
  y[middle]=0.;
}

void Projection:: FBHilbertBandlimitedFilter(Scalaire b,  
				    Scalaire alphabegin, Scalaire alphaend)
 // computes the Fan Beam Hilbert filter truncated in the Fourier domain
// h_H is such that FT(h_H)(\sigma)=-i \sign(\sigma)\chi_{[-b;b]}(\sigma)
// thus h_H(s)=1/(\pi s)  (1-\cos(2\pi b s)) = \frac{2}{\pi s} \sin^2(\pi b s)
  // see  Clacdoyle and Defrise paper  (IEEE SIGNAL PROCESSING MAGAZINE [66] JULY 2010)  
  // for the ROI reconstruction application
  // b is the bandwith, h is the angle between 2 successive projection  
  // lines for one fixed position of the source.
  // This is for equiangular FB projections
  // For FBsinnogramme use h=fanangle/(nbtrans-1)
{
  Scalaire halpha= (alphaend-alphabegin)/(Scalaire)(n-1); 
  int i;
  Scalaire pisb; //  Scalaire middle,s;
  // middle=(n-1)*.5; // "indice" corresponding to the middle...
  phi=0.0;
 for (i=0;i<n;i++)    x[i]=alphabegin+i*halpha;
  //  for (i=(int)(middle/2);i<(int)(n-middle/2);i++)
  for (i=0;i<n;i++)    {
    // //    s=b*sin((i-middle)*h);
    //s= sin(alphabegin+i*halpha);
    pisb=M_PI* sin(x[i]) *b;
    if(pisb*pisb < 0.0000000001)
      // y[i]=2* M_PI*s*b*b;  // we have approximated sin 
      y[i]=2* pisb *b;  // we have approximated sin(mpisb) by mpisb
    else
      y[i]=2* (sin(pisb)*sin(pisb))  /  (pisb) *b ;
  }
}

void Projection::FBLambdaFilter(Scalaire b,Scalaire h,Scalaire m)
// becomes the  filter in the Direct domain for Lamda Tomography 
// in the case of Fan Beam Geometry.
// see Faridani 92, SIAM J. Math Appl., Vol 52, no 2.
// see also Natterer 86, p.111. and 113
// b is the bandwith, h is the angle between 2 projections.
// m is a regularity parameter of the filter (see Faridani).
//
// REMARK "r" of Faridani is replace by "b" of Natterer 
// (detail size replaced by band limit)
{
  int i;
  Scalaire middle,s,ssquare;
  middle=(n-1)*.5; // "indice" corresponding to the middle...
  // do not forget that the indices of a
  // projection are 1,2,...n
  phi=0.0;
  for (i=0;i<n;i++)
    {
      x[i]=(double)(-1+(double)i*(double)2./(n-1));
      y[i]= 0.;
    }
  for (i=(int)(middle/2);i<(int)(n-middle/2);i++)
    {
      s=b*sin((i-middle)*h);
      ssquare=s*s;
      if(ssquare<1)
	y[i]= b*b*b*exp((m-1)*log(1-ssquare))*(1-(2*m+1)*ssquare);
      else
	y[i]=0;
    }
}

void Projection::PseudolocFilter(double m, double rho,double eps, int ninteg)
// becomes the  filter in the Direct domain for Lamda Tomography 
// in the case of parallel geometry....
// m is a regularity parameter of the filter (see Faridani).
// rho is a parameter of the pseudo-loc function 
//             (part of the support of the pseudoloc  filter ....)
// eps is 1/resolution (support of the approximation of the delta function 
//      to which the original reconstructed (pseudo-loc) function fpl 
//      is convoluted....)
// ninteg number of points for the integral computation 
//       precision of integration....
// 
{
  int i;
  Scalaire middle,s;
  middle=(n-1)*.5; // "indice" corresponding to the middle...
  // do not forget that the indices of a
  // projection are 1,2,...n
  phi=0.0;
  for (i=0;i<n;i++)
    {
      x[i]=(double)(-1+(double)i*(double)2./(n-1));
      y[i]=0.;
    }
  for (i=(int)(middle/2);i<(int)(n-middle/2);i++)
    {
      s=x[i];
      if((s>(-1)) &&  (s<1))
	y[i]= filtre_omegaprime(m, s, rho, eps, ninteg);
      else
	y[i]=0;
    }
}

void Projection::SheppLoganFilter(Scalaire b,Scalaire mu)
// we try to generalize the shepp and Logan filter 
// in the Direct domain, for the exponential Radon Transform
// see Natterer 86, p.111., p. 46 and my own research....
// b is the bandwith, mu is the attenuation coeff.
	{
	int i,l;
	Scalaire middle,coeff;
	middle=(n-1)*.5; // "indice" corresponding to the middle...
			   // do not forget that the indices of a
			   // projection are 1,2,...n
	phi=0.0;
	for (i=0;i<n;i++)
           {
	     x[i]=(double)(-1+(double)i*(double)2./(n-1));
	     l=(int)(i-middle);
	     coeff= M_PI*mu/b;
	     y[i]=(cos(coeff*l)*cos(coeff/2.)+2*l*sin(coeff*l)*sin(coeff/2.))*b*b
	       /(M_PI*M_PI*M_PI*M_PI*(1-4*l*l)) ;
           }
	}


void Projection::LinearapodFilter(Scalaire b, Scalaire epsilon,Scalaire mu)
// we try to generalize classical e-filer see Natterer 86, p.109
// in the Direct domain, for the exponential Radon Transform
// For classical Radon Transform, use the default value mu=0 (then exp(-mu t)=1)
// see Natterer 86, p.109., p. 46 and my own research....
// b is the bandwith, 
// epsilon must be set between 0 and 1 : 
//         when 0 no apodisation
//         when 1 great apodisation 
//              in the fourier domain the response(sigma) is of type
//                                    1-epsilon*sigma 
// mu is the attenuation coeff. default is zero for the classical 
//      Radon transform

	{
	int i,l;
	Scalaire middle,coeff;
	Scalaire coslpi, rho,rhopi,u,v;
	middle=(n-1)*.5; // "indice" corresponding to the middle...
			   // do not forget that the indices of a
			   // projection are 1,2,...n
	phi=0.0;
	rho = mu/b;
	for (i=0;i<n;i++)
           {
	     x[i]=(double)(-1+(double)i*(double)2./(n-1));
		l=(int)(i-middle);
		if(l==0)
		   {
		   u=(1-rho*rho)/2;
		   v=1./3.-rho/2+rho*rho*rho/6;
		   y[i]=b*b*(u-epsilon*v)/(4*M_PI*M_PI);
		   }
		else
	 	   {
		    if ((l&1)==0) // l even
			coslpi=1;
              	    else          // l odd
		      	coslpi=-1;
		    rhopi= M_PI*rho;
		    u=(coslpi-cos(rhopi*l))/(l*l*M_PI*M_PI) -
                           rho*sin(rhopi*l)/(l*M_PI); 
		    v=((2-rho)*coslpi-rho*cos(rhopi*l))/(l*l*M_PI*M_PI)+
                        2*sin(rhopi*l)/(l*l*l*M_PI*M_PI*M_PI);
	            y[i]=b*b*(u-epsilon*v)/(4*M_PI*M_PI);
		    }
           }
	}


// void Projection::Lit(char *filename)
// {
//   FILE *pfile;
//   int i;
// /* ouverture du fichier en lecture */
//   if ((pfile=fopen(filename,"r"))==NULL) { 
//     printf("fichier %s non trouve \n",filename);
//     exit(1);
//   }
// }

void Projection::Ecrit(char *filename)
{
  FILE *pfile;
  int i;
  /* ouverture du fichier en ecriture */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("ERR- file %s not found \n",filename);
      exit(1);
    }
  /* write on the first line : nbrotation   nbtranslations */
  fprintf(pfile,"# %d points\n",n);
  fprintf(pfile,"# Angle %e \n",phi);
  for (i=0;i<n;i++)
    fprintf(pfile,"%e   %e \n",x[i],y[i]);
   fclose(pfile);
}

//
// Acces structure
// ===============
//

Scalaire & Projection::operator () (int i)
// Acces a la ieme composante
	{
#if CHECK
	if (i<1 || i>n)
		fputs("Error Projection::operator (): Bad address, out of range\n",stderr);
#endif
	return *(y + (i-1));
	}

void Projection::Examine()
// examen du contenu de P
	{ 
	int i;
	cout << " Dimension : " << n << endl;
	cout << " Angle : " << phi << endl;
	cout << " X : " << endl;
	for (i=0;i<n;i++)
		cout << x[i] << " " ;
	cout  << endl;
	cout << " Y : " << endl;
	for (i=0;i<n;i++)
                cout << y[i] << " " ;
	cout  << endl;
	// cout << " X Y : " << endl;
	// for (i=0;i<n;i++)
        //         cout << x[i] << " " << y[i] << " " << endl ;
	// cout  << endl;
	cout << " ============ " << endl;
	}

int Projection::N()
// Taille
	{
	return n;
	}

Scalaire Projection::A()
// Angle
	{
	return phi;
	}

Scalaire Projection::X(int i)
// ieme composante de l'abscisse.
	{
#if CHECK
	if (i<1 || i>n )
		fputs("Err : Bad address, out of range\n",stderr);
#endif
	return  *(x + (i-1));
	}

Scalaire Projection::Y(int i)
// ieme composante de l'abscisse.
	{
#if CHECK
	if (i<1 || i>n )
		fputs("Err : Bad address, out of range\n",stderr);
#endif
	return  *(y + (i-1));
	}

void Projection::Xcopy(Projection *a)
// copy de l'abscisse
	{
	int i;
	for (i=1;i<=n;i++)
		x[i-1]=a->X(i);
	}

void Projection::Xcopy(Scalaire *a)
// copy de l'abscisse
	{
	int i;
	for (i=1;i<=n;i++)
		x[i-1]=a[i-1];
	}
void Projection::Acopy(Projection *a)
// copy de l'angle
	{
	phi=a->A();
	}
void Projection::Acopy(Scalaire s)
// copy de l'angle s 
	{
	phi=s;
	}

//
// Operateurs
// ==========
//

void Copy(Projection *a,Projection *b)
// b <- a
	{
	int i;
	b->Acopy(a);
	b->Xcopy(a);
	for (i=1;i<=a->N();i++)
		(*b)(i)=(*a)(i);
	}

void Mul(Scalaire k,Projection *u)
// u <- k.u
{
  int i;
  for (i=1;i<=u->N();i++)
    (*u)(i)=k*(*u)(i);
}
void Mul(Scalaire k,Projection *u,Projection *v)
// v <- k.u
	{
	int i;
	v->Acopy(u);
	v->Xcopy(u);
	for (i=1;i<=u->N();i++)
		(*v)(i)=k*(*u)(i);
	}

void Add(Projection *u,Projection *v,Projection *w)
// w <- u+v
	{
	int i;
#if CHECK
	if (u->N() != v->N())
	fputs("Err : Uncompatible dimenssions\n",stderr);
	for (i=1;i<=u->N();i++)
	    if (pow(u->X(i)-v->X(i),(double)2.)>.000001)
		{
		fputs("Err : Uncompatible abscisses ",stderr);
		printf ("U, V, DIFF %e %e %e \n",u->X(i),v->X(i),u->X(i)-v->X(i));
		}
#endif
	w->Acopy(u);
	w->Xcopy(u);
	for (i=1;i<=u->N();i++)
		(*w)(i)=(*u)(i)+(*v)(i);
	}

void KAdd(Projection *u,Scalaire k,Projection *v,Projection *w)
// w <- u+k.v
	{
	int i;
#if CHECK
	if (u->N() != v->N())
	fputs("Err : Uncompatible dimenssions\n",stderr);
	for (i=1;i<=u->N();i++)
		if (pow(u->X(i)-v->X(i),(double)2.)>.000001)
		fputs("Err : Uncompatible abscisses\n",stderr);
#endif
	w->Acopy(u);
	w->Xcopy(u);
	for (i=1;i<=u->N();i++)
		(*w)(i)=(*u)(i)+k*(*v)(i);
	}

void Prod(Projection *u,Projection *v,Projection *w)
// w <- u.v
	{
	int i;
#if CHECK
	if (u->N() != v->N())
	fputs("Err : Uncompatible dimenssions\n",stderr);
	for (i=1;i<=u->N();i++)
		if (pow(u->X(i)-v->X(i),(double)2.)>.000001)
		fputs("Err : Uncompatible abscisses\n",stderr);
#endif
	w->Acopy(u);
	w->Xcopy(u);
	for (i=1;i<=u->N();i++)
		(*w)(i)=(*u)(i)*(*v)(i);
	}

void Convol(Projection *u,Projection *v,Projection *w)
// w <- u*v
// circular convolution of vector u and v
	{
	int i,j,k;
#if CHECK
	if (u->N() != v->N())
	fputs("Err : Uncompatible dimenssions\n",stderr);
	for (i=1;i<=u->N();i++)
		if (pow(u->X(i)-v->X(i),(double)2.)>.000001)
		fputs("Err : Uncompatible abscisses\n",stderr);
#endif
	w->Zero();
	w->Acopy(u);
	w->Xcopy(u);
	for (i=1;i<=u->N();i++)
	for (j=1;j<=u->N();j++)
	    {
		k=(i-j);
		if(k<=0)k=k+(u->N());
		(*w)(i)=(*w)(i)+(*u)(j)*(*v)(k);
	    }
	}

Scalaire Scal(Projection *u,Projection *v)
// <- <u,v>
	{
	Scalaire s;
	int i;
#if CHECK
	if (u->N() != v->N())
	fputs("Err : Uncompatible dimenssions\n",stderr);
	for (i=1;i<=u->N();i++)
		if (pow(u->X(i)-v->X(i),(double)2.)>.000001)
		fputs("Err : Uncompatible abscisses\n",stderr);
#endif
	for (i=1,s=0.0;i<=u->N();i++)
		s+=(*u)(i)*(*v)(i);
	return s;
	}

void Normalise(Projection *u,Projection *v)
// v <- u/Norme(u)
	{
	Scalaire s;
	s=1.0/Norme(u);
	Mul(s,u,v);
	}

Scalaire Norme1(Projection *u)
// <- som|u_i|
{
  Scalaire s;
  int i;
  for (i=1,s=0.0;i<=u->N();i++)
    s+=fabs((*u)(i));
  return s;
}


Scalaire Norme2(Projection *u)
// <- <u,u>
	{
	return Scal(u,u);
	}

Scalaire Norme(Projection *u)
// <- Norme(u)
	{
	Scalaire s;
	s=Norme2(u);
	if (s<1.0e-6) return s*1.0e3;
	else return sqrt(s);
	}


