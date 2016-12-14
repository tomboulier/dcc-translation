// Operation sur les sinogrammes 
// en particulier retroprojection
//
// LD Dec 93

// (c) Copyright TIMC 1993

#ifndef __SINOGRAMME_H
#define __SINOGRAMME_H

#include <cmath>

// #include <alp.h>
// #include <alp.h> // modif 24/02/97 pour utiliser alp de Laurent
//#include <FBsinogramme.h>
#include <projection.h>
#include <realimage.h>


class Sinogramme
{
 protected:
  int p,q;	// p=Nombre de  projections (rotations de la source)
                // q=nombre d'echantillons par projection.
  Scalaire *g;	// Coefficients (mesures)
                // g[i,j] is the jth measurement of the projection no i
                //                   
  Scalaire *phi;  // angles (positions de la source)
  Scalaire *x;    // abscisses for // sampling
                 // angular position of the detector for Fan Beam sampling
  void Alloue();	// Alloue g(p.q), phi(p), x(p.q)
  public:
    //
    // Constructeurs et destructeurs
    // =============================
    //
      Sinogramme(int np,int nq);
        // Alloue un Sinogramme a np rotations et nq translations
	// for // sampling:
	// equidistant sampling on [-1 1] for x, iq=1,...,nq
	// equidistant sampling on [0 pi], for phi, ip=0,...,np-1
	// for FanBeam sampling see FBsinogramme.h
      Sinogramme(int np,int nq, Scalaire beginangle, Scalaire angle, 
		 Scalaire radius=1.); 
      // LD modif =1 22 nov 2011 // LD modif =1 17 dec 2013
      // Alloue un Sinogramme a np rotations et nq colonnes
      // for // sampling:
      // angular sampling on [beginangle angle] (default should be angle=M_pi)
      // projection sampling  [-radiu radius] (default should be radius FOV = 1
      // equidistant sampling on [-radius radius] for x, iq=1,...,nq
      // equidistant sampling on [beginangle angle], for phi, ip=0,...,np-1
      // for FanBeam sampling see FBsinogramme.h
      Sinogramme(int np,int nq, Scalaire beginangle, Scalaire angle, 
		 Scalaire *sshift, Scalaire radius=1.); 
      // the same as previous except that 
      // the sampling on [-radius radius] for x, iq=1,...,nq
      //  is  shifted by the vector sshift of length np
      Sinogramme(int np,int nq, Scalaire *vecangle,Scalaire radius=1.); 
      // angular sampling are given by vecangle (default should be angle=M_pi)
      //  the sampling on [-radius radius] for x, iq=1,...,nq
      Sinogramme(int np,int nq, Scalaire *vecangle, Scalaire *sshift,
		 Scalaire radius=1.); 
      // the same as previous except that 
      // the sampling on [-radius radius] for x, iq=1,...,nq
      //  is  shifted by the vector sshift of length np
      Sinogramme(int np,int nq, Scalaire *vecangle, 
		 Scalaire *begins, Scalaire *steps);
      // angular sampling are given by vecangle (default should be angle=M_pi)
      //  the sampling on each projection 
      //         x[ip,iq]=begins[ip]+(iq-1)steps[ip] for iq=1,...,nq
      // this allow for linogram sampling

      Sinogramme(int np,int nq,Scalaire angle,
		 char *filename,int nlig,int ncol);
       // Alloue un Sinogramme a np rotations et nq colonnes
       // read the value in the file filename (contain nlig*ncol floats)
       //             only the first p rows and q columns are used
       // equidistant sampling on [-1 1] for x, iq=1,...,nq
       // equidistant sampling on [0 angle], for phi, ip=0,...,np-1
       Sinogramme(const Sinogramme & sino);
        // constructor by copy
      ~Sinogramme();
        // free the Sinogramme
        //
	// Initialisations
	// ===============
	//
inline int addresse(int ip,int iq)
// Sinogramme addresse ta(i,j)=ip*q+iq (ip=0,...,p-1;iq=1,...,q)
{
#if CHECK
  if (ip<0 || iq<1 || ip>=p || iq>q)
    fputs("Erreur : Adressage hors Sinogramme\n",stderr);
#endif
  return (ip*q+iq-1);
}

//      void Zero();
// Initialise  g a 0, phi(i) a i*PI/p, 
// pour tout p x(p,q) equidistant sur ]-1 1[
 void Zero(Scalaire beginangle=0,Scalaire angle=M_PI,Scalaire radius=1.);
 // Initialise g a 0, phi(i) a beginangle+i*angle/p
 // pour tout p x(p,q) equidistant sur ]-radius radius[.
 void Zero(Scalaire *sshift, 
	   Scalaire beginangle=0,Scalaire angle=M_PI,Scalaire radius=1.);
 // the same as before but 
 // x(p,q) equidistant sur ]-radius radius[ + sshift(p) for all q
 void Zero(Scalaire *vecangle, Scalaire *sshift, Scalaire radius=1.);
 // Initialise g a 0, phi(i) a vecangle(i)
 // x(p,q) equidistant sur ]-radius radius[ + sshift(p) for all q
 void Un();
 // Initialise a g a 1, phi(i) a i*PI/p, 
 // pour tout p, x(p,q) equidistant sur ]-1 1[
     void Un(Scalaire angle);
       // Initialise a g a 1, phi(i) a i*PI/p, 
       // pour tout p, x(p,q) equidistant sur ]-1 1[

     void InitAVSfile(char *filename,int nlig, int ncol,
                                  Scalaire angle=M_PI);
      // Initialise  on [0,angle] with the projection read in filename 
      // containing nlig*ncol floats.
      //   only the first p rows and q columns are used
     void AddEllipse(Scalaire r, Scalaire theta, 
		     Scalaire psi,Scalaire a, Scalaire b,Scalaire density,
		     int oversampling=1);
       // Add the measurements of "density" times the indicator 
       // an allipsis in scalar tomography...
       // (r,theta) are the polar coordinates of the ellipsis
       // psi is the direction angle of the long axis
       // a and b are the repectively the length of the long 
       //             and short axis
       // density is the ellipsis density.
     // oversampling must be ODD (computing the mean of oversampling 
     //                           line integrals instead of just 1)

     void AddTranslatedEllipse(Scalaire r, Scalaire theta, Scalaire psi,
				  Scalaire a, Scalaire b,Scalaire density,
				  Scalaire* translation);
// Add the measurements of an allipsis in scalar Sinography...
  // AIJ are the four entries of the matrix A and TransScal is the 
  // scalar part of a translation such that 
  // the ellipsis indicator X is linearly transformed 
  // into \chi(x+b(\phi))
  // where $translation(\phi)$ is  $b(\phi) \cdot (\cos\phi,\sin\phi)$

     void AddMovingEllipse(Scalaire r, Scalaire theta, Scalaire psi,
			   Scalaire a, Scalaire b,Scalaire density,
			   Scalaire *A11,
			   Scalaire *A12,
			   Scalaire *A21,
			   Scalaire *A22,
			   Scalaire *TransScal
			   );
    void AddEllipse(Scalaire r,Scalaire theta,
		    Scalaire longaxis,Scalaire smallaxis,
		    Scalaire psi, Scalaire activity,
		    Scalaire mu);
      // Add the attenuated Ellipse activity to 
      // the Sinogramme.
      // (r,theta) : polar coordinats of the ellipsis
      // psi : angle of the long axis.
      // mu : attenuation parameter
    void AddTpNEllipse(Scalaire r, Scalaire theta, Scalaire psi,
		       Scalaire a, Scalaire b,Scalaire density,int expon);
// Add the measurements of an allipsis in scalar GENERALIZED tomography 
// with T as weight function, where t is the parameter of the 
// line path of integration g(\phiS)=\int_{t}f(S u + T u^\perp) T dT
// with u=(\cos \phi, \sin\phi) and u^\perp=(-\sin\phi , \cos \phi ).
    void AddAtteEllipse(Scalaire r, Scalaire theta, Scalaire psi,
			Scalaire a, Scalaire b,Scalaire density, 
			Scalaire mu);  
      // Add the measurements of "density" times the indicator 
      // an allipsis in Attenuated tomography... (ATTENUATED TRANSFORM)
      // (r,theta,psi,a,b,density see above....
      // mu is the constant attenuation coeeficient
    void MulAttHalfDisk(Scalaire mu);
      // mult by the exponential  attenuation on 
      // the unit disk of constant parameter mu
    void MulScalProj(Scalaire scal,int i);
      // mult the projection number i by scal
    void Xdef(Scalaire radius=1.);
      // Initialize the abscisse to the usual 
      //  one in analytic approach see Natterer 86
      //	 void Disque(Scalaire c1, Scalaire c2,
      //			    Scalaire rayon, Scalaire densite)
    void Xshift(Scalaire *shift);
    // shift of x vector

    void LitSino(char *filename);
      // Initialise le Sinogramme aux valeurs lues dans un fichier (read)
    void EcritSino(char *filename);
      // Copy the Sinogramme in a file  (write)

    /*
     * Acces structure
     * ===============
     */
   Scalaire & operator () (int i,int j);
      // Reference coefficient (i,j)
      // get t(i,j)
   const Scalaire & operator () (const int i,const int j) const;
      // Reference la valeur t(i,j)
      // get t(i,j)
    Scalaire & operator () (int i,Scalaire t);
// return the value of g(i,j) such that x[j] near t 
// only for equidistributed x
    void Xcopy(int i,Projection *u);
    void Xcopy(const Sinogramme *t);
    void Acopy(int i,Projection *u);
    void Acopy(const Sinogramme *t);
    void GetProj(int i,Projection *u);
      // Recopie la ligne i dans u
      // set u to t(i,.)
    void SetProj(int i,Projection *u);
      // Fixe la projection no i a partir de u
      // set t(i,.) to u
      //		void GetRot(int j,Projection *u);
      // Recopie la rotation no j dans u
      // set  u to t(.,j)
      //		void SetRot(int j,Projection *u);
      // Fixe la rotation no j a partir de projection u
      // Affect the rotation no j values to the projection u
    int NP() const;
      // Nbre de rotations
      // rotation number
    int NQ() const;
      // Nbre de translations
      // translation number
    const Scalaire A(const int i) const;
      // angle number i  for parallel projection or source parameter in Fan Beam
    void Angles(Scalaire *angle) ;
    // return all the angles in the Scalaire vector *angles (angle[i]=phi[i])
    void Detectors(Scalaire *detector) ;
    // return all the detector positions in the Scalaire vector *detector 
    // (detector[i*q+j]=x[i*q+j]
    const Scalaire X(const int i, const int j) const;
      // detector abscisse j for projection i (s_j for parallel projection,\alpha_j for FB) 
    // i=0,...,p-1 ; j=1,...,q
    // angle number i, translation j in parallel geometry and source position i and detector j in FB geometry
   void Examine();
      // represents the Sinogramme 
    void PGMWrite(char *filename);
      // Write the image in PGM ASCII format in the file filename.
    void PGMrbWrite(char *filename);
      // Write the image in PGM RAWBYTE format in the file filename.
    void PGMWriteLog(char *filename);
      // Write the image in PGM ASCII format in the file filename.
      // using a log scale on the image

    double MaxAbs();
    // return the MaxAbs of the sinogramme

    void DCC(double *dccval,int order); 
    // Compute the Moments (for each projection) of order "order" of the parallel projection (*this)
    // Compute de the Data Consitency Conditions of order order

//
// Transform
// =========
//
    
    void Deriv(); // derivation according to the scalar variable d/dq (g(p.q))
    //(*this)(i,j)=(*this(i,j+1)-*this(i,j-1))/2h
     void MulExpSqrtUnMinusSsquare(Scalaire mu);
	  //(*this)(i,j)=*this(i,j)*exp(mu*sqrt(1-x[i,j]*x[i,j]))
     void ProjPositive(Scalaire epsilon=0);
	  // *this=epsilon if *this<epsilon     
     void TransPositive();
	  // *this=*this-min (if min<0)
     void Complete(Sinogramme *a);
	  // complete the considered sinogram a
          //  from [0,pi[x[-1,1] to 
	  //  [0,2*pi[x[-1,1] with
          //  g(phi+pi,s)=g(phi,-s)
  
     //     void ModuleFourier();
     // SEE sinogrammeFourier
	  // compute the module of the 2D fourier transform
          // make the 2D DFT of *this, take the modulus of 
          // the result and store this modulus in *this
     
	};

//
// Operateurs
// ==========
//

void Copy(Sinogramme *a,Sinogramme *b);
// b <- a
void Copy(Projection *u,int i,Sinogramme *a);
// a(i,.) <- u
void Add(Sinogramme *a,Sinogramme *b,Sinogramme *c);
// c <- a+b
void Sub(Sinogramme *a,Sinogramme *b,Sinogramme *c);
// c <- a-b
void KAdd(Sinogramme *a,Scalaire k,Sinogramme *b,Sinogramme *c);
// c <- a+k.b
void Prod(Sinogramme *a,Sinogramme *b,Sinogramme *c);
// c(i,j) <- a(i,j).b(i,j)
void KMul(Scalaire k,Sinogramme *a,Sinogramme *b);
// b <- k.a
void Transpose(Sinogramme *a,Sinogramme *b);
// b <- Trans(a)
void Sin(Sinogramme *a,Sinogramme *b);
// b <- sin(a)
void Cos(Sinogramme *a,Sinogramme *b);
// b <- cos(a)
void Exp(Sinogramme *a,Sinogramme *b);
// b <- exp(a)


//
// Operateurs Projection-Sinogramme
// ==========================
//

void Convolution(Projection *u,Sinogramme *a,Sinogramme *b);
// b(i,.) <- u convol a(i,.) (for all i)
// The dimension of u must be 2*a->NQ()
// a and b contains projections whereas u is supposed to be the filter
void LocConvolution(Projection *u,Sinogramme *a,Sinogramme *b,
		    int halfwidth);
     // b(i,.) <- u convol a(i,.) (for all i)
     // The dimension of u must be 2*a->NQ()
     // do the convolution in i on indices between 
     // i-halfwidth and i+halfwidth
     // have been adapted for local and pseudo-local tomography


//void Filtre(Projection *u,Sinogramme *a,Sinogramme *b, bool real);
// SEE sinogrammeFourier
// b(i,.) <- Fourier-1(Fourier(E(a(i,.))).u) (for all i)
// the projection a(i,.) are extended to a double dimension, (zero padding)
// i.e. dim(E(a(i,.))) = 2 dim(a(i,.))=2*a->NQ()
// The dimension of u must be 2*a->NQ()
  // the BOOLEAN real indicates if the filer is real 
  // (usally even) in the Fourier domain
  // if real is false then the filter is supposed to be a pure odd imaginary 
  // filter such as the Hilbert transform  

void LocInterpol(Projection *u,Sinogramme *a, Sinogramme *b);
// b(i,.) <- Fourier-1(Fourier(a(i,.)).u) (for all i)

Scalaire InterlacedInterpol(Sinogramme *kernel,Sinogramme *sino,int pint, int qint,
			int verbose=0);
// sino(pint,qint) <- Sum_{ all i and all j possible in the corresponding INTERLACED grid}   }a(pint-i,qint-j)).kernel(i,j)) 
// Note that the sum is periodic in p (as Sinogramme but not in q)
// work only for Sinogramme defined on 2Pi
// the dimensions of *kernel must be ODD
// return the interpolated value

//
// Operateurs Sinogramme-Image
// ================================================
//

void Retroprojection(Sinogramme *t,RealImage *image, 
		     Scalaire AngularInterval=M_PI);
// image <- Retroprojection t
void Retroprojection(Sinogramme *t,RealImage *image,Scalaire mu, 
		     Scalaire AngularInterval=M_PI);
// image <- Attenuated Retroprojection t  (Dual of the exponential
//                                         Radon transform)
Scalaire Retroprojection(Sinogramme *t, Scalaire x, Scalaire y,
		     Scalaire AngularInterval=M_PI);
// return Retroprojection of *t at (x,y)
//  Scalaire x,y  :  abs and ord of the current pixel.
// see F.Natterer 
// "The mathematics of Computerized Tomography" (Wiley 86) p.47 and p.109
Scalaire RolfRetroprojection(Sinogramme *t, Scalaire xima, Scalaire yima,
			     Scalaire phimin,Scalaire phimax,
			     int n);
// return Retroprojection of *t at (xima,yima)
//  Scalaire xima,yima; // abs and ord of the current pixel.
// the back projection  angle interval is  phimin,phimax
// n is for the weight "tan(phi)^n/cos(phi)"
	
//
// tool box
//
int BadDim(Sinogramme *a,Sinogramme *b);
// ((a->NP != b->NP) || (a->NQ != b->NQ))
int BadDim2(Sinogramme *a,Projection *p);
// (2.a->NQ != b->N)
int BadDim2QmoinsUn(Sinogramme *a,Projection *p);
// (2.a->NQ-1 != b->N)

void Localize(Scalaire s,int ip,Sinogramme *t,int iq);
// result: iq such that t->X(ip,iq)<= s < t->X(ip,iq+1)

void Extend(Sinogramme *a,int i,int n,float *w);
// w(j-1)<- a(i,j-n+n/2) if n-n/2<j<2n-n/2   (w(0,...,2n) whereas a(i,1,...,n))
// w(j)<- 0.		           else
// n must be a->NQ()

void InvExtend(float *w,int i,int n,Sinogramme *a);
// a(i,j-n+n/2)<- w(j-1) if n-n/2<j<2n-n/2   (w(0,...,2n) whereas a(i,1,...,n))
// n must be a->NQ()

//
// interface subroutines 
//**********************
void SimulData(Sinogramme *p,Scalaire mu,
	       char *filename, int ecritsino, int verbose);
// simulation of a prajection stored in *p 
// mu is for exponential Radon transform (mu=0 for Radon transform)
// from a phantom defined in a file filename
// if ecrisino then save the data in a file with EcritSino
// if verbose then the soft may be verbose
void GetSimulData(Sinogramme *p,Scalaire mu=0);
// get ellipis parameter in a file
void GetSimulMovingData(Sinogramme *p, int fileformat=1);
void GetSimulDynTranslatedData(Sinogramme *p,Scalaire* translation);
// get ellipis parameter in a file and add translation for dynamic projections
// call SimulDynTranslatedData (see below)
// see AddTranslatedEllipse
void SimulDynTranslatedData(Sinogramme *p,Scalaire* translation, 
			    char *filename, char car, int verbose);
// see AddTranslatedEllipse
// filename is the name of the phantom definition file
// car is 'y' if the sinogramme is to be stored
void GetEllParam(Scalaire &er,Scalaire &etheta, Scalaire &epsi, 
		 Scalaire &ea, Scalaire &eb, Scalaire &edensity);
// to read ellipsis geometrical parameters


#endif
