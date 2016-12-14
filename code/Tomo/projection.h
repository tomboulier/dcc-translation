// Algebre en dimension n
// Partie I : scalaires, angles et vecteurs
//
// LD DEC 93 :
// (c) Copyright TIMC 1993

#ifndef __PROJECTION_H
#define __PROJECTION_H

// #include <algnvect.h>
#include <cmath>
#include <tomotypes.h>

//
// ========
// Projection
// ========
//

class Projection
{
private:
  int n;		// Dimension
  Scalaire phi;   // angle de la projection
  Scalaire *x;	// Discrete Abscissae on the direction phi 
  Scalaire *y;	// y[i] is projection(x[i])  
  void Alloue();	// Alloue x et y de dimension n.
public:
  //
  // Constructeurs et destructeurs
  // =============================
  //
  Projection(int n=32);
  // Alloue un vecteur non initialise de dimension n
  ~Projection();
  // Libere le vecteur
  //
  // Initialisations
  // ===============
  //
  void Zero();
  // Recoit (0,0,...,0), equidistant x on [-1,1]
  void Un();
  // Recoit (1,1,...,1), equidistant x on [-1,1]
  void Gaussian(Scalaire sigma, Scalaire mean,
		Scalaire sbegin=-1, Scalaire send=1);
  // Recoit 1/(sigma * sqrt(2 M_PI)) * exp(-1/2 (x(i)-mean)^2/sigma^2), 
  void FourierSigneFilter();  
  // becomes the Fourier SIGNE filter in the fourier domain, i.e.,
  // HILBERT transform in the direct domain when used with the method Filtre
  void FourierIdealFilter(Scalaire b);
  // becomes the Fourier ideal filter in the fourier domain, i.e.,
  // |sigma| if |sigma| < b, 0 else.... See Natterer 86 p. 18
  void FourierExponentialIdealFilter(Scalaire b,Scalaire mu);
  // becomes the Fourier ideal filter in the fourier domain, i.e.,
  // |sigma| if abs(mu) < |sigma| < b, 0 else.... See Natterer 86 p. 18
  //

  void SheppLoganFilter(Scalaire b,Scalaire sbegin=-2, Scalaire send=2);
// becomes the shepp and Logan filter 
  // in the Direct domain see Natterer 86, p.111.
  void FBSheppLoganFilter(Scalaire b,Scalaire h);
  // becomes the shepp and Logan filter in the Direct domain for Fan Beam
  // see Natterer 86, p.111. and 113
  // b is the bandwith, h is the angle between 2 successive measurements  
  // lines for one fixed position of the source.
  // For FBsinnogramme use h=fanangle/(nbtrans-1)

  void LowPassHilbertFilter(int q, Scalaire h);
// compute the LowPassHilbertFilter projection of length 2n-1 for a parallel proj of length q
// of the form $(1-cos(b s)) / (\pi s)$ or equiv $(1-cos(2 \pi c s)) / (\pi s) $ with $b=2\pi c$
// with the cutoff frequency c=1/(2h) or $b=\pi/h$ !!!!
// at ih i=-n+1:n-1, LPHF(i)=(1-\cos(\pi i)) / (\pi i h) and is 0 at i=0.
// and the indices are translated for starting at 0 (+n-1)

  void LowPassHilbertFilter(Scalaire b, 
			    Scalaire sbegin=-2., Scalaire send=2.);
  // build low pass hilbert  filter in the Direct domain
  // see Natterer 2001 inverse problems 17,pp. 113-119
  // b is the cut off bandwith  (with the Natterer FT def... WARNING)
  // return        x[i]=sbegin + i*h;      bxi=b*x[i];
  //                  y[i]= ((1- cos( b*xi  ))/ x[i])/M_PI;
  //      instead of ( (1-cos(2M_PI c x[i]) / x[i] ) / M_PI  for the classical FT (with 2\pi factor)
  // in fact, is the same with 2M_PI c=b
  // sbegin and send define the interval (support) of the filter
  // warning : sbegin and send are sampled 

  void FBHilbertBandlimitedFilter(Scalaire b,
				  Scalaire alphabegin, Scalaire alphaend);
  // computes the Fan Beam Hilbert filter truncated in the Fourier domain
// h_H is such that FT(h_H)(\sigma)=-i \sign(\sigma)\chi_{[-b;b]}(\sigma)
// thus h_H(s)=1/(\pi s)  (1-\cos(2\pi b s)) = \frac{2}{\pi s} \sin^2(\pi b s)
  // see  Clacdoyle and Defrise paper  (IEEE SIGNAL PROCESSING MAGAZINE [66] JULY 2010)  
  // for the ROI reconstruction application
  // b is the bandwith, [alphabegin, alphaend] is the covered angular fan 
  //alphabegin and alphaend are the first and last sampled angle of 
  //the equiangular projection of angular step halpha=(alphaend-alphabegin)/(Scalaire)(n-1) 
  // (the angle between 2 successive projection lines for one fixed position of the source).
  // This is for equiangular FB projections
  // For FBsinnogramme use h=fanangle/(nbtrans-1)

  void FBLambdaFilter(Scalaire b,Scalaire h,Scalaire m);
  // becomes the  filter in the Direct domain for Lamda Tomography 
  // in the case of Fan Beam Geometry.
  // see Faridani 92, SIAM J. Math Appl., Vol 52, no 2.
  // see also Natterer 86, p.111. and 113
  // b is the bandwith, h is the angle between 2 successive measurement  
  // lines for one fixed position of the source.
  // m is a regularity parameter of the filter (see Faridani).
  //
  // REMARK "r" of Faridani is replace by "b" of Natterer 
  // (detail size replaced by band limit)
  //
  void PseudolocFilter(double m,double rho,double eps,int ninteg);
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
  void SheppLoganFilter(Scalaire b,Scalaire mu);
  // we try to generalize the shepp and Logan filter 
  // in the Direct domain, for the exponential Radon Transform
  // see Natterer 86, p.111., p. 46 and my own research....
  // b is the bandwith, mu is the attenuation coeff.
  void LinearapodFilter(Scalaire b, Scalaire epsilon,Scalaire mu=0);
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
  //
  //		void LitP();
  // lit P dans un fichier

  //  void Lit(char* filename);
  void Ecrit(char* filename);
  // store a projection in a file

  //
  // Acces structure
  // ===============
  //
  Scalaire & operator () (int i);
  // Acces a la ieme composante de la projection.
  Scalaire X (int i);
  // ieme composante de l'abscisse.
  Scalaire Y (int i);
  // ieme composante de la projection.
  void Examine();
  // examen du contenu de P.
  int N();
  // Taille.
  Scalaire A();
  // Angle
  void Xcopy(Projection *a);
  // copy de l'abscisse
  void Xcopy(Scalaire *a);
                        // copy de l'abscisse
  void Acopy(Projection *a);
  // copy de l'angle de la projection *a
  void Acopy(Scalaire s);
  // copy de l'angle s
};
//
// Operateurs
// ==========
//

void Copy(Projection *a,Projection *b);
// b <- a
void Mul(Scalaire k,Projection *u);
// u <- k.u
void Mul(Scalaire k,Projection *u,Projection *v);
// v <- k.u
void Add(Projection *u,Projection *v,Projection *w);
// w <- u+v
void KAdd(Projection *u,Scalaire k,Projection *v,Projection *w);
// w <- u+k.v
void Prod(Projection *u,Projection *v,Projection *w);
// w <- u.v
void Convol(Projection *u,Projection *v,Projection *w);
// w <- u*v
// circular convolution of vector u and v
Scalaire Scal(Projection *u,Projection *v);
// <- <u,v>
void Normalise(Projection *u,Projection *v);
// v <- u/Norme(u)
Scalaire Norme1(Projection *u);
// <- som|u_i|
Scalaire Norme2(Projection *u);
// <- <u,u>
Scalaire Norme(Projection *u);
// <- Norme(u)


#endif
