// FanBeam sinogramme
//

#ifndef __FBSINOGRAMME_H
#define __FBSINOGRAMME_H

#include <sinogramme.h>
#include <math.h>

class FBSinogramme : public Sinogramme //
{
protected:
  Scalaire radius; // ratio: radius of the source circle/radius of the support
  Scalaire fanangle; // connic angle considered for the data
public:
  //
  // Constructors et destructors
  // =============================
  //
  FBSinogramme(int np,int nq,Scalaire sr,Scalaire fa=M_PI);
  // Allocation of Sinogramme of np rotations and nq translations
     // ratio = sr is the source radius (> 1 because the measured and reconstructed  
     //         region is the the disk of radius 1
     // equidistant sampling on [0 2*pi], for phi, ip=0,...,np-1
     //     the vertex path if sr*(cos(phi),sin(phi))
     // equidistant sampling on [-fa fa] for x, iq=1,...,nq 
     //      (fan angles between the mesured line and the line throught the source vertex and the origine)  
     // by default, the sources are on the whole circle (angle=2*M_PI)

  FBSinogramme(int np,int nq,Scalaire angle,Scalaire sr,Scalaire fa);
  // Alloue un Sinogramme a np rotations et nq colonnes
  // equidistant sampling on [0 angle], for phi, ip=0,...,np-1
  //     the vertex path if sr*(cos(phi),sin(phi))
  // equidistant sampling on [-fa/2 fa/2] for x, iq=1,...,nq (by def fa=pi);
  //      (fan angles between the mesured line and the line throught the source vertex and the origine)  

  //
  // structure initialisation
  //
  void FBXdef();
  // Initialize the abscisse,i.e., the angle alpha, to the usual one in 
  // analytic approach, see Natterer 86 fan-beam approach
  //
  //
  // structure access
  //
  Scalaire FA();
  // Fan angle (general case M_PI: practically less...)
  Scalaire R(); 
  // ratio of the radius of the source/radius of the support disk
  // by convention  the radius of the support disk is normalized to 1.
  //
  // transformation
  //
  void FBAddEllipse(Scalaire r, Scalaire theta, Scalaire psi,
		    Scalaire a, Scalaire b,Scalaire density,
		    int oversampling=1);
  void ConvolPorte(int porte);
// convolution par une porte de largeur porte
// porte MUST BE ODD

};

//
// Operators
// ==========
//
void Convolution(Projection *u,FBSinogramme *a,FBSinogramme *b);
// For each vertex position (i.e. each Fan Beam projection)
// each Fan Beam projection is convolved by the same vector u
//  (*b(i,l)) = \int (*a(i,l)) * (*u)(j-l+n)
// the dimension of (*u) must be 2n-1 ;  (*u)(n ) is "u(0)"  (*u)(l)=u( (l-n)*halpha)
//      with halpha=fangle/(n-1) 
void FBConvolution(Projection *u,FBSinogramme *a,FBSinogramme *b);
// b(i,.) <- u convol a(i,.) (for all i)
// The dimension of u must be 2*a->NQ()-1
// FBConvolution is for the FanBeam FBP (acosinus factor appears in the filter) 
void FBRetroprojection(FBSinogramme *t,RealImage *image);
// image <- Retroprojection t
// see F.Natterer 
// "The mathematics of Computerized Sinography" (Wiley 86) p.113
void FBLambdaRetroprojection(FBSinogramme *t,RealImage *image);
// image <- Retroprojection t
// in the case of Lambda Tomography
// see F.Natterer 
// "The mathematics of Computerized Sinography" (Wiley 86) p.113
// see Faridani 92 SIAM J. Appl. Math.

void ConvertIeq(FBSinogramme & div, Sinogramme *para);
// convert the Fan Beam sinogramme FBSinogramme *div 
// into *para a parallel beam Sinogramme
// the projections are supposed to be equi sampled


//
void GetSimulData(FBSinogramme *p);

#endif

