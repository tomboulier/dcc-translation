// General FanBeam sinogramme
//
// Operation sur les Fan Beam sinogrammes 
//
// Created : LD Dec 2013
// last modification: 
// (c) Copyright TIMC 2013
 
#ifndef __GFBSINOGRAMME_H
#define __GFBSINOGRAMME_H

#include <cmath>
//
#include <tomotypes.h>
#include <sinogramme.h>
#include <realimage.h>

// #########################################

class GFBSinogramme //
{
 protected:
  int p,q; // p is the projection number = number of source position 
  // q is the fan number per projection = number of fans (must be >1);
  Scalaire *xsource; // abscisses of the source positions
  Scalaire *ysource; // ordonnées of the source positions
  Scalaire *fananglewidths; // fan angle (angle of view from the source
  Scalaire *fananglebegins; // the projection direction are supposed to be 
  // equi-sampled 
  // angle_{ip,j)=fananglebegin(ip)+j*hangle_ip; 
  // with hfanangle_ip= fanangle(ip)/(q-1); (q>=2)
  // the fan direction is then \vzeta(phi)=(-sin(phi),cos(phi))
  //     with phi=fananglebegins+M_PI/2+circleangle;
  //        where circleangle=acos(xsource[i]/sqrt(xsource[i]*xsource[i]+ysource[i]*ysource[i])); if(ysource[i]<0) circleangle=2*M_PI- circleangle;
  Scalaire *data; // data(q*ip+j) is the j^th fan integral in the ip^th projection
  void Alloue();    // Alloue xsource(p), ysource(p), 
                    // fananglewidths(p),fananglebegins(p),
                    // data(p.q) 
 public:
  //
  // Constructors et destructors
  // =============================
  //
  GFBSinogramme(int np,int nq, // np sources position of nq fans
		Scalaire *xs,Scalaire *ys, // np source positions
		Scalaire *FAIwidths, Scalaire *FAbegs // fan angles widths; fan angle begins
		);
  // Allocation of GFBSinogramme of np source positions each of nq fans
  // equidistant sampling on [fa fa] for x, iq=1,...,nq 
     //      (fan angles between the mesured line and the line throught the source vertex and the origine)  

  GFBSinogramme(int np,int nq, Scalaire circulartrajradius, 
		Scalaire circleanglewidth=2*M_PI,Scalaire circleanglebegin=0.,
		Scalaire FOVradius=1.);
  // by default a circular trajectory of the source on a circle 
  // of radius  circulartrajradius
  // from circleanglebegin to circleanglebegin+circleanglewidth-hangle
  // FOVradius is the radius of the Field Of View
  // the Fan Angle is constant FAconstwidth= 2*asin(FOVradius/circulartrajradius), 
  // sampled around the origine symetricaly from -FAconstwidth/2 to FAconstwidth/2
  // fananglebegins[i]=circleangle+M_PI/2-FAconstwidth/2 

  GFBSinogramme(const GFBSinogramme & gfbsino);
  // constructor by copy
   ~GFBSinogramme();
   // free the GFBSinogramme

  //
  // structure initialisation
  //
  void Zero();
// set data[i*q+j] i=0...p-1, j=0...q-1, to 0

  //
  // structure access
  //
  Scalaire & operator () (int i,int j);
      // Reference coefficient data(i,j)
      // get data(i,j)
  const Scalaire & operator () (const int i,const int j) const;
      // Reference la valeur data(i,j)
      // get data(i,j)
  int NP() const;
  // Nbre de rotations
  // rotation number
  int NQ() const;
  // Nbre de translations
  // translation number
  Scalaire Xsource(const int i) const;
  //return abscisse of the i^{th} source positions
  Scalaire Ysource(const int i) const;
  //return ordonnee of the i^{th} source positions
  Scalaire FAbegin(const int i) const;
  //return the i^{th} fan angle begin = fananglebegins[i]
  Scalaire FAwidth(const int i) const;
  //return the i^{th} fan angle width = fananglewidths[i]

  void EcritSino(char *filename);
  // Copy the GFBSinogramme data in a file
  void LitSino(char *filename);
  // read the GFBSinogramme data in a file
  void PGMWrite(char *filename);
  // Save the GFBSinogramme in a PGM file (raw byte format)
  void Examine();
  // cout all the content of a GFBSinogramme (except the data)

  //
  // transformation
  //
  void AddEllipse(Scalaire r, Scalaire theta, Scalaire psi,
		  Scalaire a, Scalaire b,Scalaire density);
		  //		    ,int oversampling=1);
  void AddTranslatedEllipse(Scalaire r, Scalaire theta, Scalaire psi,
			    Scalaire a, Scalaire b,Scalaire density,
			    Scalaire* translation);
  //WARNING here translation contain the Xtranslation of the phantom
  Scalaire RolfRetroprojection(Scalaire xima, Scalaire yima,
			       Scalaire FOVradius, int n,
			       int UP=1);
  // return Retroprojection of *this at (xima,yima)
  // The GFBSinogramme *this is supposed to be a FBacquisition 
  //    with a given centered FOV (of radius FOVradius)
  // n is for the weight "tan(phi)^n/cos(phi)"
  // // NOT IMPLEMENTED yet
  // the back projection  angle interval is  anglemin,anglemax
  //			       Scalaire anglemin,Scalaire anglemax,
  // UP==1 means that we compute the RolfBackProj for [muR, muL] (default)
  //    else UP==0 we compute for [muR, 2\pi+muL] with the sign (-1)
  //       or else UP==2 for [0, 2\pi[ traj : 
  //                       for [muR, muL] sign +
  //                       for [muR, 2\pi+muL (mod 2\pi)] sign -
  void RolfRetroprojection(RealImage *image,Scalaire FOVradius,
			   int n);
  // same as the previous for all pixel points (xima,yima) of the image
  // if n<0 simple backprojection else Rolf weighted backprojection (for DCC)
  // // NOT IMPLEMENTED yet
  //		Scalaire anglemin,Scalaire anglemax,

};

//
// Operators
// ==========
//

void Add(GFBSinogramme *a,GFBSinogramme *b,GFBSinogramme *c);
// c <- a+b
void Sub(GFBSinogramme *a,GFBSinogramme *b,GFBSinogramme *c);
// c <- a-b
// std::cout<< " WARNING :Add and Sub(GFBSinogramme *a,GFBSinogramme *b,GFBSinogramme *c) must BE TERMINATED !!! " << std::endl; 
//
void GetSimulData(GFBSinogramme *gfbsino);
// simulation of data (asking questions to users)
void SimulData(GFBSinogramme *gfbsino,// Scalaire mu, pourrait être prévu
	       char *filename, int ecritsino, int verbose);
// silent simulation of data from the phantom definition in file "filename"
// if ecritsino == 1 if the GFBSinogramme must be saved (no more silent)
// if verbose ==1 then some trace of what is done  
void GetSimulDynTranslatedData(GFBSinogramme *p,Scalaire *Xtranslation);
// Warning here translation contains the X translation of the center 
// dor each source position 
void SimulDynTranslatedData(GFBSinogramme *p,Scalaire* Xtranslation, 
			    char *filename, char car, int verbose);
// see AddTranslatedEllipse
// Warning here translation contains the X translation of the center 
// dor each source position 

#endif

