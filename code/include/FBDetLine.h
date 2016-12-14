// General FanBeam sinogramme
//
// Operation sur les Fan Beam sinogrammes 
//
// Created : LD Dec 2013
// last modification: 
// (c) Copyright TIMC 2013
 
#ifndef __FBDETLINE_H
#define __FBDETLINE_H

#include <cmath>
//
#include <tomotypes.h>
#include <sinogramme.h>
#include <realimage.h>

// #########################################

class FBDetLine //
{
 protected:
  int p,q; // p is the projection number = number of source position 
  // q is the fan number per projection = number of fans (must be >1);
  Scalaire *xsource; // abscisses of the source positions
  Scalaire *ysource; // ordonnées of the source positions
  Scalaire *xdetector; // abscisses of the original detector positions
  Scalaire *ydetector; // ordonnées of the original detector positions
  Scalaire *xvector; // (xvector, yvector) vector directions of the detector line
  Scalaire *yvector; // 
  // equi-sampled detector : detector_{i,j}= detector_i + j*vector_i, j=0,...,q-1

  Scalaire *data; // data(q*ip+j) is the j^th fan integral in the ip^th projection
  // i.e. on the line joining source_ip and detector_{ip,j}
  void Alloue();    // Alloue xsource(p), ysource(p), 
                    // xdetector(p),ydetector(p)
                    // xvector(p),yvector(p)
                    // data(p.q) 
 public:
  //
  // Constructors et destructors
  // =============================
  //
  FBDetLine(int np,int nq, // np source positions of nq detectors per source position 
		Scalaire *xs,Scalaire *ys, // np source positions
		Scalaire *xd,Scalaire *yd, // np detector positions
		Scalaire *xv,Scalaire *yv // np detector positions
		);
  // Allocation of FBDetLine of np source positions each of nq detectors on a line 
  // equidistant sampling detector_{ip,jq}= detector_ip + jq*vector_ip, jq=0,...,q-1
     //      (fan angles between the mesured line and the line throught the source vertex and the origine)  

  FBDetLine(const FBDetLine & fbdl);
  // constructor by copy
   ~FBDetLine();
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
  // Nb of source positions (return p)
  int NQ() const;
  // Nb of detectors per projection (return q)
  // translation number
  Scalaire Xsource(const int i) const;
  //return abscisse of the i^{th} source position
  Scalaire Ysource(const int i) const;
  //return ordonnee of the i^{th} source position
  Scalaire Xdetector(const int i) const;
  //return abscisse of the i^{th} ORIGINE detector position
  Scalaire Ydetector(const int i) const;
  //return ordonnee of the i^{th} ORIGINE detector position
  Scalaire Xvector(const int i) const;
  //return abscisse of the i^{th} vector direction 
  Scalaire Yvector(const int i) const;
  //return ordonnee of the i^{th} vector direction 

  void EcritSino(char *filename);
  // Copy the FBDetLine data in a file
  void LitSino(char *filename);
  // read the FBDetLine data in a file
  void PGMWrite(char *filename);
  // Save the FBDetLine in a PGM file (raw byte format)
  void Examine();
  // cout all the content of a GFBSinogramme (except the data)

  //
  // transformation
  //
  void AddEllipse(Scalaire r, Scalaire theta, Scalaire psi,
		  Scalaire a, Scalaire b,Scalaire density);
		  //		    ,int oversampling=1);
//
// Operators
// ==========
//

void Add(FBDetLine  & a,FBDetLine & b);
// (*this) <- a+b
void Sub(FBDetLine  & a,FBDetLine & b);
// (*this) <- a-b
void GetSimulData();
// simulation of data (asking questions to users) Scalaire mu, pourrait être prévu
void SimulData(char *filename, int ecritsino, int verbose);// Scalaire mu, pourrait être prévu
// silent simulation of data from the phantom definition in file "filename"
// if ecritsino == 1 if the FBDetLine must be saved (no more silent)
// if verbose ==1 then some trace of what is done  
//void GetSimulDynTranslatedData(FBDetLine *p,Scalaire *Xtranslation);
// Warning here translation contains the X translation of the center 
// dor each source position 
//void SimulDynTranslatedData(FBDetLine *p,Scalaire* Xtranslation, 
//			    char *filename, char car, int verbose);
// see AddTranslatedEllipse
// Warning here translation contains the X translation of the center 
// dor each source position 
};


#endif

