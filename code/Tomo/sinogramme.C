// Operation sur les sinogrammes 
// en particulier retroprojection
//
// Created : LD Dec 93
// last modification: LD Nov 95
// (c) Copyright TIMC 1993

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <iostream>
#include <fstream>

#include <sinogramme.h>
#include <projection.h>
#include <realimage.h>

// #include <fourier.h>

using namespace std;


void Sinogramme::Alloue()
     // Alloue g,phi,x (prive)
{
  g=new Scalaire [p*q];
  phi=new Scalaire [p];
  x=new Scalaire [p*q];
}

//
// Constructeurs et destructeurs
// =============================
//

Sinogramme::Sinogramme(int np,int nq)
// Alloue un Sinogramme a np rotations et nq colonnes
// equidistant sampling on [-1 1] for x, iq=1,...,nq
// equidistant sampling on [0 pi], for phi, ip=0,...,np-1
{
  p=np;q=nq;
  Alloue();
  Zero();
}

Sinogramme::Sinogramme(int np,int nq,
		       Scalaire beginangle, Scalaire angle, Scalaire radius)
      // Alloue un Sinogramme a np rotations et nq colonnes
      // for // sampling:
      // angular sampling on [beginangle angle] (default should be angle=M_pi)
      // projection sampling  [-radiu radius] (default should be radius FOV = 1
      // equidistant sampling on [-radius radius] for x, iq=1,...,nq
      // equidistant sampling on [beginangle angle], for phi, ip=0,...,np-1
      // for FanBeam sampling see FBsinogramme.h
{
  p=np;q=nq;
  Alloue();
  Zero(beginangle,angle,radius);
}
Sinogramme::Sinogramme(int np,int nq,
		       Scalaire beginangle, Scalaire angle, 
		       Scalaire *sshift, Scalaire radius)
// the same as previous except that the sampling on [-radius radius] 
//for x, iq=1,...,nq   is  shifted by the vector sshift of length np
{
  p=np;q=nq;
  Alloue();
  Zero(sshift,beginangle,angle,radius);
}

Sinogramme::Sinogramme(int np,int nq, Scalaire *vecangle,Scalaire radius)
// angular sampling are given by vecangle (default should be angle=M_pi)
//  the sampling on [-radius radius] for x, iq=1,...,nq
{
  p=np;q=nq;
  Alloue();
  Xdef(radius);
  for (int i=0;i<p;i++){
    phi[i]=vecangle[i];
    for (int j=1;j<=q;j++)
      g[addresse(i,j)]=0.0;
  }
}
Sinogramme::Sinogramme(int np,int nq, Scalaire *vecangle, 
		       Scalaire *sshift,Scalaire radius) 
// the same as previous except that 
// the sampling on [-radius radius] for x, iq=1,...,nq
//  is  shifted by the vector sshift of length np
{
  p=np;q=nq;
  Alloue();
  Zero(vecangle,sshift,radius);
}
Sinogramme::Sinogramme(int np,int nq, Scalaire *vecangle, 
		       Scalaire *begins, Scalaire *steps)
// angular sampling are given by vecangle (default should be angle=M_pi)
//  the sampling on each projection 
//         x[ip,iq]=begins[ip]+(iq-1)steps[ip] for iq=1,...,nq
// this allow for linogram sampling
{
  p=np;q=nq;
  Alloue();
  Scalaire s,hs;
  for (int i=0;i<p;i++){
    phi[i]=vecangle[i];
    s=begins[i];hs=steps[i];
    for (int j=1;j<=q;j++){
      x[addresse(i,j)]=s;s+=hs;
      g[addresse(i,j)]=0.0;
    }
  }
}

Sinogramme::Sinogramme(int np,int nq,Scalaire angle,
		       char *filename,int nlig,int ncol)
     // Alloue un Sinogramme a np rotations et nq colonnes
     // read the value in the file filename (contain nlig*ncol floats)
     //             only the first p rows and q columns are used
     // equidistant sampling on [-1 1] for x, iq=1,...,nq
     // equidistant sampling on [0 angle], for phi, ip=0,...,np-1
{
  p=np;q=nq;
  Alloue();
  InitAVSfile(filename,nlig,ncol,angle);
}

Sinogramme::Sinogramme(const Sinogramme & sino)
     // constructor by copy
{
  int i,j;
  p=sino.NP();
  q=sino.NQ();
  Alloue();
  for (i=0;i<p;i++)
    {
      this->Xcopy(&sino);
      this->Acopy(&sino);
      for (j=1;j<=q;j++)
	  g[addresse(i,j)]=sino(i,j);
    }
}

Sinogramme::~Sinogramme()
     // Libere la Sinogramme
{
  delete [] g;
  delete [] phi;
  delete [] x;
}

//
// Initialisations
// ===============
//

// void Sinogramme::Zero()
// // Initialise a 0
// 	{
// 	int i,j;
// 	Xdef();
// 	for (i=0;i<p;i++)
// 	  {
// 	  phi[i]=(double)i*M_PI/(double)p;
// 	  for (j=1;j<=q;j++)
// 	    {
// 		g[addresse(i,j)]=0.0;
// 	    }
// 	  }
// 	}

void Sinogramme::Zero(Scalaire beginangle,Scalaire angle,Scalaire radius)
// Initialise a 0 on [beginangle,beginangle+angle]
{
  int i,j;
  Xdef(radius);
  Scalaire dphi=angle/(double)p;
  Scalaire anglephi=beginangle;
  for (i=0;i<p;i++,anglephi += dphi){
    phi[i]=anglephi;
    for (j=1;j<=q;j++)
      g[addresse(i,j)]=0.0;
  }
}
void Sinogramme::Zero(Scalaire *sshift, Scalaire beginangle,Scalaire angle,
		      Scalaire radius)
// Initialise a 0 on [beginangle,beginangle+angle]
// introduce shift in s
{
  Zero(beginangle,angle,radius);
  for (int i=0;i<p;i++)
    for (int j=1;j<=q;j++)
      x[addresse(i,j)]+=sshift[i];
}

void Sinogramme::Zero(Scalaire *vecangle, Scalaire *sshift, Scalaire radius)
 // Initialise g a 0, phi(i) a vecangle(i)
 // x(p,q) equidistant sur ]-radius radius[ + sshift(p) for all q
{
  Xdef(radius);
  for (int i=0;i<p;i++){
    phi[i]=vecangle[i];
    for (int j=1;j<=q;j++){
      x[addresse(i,j)]+=sshift[i];
      g[addresse(i,j)]=0.0;
    }
  }
}

void Sinogramme::Un()
// Initialise a 1
	{
	int i,j;
	Xdef();
	for (i=0;i<p;i++)
	  {
	  phi[i]=(double)i*M_PI/(double)p;
	  for (j=1;j<=q;j++)
	    {
		g[addresse(i,j)]=1.0;
	    }
	  }
	}

void Sinogramme::Un(Scalaire angle)
// Initialise a 0 on [0,angle]
	{
	int i,j;
	Xdef();
	for (i=0;i<p;i++)
	  {
	  phi[i]=(double)i*angle/(double)p;
	  for (j=1;j<=q;j++)
	    {
		g[addresse(i,j)]=1.0;
	    }
	  }
	}


void Sinogramme::InitAVSfile(char *filename,int nlig,int ncol,Scalaire angle)
// Initialise  on [0,angle] with the projection read in filename 
// containing nlig*ncol floats.
//             only the first p rows and q columns are used
{
  FILE *pfile;
  int i,j;
  float *uc;

  Xdef();

  if (nlig<p)
    {
      cout << "ERR InitAVSfile: to few lines in your file " <<filename<<endl;
      exit(1);
    }
  if (ncol<q)
    {
      cout << "ERR InitAVSfile: to few columns in your file " <<filename<<endl;
      exit(1);
    }
/* ouverture du fichier en lecture */
  if ((pfile=fopen(filename,"r"))==NULL)
    {
      printf("ERR: (InitAVSfile) fichier %s non trouve \n",filename);
      exit(1);
    }
  uc = new float[nlig*ncol];
/* lecture des lignes suivantes */
  int retour=fread(uc,sizeof(float),nlig*ncol,pfile);
  fclose(pfile);

  /* only the first q column are taken into account!!! */
  for (i=0;i<p;i++)
    {
      phi[i]=(double)i*angle/(double)p;
      for (j=1;j<=q;j++)
	{
	  g[addresse(i,j)]=uc[i*ncol+j-1];
	}
    }
  delete [] uc;
}

void Sinogramme::AddEllipse(Scalaire r, Scalaire theta, Scalaire psi,
			    Scalaire a, Scalaire b,Scalaire density, 
			    int oversampling)
// Add the measurements of an allipsis in scalar Sinography...
{
  int i,j;
  Scalaire angle,deca_s,decalage,c,d,e;
  Scalaire asquare,bsquare;
  Scalaire co,si,cosquare,sisquare;
  Scalaire tone,ttwo;

  Scalaire dh=(x[addresse(0,q)]-x[addresse(0,1)])/(Scalaire)(q-1);
  //  Scalaire dhover=((Scalaire)2./(Scalaire)(q-1))/((Scalaire)oversampling);
  Scalaire dhover=dh/((Scalaire)oversampling);
  Scalaire somme;
  Scalaire *weight;
  Scalaire totalweight;
  weight=new Scalaire[oversampling];
  for(int k=0;k<oversampling;k++)
    weight[k]=1.;
  totalweight=0.;
  for(int k=0;k<oversampling;k++)
    totalweight+=weight[k];
  
  asquare=a*a;
  bsquare=b*b;
  for (i=0;i<p;i++) {
    angle=psi-phi[i];
    decalage=r*cos(theta-phi[i]);
    co=cos(angle);
    si=sin(angle);
    cosquare=co*co;
    sisquare=si*si;
    c=asquare*cosquare+bsquare*sisquare;
    for (j=1;j<=q;j++) {
      deca_s=x[addresse(i,j)]-dhover*(oversampling/2)-decalage; 
      //deca_s=x[addresse(i,j)]-decalage; 
      // if oversampling =  1   then  oversampling/2 = 0 integer !!!
      somme=0.;
      for(int k=0;k<oversampling;k++)
	if((deca_s*deca_s)<c)  {
	  d=a*b*sqrt(c-(deca_s*deca_s));
	  e=deca_s*si*co*(bsquare-asquare);
	  tone=(d+e)/c;
	  ttwo=(-d+e)/c;
	  somme+=density*(tone-ttwo)*weight[k];
	  deca_s+=dhover;
	} 
      g[addresse(i,j)]+=somme/totalweight;
	// we know that e has not to be computed bute we want to keep
	// tone and ttwo for the attenuated transformation....
    }
  }
  delete []  weight;
}

void Sinogramme::AddTranslatedEllipse(Scalaire r, Scalaire theta, Scalaire psi,
				  Scalaire a, Scalaire b,Scalaire density,
				  Scalaire* translation)
// Add the measurements of an allipsis in scalar Sinography...
// where $translation(\phi)$ is  $b(\phi) \cdot (\cos\phi,\sin\phi)$
{
  int i,j;
  Scalaire angle,deca_s,decalage,c,d,e;
  Scalaire asquare,bsquare;
  Scalaire co,si,cosquare,sisquare;
  Scalaire tone,ttwo;

  Scalaire dh=(x[addresse(0,q)]-x[addresse(0,1)])/(Scalaire)(q-1);
  
  asquare=a*a;
  bsquare=b*b;
  for (i=0;i<p;i++) {
    angle=psi-phi[i];
    decalage=r*cos(theta-phi[i])-translation[i];
    co=cos(angle);
    si=sin(angle);
    cosquare=co*co;
    sisquare=si*si;
    c=asquare*cosquare+bsquare*sisquare;
    for (j=1;j<=q;j++) {
      deca_s=x[addresse(i,j)]-decalage; 
      if((deca_s*deca_s)<c)  {
	d=a*b*sqrt(c-(deca_s*deca_s));
	e=deca_s*si*co*(bsquare-asquare);
	tone=(d+e)/c;
	ttwo=(-d+e)/c;
	g[addresse(i,j)]+=density*(tone-ttwo);
      } 
      // we know that e has not to be computed bute we want to keep
      // tone and ttwo for the attenuated transformation....
    }
  }
}


void Sinogramme::AddMovingEllipse(Scalaire r, Scalaire theta, Scalaire psi,
				  Scalaire a, Scalaire b,Scalaire density,
				  Scalaire* A11,
				  Scalaire* A12,
				  Scalaire* A21,
				  Scalaire* A22,
				  Scalaire* TransScal
 )
// Add the measurements of an allipsis in scalar Sinography...
  // AIJ are the four entries of the matrix A and TransScal is the 
  // scalar part of a translation such that 
  // the ellipsis indicator X is linearly transformed 
  // into \chi(A(\phi)x+b(\phi))
  // where $TransScal(\phi)$ is 
  // $b(\phi) \cdot A^{-t}(\cos\phi,\sin\phi)^t$
{
  int i,j;
  Scalaire angle,deca_s,decalage,c,d,e;
  Scalaire asquare,bsquare;
  Scalaire co,si,cosquare,sisquare;
  Scalaire tone,ttwo;

  Scalaire newphi; // new angles after transformation by A(\phi)
  Scalaire cosphi, sinphi; //
  Scalaire epsilon=0.00001;
  Scalaire detAphi;
  Scalaire invdetAphi;
  Scalaire AmtThetaUn,AmtThetaDeux;
  Scalaire NormAmtTheta, InvNormAmtTheta;
  Scalaire news;

  asquare=a*a;
  bsquare=b*b;
  for (i=0;i<p;i++) {
    cosphi=cos(phi[i]);
    sinphi=sin(phi[i]);
    detAphi=A11[i]*A22[i]-A12[i]*A21[i];
    if(detAphi<epsilon){
      cout << " error in Sinogramme::AddMovingEllipse" << endl;
      cout << " A[ "<<i<<"] is almost singular: det A="<< detAphi<< endl;
      exit(0);
    }
    invdetAphi=1/detAphi;
    AmtThetaUn=(A22[i]*cosphi-A21[i]*sinphi)*invdetAphi;
    AmtThetaDeux=(-A12[i]*cosphi+A22[i]*sinphi)*invdetAphi;
    NormAmtTheta=sqrt(AmtThetaUn*AmtThetaUn+AmtThetaDeux*AmtThetaDeux);
    InvNormAmtTheta=1/NormAmtTheta;
    AmtThetaUn=AmtThetaUn*InvNormAmtTheta;
    AmtThetaDeux=AmtThetaDeux*InvNormAmtTheta;
    newphi=atan2(AmtThetaDeux,AmtThetaUn);
    
    angle=psi-newphi;
    decalage=r*cos(theta-newphi);
    co=cos(angle);
    si=sin(angle);
    cosquare=co*co;
    sisquare=si*si;
    for (j=1;j<=q;j++) {
      news=(x[addresse(i,j)]+TransScal[i])*InvNormAmtTheta;
      deca_s=news-decalage;
      c=asquare*cosquare+bsquare*sisquare;
      if((deca_s*deca_s)<c)  {
	d=a*b*sqrt(c-(deca_s*deca_s));
	e=deca_s*si*co*(bsquare-asquare);
	tone=(d+e)/c;
	ttwo=(-d+e)/c;
	g[addresse(i,j)]=g[addresse(i,j)]+density*(tone-ttwo); 
	// we know that e has not to be computed bute we want to keep
	// tone and ttwo for the attenuated transformation....
      }
    }
  }
}


void Sinogramme::AddTpNEllipse(Scalaire r, Scalaire theta, Scalaire psi,
			       Scalaire a, Scalaire b,Scalaire density,int expon)
// Add the measurements of an allipsis in scalar GENERALIZED Sinography with T as weight function...
{
  int i,j,k;
  Scalaire angle,deca_s,decalages,decalaget,c,d,e;
  Scalaire asquare,bsquare;
  Scalaire co,si,cosquare,sisquare;
  Scalaire tone,ttwo;
  double dexponp1;
  
  dexponp1=(double)(expon+1);
  asquare=a*a;
  bsquare=b*b;
  for (i=0;i<p;i++){
    angle=psi-phi[i];
    decalages=r*cos(theta-phi[i]);
    decalaget=r*sin(theta-phi[i]);
    co=cos(angle);
    si=sin(angle);
    cosquare=co*co;
    sisquare=si*si;
    for (j=1;j<=q;j++){
      deca_s=x[addresse(i,j)]-decalages;
      c=asquare*cosquare+bsquare*sisquare;
      if((deca_s*deca_s)<c){
	d=a*b*sqrt(c-(deca_s*deca_s));
	e=deca_s*si*co*(bsquare-asquare);
	tone=decalaget+(d+e)/c;
	ttwo=decalaget+(-d+e)/c;
	for(k=1;k<=q;k++)
	  g[addresse(i,j)]=g[addresse(i,j)]+density*(pow(tone,dexponp1)-pow(ttwo,dexponp1))/dexponp1;  
	// g[addresse(i,j)]=g[addresse(i,j)]+density*(tone-ttwo)/q*(ttwo+k*((tone-ttwo)/q)); 
// we know that e has not to be computed bute we want to keep
// tone and ttwo for the attenuated transformation....
      }
    }
  }
}

void Sinogramme::AddAtteEllipse(Scalaire r, Scalaire theta, Scalaire psi,
Scalaire a, Scalaire b,Scalaire density, Scalaire mu)
// Add the measurements of "density" times the indicator 
// an allipsis in Attenuated Sinography... (ATTENUATED TRANSFORM)
// (r,theta) are the polar coordinates of the ellipsis
// psi is the direction angle of the long axis
// a and b are the repectively the length of the long 
//             and short axis
// density is the ellipsis density.
// mu is the constant attenuation coefficient
{
  int i,j;
  Scalaire angle,deca_s,decalages,decalaget,c,d,e;
  Scalaire asquare,bsquare;
  Scalaire co,si,cosquare,sisquare;
  Scalaire tone,ttwo;
  
  asquare=a*a;
  bsquare=b*b;
  for (i=0;i<p;i++)  {
    angle=psi-phi[i];
    decalages=r*cos(theta-phi[i]);
    decalaget=r*sin(theta-phi[i]);
    co=cos(angle);
    si=sin(angle);
    cosquare=co*co;
    sisquare=si*si;
    c=asquare*cosquare+bsquare*sisquare;
    for (j=1;j<=q;j++){
      deca_s=x[addresse(i,j)]-decalages;
      if((deca_s*deca_s)<c){
	d=a*b*sqrt(c-(deca_s*deca_s));
	e=deca_s*si*co*(bsquare-asquare);
	tone=decalaget+(d+e)/c;
	ttwo=decalaget+(-d+e)/c;
	if(mu!=0)
	  g[addresse(i,j)]=g[addresse(i,j)]+density*
	    (exp(mu*tone)-exp(mu*ttwo))/mu; 
	else
	  g[addresse(i,j)]=g[addresse(i,j)]+density*(tone-ttwo); 
	// we know that e has not to be computed but we want to keep
	// tone and ttwo for the attenuated transformation....
      }
    }
  }
}


void Sinogramme::Xdef(Scalaire radius)
// Initialize the abscisse to the usual one in analytic approach
// see Natterer 86
//
{
  int i,j;
  for (i=0;i<p;i++)
    for (j=1;j<=q;j++)
      x[addresse(i,j)]=(Scalaire)(-radius+(Scalaire)(j-1)*
				  (Scalaire)(2*radius)/(Scalaire)(q-1));
  // for slices (alebraic approach)
  //		x[addresse(i,j)]=(Scalaire)(-1+(Scalaire)1/
  //				(Scalaire)q+(Scalaire)(j-1)*(Scalaire)2./q);
}

void Sinogramme::Xshift(Scalaire *shift)
// shift of x vector
{ 
  int i,j;
  for (i=0;i<p;i++)
    for (j=1;j<=q;j++)
      x[addresse(i,j)]+=shift[i];
}

void Sinogramme::LitSino(char *filename)
// Initialise le Sinogramme aux valeurs lues dans un fichier
{
  FILE *pfile;
  int i,j;
  int nbrot, nbtrans; /* nb de rotations et de translations */
  float val;
  int retour;

/* ouverture du fichier en lecture */
  if ((pfile=fopen(filename,"r"))==NULL)
    { 
      printf("fichier %s non trouve \n",filename);
      exit(1);
    }

/* lecture de la premiere ligne : nbrotation   nbtranslations */
  retour=fscanf(pfile,"%d %d\n",&nbrot,&nbtrans);
  if (nbrot!=p)
    {printf("bad rotation number: %d different de %d\n",nbrot,p);
     exit(1);
   }
  if (nbtrans!=q)
    {printf("bad translation number: %d != %d\n",
	    nbtrans,q); 
     exit(1);
   }
  //	p=nbrot;
  //	q=nbtrans;
/* lecture des lignes suivantes */
  for (i=0;i<p;i++)
    for (j=1; j<=q;j++){
      retour=fscanf(pfile,"%e",&val);
      g[addresse(i,j)] = val;
    }
  fclose(pfile);
// now we use defaults for x and phi
// i.e. equi-sampling for x on [-1,1] and
// for phi on [0 PI]
//  Xdef();
//  for (i=0;i<p;i++)
//    phi[i]=(Scalaire)i*M_PI/(Scalaire)p;
//  std::cout<<" LitSino is no more doing Xdef and initializing phi"<<std::endl;
}

void Sinogramme::EcritSino(char *filename)
// Copy the Sinogramme in a file
{
  FILE *pfile;
  int i,j;
  //  int nbrot,nbtrans; /* nb de rotations et de translations */
  //float val;
  
  /* ouverture du fichier en ecriture */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("ERR- file %s not found \n",filename);
      exit(1);
    }

  /* write on the first line : nbrotation   nbtranslations */
  fprintf(pfile,"%d %d\n",p,q);
  /* write on the following lines */
  for (i=0;i<p;i++)
    for (j=1; j<=q;j++)
      fprintf(pfile,"%e ",g[addresse(i,j)]);
  fclose(pfile);
}



//
// Acces structure
// ===============
//

Scalaire & Sinogramme::operator () (int i,int j)
// Reference coefficient (i,j)
	{
#if CHECK
	if (i<0 || j<1 || i>=p || j>q)
		fputs("Erreur : Adressage hors Sinogramme\n",stderr);
#endif
	return *( g + (q*i+(j-1)) );
	}

const Scalaire & Sinogramme::operator () (const int i,const int j) const 
// Reference coefficient (i,j)
	{
#if CHECK
	if (i<0 || j<1 || i>=p || j>q)
		fputs("Erreur : Adressage hors Sinogramme\n",stderr);
#endif
	return *( g + (q*i+(j-1)) );
	}

Scalaire & Sinogramme::operator () (int i,Scalaire t)
// return the value of g(i,j) such that x[j] near t 
// only for equidistributed x
{
  int j;
  j=(int) (1 + (q-1) * (1.+ t)/2.);
  if(x[addresse(i,j)]>t)
    {
      fputs("Error in operator () : Not equaly spaced abscisse x \n",stderr);
      fputs("                     : Not supported \n",stderr);
      exit(1);
    }
  if(j<q)
    if(x[addresse(i,j+1)]<t)
      {
	fputs("Error in operator () : Not equaly spaced abscisse x \n",stderr);
	fputs("                     : Not supported \n",stderr);
	exit(1);
      }
#if CHECK
  if (i<0 || j<1 || i>=p || j>q)
    {
      fputs("Erreur : Adressage hors Sinogramme\n",stderr);
      exit(1);
    }
#endif
  return *( g + (q*i+(j-1)) );
}

void Sinogramme::Xcopy(int i,Projection *u)
	{
	int j;
#if CHECK
	if (i<0 || i>=p )
		fputs("Err : Adressage hors Sinogramme\n",stderr);
	if (u->N() != nq)
		fputs("Err : Bad projection length\n",stderr);
#endif
	for (j=1;j<=q;j++)
		x[addresse(i,j)]=u->X(j);	
	}

void Sinogramme::Xcopy(const Sinogramme *t)
	{
	int i,j;
#if CHECK
	if ((t->NQ() != q)||(t->NP() != p))
		fputs("Err : Bad Dimmensions\n",stderr);
#endif
	for (i=0;i<p;i++)
	for (j=1;j<=q;j++)
		x[addresse(i,j)]=t->X(i,j);	
	}

void Sinogramme::Acopy(int i,Projection *u)
	{
#if CHECK
	if (i<0 || i>=p )
		fputs("Err : Adressage hors Sinogramme\n",stderr);
	if (u->N() != nq)
		fputs("Err : Bad projection length\n",stderr);
#endif
	phi[i]=u->A();
	}
void Sinogramme::Acopy(const Sinogramme *t)
	{
	int i;
#if CHECK
	if ((t->NP() != p))
		fputs("Err : Bad Dimmensions\n",stderr);
#endif
	for (i=0;i<=p;i++)
		phi[i]=t->A(i);
      }

void Sinogramme::MulAttHalfDisk(Scalaire mu)
// mult by the exponential  attenuation on 
// the unit disk of constant parameter mu
	{
	int i,j;
	float s;
	for (j=1;j<=q;j++)
	   for (i=0;i<p;i++)
	   {
		   s= (float) (x[addresse(i,j)]);
//		   printf(" s= %f, x[addresse(i,j)=%e, fabs (x[addresse(i,j)])=%f \n",s,x[addresse(i,j)],fabs (x[addresse(i,j)]));
		g[addresse(i,j)]=g[addresse(i,j)]*(Scalaire) 
				exp((Scalaire) (mu*sqrt(1.-s*s)));
		printf(" g= %e\n",g[addresse(i,j)]);
	   }
        }

void Sinogramme::MulScalProj(Scalaire scal,int i)
     // mult the projection number i by scal
{
  int j;
#if CHECK
  if (i<0 || i>=p )
    fputs("Err : Adressage hors Sinogramme\n",stderr);
#endif
  for (j=1;j<=q;j++)
    g[addresse(i,j)]=g[addresse(i,j)]*scal;
}

void Sinogramme::GetProj(int i,Projection *u)
// copy the projection number i in  u
	{
	int j;
#if CHECK
	if (i<0 || i>=p || u->N()!=q)
		fputs("Erreur : Adressage hors Sinogramme\n",stderr);
#endif
	for (j=1;j<=u->N();j++)
		(*u)(j)=(*this)(i,j);
	u->Xcopy(x+q*i);
	u->Acopy(phi[i]);
	}

void Sinogramme::SetProj(int i,Projection *u)
// Fixe la projection i a partir de u
	{
	int j;
#if CHECK
	if (i<0 || i>=p || u->N()!=q)
		fputs("Erreur : Adressage hors Sinogramme\n",stderr);
#endif
	for (j=1;j<=q;j++)
		(*this)(i,j)=(*u)(j);	
	for (j=1;j<=q;j++)
		x[addresse(i,j)]=u->X(j);	
	phi[i]=u->A();
	}


int Sinogramme::NP() const
// Nbre de lignes
	{
	return p;
	}

int Sinogramme::NQ() const
// Nbre de colonnes
	{
	return q;
	}

const Scalaire Sinogramme::A(const int i) const
// angle number i
	{
#if CHECK
	if (i<0 || i>=p)
		fputs("Erreur : Adressage hors Sinogramme\n",stderr);
#endif
	return phi[i];
	}
const Scalaire Sinogramme::X(const int i,const int j) const
// detector abscisse j for projection i (s_j for parallel projection,\alpha_j for FB) 
// i=0,...,p-1 ; j=1,...,q
// angle number i, translation j in parallel geometry and source position i and detector j in FB geometry
	{
#if CHECK
	if (i<0 || i>=p || j<1 || j>q)
	    {
		fputs("Erreur : Adressage hors Sinogramme\n",stderr);
		return 0.;
	    }	
#endif
	return x[i*q+j-1];
	}
void Sinogramme::Angles(Scalaire *angle) 
   // return all the angles in the Scalaire vector *angles (angle[i]=phi[i])
{ 
  for (int i=0;i<p;i++) angle[i]=phi[i];
}
   void Sinogramme::Detectors(Scalaire *detector) 
    // return all the detector positions in the Scalaire vector *detector 
    // (detector[i*q+j]=x[i*q+j]
{ 
  for (int i=0;i<p;i++) 
    for (int j=0;j<q;j++)
      detector[i*q+j]=x[i*q+j];
}

 void Sinogramme::Examine()
// show the contents of the Sinogramme
        { 
        int i,j;
        cout << " Projection : " << p << endl;
        cout << " Translation : " << q << endl;
        for (i=0;i<p;i++)
            {	
		cout << " Projection no : " << i <<":" << endl;
	        cout << " Angle : " << phi[i] << endl;
        	cout << " X : " << endl;
        	for (j=1;j<=q;j++)
                	cout << x[addresse(i,j)] << " " ;
        	cout  << endl;
        	cout << " G : " << endl;
        	for (j=1;j<=q;j++)
           	     cout << g[addresse(i,j)] << " " ;
        	cout  << endl;
        	cout << " ============ " << endl;
	    }
	cout << " ======================== " << endl;
       	cout << " ======================== " << endl;
       	cout << " ======================== " << endl;
        }

void Sinogramme::PGMWrite(char *filename)
// Write the image in vmi format in the file filename.

	{
	FILE *pfile;
	int i,j,row,col,rows,cols;
	float val;
	int ival;
	float min,max,scale;

	rows=p;
	cols=q;
/* open the file for writing in */
	if ((pfile=fopen(filename,"w"))==NULL)
	     { printf("file %s non found \n",filename);
	       exit(1);
	     }
	max=g[0];
	min=g[0];
	for (i=0;i<rows;i++)
	 for (j=0; j<cols;j++)
	   {
		val=g[i*cols+j];
		if(max<val) max = val;
		if(min>val) min = val;
	   }
	if(max!=min) 
	   scale=255./(max-min);
	else
	   scale=1.;

	fprintf(pfile,"P2\n");
	fprintf(pfile,"# done by Sinogramme::PGMWrite \n");
	fprintf(pfile,"%d %d \n",cols,rows);
	fprintf(pfile,"%d \n",255);
        for ( row = 0; row < rows; ++row)
            {
            for ( col = 0; col < cols; ++col)
                {
                ival = (int) ( scale * ( g[row*cols+col] - min) );
		fprintf(pfile,"%d ",ival);
		}
	    fprintf(pfile," \n");
	    }
	fclose(pfile);
	}

void Sinogramme::PGMrbWrite(char *filename)
     // Write the image in PGM RawByte format in the file filename.

{
  FILE *pfile;
  int i,j,row,col,rows,cols;
  float val;
  int ival;
  float min,max,scale;

  rows=p;
  cols=q;
  /* open the file for writing in */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
  max=g[0];
  min=g[0];
  for (i=0;i<rows;i++)
    for (j=0; j<cols;j++)
      {
	val=g[i*cols+j];
	if(max<val) max = val;
	if(min>val) min = val;
      }
  cout << "-----------------------------------"<< endl;
  cout << " the minimum is : " << min << endl;
  cout << " the maximum is : " << max << endl;
  cout << "-----------------------------------"<< endl;
  if(max!=min) 
    scale=255./(max-min);
  else
    scale=1.;

  unsigned char *uc;

  uc = new unsigned char[rows*cols];
  
  fprintf(pfile,"P5\n");
  fprintf(pfile,"# done by Sinogramme::PGMrbWrite \n");
  fprintf(pfile,"# Min = %f\n",min);
  fprintf(pfile,"# Max = %f\n",max);
  fprintf(pfile,"%d %d \n",cols,rows);
  fprintf(pfile,"%d\n",255);   /* NO WHITE SPACE !!! */
  for ( row = 0; row < rows; ++row)
    for ( col = 0; col < cols; ++col)
      uc[row*cols+col] = 
	(unsigned char) ( scale * ( g[row*cols+col] - min) );
  fwrite(uc,1,rows*cols,pfile);
  fclose(pfile);
  delete [] uc;
}

void Sinogramme::PGMWriteLog(char *filename)
// Write the image in vmi format in the file filename.
// using a log scale on the image
{
  FILE *pfile;
  int i,j,row,col,rows,cols;
  float val;
  int ival;
  float min,max,scale;
  
  rows=p;
  cols=q;
  /* open the file for writing in */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
  max=g[0];
  min=g[0];
  for (i=0;i<rows;i++)
    for (j=0; j<cols;j++)
      {
	val=g[i*cols+j];
	if(max<val) max = val;
	if(min>val) min = val;
      }
  if(min<=0)
    {
      printf(" error: Negative value, Log scale cannot be used");
      return;
    }
  
  if(max!=min) 
    scale=255./(log(max)-log(min));
  else
    scale=1.;
  
  fprintf(pfile,"P2\n");
  fprintf(pfile,"# done by Sinogramme::PGMWrite \n");
  fprintf(pfile,"%d %d \n",cols,rows);
  fprintf(pfile,"%d \n",255);
  for ( row = 0; row < rows; ++row)
    {
      for ( col = 0; col < cols; ++col)
	{
	  ival = (int) ( scale * ( log(g[row*cols+col]) - log(min) ) );
	  fprintf(pfile,"%d ",ival);
	}
      fprintf(pfile," \n");
    }
  fclose(pfile);
}

double Sinogramme::MaxAbs()
// return the MaxAbs of the sinogramme
{ 
  double maximum=0;
  for(int i=0;i<p;i++)
    for(int j=0;j<q-1;j++)
      if(fabs(g[i*q+j])>maximum)maximum=fabs(g[i*q+j]);
  return(maximum);
}

void Sinogramme::DCC(double *dccval,int order)
// Compute the Moments (for each projection) of order "order" of the parallel projection (*this)
// Compute de the Data Consitency Conditions of order order
{
  double dcci;
  for(int i=0;i<p;i++){
    dcci=0;
    //formule des trapèzes si échantillonnage irrégulier
    for(int j=1;j<q-1;j++)
      dcci+=g[i*q+j]*pow(x[i*q+j],(double)order)*(x[i*q+j+1]-x[i*q+j-1])/2;
    dcci+=g[i*q]*pow(x[i*q],(double)order)*(x[i*q+1]-x[i*q])/2;
    dcci+=g[i*q+q-1]*pow(x[i*q+q-1],(double)order)*(x[i*q+q-1]-x[i*q+q-2])/2;
    dccval[i]=dcci;
  }
}

//
// Transform
// =========
//
void Sinogramme::Deriv() // derivation according to the scalar variable d/dq (g(p.q))
    //(*this)(i,j)=(*this(i,j+1)-*this(i,j-1))/2h
{
  int i,j;
  Scalaire gprime[q];
  for (i=0;i<p;i++){
    //    gprime[0]=(g[i*q+1]-g[i*q])/(x[i*q+1]*x[i*q]);
    gprime[0]=(g[i*q+1]-g[i*q]);
    for (j=1; j<q-1;j++)
      //      gprime[j]=(g[i*q+j+1]-g[i*q+j-1])/(x[i*q+j+1]*x[i*q+j-1]);
      gprime[j]=(g[i*q+j+1]-g[i*q+j-1]);
    // gprime[q-1]=(g[i*q+q-1]-g[i*q+q-2])/(x[i*q+q-1]*x[i*q+q-2]);
    gprime[q-1]=(g[i*q+q-1]-g[i*q+q-2]);
    for (j=0; j<q;j++)
      g[i*q+j]=gprime[j];
  }
}

void Sinogramme::MulExpSqrtUnMinusSsquare(Scalaire mu)
     //(*this)(i,j)=*this(i,j)*exp(mu*sqrt(1-x[i,j]*x[i,j]))
{
  int i,j;
  for (i=0;i<p;i++)
    for (j=0; j<q;j++)
	if(x[i*q+j]*x[i*q+j]<1)
	  g[i*q+j]=g[i*q+j]*exp(mu*sqrt(1-x[i*q+j]*x[i*q+j]));
}

void Sinogramme::ProjPositive(Scalaire epsilon)
// *this=epsilon if *this<epsilon
{
  int i,j;
  for (i=0;i<p;i++)
    for (j=0; j<q;j++)
	if(g[i*q+j]<epsilon)g[i*q+j]=epsilon;
}
void Sinogramme::TransPositive()
// *this=*this-min (if min<0)
{
  int i,j;
  Scalaire min,val;
  min=g[0];
  for (i=0;i<p;i++)
    for (j=0; j<q;j++)
      {
	val=g[i*q+j];
	if(min>val) min = val;
      }
  if(min<0)
    {
      for (i=0;i<p;i++)
	for (j=0; j<q;j++)
	  g[i*q+j]=g[i*q+j]-min;
    }
}
      
  
void Sinogramme::Complete(Sinogramme *a)
/* complete the considered Sinogram a
    from [0,pi[x[-1,1] to 
    [0,2*pi[x[-1,1] with
    g(phi+pi,s)=g(phi,-s)
*/
{
  int i,j;
#if CHECK
  if ( (2*a->NP() != p) || (a->NQ() != q) )
    fputs("Err : Bad dimmenssions in Complete\n",stderr);
#endif
  for (i=0;i<p/2;i++)
    for (j=1;j<=q;j++)
      g[addresse(i,j)]=(*a)(i,j);
  for (i=p/2;i<p;i++)
    for (j=1;j<=q;j++)
      g[addresse(i,j)]=(*a)(i-p/2,q-j+1);
  for (i=0;i<=p/2;i++)
    {
      phi[i]=a->A(i);
      phi[i+p/2]=a->A(i)+M_PI;
      for (j=1;j<=q;j++)
	{
	  x[addresse(i,j)]=a->X(i,j);	
	  x[addresse(i+p/2,j)]=a->X(i,j);	
	}
    }
}


//
// Operateurs
// ==========
//

	     
void Copy(Sinogramme *a,Sinogramme *b)
// b <- a
	{
	int i,j;
#if CHECK
	if (BadDim(a,b))
	  {
		fputs("Err : Bad dimmenssions in Copying\n",stderr);
		printf("Err : Bad dimmenssions in Copying\n");
	  }
#endif
	for (i=0;i<b->NP();i++)
	for (j=1;j<=b->NQ();j++)
		(*b)(i,j)=(*a)(i,j);
	b->Xcopy(a);
	b->Acopy(a);
	}

void Copy(Projection *u,int i,Sinogramme *a)
// a(i,.) <- u
	{
	int j;
#if CHECK
	if (i<0 || i>=p || u->N()!=q)
		fputs("Err : Adress out of range of Sinogramme\n",stderr);
	if (a->NQ()!=u->N())
		fputs("Err : Bad dimmenssions in Copying\n",stderr);
#endif
	for (j=1;j<=a->NQ();j++)
		(*a)(i,j)=(*u)(j);
	a->Xcopy(i,u);
	a->Acopy(i,u);
	}

void Add(Sinogramme *a,Sinogramme *b,Sinogramme *c)
// c <- a+b
	{
	int i,j;
#if CHECK
	if ((BadDim(a,b)||BadDim(a,c)))
		fputs("Err : Bad dimmenssions in Copying\n",stderr);
#endif
	for (i=0;i<c->NP();i++)
	for (j=1;j<=c->NQ();j++)
		(*c)(i,j)=(*a)(i,j)+(*b)(i,j);
	c->Xcopy(a);
	c->Acopy(a);
	}

void Sub(Sinogramme *a,Sinogramme *b,Sinogramme *c)
// c <- a-b
{
  int i,j;
#if CHECK
  if ((BadDim(a,b)||BadDim(a,c)))
    fputs("Err : Bad dimmenssions in Copying\n",stderr);
#endif
  for (i=0;i<c->NP();i++)
    for (j=1;j<=c->NQ();j++)
      (*c)(i,j)=(*a)(i,j) - (*b)(i,j);
  c->Xcopy(a);
  c->Acopy(a);
}

void KAdd(Sinogramme *a,Scalaire k,Sinogramme *b,Sinogramme *c)
// c <- a+k.b
	{
	int i,j;
#if CHECK
	if ((BadDim(a,b)||BadDim(a,c)))
		fputs("Err : Bad dimmenssions in KAdd\n",stderr);
#endif
	for (i=0;i<c->NP();i++)
	for (j=1;j<=c->NQ();j++)
		(*c)(i,j)=(*a)(i,j)+k*(*b)(i,j);
	c->Xcopy(a);
	c->Acopy(a);
	}

void Prod(Sinogramme *a,Sinogramme *b,Sinogramme *c)
// c <- a.b
	{
	int i,j;

#if CHECK
	if ((BadDim(a,b)||BadDim(a,c)))
		fputs("Err : Bad dimmenssions in Prod\n",stderr);
#endif
	for (i=0;i<a->NP();i++)
	for (j=1;j<=a->NQ();j++)
		(*c)(i,j)=(*a)(i,j) * (*b)(i,j);
	c->Xcopy(a);
	c->Acopy(a);
	}

void KMul(Scalaire k,Sinogramme *a,Sinogramme *b)
// b <- k.a
	{
	int i,j;
#if CHECK
	if (BadDim(a,b))
		fputs("Err : Bad dimmenssions in KMul\n",stderr);
#endif
	for (i=0;i<b->NP();i++)
	for (j=1;j<=b->NQ();j++)
		(*b)(i,j)=k*(*a)(i,j);
	b->Xcopy(a);
	b->Acopy(a);
	}

void Transpose(Sinogramme *a,Sinogramme *b)
// b <- Trans(a)
{
  // On cree une copie du resultat pour permettre a=b
  Sinogramme *c;
  int i,j;
#if CHECK
  if (BadDim(a,b))
    fputs("Err : Bad dimmenssions in Transpose\n",stderr);
#endif

  c=new Sinogramme(b->NP(),b->NQ());
  for (i=0;i<c->NP();i++)
    for (j=1;j<=c->NQ();j++)
      (*c)(i,j)=(*a)(j,i);
  c->Xcopy(a);
  c->Acopy(a);
  Copy(c,b);
  delete c;
}

void Sin(Sinogramme *a,Sinogramme *b)
// b <- sin(a)
{
  int i,j;
#if CHECK
  if (BadDim(a,b))
    fputs("Err : Bad dimmenssions in Transpose\n",stderr);
#endif
  for (i=0;i<b->NP();i++)
    for (j=1;j<=b->NQ();j++)
      (*b)(i,j)=sin((*a)(i,j));
}

void Cos(Sinogramme *a,Sinogramme *b)
// b <- cos(a)
{
  int i,j;
#if CHECK
  if (BadDim(a,b))
    fputs("Err : Bad dimmenssions in Transpose\n",stderr);
#endif
  for (i=0;i<b->NP();i++)
    for (j=1;j<=b->NQ();j++)
      (*b)(i,j)=cos((*a)(i,j));
}

void Exp(Sinogramme *a,Sinogramme *b)
// b <- exp(a)
{
  int i,j;
#if CHECK
  if (BadDim(a,b))
    fputs("Err : Bad dimmenssions in Transpose\n",stderr);
#endif
  for (i=0;i<b->NP();i++)
    for (j=1;j<=b->NQ();j++)
      (*b)(i,j)=exp((*a)(i,j));
}

//
// Operateurs Projection-Sinogramme
// ==========================
//

void Convolution(Projection *u,Sinogramme *a,Sinogramme *b)
// b(i,.) <- u convol a(i,.) (for all i)
// The dimension of u must be 2*a->NQ()
// a and b contains projections whereas u is supposed to be the filter
{
  int i,j,l,n;
  Scalaire h,somme;
#if CHECK
  if ( (BadDim(a,b)) || (BadDim2QmoinsUn(a,u)) )
    fputs("Err : Bad dimmenssions in Convolution\n",stderr);
#endif
  n=a->NQ();
  h=2./(Scalaire)(n-1);
  for (i=0;i<a->NP();i++) {
    for (j=1;j<=n;j++){
      somme=0.;
      for (l=1;l<=n;l++)
	somme=somme+(*u)(j-l+n)*(*a)(i,l);
      (*b)(i,j)=somme*h;
    }
  }	
}


void LocConvolution(Projection *u,Sinogramme *a,Sinogramme *b,
		    int halfwidth)
     // b(i,.) <- u convol a(i,.) (for all i)
     // The dimension of u must be 2*a->NQ()
     // do the convolution in i on indices between 
     // i-halfwidth and i+halfwidth
     // have been adapted for local and pseudo-local tomography
{
  int i,j,l,n,nproj;
  Scalaire h,somme;
  //  Scalaire fangle;
  n=a->NQ();
  nproj=u->N();
  //  fangle=a->FA();
#if CHECK
  if ( (BadDim(a,b)) || (BadDim2QmoinsUn(a,u)) )
    {
      fputs("Err : Bad dimensions in Convolution\n",stderr);
      return;
    }
  if ( nproj<(2*halfwifth+1) )
    {
      fputs("Err : Bad dimensions in Convolution\n",stderr);
      fputs("Err : convolutor array dim must be > 2*halfwifth+1 \n",stderr);
      return;
    }
#endif
  if(n==1)
    {
	 puts("Error : Cannot filter with only 1 translation!");
      return;
    }
  else
    {
      h=2./(Scalaire)(n-1);
      //      h=fangle/(n-1);
    }

  for (i=0;i<a->NP();i++)
    {
      for (j=halfwidth+1;j<=n-halfwidth-1;j++)
	{
	  somme=0.;
	  for (l=j-halfwidth;l<=j+halfwidth;l++)
	    somme=somme+(*u)(j-l+(nproj+1)/2)*(*a)(i,l);
	  (*b)(i,j)=somme*h;
	}
      for (j=1;j<halfwidth+1;j++)
	(*b)(i,j)=(*b)(i,halfwidth+1);
      for (j=n-halfwidth;j<=n;j++)
	(*b)(i,j)=(*b)(i,n-halfwidth-1);
    }	
}


// ICI
Scalaire InterlacedInterpol(Sinogramme *f,Sinogramme *a,int pint, int qint,int verbose)
// a(pint,qint) <- Sum_{ all i and all j possible in the corresponding INTERLACED grid}   }a(pint-i,qint-j)).f(i,j)) 
// Note that the sum is periodic in p (as Sinogramme but not in q)
// work only for Sinogramme defined on 2Pi
// the dimensions of *f must be ODD
{
  int ip,iq;
  int decap, decaq;
  int fp=f->NP();
  int fq=f->NQ();
  int fpsd=fp/2;
  int fqsd=fq/2;
  if ( ((fpsd*2) == fp) || ((fqsd*2) == fq) ) {
    cout << " Erreur dans InterlacedInterpol :: " << endl ; 
    cout << " les dimensions du Sinogramme de filtre doivent etre impaires" << endl ; 
    exit(1);
  }
  int modecap;
  int ap=a->NP();
  int aq=a->NQ();
  Scalaire normal=0;
  Scalaire interp=0;
  for (ip=0;ip<fp;ip++)
    for (iq=1;iq<(fq+1);iq++){
      decap=(ip-fpsd);
      decaq=(iq-1-fqsd);
	if(abs(decap%2)!=abs(decaq%2)) { // entrelace
	if(((qint-decaq)>0) && ((qint-decaq)<=aq)){
	  modecap=(pint-decap)%ap;
	  if(modecap<0)modecap=modecap+ap;
	  interp+=(*f)(ip,iq) * (*a)((pint-decap)%ap,qint-decaq);
	}
	normal+=(*f)(ip,iq);
      }
    }
  if(verbose){
    Scalaire interpclassique;
    interpclassique=(*a)(pint,qint+1) + (*a)(pint,qint-1);
    interpclassique += ( (*a)(pint+1,qint+1) + (*a)(pint+1,qint-1) +
			 (*a)(pint-1,qint+1) + (*a)(pint-1,qint-1))/sqrt(2.);
    interpclassique = interpclassique / (4/sqrt(2.)+2);
    cout << " valeur `a interpoler " << (*a)(pint,qint) << endl;
    cout << " Normalisation " << normal <<endl ;
    cout << " valeur `a interpoler normalis'ee " << interp/normal << endl ;
    cout << " Valeur classiquement utilis'ee : " << interpclassique << endl ;
  }
  return interp/normal ;
}




//
// Operateurs Sinogramme-Image
// ================================================
//

void Retroprojection(Sinogramme *t,RealImage *image, Scalaire AngularInterval)
// image <- Retroprojection t
// see F.Natterer 
// "The mathematics of Computerized Sinography" (Wiley 86) p.109
{
  int i,j,nx,ny; // for image(i,j),i=0,...,nx-1;j=0,...,ny-1
  int ip,iq,p,q; // for t(ip,iq),i=0,...,p-1;j=1,...,q
  Scalaire co,si,s; // cosinus, sinus, and s=<x,theta> (see Natterer).
  Scalaire xima,yima; // abs and ord of the current pixel.
  Scalaire xbegin,xend,ybegin,yend; // describe the 4 image corners. 
  Scalaire h;// in case of regularly sampled projections, h is
  //the distance between 2 "s" points.
  Scalaire u;// proportion 
  
  nx=image->NX();
  ny=image->NY();
  p=t->NP();
  q=t->NQ();
  Scalaire dphi [p] ;// integral weight for each proje in the backproj (jacobian     
  std::cout <<" Retroprojection angular Interval "<< AngularInterval<< std::endl;
  if((t->A(0)+AngularInterval)<t->A(p-1)){
    std::cout <<" ERR  Retroprojection(Sinogramme *t,RealImage *image): (t->A(0)+AngularInt)<t->A(p-1) " << std::endl;
    exit(0);
  }

  dphi[0]=(t->A(1)-t->A(0)+t->A(0)+AngularInterval-t->A(p-1))/2;
  dphi[p-1]=(t->A(p-1)-t->A(p-2)+t->A(0)+AngularInterval-t->A(p-1))/2;
  for (ip=1;ip<p-1;ip++)
    dphi[ip]=(t->A(ip+1)-t->A(ip-1))/2;
  //printf("avant GetCorners\n");
  image->GetCorners(&xbegin,&xend,&ybegin,&yend);
//printf("xbegin,xend,ybegin,yend %e %e %e %e \n",xbegin,xend,ybegin,yend);
//printf("avant for i =0 to nx=%d\n",nx);
  for (i=0;i<nx;i++){
    //		  printf("i = %d \n", i);
    if((i%10)==0)std::cout<<"*"<<std::flush;    
    xima=xbegin+(i+.5)*((xend-xbegin)/nx);
    for (j=0;j<ny;j++)	{
      (*image)(i,j)=0.;
//printf("j = %d \n", j);
      yima=ybegin+(j+.5)*((yend-ybegin)/ny);
      if ((xima*xima+yima*yima)<1)
	for (ip=0;ip<p;ip++){
	  co=cos(t->A(ip));
	  si=sin(t->A(ip));
	  s=co*xima+si*yima;
	  h=t->X(ip,2)-t->X(ip,1);
	  iq=int((s-t->X(ip,1))/h)+1;
	  // LOCALIZE
	  //				if((s<t->X(ip,iq))||(s>t->X(ip,iq+1)))
	  //				   {
	  //				   Localize(s,ip,t,iq);
	  //printf("x= %f,y=%f\n",xima,yima); 
	  //printf("iq= %d, t->X(ip,iq)= %f,s=%f, t->X(ip,iq+1)= %f \n",iq,t->X(ip,iq),s,t->X(ip,iq+1)); 
	  //				   }
	  u=(s-t->X(ip,iq))/h;
	  if (iq!=q)
	    (*image)(i,j)=(*image)(i,j)
	      +((1.-u)*(*t)(ip,iq)+u*(*t)(ip,iq+1))*dphi[ip];
	  else
	    (*image)(i,j)=(*image)(i,j)+((*t)(ip,iq)*dphi[ip]);
	}
      //      (*image)(i,j)=2.*M_PI*(*image)(i,j)/t->NP();
    }
  }
  std::cout<<std::endl;    
}


void Retroprojection(Sinogramme *t,RealImage *image,Scalaire mu,
		     Scalaire AngularInterval)
// image <- Attenuated Retroprojection t  (Dual of the exponential
//                                         Radon transform)
// see F.Natterer 
// "The mathematics of Computerized Sinography" (Wiley 86) p.47 and p.109
{
  int i,j,nx,ny; // for image(i,j),i=0,...,nx-1;j=0,...,ny-1
  int ip,iq,p,q; // for t(ip,iq),i=0,...,p-1;j=1,...,q
  Scalaire co,si,s,sortho; // cosinus, sinus, 
  // s=<x,theta> , sortho=<x,thetaortho> (see Natterer).
  Scalaire xima,yima; // abs and ord of the current pixel.
  Scalaire xbegin,xend,ybegin,yend; // describe the 4 image corners. 
  Scalaire h;// in case of regularly sampled projections, h is
  //the distance between 2 points.
  Scalaire u;// proportion 
  
  nx=image->NX();
  ny=image->NY();
  p=t->NP();
  q=t->NQ();
  Scalaire dphi [p] ;// integral weight for each proj in the backproj jacobian 
  std::cout <<" Retroprojection angular Interval "<< AngularInterval<< std::endl;
  if((t->A(0)+AngularInterval)<t->A(p-1)){
    std::cout <<" ERR  Retroprojection(Sinogramme *t,RealImage *image): (t->A(0)+AngularInt)<t->A(p-1) " << std::endl;
    exit(0);
  }
  dphi[0]=(t->A(1)-t->A(0)+t->A(0)+AngularInterval-t->A(p-1))/2;
  dphi[p-1]=(t->A(p-1)-t->A(p-2)+t->A(0)+AngularInterval-t->A(p-1))/2;
  for (ip=1;ip<p-1;ip++)
    dphi[ip]=(t->A(ip+1)-t->A(ip-1))/2;
  //printf("avant GetCorners\n");
  image->GetCorners(&xbegin,&xend,&ybegin,&yend);
  //printf("xbegin,xend,ybegin,yend %e %e %e %e \n",xbegin,xend,ybegin,yend);
  //printf("avant for i =0 to nx=%d\n",nx);
  for (i=0;i<nx;i++)
    {
      printf("i =  %d \n", i);
      xima=xbegin+(i+.5)*((xend-xbegin)/nx);
      for (j=0;j<ny;j++)
	{
	  (*image)(i,j)=0.;
//printf("j = %d \n", j);
	  yima=ybegin+(j+.5)*((yend-ybegin)/ny);
	  if ((xima*xima+yima*yima)<1)
	    for (ip=0;ip<t->NP();ip++) {
	      co=cos(t->A(ip)); 
	      si=sin(t->A(ip));
	      s=co*xima+si*yima;
	      sortho=-si*xima+co*yima;
	      h=t->X(ip,2)-t->X(ip,1);
	      iq=int((s-t->X(ip,1))/h)+1;
	      // LOCALIZE
	      //				if((s<t->X(ip,iq))||(s>t->X(ip,iq+1)))
	      //				   {
	      //				   Localize(s,ip,t,iq);
	      //printf("x= %f,y=%f\n",xima,yima); 
	      //printf("iq= %d, t->X(ip,iq)= %f,s=%f, t->X(ip,iq+1)= %f \n",iq,t->X(ip,iq),s,t->X(ip,iq+1)); 
	      //				   }
	      u=(s-t->X(ip,iq))/h;
	      if (iq!=q)
		(*image)(i,j)=(*image)(i,j) +
		  ((1.-u)*(*t)(ip,iq)+u*(*t)(ip,iq+1))
		  * exp(mu*sortho)*dphi[ip];
	      else
		(*image)(i,j)=(*image)(i,j)+
		  (*t)(ip,iq) * exp(mu*sortho)*dphi[ip];
	      }
	  //	  (*image)(i,j)=2.*M_PI*(*image)(i,j)/t->NP();
	}
    }
}

Scalaire Retroprojection(Sinogramme *t, Scalaire xima, Scalaire yima, 
			 Scalaire AngularInterval)
// return Retroprojection of *t at (xima,yima)
//  Scalaire xima,yima; // abs and ord of the current pixel.

// see F.Natterer 
// "The mathematics of Computerized Sinography" (Wiley 86) p.109
{
  int ip,iq,p,q; // for t(ip,iq),i=0,...,p-1;j=1,...,q
  Scalaire co,si,s; // cosinus, sinus, and s=<(x,y),theta> (see Natterer).
  Scalaire h;// in case of regularly sampled projections, h is
  //the distance between 2 points.
  Scalaire u;// proportion 
  
  p=t->NP();
  q=t->NQ();
  Scalaire dphi [p] ;// integral weight for each proj in the backproj jacobian 
  // std::cout <<" Retroprojection angular Interval "<< AngularInterval<< std::endl;
  if((t->A(0)+AngularInterval)<t->A(p-1)){
    std::cout <<" ERR  Retroprojection(Sinogramme *t,RealImage *image): (t->A(0)+AngularInt)<t->A(p-1) " << std::endl;
    exit(0);
  }

  dphi[0]=(t->A(1)-t->A(0)+t->A(0)+AngularInterval-t->A(p-1))/2;
  dphi[p-1]=(t->A(p-1)-t->A(p-2)+t->A(0)+AngularInterval-t->A(p-1))/2;
  for (ip=1;ip<p-1;ip++)
    dphi[ip]=(t->A(ip+1)-t->A(ip-1))/2;
  Scalaire backproj=0;
  int OK=1; 
  for (ip=0;ip<p;ip++){
    co=cos(t->A(ip));
    si=sin(t->A(ip));
    s=co*xima+si*yima;
    if((s>t->X(ip,1))&& (s<t->X(ip,q))){
      h=t->X(ip,2)-t->X(ip,1);
      //iq=int((s-t->X(ip,1))/h)+1;
      iq=(int)floor((s-t->X(ip,1))/h)+1;
      // LOCALIZE
      //				if((s<t->X(ip,iq))||(s>t->X(ip,iq+1)))
      //				   {
      //				   Localize(s,ip,t,iq);
      //printf("x= %f,y=%f\n",xima,yima); 
      //printf("iq= %d, t->X(ip,iq)= %f,s=%f, t->X(ip,iq+1)= %f \n",iq,t->X(ip,iq),s,t->X(ip,iq+1)); 
	  //				   }
      u=(s-t->X(ip,iq))/h;
      if (iq!=q)
	backproj += ((1.-u)*(*t)(ip,iq)+u*(*t)(ip,iq+1))*dphi[ip];
      else
	    backproj += (*t)(ip,iq)*dphi[ip];
    }
    else{
      OK=0;
      //      std::cout<<" s smin smax " << s<<" , "<<t->X(ip,1)<<" , "<<t->X(ip,q)<<std::endl;
    }
  }
  if(!OK)    backproj=0;
    //    backproj*=2.*M_PI/t->NP();
    //  else{
    //    std::cout<<" Retroprojection(Sinogramme *t, Scalaire xima, Scalaire yima) : pixel is out of [smin smax]" <<t->X(ip,1)<<" , "<<t->X(ip,q)<<std::endl;    
  //}
  return(backproj);
}

Scalaire RolfRetroprojection(Sinogramme *t, Scalaire xima, Scalaire yima,
			     Scalaire phimin,Scalaire phimax,
			     int n)
// return Retroprojection of *t at (xima,yima)
//  Scalaire xima,yima; // abs and ord of the current pixel.
// the back projection  angle interval is  phimin,phimax
// n is for the weight "tan(phi)^n/cos(phi)"


// see F.Natterer 
// "The mathematics of Computerized Sinography" (Wiley 86) p.109
{
  int ip,iq,p,q; // for t(ip,iq),i=0,...,p-1;j=1,...,q
  Scalaire co,si,s; // cosinus, sinus, and s=<(x,y),theta> (see Natterer).
  Scalaire h;// in case of regularly sampled projections, h is
  //the distance between 2 points.
  Scalaire u;// proportion 
  
  p=t->NP();
  q=t->NQ();
  Scalaire backproj=0;
  int OK=1; 
  int ipmin,ipmax;
  ipmin=0;ipmax=p;
  for (ip=0;ip<p;ip++){
    if(phimin>t->A(ip))ipmin=ip;
    if(phimax<t->A(ip))ipmax=ip;
  }
  if(phimin>t->A(ipmin))ipmin++;
  //  if(phimax<t->A(ipmax))ipmax--;

  for (ip=ipmin;ip<ipmax;ip++){
    co=cos(t->A(ip));
    si=sin(t->A(ip));
    s=co*xima+si*yima;
    if((s>t->X(ip,1))&& (s<t->X(ip,q))){
      h=t->X(ip,2)-t->X(ip,1);
      //iq=int((s-t->X(ip,1))/h)+1;
      iq=(int)floor((s-t->X(ip,1))/h)+1;
      // LOCALIZE
      //				if((s<t->X(ip,iq))||(s>t->X(ip,iq+1)))
      //				   {
      //				   Localize(s,ip,t,iq);
      //printf("x= %f,y=%f\n",xima,yima); 
      //printf("iq= %d, t->X(ip,iq)= %f,s=%f, t->X(ip,iq+1)= %f \n",iq,t->X(ip,iq),s,t->X(ip,iq+1)); 
	  //				   }
      u=(s-t->X(ip,iq))/h;
      if (iq!=q)
	backproj += ((1.-u)*(*t)(ip,iq)+u*(*t)(ip,iq+1))*pow((si/co),(double)n)/co;
      else
	backproj += (*t)(ip,iq)*pow((si/co),(double)n)/co;
    }
    else{
      OK=0;
      //      std::cout<<" s smin smax " << s<<" , "<<t->X(ip,1)<<" , "<<t->X(ip,q)<<std::endl;
    }
  }
  if(OK)
    backproj*=(t->A(ipmax-1)-t->A(ipmin))/(ipmax-ipmin-1);//2.*M_PI/(ipmax-ipmin+1); NON backproj on PI-epslison, i.e. on  [t->A(ipmin) t->A(ipmax-1]
  else{
    backproj=0;
    //    std::cout<<" Retroprojection(Sinogramme *t, Scalaire xima, Scalaire yima) : pixel is out of [smin smax]" <<t->X(ip,1)<<" , "<<t->X(ip,q)<<std::endl;    
  }
  return(backproj);
}

//void FanBeamRetroprojection(Sinogramme *t,RealImage *image)
// image <- Retroprojection t
// see F.Natterer 
// "The mathematics of Computerized Sinography" (Wiley 86) p.113


//
// tool box
//
int BadDim(Sinogramme *a,Sinogramme *b)
// ((a->NP != b->NP) || (a->NQ != b->NQ))
{
return ( (a->NP() != b->NP()) || (a->NQ() != b->NQ()) );
}

int BadDim2(Sinogramme *a,Projection *p)
// (2.a->NQ != b->N)
{
return ( 2.*a->NQ() != p->N() );
}

int BadDim2QmoinsUn(Sinogramme *a,Projection *p)
// (2.a->NQ-1 != b->N)
{
return ( (2.*a->NQ()-1) != p->N() );
}

void Localize(Scalaire s,int ip,Sinogramme *t,int iq)
// result: iq such that t->X(ip,iq)<= s < t->X(ip,iq+1)
{
printf("Err : unregular sampling\n");
//fputs("Err : I still can not localize s for unregular sampling\n",stderr);
// fputs(" but I will soon do ....\n",stderr);
}

void Extend(Sinogramme *a,int i,int n,float *w)
// w(j-1)<- a(i,j-n+n/2) if n-n/2<j<2n-n/2   (w(0,...,2n) whereas a(i,1,...,n))
// w(j)<- 0.		           else
// n must be a->NQ()
{
  int j;
#if CHECK
  if (n!=a->NQ())
    fputs("Err : Bad dimmenssions in Extend\n",stderr);
#endif
  for (j=0;j<2*n;j++)
    w[j]=0.;
  for (j=n-n/2+1;j<2*n-n/2+1;j++)
    w[j-1]=(float)(*a)(i,j-(n-n/2));
}

void InvExtend(float *w,int i,int n,Sinogramme *a)
// a(i,j-n+n/2)<- w(j-1) if n-n/2<j<2n-n/2   (w(0,...,2n) whereas a(i,1,...,n))
// n must be a->NQ()
{
  int j;
#if CHECK
  if (n!=a->NQ())
    fputs("Err : Bad dimmenssions in Extend\n",stderr);
#endif
  for (j=n-n/2+1;j<2*n-n/2+1;j++)
    (*a)(i,j-(n-n/2))=(double)w[j-1];
}


// interface subroutines 
void GetSimulMovingData(Sinogramme *p, int fileformat)
  /* simulation of data from a moving object
 the movement is coded by a matrix A(\phi) and a vector b(\phi)
   We simulate Rf_{A,b}(\phi,s)
     where f_{A,b}(x) is f(Ax+b) 
        note that A is a function of \phi as b
      note that only b\dot A^{-t}\theta has to be given
 *p is the result
 format if 1 then the number of projection is not present in the movement file
  */
{
  int i;
  char *response,*filename,car;
  FILE *pfile;
  int nbell,nbdisk,nbtri; // number of ellipsis, number of disks,
  //number of triangles (trinagles are not considered)
  response=new char[80];
  filename=new char[80];
  Scalaire er, etheta, epsi, ea, eb, edensity; // ellipse parameters
  float fr, ftheta, fpsi, fa, fb, fdensity; // ellipse parameters
  int nbangle;
  int retour;

  printf(" Movement file name ? : ");
  retour=scanf("%s",filename);
/* open the file for reading in */
  if ((pfile=fopen(filename,"r"))==NULL)
    { printf("file %s non found \n",filename);
    exit(1);
    }
  if(fileformat){
    retour=fscanf(pfile,"%d",&nbangle);
    if(p->NP()!=nbangle){
      cout << " ERR in GetSimulMovingData : Unconsistent projection number "<< 
	nbangle << " in file "<<filename<<endl;
      exit(1);
    }
  }
  else
    nbangle=p->NP();

  Scalaire *A11,*A12,*A21,*A22,*TransScal;
  A11 = new Scalaire [nbangle];
  A12 = new Scalaire [nbangle];
  A21 = new Scalaire [nbangle];
  A22 = new Scalaire [nbangle];
  TransScal = new Scalaire [nbangle];
  for(i=0;i<nbangle;i++)
    retour=fscanf(pfile,"%lf %lf %lf %lf %lf",
	   &A11[i],&A12[i],&A21[i],&A22[i],&TransScal[i]);
  fclose(pfile);

  printf(" phantom definition file name ? : ");
  retour=scanf("%s",filename);
/* open the file for reading in */
  if ((pfile=fopen(filename,"r"))==NULL)
    { printf("file %s non found \n",filename);
    exit(1);
    }
  retour=fscanf(pfile,"%d",&nbell);
  retour=fscanf(pfile,"%d",&nbdisk);
  retour=fscanf(pfile,"%d",&nbtri);
  if(nbtri!=0)
    printf(" triangles are not considered ...\n");
  for(i=0;i<nbell;i++){
    retour=fscanf(pfile,"%f %f %f %f %f %f ",
	   &fr,&ftheta,&fa,&fb,&fpsi,&fdensity);
    er=fr; etheta=ftheta; epsi=fpsi; ea=fa; eb=fb; edensity=fdensity;
    printf("r, theta, psi, a, b, density = :%e %e %e %e %e %e\n ", 
	   er, etheta, epsi, ea, eb,edensity);
    p->AddMovingEllipse(er, etheta, epsi, ea, eb, edensity,
			A11,A12,A21,A22,TransScal);
  }
  for(i=0;i<nbdisk;i++){
    retour=fscanf(pfile,"%f %f %f %f ",
	   &fr,&ftheta,&fa,&fdensity);
    er=fr; etheta=ftheta; epsi=0.; ea=fa; eb=fa; edensity=fdensity;
    printf("r, theta, psi, a, b, density = :%e %e %e %e %e %e\n ", 
	   er, etheta, epsi, ea, eb,edensity);
    p->AddMovingEllipse(er, etheta, epsi, ea, eb,edensity,
			A11,A12,A21,A22,TransScal);
  }
  fclose(pfile);

  printf(" Save the data ? (y/n): ");
  retour=scanf("%s",response);
  car=response[0];
  if((car=='O')||(car=='o')||(car=='y')||(car=='Y'))
    {
      printf(" data file name ? : ");
      retour=scanf("%s",filename);
      p->EcritSino(filename);
    }
  
  delete  [] filename;
  delete  [] response;
  delete  [] A11;
  delete  [] A12;
  delete  [] A21;
  delete  [] A22;
  delete  [] TransScal;
}

void GetSimulDynTranslatedData(Sinogramme *p,Scalaire *translation)
{
  char *filename;
  filename=new char[80];
  char *response,car;
  response=new char[80];
  printf(" phantom definition file name ? : ");
  int retour=scanf("%s",filename);
  printf(" Save the data with EcritSino ? (y/n): ");
  retour=scanf("%s",response);
  car=response[0];
  SimulDynTranslatedData(p, translation, filename,car,1); // 1 is for verbose
  delete  [] filename;
  delete  [] response;
}

void SimulDynTranslatedData(Sinogramme *p,Scalaire* translation, 
			    char *filename, char car, int verbose)
// see AddTranslatedEllipse
{
  int onemore; // logical 
  int i;
  FILE *pfile;
  int nbell,nbdisk,nbtri; // number of ellipsis, number of disks,
  //number of triangles (trinagles are not considered)
  Scalaire er, etheta, epsi, ea, eb, edensity; // ellipse parameters
  float fr, ftheta, fpsi, fa, fb, fdensity; // ellipse parameters
  int retour;
  // printf(" Manual or File (*/f) ? ");
  // retour=scanf("%s",response);
  // printf(" votre reponse %s \n",response);
  // if(response[0]=='f')
  // //response is file HERE (only file def are supported no manual)
  //    {
/* open the file for reading in */
  if ((pfile=fopen(filename,"r"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
  retour=fscanf(pfile,"%d",&nbell);
  retour=fscanf(pfile,"%d",&nbdisk);
  retour=fscanf(pfile,"%d",&nbtri);
  if(nbtri!=0)
    printf(" triangles are not considered ...\n");
  for(i=0;i<nbell;i++){
    retour=fscanf(pfile,"%f %f %f %f %f %f ",
		  &fr,&ftheta,&fa,&fb,&fpsi,&fdensity);
    er=fr; etheta=ftheta; epsi=fpsi; ea=fa; eb=fb; edensity=fdensity;
    if(verbose)
      printf("r, theta, psi, a, b, density = :%e %e %e %e %e %e\n ", 
	     er, etheta, epsi, ea, eb,edensity);
    //    if(mu!=0)
    //p->AddAtteEllipse(er, etheta, epsi, ea, eb, edensity,mu);
    //else
      //	    p->AddEllipse(er, etheta, epsi, ea, eb, edensity,3);
    p->AddTranslatedEllipse(er, etheta, epsi, ea, eb, edensity,translation);
  }
  for(i=0;i<nbdisk;i++){
    retour=fscanf(pfile,"%f %f %f %f ",
		  &fr,&ftheta,&fa,&fdensity);
    er=fr; etheta=ftheta; epsi=0.; ea=fa; eb=fa; edensity=fdensity;
    if(verbose)
      printf("r, theta, psi, a, b, density = :%e %e %e %e %e %e\n ", 
	     er, etheta, epsi, ea, eb,edensity);
    // if(mu!=0)
    //   p->AddAtteEllipse(er, etheta, epsi, ea, eb, edensity,mu);
    // 	  else
	    //	    p->AddEllipse(er, etheta, epsi, ea, eb, edensity,3);
    p->AddTranslatedEllipse(er, etheta, epsi, ea, eb, edensity,translation);
  }
  fclose(pfile);
  
  if((car=='O')||(car=='o')||(car=='y')||(car=='Y')){
    printf(" data file name ? : ");
    retour=scanf("%s",filename);
    p->EcritSino(filename);
  }
  
}

void SimulData(Sinogramme *p,Scalaire mu,
	       char *filename, int ecritsino, int verbose)
{
  // if ecrisino then save the data in a file with ecritsino
  
  int i;
  FILE *pfile;
  int retour;
  Scalaire er, etheta, epsi, ea, eb, edensity; // ellipse parameters
  float fr, ftheta, fpsi, fa, fb, fdensity; // ellipse parameters
  int nbell,nbdisk,nbtri; // number of ellipsis, number of disks,
  //number of triangles (trinagles are not considered)

/* open the file for reading in */
  if ((pfile=fopen(filename,"r"))==NULL)
    { printf("file %s non found \n",filename);
      exit(1);
    }
  retour=fscanf(pfile,"%d",&nbell);
  retour=fscanf(pfile,"%d",&nbdisk);
  retour=fscanf(pfile,"%d",&nbtri);
  if(nbtri!=0)
    printf(" triangles are not considered ...\n");
  for(i=0;i<nbell;i++)   {
    retour=fscanf(pfile,"%f %f %f %f %f %f ",
		  &fr,&ftheta,&fa,&fb,&fpsi,&fdensity);
    er=fr; etheta=ftheta; epsi=fpsi; ea=fa; eb=fb; edensity=fdensity;
    if(verbose==1)
      printf("r, theta, psi, a, b, density = :%e %e %e %e %e %e\n ", 
	     er, etheta, epsi, ea, eb,edensity);
    if(mu!=0)
      p->AddAtteEllipse(er, etheta, epsi, ea, eb, edensity,mu);
    else
      //	    p->AddEllipse(er, etheta, epsi, ea, eb, edensity,3);
      p->AddEllipse(er, etheta, epsi, ea, eb, edensity,1);//last parameter is oversampling
  }
  for(i=0;i<nbdisk;i++) {
    retour=fscanf(pfile,"%f %f %f %f ",
		  &fr,&ftheta,&fa,&fdensity);
    er=fr; etheta=ftheta; epsi=0.; ea=fa; eb=fa; edensity=fdensity;
    if(verbose==1)
      printf("r, theta, psi, a, b, density = :%e %e %e %e %e %e\n ", 
	     er, etheta, epsi, ea, eb,edensity);
    if(mu!=0)
      p->AddAtteEllipse(er, etheta, epsi, ea, eb, edensity,mu);
    else
      //	    p->AddEllipse(er, etheta, epsi, ea, eb, edensity,3);
      p->AddEllipse(er, etheta, epsi, ea, eb, edensity,1);//last parameter is oversampling
  }
  fclose(pfile);
  if(ecritsino==1){
    printf(" data file name ? : ");
    retour=scanf("%s",filename);
    p->EcritSino(filename);
  }
}      

void GetSimulData(Sinogramme *p,Scalaire mu)
{
  int onemore; // logical 
  int i;
  char *response,*filename,car;
  response=new char[80];
  filename=new char[80];
  Scalaire er, etheta, epsi, ea, eb, edensity; // ellipse parameters
  float fr, ftheta, fpsi, fa, fb, fdensity; // ellipse parameters
  int retour;

  printf(" Manual or File (*/f) ? ");
  retour=scanf("%s",response);
  printf(" votre reponse %s \n",response);
  if(response[0]=='f')
    {
      printf(" phantom definition file name ? : ");
      retour=scanf("%s",filename);
      // ICI ICI ICI
      SimulData(p, mu, filename, 0, 1);
    }
  else
    {
      onemore=1;
      while(onemore)
	{
	  GetEllParam(er,etheta, epsi, ea, eb, edensity);
	  if(mu!=0)
	    p->AddAtteEllipse(er, etheta, epsi, ea, eb, edensity,mu);
	  else
	    //	    p->AddEllipse(er, etheta, epsi, ea, eb, edensity,3);
	    p->AddEllipse(er, etheta, epsi, ea, eb, edensity,1);//last parameter is oversampling
	  //	  p->AddAtteEllipse(er, etheta, epsi, ea, eb, edensity,mu);
	  printf(" One more ellipsis (*/o) ? ");
	  retour=scanf("%s",response);
	  printf(" votre reponse %s \n",response);
	  if(response[0]!='o')
	    onemore=0;
	  
	}
    }
  printf(" Save the data ? (y/n): ");
  retour=scanf("%s",response);
  car=response[0];
  if((car=='O')||(car=='o')||(car=='y')||(car=='Y'))
    {
      printf(" data file name ? : ");
      retour=scanf("%s",filename);
      p->EcritSino(filename);
    }
  
  delete  [] filename;
  delete  [] response;
}

void GetEllParam(Scalaire &er,Scalaire &etheta, Scalaire &epsi, 
		 Scalaire &ea, Scalaire &eb, Scalaire &edensity)
{
  float dummy; // used for reading scalaire values...
  int retour;

  printf(" ellipsis parameters : ");
  printf(" r ?:");
  retour=scanf("%f", &dummy);
  er=(Scalaire) dummy;
  printf(" theta ?:");
  retour=scanf("%f", &dummy);
  etheta=(Scalaire) dummy;
  printf(" psi ?:");
  retour=scanf("%f", &dummy);
  epsi=(Scalaire) dummy;
  printf("  a ?:");
  retour=scanf("%f", &dummy);
  ea=(Scalaire) dummy;
  printf("  b ?:");
  retour=scanf("%f", &dummy);
  eb=(Scalaire) dummy;
  printf(" density  ?:");
  retour=scanf("%f", &dummy);
  edensity=(Scalaire) dummy;
  printf("r, theta, psi, a, b, density = :%e %e %e %e %e %e\n ", 
	 er, etheta, epsi, ea, eb,edensity);
}


