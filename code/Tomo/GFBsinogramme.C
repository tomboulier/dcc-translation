#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <iostream>
#include <fstream>

//#include <tomotypes.h>
//#include <sinogramme.h>
#include <GFBsinogramme.h>

using namespace std;

void GFBSinogramme::Alloue()
     // Alloue g,phi,x (prive)
{
  xsource= new Scalaire [p];
  ysource= new Scalaire [p];
  fananglewidths= new Scalaire [p];
  fananglebegins= new Scalaire [p];
  data=new Scalaire [p*q];
}

  //
  // Constructors et destructors
  // =============================
  //
GFBSinogramme::GFBSinogramme(int np,int nq, // np sources position of nq fans
		Scalaire *xs,Scalaire *ys, // np source positions
			     Scalaire *FAIwidths, Scalaire *FAbegs) // fan angles widths; fan angle begins
  // Allocation of GFBSinogramme of np source positions each of nq fans
  // equidistant sampling on FAbegs(i)+[0 ; FAIwidths(i)] for x, i=0,...,np-1 
  // i.e. fanangle_{i,j}= FAbegs(i)+j* FAIwidths(i)/(nq-1), j=0,..., nq-1 ; i=0,...,np-1
   //   fan angles = (angles between the mesured line and 
   //       the line throught the source vertex and the origine )
// belongs to the interval (-M_PI/2 ; M_PI/2)
//  NO NO NO +  circleangle+M_PI/2
   //      it is the angle of the corresponding parallel projection
{
  p=np;q=nq;
  Alloue();
  for(int i=0;i<p;i++){
    xsource[i]=xs[i];
    ysource[i]=ys[i];
    fananglewidths[i]=FAIwidths[i]; 
    fananglebegins[i]=FAbegs[i]; 
    for(int j=0;j<q;j++)
      data[i*q+j]=0.;
  }
}

GFBSinogramme::GFBSinogramme(int np,int nq, Scalaire circulartrajradius, 
			     Scalaire circleanglewidth, Scalaire circleanglebegin,
			     Scalaire FOVradius)
// by default a circular trajectory of the source on a circle of radius 
// circulartrajradius
// from circleanglebegin to circleanglebegin+circleanglewidth-hangle
// FOVradius is the radius of the Field Of View
// the Fan Angle is constant FAconstwidth= 2*asin(FOVradius/circulartrajradius), 
// sampled around the origine symetricaly from -FAconstwidth/2 to FAconstwidth/2
// fananglebegins[i]=-FAconstwidth/2 
// NO NO NO fananglebegins[i]=circleangle+M_PI/2-FAconstwidth/2 
{
  p=np;q=nq;
  Alloue();
  //  std::cout<<"after Alloue"<<std::endl;
  Scalaire hcircleangle=circleanglewidth/np;
  Scalaire circleangle=circleanglebegin;
  Scalaire FAconstwidth=2*asin(FOVradius/circulartrajradius);
  for(int i=0;i<p;i++,circleangle+=hcircleangle){
    //std::cout<<i<<" "<<std::flush;
    xsource[i]=circulartrajradius*cos(circleangle);
    ysource[i]=circulartrajradius*sin(circleangle);;
    fananglewidths[i]=FAconstwidth; 
    //    fananglebegins[i]=circleangle+M_PI/2-FAconstwidth/2; 
    fananglebegins[i]=-FAconstwidth/2; 
    for(int j=0;j<q;j++)
      data[i*q+j]=0.;
  }
}

GFBSinogramme::GFBSinogramme(const GFBSinogramme & gfbsino)
// constructor by copy
{
  int i,j;
  p=gfbsino.NP();
  q=gfbsino.NQ();
  Alloue();
  for (i=0;i<p;i++){
      xsource[i]= gfbsino.Xsource(i);
      ysource[i]= gfbsino.Ysource(i);
      fananglewidths[i]= gfbsino.FAwidth(i);
      fananglebegins[i]= gfbsino.FAbegin(i);
    for (j=0;j<q;j++)
      data[i*q+j]=gfbsino(i,j);
    }
}
GFBSinogramme::~GFBSinogramme()
// free the GFBSinogramme
{
  delete [] xsource;
  delete [] ysource;
  delete [] fananglewidths;
  delete [] fananglebegins;
  delete [] data;
}
//
// structure initialisation
//
void GFBSinogramme::Zero()
// set data[i*q+j] i=0...p-1, j=0...q-1, to 0
{
  for(int i=0;i<p;i++)
    for(int j=0;j<q;j++)
      data[i*q+j]=0.;
}

//
// structure access
//
Scalaire & GFBSinogramme::operator () (int i,int j)
// Reference coefficient (i,j)
{
  //#if CHECK
  if (i<0 || j<0 || i>=p || j>=q){
    //      fputs("Erreur : Adressage hors GFBSinogramme\n",stderr);
    std::cout<<"Erreur : operator () Adressage hors GFBSinogramme "<<std::flush;
    std::cout<<"i= "<<i<<" ; j="<<j<<std::endl;
    }
    //#endif
    return *( data + (q*i+j) );
}

const Scalaire & GFBSinogramme::operator () (const int i,const int j) const 
// Reference coefficient (i,j)
{
  //#if CHECK
  if (i<0 || j<0 || i>=p || j>=q){
      fputs("Erreur : Adressage hors GFBSinogramme\n",stderr);
    }
    //#endif
    return *( data + (q*i+j) );
}

int GFBSinogramme::NP() const
// Nbre de lignes
{
  return p;
}

int GFBSinogramme::NQ() const
// Nbre de colonnes
{
  return q;
}

Scalaire GFBSinogramme::Xsource(const int i) const
//return abscisse of the i^{th} source positions
{ 
  if (i<0 || i>=p) fputs("Err GFBSinogramme::Xsource : Adressage hors GFBSinogramme\n",stderr);
  return xsource[i];
}
Scalaire GFBSinogramme::Ysource(const int i) const
//return ordonnee of the i^{th} source positions
{  if (i<0 || i>=p) fputs("Err GFBSinogramme::Ysource : Adressage hors GFBSinogramme\n",stderr);
return ysource[i];
}

Scalaire GFBSinogramme::FAbegin(const int i) const
//return the i^{th} fan angle begin = fananglebegins[i]
{  if (i<0 || i>=p) fputs("Err GFBSinogramme::FAbegin : Adressage hors GFBSinogramme\n",stderr);
return fananglebegins[i];
}
Scalaire GFBSinogramme::FAwidth(const int i) const
//return the i^{th} fan angle width = fananglewidths[i]
{  if (i<0 || i>=p) fputs("Err GFBSinogramme::FAwidth : Adressage hors GFBSinogramme\n",stderr);
return fananglewidths[i];
}

void  GFBSinogramme::LitSino(char *filename)
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
    for (j=0; j<q;j++){
      retour=fscanf(pfile,"%e",&val);
      data[i*q+j] = val;
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

void GFBSinogramme::EcritSino(char *filename)
// Copy the Sinogramme in a file
{
  FILE *pfile;
  int i,j;
  //  int nbrot,nbtrans; /* nb de rotations et de translations */
  //  float val;
  
  /* ouverture du fichier en ecriture */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("ERR- file %s not found \n",filename);
      exit(1);
    }

  /* write on the first line : nbrotation   nbtranslations */
  fprintf(pfile,"%d %d\n",p,q);
  /* write on the following lines */
  for (i=0;i<p;i++)
    for (j=0; j<q;j++)
      fprintf(pfile,"%e ",data[i*q+j]);
  fclose(pfile);
}

void GFBSinogramme::PGMWrite(char *filename)
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
  max=data[0];
  min=data[0];
  for (i=0;i<rows;i++)
    for (j=0; j<cols;j++){
      val=data[i*cols+j];
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
	(unsigned char) ( scale * ( data[row*cols+col] - min) );
  fwrite(uc,1,rows*cols,pfile);
  fclose(pfile);
  delete [] uc;
}

void GFBSinogramme::Examine()
// cout all the content of a GFBSinogramme (except the data)
{
  std::cout<< " projection number = " << p << std::endl;
  std::cout<< " fan number per projection = " << q << std::endl;
  std::cout<< " sources abscisses: " << std::endl;
  for(int i=0;i<p;i++)
    std::cout<< xsource[i]<<" " << std::flush;
  std::cout<< std::endl << " sources ordo: " << std::endl;
  for(int i=0;i<p;i++)
    std::cout<< ysource[i]<<" " << std::flush;
  std::cout<< std::endl << " Fan Angle Widths: " << std::endl;
  for(int i=0;i<p;i++)
    std::cout<< fananglewidths[i]<<" " << std::flush;
  std::cout<< std::endl << " Fan Angle Begins: " << std::endl;
  for(int i=0;i<p;i++)
    std::cout<< fananglebegins[i]<<" " << std::flush;
  std::cout<< std::endl ;
}
//
// transformation
//
void GFBSinogramme::AddEllipse(Scalaire r, Scalaire theta, Scalaire psi,
			       Scalaire a, Scalaire b,Scalaire density)
//, int oversampling)
{
  // if(oversampling!=1){
  //   std::cout<<"Err FBAddEllipse :oversampling must be =1 !"<< std::endl;
  //   exit(0);
  // }
  Scalaire epsilon=0.0000001;
  int i,j;
  Scalaire phi,angle,circleangle,hfanangle,deca_s,decalage,c,d,e;
  Scalaire asquare,bsquare;
  Scalaire co,si,cosquare,sisquare;
  Scalaire tone,ttwo;

  // this is the same algorithm as for parallel beam
  // we just change the fan beam variable to parallel beam variable
  // with s= (xsource,ysource) \cdot \vtheta(\fanangle)
  //      s = xsource \cos(\fanangle) + ysource \sin(\fanangle)
  // and \phi=\fanangle
  asquare=a*a;
  bsquare=b*b;
  for (i=0;i<p;i++) {
    circleangle=acos(xsource[i]/sqrt(xsource[i]*xsource[i]+ysource[i]*ysource[i]));
    if(ysource[i]<0) circleangle=2*M_PI- circleangle;
    phi=fananglebegins[i]+M_PI/2+circleangle;
    // phi=fananglebegins[i] avant car fananglebegins[i] était circleangle+M_PI/2-FAconstwidth/2; 
    angle=psi-phi;
    hfanangle=fananglewidths[i]/(q-1);
    //     for (j=0;j<q;j++,phi+=hfanangle,angle-=hfanangle) {
    for (j=0;j<q;j++,phi+=hfanangle,angle=psi-phi) {
      decalage=r*cos(theta-phi);
      co=cos(angle);
      si=sin(angle);
      cosquare=co*co;
      sisquare=si*si;
      c=asquare*cosquare+bsquare*sisquare;
      // here phi is fanangle
      deca_s=xsource[i]*cos(phi)+ysource[i]*sin(phi)-decalage;
      if((deca_s*deca_s)<c)  {
	d=a*b*sqrt(c-(deca_s*deca_s));
	e=deca_s*si*co*(bsquare-asquare);
	tone=(d+e)/c;
	ttwo=(-d+e)/c;
	data[i*q+j]+=density*(tone-ttwo);
      } 
	// we know that e has not to be computed bute we want to keep
	// tone and ttwo for the attenuated transformation....
    }
  }
}

void GFBSinogramme::AddTranslatedEllipse(Scalaire r, Scalaire theta, Scalaire psi,
					 Scalaire a, Scalaire b,Scalaire density,
					 Scalaire* Xtranslation)
//, int oversampling)
{
  // if(oversampling!=1){
  //   std::cout<<"Err FBAddEllipse :oversampling must be =1 !"<< std::endl;
  //   exit(0);
  // }
  Scalaire epsilon=0.0000001;
  int i,j;
  Scalaire phi,angle,circleangle,hfanangle,deca_s,decalage,c,d,e;
  Scalaire asquare,bsquare;
  Scalaire co,si,cosquare,sisquare;
  Scalaire tone,ttwo;

  // this is the same algorithm as for parallel beam
  // we just change the fan beam variable to parallel beam variable
  // with s= (xsource,ysource) \cdot \vtheta(\fanangle)
  //      s = xsource \cos(\fanangle) + ysource \sin(\fanangle)
  // and \phi=\fanangle
  asquare=a*a;
  bsquare=b*b;
  for (i=0;i<p;i++) {
    circleangle=acos(xsource[i]/sqrt(xsource[i]*xsource[i]+ysource[i]*ysource[i]));
    if(ysource[i]<0) circleangle=2*M_PI- circleangle;
    phi=fananglebegins[i]+M_PI/2+circleangle;
    // phi=fananglebegins[i] avant car fananglebegins[i] était circleangle+M_PI/2-FAconstwidth/2; 
    angle=psi-phi;
    hfanangle=fananglewidths[i]/(q-1);
    //     for (j=0;j<q;j++,phi+=hfanangle,angle-=hfanangle) {
    for (j=0;j<q;j++,phi+=hfanangle,angle=psi-phi) {
      decalage=r*cos(theta-phi)+Xtranslation[i]*cos(phi); // WARNING in FB here 
      // the translation can not be pre-computed as in parallel as xtrans*cos(phi)
      // because cos(phi) depend on the fan j ans source position i...
      co=cos(angle);
      si=sin(angle);
      cosquare=co*co;
      sisquare=si*si;
      c=asquare*cosquare+bsquare*sisquare;
      // here phi is fanangle
      deca_s=xsource[i]*cos(phi)+ysource[i]*sin(phi)-decalage;
      if((deca_s*deca_s)<c)  {
	d=a*b*sqrt(c-(deca_s*deca_s));
	e=deca_s*si*co*(bsquare-asquare);
	tone=(d+e)/c;
	ttwo=(-d+e)/c;
	data[i*q+j]+=density*(tone-ttwo);
      } 
	// we know that e has not to be computed bute we want to keep
	// tone and ttwo for the attenuated transformation....
    }
  }
}

Scalaire GFBSinogramme::RolfRetroprojection(Scalaire xima, Scalaire yima,
					    Scalaire FOVradius,int n,
					    int UP)
// return Retroprojection of *this at (xima,yima)
// The GFBSinogramme *this is supposed to be a FBacquisition 
//    with a given centered FOV (of radius FOVradius)
// n is for the weight "tan(phi)^n/cos(phi)"
// // NOT IMPLEMENTED the back projection  angle interval is  anglemin,anglemax
// // Scalaire anglemin,Scalaire anglemax,
{
  int ip, iq;
  Scalaire hgamma,begingamma;
  // 
  //  exit(0);
  Scalaire val=0.;
  int out=0; // boolean : true if some pixels are outside the projection ranges
  Scalaire smpx, smpy, scal, scalortho, gamma, u;

  Scalaire  circleangle, dangle, mu, weight;

  Scalaire R=sqrt(xsource[0]*xsource[0]+ysource[0]*ysource[0]);
  Scalaire muL=acos(yima/R), muR=-muL;

  if ((xima*xima+yima*yima)<(FOVradius*FOVradius)){
    for (ip=0;ip<p;ip++){
      hgamma=fananglewidths[ip]/(q-1); 
      begingamma=fananglebegins[ip];
      smpx=xsource[ip]-xima;
      smpy=ysource[ip]-yima;
      scal=smpx*xsource[ip]+smpy*ysource[ip];
      scalortho=-smpx*ysource[ip]+smpy*xsource[ip];
      gamma = atan(scalortho/scal);
      //      iq=floor((gamma-begingamma)/hgamma);
      //      u=(gamma-begingamma)/hgamma-iq;
      u=(gamma-begingamma)/hgamma;
      iq=floor(u);
      u-=iq;

      // weight computation
      if((iq<0)||(iq>(q-2)))
	out=1;
      else{
	if(fabs(sqrt(xsource[ip]*xsource[ip]+ysource[ip]*ysource[ip])-R)>.01){
	  std::cout<<" ERR RolfRetroprojection is designed for circular trajectory ONLY " << std::endl;
	  exit(0);
	}
	if(ip<(p-1))
	  dangle=sqrt((xsource[ip+1]-xsource[ip])*(xsource[ip+1]-xsource[ip])
		      +(ysource[ip+1]-ysource[ip])*(ysource[ip+1]-ysource[ip]))
	    /R;
	else
	  dangle=sqrt((xsource[p-1]-xsource[p-2])*(xsource[p-1]-xsource[p-2])
		      +(ysource[p-1]-ysource[p-2])*(ysource[p-1]-ysource[p-2]))
	    /R;
	  
	if(n<0)
	  val+=((1.-u)*data[ip*q+iq]+u*data[ip*q+iq+1])*dangle;
	else{
	  circleangle=acos(xsource[ip]/R);
	  if(ysource[ip]<0) circleangle=2*M_PI- circleangle;
	  mu=circleangle-M_PI/2;
	  weight=pow(R*sin(mu)+xima,(double)n)/
	    pow(R*cos(mu)-yima,(double)(n+1));
	  weight*=R*(R+xima*sin(mu)-yima*cos(mu))/
	    sqrt(R*R+xima*xima+yima*yima+2*R*(xima*sin(mu)-yima*cos(mu)) );
	  // NEW
	  if( (mu>muR)&&(mu<muL) ){
	    if((UP==1)||(UP==2)) // LD 2 aout 2014 ajout de ||(UP==2)	
	  // this if is NEW
	      val+=(((1.-u)*data[ip*q+iq]+u*data[ip*q+iq+1])*weight)*dangle; 
	  }
	  else
	    //	    if(UP!=1) LD 2 aout 2014 modif	
	    if((UP==0)||(UP==2))	
	 // this if is NEW
	      val-=(((1.-u)*data[ip*q+iq]+u*data[ip*q+iq+1])*weight)*dangle; 
		       // ICI ICI ICI ICICI
	}
      }
    }
  }
  if(out)std::cout<<"GFBSinogramme::RolfRetroprojection : WARNING some pixels out of projection range "<<std::endl;
  return(val);
}

void GFBSinogramme::RolfRetroprojection(RealImage *image,Scalaire FOVradius,
					int n)
// same as the previous for all pixel points (xima,yima) of the image
// if n<0 simple backprojection else Rolf weighted backprojection (for DCC)
// Salaire anglemin,Scalaire anglemax,  NOT IMPLEMENTED
{
  int i,j,nx,ny; // for image(i,j),i=0,...,nx-1;j=0,...,ny-1
  Scalaire xbegin,xend,ybegin,yend; // describe the 4 image corners. 
  Scalaire xima,yima;
  nx=image->NX();
  ny=image->NY();
  image->GetCorners(&xbegin,&xend,&ybegin,&yend);
  Scalaire hx=(xend-xbegin)/nx,hy=(yend-ybegin)/ny;
  for (i=0;i<nx;i++) {
    std::cout<<"." <<std::flush;
    xima=xbegin+(i+.5)*hx;
    for (j=0;j<ny;j++){
      yima=ybegin+(j+.5)*hy;
      (*image)(i,j)+=this->RolfRetroprojection(xima,yima,FOVradius,n);
    }
  }
}
 
// Operators 
// ==========
//
void Add(GFBSinogramme *a,GFBSinogramme *b,GFBSinogramme *c)
  // c <- a+b
{
  for (int i=0;i<c->NP();i++)
    for (int j=0;j<c->NQ();j++)
      (*c)(i,j)=(*a)(i,j)+(*b)(i,j);
  // std::cout<< " WARNING WARNING WARNING WARNING WARNING WARNING "<<std::endl; 
  // std::cout<< " WARNING :Add(GFBSinogramme *a,GFBSinogramme *b,GFBSinogramme *c) must BE TERMINATED !!! " << std::endl; 
  // std::cout<< " WARNING WARNING WARNING WARNING WARNING WARNING "<<std::endl;
 }

void Sub(GFBSinogramme *a,GFBSinogramme *b,GFBSinogramme *c)
  // c <- a-b
{
  for (int i=0;i<c->NP();i++)
    for (int j=0;j<c->NQ();j++)
      (*c)(i,j)=(*a)(i,j)-(*b)(i,j);
  // std::cout<< " WARNING WARNING WARNING WARNING WARNING WARNING "<<std::endl; 
  // std::cout<< " WARNING :Add and Sub(GFBSinogramme *a,GFBSinogramme *b,GFBSinogramme *c) must BE TERMINATED !!! " << std::endl; 
  // std::cout<< " WARNING WARNING WARNING WARNING WARNING WARNING "<<std::endl; 

 }

void SimulData(GFBSinogramme *gfbsino,// Scalaire mu, pourrait être prévu
	       char *filename, int ecritsino, int verbose)
// if ecrisino then save the data in a file with ecritsino
{
  int i, retour;
  FILE *pfile;
  Scalaire er, etheta, epsi, ea, eb, edensity; // ellipse parameters
  float fr, ftheta, fpsi, fa, fb, fdensity; // ellipse parameters
  int nbell,nbdisk,nbtri; // number of ellipsis, number of disks,
  //number of triangles (trinagles are not considered)

  /* open the file for reading in */
  if ((pfile=fopen(filename,"r"))==NULL){ 
    printf("file %s non found \n",filename);
    exit(1);
  }
  retour=fscanf(pfile,"%d",&nbell);
  retour=fscanf(pfile,"%d",&nbdisk);
  retour=fscanf(pfile,"%d",&nbtri);
  if(nbtri!=0)
    printf(" triangles are not considered ...\n");
  for(i=0;i<nbell;i++) {
    retour=fscanf(pfile,"%f %f %f %f %f %f ",
		  &fr,&ftheta,&fa,&fb,&fpsi,&fdensity);
    er=fr; etheta=ftheta; epsi=fpsi; ea=fa; eb=fb; edensity=fdensity;
    if(verbose==1)
      printf("r, theta, psi, a, b, density = :%e %e %e %e %e %e\n ", 
	   er, etheta, epsi, ea, eb,edensity);
    gfbsino->AddEllipse(er, etheta, epsi, ea, eb, edensity);
  }
  for(i=0;i<nbdisk;i++){
    retour=fscanf(pfile,"%f %f %f %f ",
		  &fr,&ftheta,&fa,&fdensity);
    er=fr; etheta=ftheta; epsi=0.; ea=fa; eb=fa; edensity=fdensity;
    if(verbose==1)
      printf("r, theta, psi, a, b, density = :%e %e %e %e %e %e\n ", 
	   er, etheta, epsi, ea, eb,edensity);
    gfbsino->AddEllipse(er, etheta, epsi, ea, eb,edensity);
  }
  fclose(pfile);
  if(ecritsino==1){
    printf(" data file name ? : ");
    retour=scanf("%s",filename);
    gfbsino->EcritSino(filename);
  }
}

//
void GetSimulData(GFBSinogramme *gfbsino)
{
  int onemore; // logical 
  int i;
  char *response,*filename,car;
  response=new char[80];
  filename=new char[80];
  Scalaire er, etheta, epsi, ea, eb, edensity; // ellipse parameters
  float fr, ftheta, fpsi, fa, fb, fdensity; // ellipse parameters
  int retour ;

  printf(" Manual or File (*/f) ? ");
  retour=scanf("%s",response);
  printf(" votre reponse %s \n",response);
  if(response[0]=='f')    {
    printf(" phantom definition file name ? : ");
    retour=scanf("%s",filename);
      SimulData(gfbsino, filename, 0, 1); 
  }
  else  {
    onemore=1;
    while(onemore)	{
      GetEllParam(er,etheta, epsi, ea, eb, edensity);
      gfbsino->AddEllipse(er, etheta, epsi, ea, eb, edensity);
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
  if((car=='O')||(car=='o')||(car=='y')||(car=='Y')){
    printf(" data file name ? : ");
    retour=scanf("%s",filename);
    gfbsino->EcritSino(filename);
  }
  delete  [] filename;
  delete  [] response;
}

void GetSimulDynTranslatedData(GFBSinogramme *p,Scalaire *Xtranslation)
// Warning here Xtranslation contains the X translation of the center 
// for each source position 
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
  SimulDynTranslatedData(p, Xtranslation, filename,car,1); // 1 is for verbose
  delete  [] filename;
  delete  [] response;
}

void SimulDynTranslatedData(GFBSinogramme *p,Scalaire* Xtranslation, 
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
    p->AddTranslatedEllipse(er, etheta, epsi, ea, eb, edensity,Xtranslation);
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
    p->AddTranslatedEllipse(er, etheta, epsi, ea, eb, edensity,Xtranslation);
  }
  fclose(pfile);
  
  if((car=='O')||(car=='o')||(car=='y')||(car=='Y')){
    printf(" data file name ? : ");
    retour=scanf("%s",filename);
    p->EcritSino(filename);
  }
  
}
