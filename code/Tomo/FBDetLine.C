#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <iostream>
#include <fstream>
#include <FBDetLine.h>

using namespace std;

void FBDetLine::Alloue()
     // Alloue xsource,xdetector,xvector,ysource,ydetector,yvector
{
  xsource= new Scalaire [p];
  ysource= new Scalaire [p];
  xdetector= new Scalaire [p];
  ydetector= new Scalaire [p];
  xvector= new Scalaire [p];
  yvector= new Scalaire [p];
  data=new Scalaire [p*q];
}

  //
  // Constructors et destructors
  // =============================
  //
FBDetLine::FBDetLine(int np,int nq, // np source positions of nq detectors per source position 
		Scalaire *xs,Scalaire *ys, // np source positions
		Scalaire *xd,Scalaire *yd, // np detector positions
		Scalaire *xv,Scalaire *yv // np detector positions
		)
  // Allocation of FBDetLine of np source positions each of nq detectors on a line 
  // equidistant sampling detector_{ip,jq}= detector_ip + jq*vector_ip, jq=0,...,q-1
     //      (fan angles between the mesured line and the line throught the source vertex and the origine)  
{
  p=np;q=nq;
  Alloue();
  for(int i=0;i<p;i++){
    xsource[i]=xs[i];
    ysource[i]=ys[i];
    xdetector[i]=xd[i];
    ydetector[i]=yd[i];
    xvector[i]=xv[i];
    yvector[i]=yv[i];
    for(int j=0;j<q;j++)
      data[i*q+j]=0.;
  }
}

FBDetLine::FBDetLine(const FBDetLine & fbdl)
// constructor by copy
{
  int i,j;
  p=fbdl.NP();
  q=fbdl.NQ();
  Alloue();
  for (i=0;i<p;i++){
    xsource[i]= fbdl.Xsource(i);
    ysource[i]= fbdl.Ysource(i);
    xdetector[i]= fbdl.Xdetector(i);
    ydetector[i]= fbdl.Ydetector(i);
    xvector[i]= fbdl.Xvector(i);
    yvector[i]= fbdl.Yvector(i);
    for (j=0;j<q;j++)
      data[i*q+j]=fbdl(i,j);
  }
}

FBDetLine::~FBDetLine()   // free the FBDetLine
{
  delete [] xsource;
  delete [] ysource;
  delete [] xdetector;
  delete [] ydetector;
  delete [] xvector;
  delete [] yvector;
  delete [] data;
}

//
// structure initialisation
//
void FBDetLine::Zero()
// set data[i*q+j] i=0...p-1, j=0...q-1, to 0
{
  for(int i=0;i<p;i++)
    for(int j=0;j<q;j++)
      data[i*q+j]=0.;
}
//
// structure access
//
Scalaire & FBDetLine::operator () (int i,int j)
// Reference coefficient (i,j)
{
  //#if CHECK
  if (i<0 || j<0 || i>=p || j>=q){
    std::cout<<"Erreur : operator () Adressage hors FBDetLine "<<std::flush;
    std::cout<<"i= "<<i<<" ; j="<<j<<std::endl;
    }
    //#endif
    return *( data + (q*i+j) );
}

const Scalaire & FBDetLine::operator () (const int i,const int j) const 
// Reference coefficient (i,j)
{
  //#if CHECK
  if (i<0 || j<0 || i>=p || j>=q){
    std::cout<<"Err: operator () const => Adressage hors FBDetLine "<<std::flush;
    }
    //#endif
    return *( data + (q*i+j) );
}

int FBDetLine::NP() const
// Nbre de lignes
{return p;}

int FBDetLine::NQ() const
// Nbre de colonnes
{return q;}

Scalaire FBDetLine::Xsource(const int i) const
//return abscisse of the i^{th} source position
{ 
  if (i<0 || i>=p) fputs("Err FBDetLine::Xsource : Adressage hors FBDetLine\n",stderr);
  return xsource[i];
}
Scalaire FBDetLine::Ysource(const int i) const
//return ordonnee of the i^{th} source position
{  if (i<0 || i>=p) fputs("Err FBDetLine::Ysource : Adressage hors FBDetLine\n",stderr);
return ysource[i];
}

Scalaire FBDetLine::Xdetector(const int i) const
//return abscisse of the i^{th} ORIGINE detector position
{ 
  if (i<0 || i>=p) fputs("Err FBDetLine::Xdetector : Adressage hors FBDetLine\n",stderr);
  return xdetector[i];
}
Scalaire FBDetLine::Ydetector(const int i) const
//return ordonnee of the i^{th} ORIGINE detector position
{  if (i<0 || i>=p) fputs("Err FBDetLine::Ydetector : Adressage hors FBDetLine\n",stderr);
return ydetector[i];
}

Scalaire FBDetLine::Xvector(const int i) const
//return abscisse of the i^{th} vector direction
{ 
  if (i<0 || i>=p) fputs("Err FBDetLine::Xvector : Adressage hors FBDetLine\n",stderr);
  return xvector[i];
}
Scalaire FBDetLine::Yvector(const int i) const
//return ordonnee of the i^{th} vector direction
{  if (i<0 || i>=p) fputs("Err FBDetLine::Yvector : Adressage hors FBDetLine\n",stderr);
return yvector[i];
}

void  FBDetLine::LitSino(char *filename)
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
/* lecture des lignes suivantes */
  for (i=0;i<p;i++)
    for (j=0; j<q;j++){
      retour=fscanf(pfile,"%e",&val);
      data[i*q+j] = val;
    }
  fclose(pfile);
}

void FBDetLine::EcritSino(char *filename)
// Copy the Sinogramme in a file
{
  FILE *pfile;
  int i,j;
  /* ouverture du fichier en ecriture */
  if ((pfile=fopen(filename,"w"))==NULL)
    { printf("ERR- file %s not found \n",filename);
      exit(1);
    }
  /* write on the first line : p and q */
  fprintf(pfile,"%d %d\n",p,q);
  /* write data on the following lines */
  for (i=0;i<p;i++)
    for (j=0; j<q;j++)
      fprintf(pfile,"%e ",data[i*q+j]);
  fclose(pfile);
}

void FBDetLine::PGMWrite(char *filename)
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

void FBDetLine::Examine()
// cout all the content of a FBDetLine (except the data)
{
  std::cout<< " projection number = " << p << std::endl;
  std::cout<< " detector number per projection = " << q << std::endl;
  std::cout<< " sources abscisses: " << std::endl;
  for(int i=0;i<p;i++)
    std::cout<< xsource[i]<<" " << std::flush;
  std::cout<< std::endl << " sources ordo: " << std::endl;
  for(int i=0;i<p;i++)
    std::cout<< ysource[i]<<" " << std::flush;
  std::cout<< " origin detectors abscisses: " << std::endl;
  for(int i=0;i<p;i++)
    std::cout<< xdetector[i]<<" " << std::flush;
  std::cout<< std::endl << " origin detectors ordo: " << std::endl;
  for(int i=0;i<p;i++)
    std::cout<< ydetector[i]<<" " << std::flush;
  std::cout<< " direction vectors abscisses: " << std::endl;
  for(int i=0;i<p;i++)
    std::cout<< xvector[i]<<" " << std::flush;
  std::cout<< std::endl << " direction vectors ordo: " << std::endl;
  for(int i=0;i<p;i++)
    std::cout<< yvector[i]<<" " << std::flush;

  std::cout<< std::endl ;
}
//
// transformation
//
void  FBDetLine::AddEllipse(Scalaire r, Scalaire theta, Scalaire psi,
			    Scalaire a, Scalaire b,Scalaire density)
{
  int i,j;
  double xdet,ydet;// current detector position 
  double xsou,ysou;// current source position
  //double xzeta,yzeta;// current unit direction from the source to the detector
  double cosphi, sinphi;// cosphi=yzeta; sinphi=-xzeta
  double normezeta;
  Scalaire phi,angle,deca_s,decalage,c,d,e;
  Scalaire asquare,bsquare;
  Scalaire co,si,cosquare,sisquare;
  Scalaire tone,ttwo;
  // this is the same algorithm as for parallel beam
  // we just change the fan beam variables to parallel beam variable
  // with s= (xsource,ysource) \cdot \vtheta(\fanangle)
  //      s = xsource \cos(\fanangle) + ysource \sin(\fanangle)
  // and \phi=
  asquare=a*a;
  bsquare=b*b;
  for (i=0;i<p;i++) {
    xsou=xsource[i];
    ysou=ysource[i];
    xdet=xdetector[i];
    ydet=ydetector[i];
    for (j=0;j<q;j++,xdet+=xvector[i],ydet+=yvector[i]) {
      // phi depends on j (because depends on both source and detector pose
      sinphi=-xdet+xsou;//xzeta=xdet-xsou;
      cosphi=ydet-ysou;//yzeta=ydet-ysou;
      normezeta=sqrt(cosphi*cosphi+sinphi*sinphi);//=sqrt(xzeta*xzeta+yzeta*yzeta);
      sinphi/=normezeta;// xzeta/=normezeta;
      cosphi/=normezeta;// yzeta/=normezeta;
      phi=asin(sinphi);
      if(cosphi<0)phi=M_PI-phi;
      //      std::cout<<"phi,cosphi,sinphi: " <<phi<< ", "<<cosphi<< ", "<<sinphi<< " ;  "<<std::flush;
      angle=psi-phi;
      decalage=r*cos(theta-phi);
      co=cos(angle);
      si=sin(angle);
      cosquare=co*co;
      sisquare=si*si;
      c=asquare*cosquare+bsquare*sisquare;
      deca_s=xsource[i]*cosphi+ysource[i]*sinphi-decalage; 
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
    //    std::cout<< std::endl;
  }
  //  std::cout<< std::endl;
}

void FBDetLine::Add(FBDetLine & a, FBDetLine & b)
// (*this) <- a+b
{
  if( (p!=a.NP()) ||(p!=b.NP()) ||(q!=a.NQ()) ||(q!=b.NQ()) ){
    std::cout<<"ERR(FBDetLine::Add): dimension problems"<<std::endl;
    exit(0);
  }
  for (int i=0;i<p;i++)
    for (int j=0;j<q;j++)
	   data[i*q+j]=a(i,j)+b(i,j);
 }
void FBDetLine::Sub(FBDetLine & a, FBDetLine & b)
// (*this) <- a-b
{
  if( (p!=a.NP()) ||(p!=b.NP()) ||(q!=a.NQ()) ||(q!=b.NQ()) ){
    std::cout<<"ERR(FBDetLine::Add): dimension problems"<<std::endl;
    exit(0);
  }
  for (int i=0;i<p;i++)
    for (int j=0;j<q;j++)
	   data[i*q+j]=a(i,j)-b(i,j);
 }

void FBDetLine::GetSimulData()
// simulation of data (asking questions to users) 
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
      this->SimulData(filename, 0, 1); 
  }
  else  {
    onemore=1;
    while(onemore)	{
      GetEllParam(er,etheta, epsi, ea, eb, edensity);
      this->AddEllipse(er, etheta, epsi, ea, eb, edensity);
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
    this->EcritSino(filename);
  }
  delete  [] filename;
  delete  [] response;
}

void FBDetLine::SimulData(char *filename, int ecritsino, int verbose)
// silent simulation of data from the phantom definition in file "filename"
// if ecritsino == 1 if the FBDetLine must be saved (no more silent)
// if verbose ==1 then some trace of what is done  
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
    this->AddEllipse(er, etheta, epsi, ea, eb, edensity);
  }
  for(i=0;i<nbdisk;i++){
    retour=fscanf(pfile,"%f %f %f %f ",
		  &fr,&ftheta,&fa,&fdensity);
    er=fr; etheta=ftheta; epsi=0.; ea=fa; eb=fa; edensity=fdensity;
    if(verbose==1)
      printf("r, theta, psi, a, b, density = :%e %e %e %e %e %e\n ", 
	   er, etheta, epsi, ea, eb,edensity);
    this->AddEllipse(er, etheta, epsi, ea, eb,edensity);
  }
  fclose(pfile);
  if(ecritsino==1){
    printf(" data file name ? : ");
    retour=scanf("%s",filename);
    this->EcritSino(filename);
  }
}
