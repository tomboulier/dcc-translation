// Operation sur les FBsinogrammes 
// en particulier retroprojection
//
// LD Dec 93

// (c) Copyright TIMC 1993


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sinogramme.h>
#include <FBsinogramme.h>
#include <iostream>
#include <fstream>

using namespace std;

// Constructeurs et destructeurs
// =============================
//
FBSinogramme::FBSinogramme(int np,int nq,Scalaire sr,Scalaire fa) : Sinogramme(np,nq,0.,2*M_PI,1.)
     // Alloue un Sinogramme a np rotations et nq translations
     // ratio = sr is the source radius (> 1 because the measured and reconstructed  
     //         region is the the disk of radius 1
     // equidistant sampling on [0 2*pi], for phi, ip=0,...,np-1
     //     the vertex path if sr*(cos(phi),sin(phi))
     // equidistant sampling on [-fa fa] for x, iq=1,...,nq 
     //      (fan angles between the mesured line and the line throught the source vertex and the origine)  
     // by default, the sources are on the whole circle (angle=2*M_PI)
{
  radius=sr;
  fanangle=fa;
  FBXdef();
}

FBSinogramme::FBSinogramme(int np,int nq,Scalaire angle,Scalaire sr,Scalaire fa) : Sinogramme(np,nq,0.,angle,1.)
     // Alloue un Sinogramme a np rotations et nq colonnes
     // equidistant sampling on [0 angle], for phi, ip=0,...,np-1
     //     the vertex path if sr*(cos(phi),sin(phi))
     // equidistant sampling on [-fa fa] for x, iq=1,...,nq 
     //      (fan angles between the mesured line and the line throught the source vertex and the origine)  
{
  radius=sr;
  fanangle=fa;
  FBXdef();
}

//
// structure initialisation
//
void FBSinogramme::FBXdef()
// Initialize the abscisse,i.e., the angle alpha, to the usual one in 
// analytic approach, see Natterer 86 fan-beam approach
//
{
  int i,j;
  for (i=0;i<p;i++)
    for (j=1;j<=q;j++)
      x[addresse(i,j)]=(Scalaire)(-fanangle/2+(Scalaire)(j-1)*
				  (Scalaire)fanangle/(Scalaire)(q-1));
}

//
// structure access
//
Scalaire FBSinogramme::FA()
// return the fan angle
{
return(fanangle);
}

Scalaire FBSinogramme::R()
// return the radius
{
return(radius);
}

void FBConvolution(Projection *u,FBSinogramme *a,FBSinogramme *b)
// b(i,.) <- u convol a(i,.) (for all i)
// The dimension of u must be 2*a->NQ()
// FBConvolution is for the FanBeam FBP (acosinus factor appears in the filter) 
{
  int i,j,l,n;
  Scalaire h,somme;
  Scalaire fangle;
#if CHECK
  if ( (BadDim(a,b)) || (BadDim2QmoinsUn(a,u)) )
    fputs("Err : Bad dimmenssions in Convolution\n",stderr);
#endif
  n=a->NQ();
  fangle=a->FA();
  if(n==1){
    puts("Error : Cannot filter with only 1 fan!");
    return;
  }
  else {
    h=fangle/(n-1);
  }

  for (i=0;i<a->NP();i++)  {
    for (j=1;j<=n;j++){
      somme=0.;
      for (l=1;l<=n;l++)
	somme+=(*u)(j-l+n)*(*a)(i,l)*(Scalaire)cos(-fangle/2+
						   (l-1)*(fangle/(n-1)));
      (*b)(i,j)=somme*h;
    }
  }	
}

void Convolution(Projection *u,FBSinogramme *a,FBSinogramme *b)
// For each vertex position (i.e. each Fan Beam projection)
// each Fan Beam projection is convolved by the same vector u
//  (*b(i,l)) = \int (*a(i,l)) * (*u)(j-l+n)
// the dimension of (*u) must be 2n-1 ;  (*u)(n ) is "u(0)"  (*u)(l)=u( (l-n)*halpha)
//      with halpha=fangle/(n-1) 
{
  int i,j,l;
  Scalaire somme;
  int n=a->NQ();
  if(n==1){puts("Error : Cannot filter with only 1 fan!");return;}  
  Scalaire fangle=a->FA();
  Scalaire h=fangle/(n-1);

#if CHECK
  if ( (BadDim(a,b)) || (BadDim2QmoinsUn(a,u)) )
    fputs("Err : Bad dimmenssions in Convolution\n",stderr);
#endif
  for (i=0;i<a->NP();i++)  
    for (j=1;j<=n;j++){
      somme=0.;
      for (l=1;l<=n;l++)
	somme+=(*u)(j-l+n)*(*a)(i,l);
      (*b)(i,j)=somme*h;
    }  
}
void FBSinogramme::FBAddEllipse(Scalaire r, Scalaire theta, Scalaire psi,
				Scalaire a, Scalaire b,Scalaire density, 
				int oversampling)
// Add the measurements of an allipsis in scalar Sinography...
{
  int i,j;
  Scalaire angle,deca_s,decalage,c,d,e;
  Scalaire asquare,bsquare;
  Scalaire co,si,cosquare,sisquare;
  Scalaire tone,ttwo;
  
  Scalaire deltaalpha=fanangle/((oversampling+1)*(q-1));
  Scalaire alpha, somme;
  Scalaire *weight;
  Scalaire totalweight;

  weight=new Scalaire[oversampling];
  /*  for(int k=0;k<(oversampling/2+1);k++)
    weight[k]=k;
  for(int k=(oversampling/2+1);k<oversampling;k++)
    weight[k]=oversampling-k;
  */
  for(int k=0;k<oversampling;k++){
    weight[k]=1.;
    //    std::cout<<"weihght["<<k<<"]"<<  weight[k]<<std::endl;
  }
  totalweight=0;
  for(int k=0;k<oversampling;k++)
    totalweight+=weight[k];
  //  std::cout<<"totalweihght "<< totalweight<<std::endl;

  asquare=a*a;
  bsquare=b*b;
  for (i=0;i<p;i++) {
    for (j=1;j<=q;j++) {
      alpha=x[addresse(i,j)]-deltaalpha*(oversampling/2);
      somme=0.;
      for(int k=0;k<oversampling;k++){
	angle=phi[i]+alpha-M_PI/2-psi;
	//phi=beta+alpha-pi/2 (but beta is in phi[],alpha is in x[])
	decalage=r*cos(phi[i]+alpha-M_PI/2-theta);
	co=cos(angle);
	si=sin(angle);
	cosquare=co*co;
	sisquare=si*si;
	deca_s=radius*sin(alpha)-decalage;
	//s=r sin(alpha) (alpha is in x[])
	c=asquare*cosquare+bsquare*sisquare;
	if((deca_s*deca_s)<c)    {
	  d=a*b*sqrt(c-(deca_s*deca_s));
	  e=deca_s*si*co*(bsquare-asquare);
	  tone=(d+e)/c;
	  ttwo=(-d+e)/c;
	  somme+=density*(tone-ttwo)*weight[k];
	  // we know that e has not to be computed bute we want to keep
	  // tone and ttwo for the attenuated transformation....
	  alpha+=deltaalpha;
	}
      }
      //	g[addresse(i,j)]+=(somme/oversampling);
      g[addresse(i,j)]+=(somme/totalweight);
    }
  }
}

void FBSinogramme::ConvolPorte(int porte)
// convolution par une porte de largeur porte
// porte MUST BE ODD
{
  int i,j,k;
  Scalaire somme;
  if(porte%2==0){
    std::cout<<"  Err in ConvolPorte (FBsinogramme.C) : porte muste be odd" << std::endl;
    exit(0);
  }
  for (i=0;i<p;i++) {
    for (j=1+porte/2;j<=q-porte/2;j++) {
      somme=0.;
      for(int k=0;k<porte;k++){
	somme+=g[addresse(i,j-porte/2+k)];
      }
      g[addresse(i,j)]=somme/porte;
    }
  }
}
void FBRetroprojection(FBSinogramme *t,RealImage *image)
// image <- Retroprojection t
// see F.Natterer 
// "The mathematics of Computerized Sinography" (Wiley 86) p.113
{
  int i,j,nx,ny; // for image(i,j),i=0,...,nx-1;j=0,...,ny-1
  int ip,iq,p,q; // for t(ip,iq),i=0,...,p-1;j=1,...,q
  Scalaire co,si,s; // cosinus, sinus, and s=<x,theta> (see Natterer).
  Scalaire xima,yima; // abs and ord of the current pixel.
  Scalaire sx,sy; // abs and ord of the source
  Scalaire smpx,smpy; // abs and ord of (the source - the current pixel)
  Scalaire nsmp,scal; // Norm**2 of source-pixel, dummy
  Scalaire xbegin,xend,ybegin,yend; // describe the 4 image corners. 
  Scalaire h;// in case of regularly sampled projections, h is
  //the angle distance between 2 successive fans.
  Scalaire radius,fangle;
  Scalaire gamma;
  Scalaire u;// proportion 
  int out; // boolean : true if some pixels are outside the projection ranges
  
  out=0;
  nx=image->NX();
  ny=image->NY();
  p=t->NP();
  q=t->NQ();
  radius=t->R();
  fangle=t->FA();
  h=fangle/(q-1);
  //printf("avant GetCorners\n");
  image->GetCorners(&xbegin,&xend,&ybegin,&yend);
  //printf("xbegin,xend,ybegin,yend %e %e %e %e \n",xbegin,xend,ybegin,yend);
  //printf("avant for i =0 to nx=%d\n",nx);
  for (i=0;i<nx;i++)
    {
      std::cout<<"." <<std::flush;
      xima=xbegin+(i+.5)*((xend-xbegin)/nx);
      for (j=0;j<ny;j++)
	{
	  (*image)(i,j)=0.;
	  //printf("j = %d \n", j);
	  yima=ybegin+(j+.5)*((yend-ybegin)/ny);
	  if ((xima*xima+yima*yima)<1)
	    for (ip=0;ip<t->NP();ip++)
	      {
		sx=radius*cos(t->A(ip));
		sy=radius*sin(t->A(ip));
		smpx=sx-xima;
		smpy=sy-yima;
		scal=smpx*sx+smpy*sy;
		nsmp=smpx*smpx+smpy*smpy;
		gamma = acos(scal/sqrt(nsmp*(sx*sx+sy*sy)));
		sx=radius*cos(t->A(ip)+M_PI/2);
		sy=radius*sin(t->A(ip)+M_PI/2);
		if((sx*xima+sy*yima)>0)
		  gamma=-gamma;
		//		iq=(int) (gamma/h);// LD 1er avril 2014
		iq=floor(gamma/h);// LD 1er avril 2014
		u=gamma/h-iq;
		iq=iq+q/2;
		if((iq<1)||(iq>q))
                   out=1;
		else
                  {
		    if (iq!=q)
		      (*image)(i,j)=(*image)(i,j)
			+((1.-u)*(*t)(ip,iq)+u*(*t)(ip,iq+1))/nsmp;
		    else
		      (*image)(i,j)=(*image)(i,j)+(*t)(ip,iq)/nsmp;
		  }
	      }
	  (*image)(i,j)=2.*radius*M_PI*(*image)(i,j)/t->NP();
	}
    }
  std::cout<<std::endl;

  if(out)
      printf("WARNING some pixels out of projection range ");
}

void FBLambdaRetroprojection(FBSinogramme *t,RealImage *image)
// image <- Retroprojection t 
// in the case of Lambda Tomography
// see F.Natterer 
// "The mathematics of Computerized Sinography" (Wiley 86) p.113
// see Faridani 92 SIAM J. Appl. Math.
{
  int i,j,nx,ny; // for image(i,j),i=0,...,nx-1;j=0,...,ny-1
  int ip,iq,p,q; // for t(ip,iq),i=0,...,p-1;j=1,...,q
  Scalaire co,si,s; // cosinus, sinus, and s=<x,theta> (see Natterer).
  Scalaire xima,yima; // abs and ord of the current pixel.
  Scalaire sx,sy; // abs and ord of the source
  Scalaire smpx,smpy; // abs and ord of (the source - the current pixel)
  Scalaire nsmp,nsmptroisdemi,scal; // Norm**2 of source-pixel, dummy
  Scalaire xbegin,xend,ybegin,yend; // describe the 4 image corners. 
  Scalaire h;// in case of regularly sampled projections, h is
  //the angle distance between 2 successive fans.
  Scalaire radius,fangle;
  Scalaire gamma;
  Scalaire u;// proportion 
  int out; // boolean : true if some pixels are outside the projection ranges
  
  out=0;  
  nx=image->NX();
  ny=image->NY();
  p=t->NP();
  q=t->NQ();
  radius=t->R();
  fangle=t->FA();
  h=fangle/(q-1);
  //printf("avant GetCorners\n");
  image->GetCorners(&xbegin,&xend,&ybegin,&yend);
  //printf("xbegin,xend,ybegin,yend %e %e %e %e \n",xbegin,xend,ybegin,yend);
  //printf("avant for i =0 to nx=%d\n",nx);
  for (i=0;i<nx;i++)
    {
      printf("i = %d \n", i);
      xima=xbegin+(i+.5)*((xend-xbegin)/nx);
      for (j=0;j<ny;j++)
	{
	  (*image)(i,j)=0.;
	  //printf("j = %d \n", j);
	  yima=ybegin+(j+.5)*((yend-ybegin)/ny);
	  if ((xima*xima+yima*yima)<1)
	    for (ip=0;ip<t->NP();ip++)
	      {
		sx=radius*cos(t->A(ip));
		sy=radius*sin(t->A(ip));
		smpx=sx-xima;
		smpy=sy-yima;
		scal=smpx*sx+smpy*sy;
		nsmp=smpx*smpx+smpy*smpy;
		nsmptroisdemi=sqrt(nsmp*nsmp*nsmp); 
		// cf |x-a|**3 for LambdaCT, Faridani 92, SIMA J. Appl Math
		gamma = acos(scal/sqrt(nsmp*(sx*sx+sy*sy)));
		sx=radius*cos(t->A(ip)+M_PI/2);
		sy=radius*sin(t->A(ip)+M_PI/2);
		if((sx*xima+sy*yima)>0)
		  gamma=-gamma;
		iq=(int) (gamma/h);
		u=gamma/h-iq;
		iq=iq+q/2;
		if((iq<1)||(iq>q))
		  out=1;
		else
		  {
		    if (iq!=q)
		      (*image)(i,j)=(*image)(i,j)
			+((1.-u)*(*t)(ip,iq)+u*(*t)(ip,iq+1))/nsmptroisdemi;
		    else
		      (*image)(i,j)=(*image)(i,j)+(*t)(ip,iq)/nsmptroisdemi;
		  }
	      }
	  (*image)(i,j)=radius/2.*(*image)(i,j)/t->NP();
	}
    }
  if(out)
      printf("WARNING some pixels out of projection range ");
}

void ConvertIeq(FBSinogramme & div, Sinogramme *para)
// convert the Fan Beam sinogramme FBSinogramme *div 
// into *para a parallel beam Sinogramme
// the projections are supposed to be equi sampled
{
  int i,j;
  Scalaire rayon=div.R();
  int nsource=div.NP() ;
  int ndetFB=div.NQ() ;
  //Scalaire angleslambda[nsource];
  //Scalaire anglesdetector[nsource*ndet];
  //div.Angles(angleslambda);
  //  Scalaire beginlambda=  angleslambda(0);
  // Scalaire endlambda=angleslambda(nsource-1);
  Scalaire beginlambda=  div.A(0);
 Scalaire endlambda=div.A(nsource-1);
  Scalaire hlambda=(endlambda-beginlambda)/(nsource-1);// WARNING we suppose equidistant sampling !!!
  Scalaire  lambdaone,lambdatwo,deltalambda;
  int ilambdaone,ilambdatwo;
  //  div.Detectors(anglesdetector);
  Scalaire alphamax=div.FA()/2;
  Scalaire alphamin=-alphamax;
  if( ((alphamax+alphamin)* (alphamax+alphamin ) )>0.000001){
    std::cout<<"Err(Sinogramme::ConvertIeq) alphamax!=alphamin ";
    exit(0); //  we suppose that alphamax=-\alphamin
  } 
  Scalaire halpha=(alphamax- alphamin)/(ndetFB-1);

  int p=para->NP() ;
  int q=para->NQ() ;
  Scalaire ssurr, acosssurr[q];// is s/rayon  and acos(s_j/rayon), j=0,...,q-1
  Scalaire s; Scalaire smin=para->X(0,1) ;  
  Scalaire ds=(para->X(0,q) - smin)/(q-1); 
  Scalaire alpha[q]; // contains \alpha(s_j) :
  Scalaire deltaalpha[q];  
  //  Scalaire  alphaone, alphatwo;
  int indexalphamin, indexalpha[q];// we suppose that 
  //the alpha can be pre computed here because we assume the same equidistant sampling for all projections
  s=smin;  indexalphamin=0;
  for (j=0;j<q;j++,s+=ds){
    ssurr=s/rayon;
    acosssurr[j]=acos(ssurr); // acos is in [0 \pi[
    alpha[j]=asin(ssurr); //asin return the principal value of the arc sine of x in radians; the return value is in the range [-pi/2, pi/2].
    indexalpha[j]= (int) floor((alpha[j]-alphamin)/halpha); 
    //    deltaalpha[j]=alpha[j] - (alphamin+ indexalpha[j]*halpha); on divise par halpha
    deltaalpha[j]=(alpha[j] - alphamin)/halpha- indexalpha[j];
    if(indexalpha[j]<0) indexalphamin++;
    //   if(ilambdaone>(-1) && ilambdaone< (nsource-1))
  }  
  int ndroite=0;Scalaire valone,valtwo;
  Scalaire phii;
  for(i=0;i<p;i++){
    phii=para->A(i);
    for (j= indexalphamin;j<q- indexalphamin ; j++){
      lambdaone= phii- acosssurr[j];
      ilambdaone= (int) (lambdaone-beginlambda)/hlambda;
      ndroite=0;valone=0;valtwo=0;
      if(ilambdaone>(-1) && ilambdaone< (nsource-1)){
	// interpolation bi-lineaire  [ilambdaone  ilambdaone+1] x  [indexalpha[j] indexalpha[j]+1] 
	deltalambda=(lambdaone-beginlambda)/hlambda-ilambdaone;
	valone=(1-deltalambda) * (
				deltaalpha[j] * div(ilambdaone,indexalpha[j]+1+1)
				+ (1-deltaalpha[j]) * div(ilambdaone,indexalpha[j]+1) )
	  + deltalambda * (
			   deltaalpha[j] * div(ilambdaone+1,indexalpha[j]+1+1)
			   + (1-deltaalpha[j]) * div(ilambdaone,indexalpha[j]+1)  );
	ndroite++;
      }
      lambdatwo= phii + acosssurr[j];
      ilambdatwo=  (int) (lambdatwo-beginlambda)/hlambda;
      if(ilambdatwo>(-1) && ilambdatwo< (nsource-1)){
	deltalambda=(lambdatwo-beginlambda)/hlambda-ilambdatwo; // below we use the symetry...
	valtwo=(1-deltalambda) * (  
				(1-deltaalpha[j])* div(ilambdatwo,(q-indexalpha[j])+1) //(q-indexalpha[j]-1)+1+1
				+ deltaalpha[j] * div(ilambdatwo,q-indexalpha[j]) )
	  + deltalambda * (
			   (1-deltaalpha[j])* div(ilambdatwo+1,q-indexalpha[j]+1)
			   + deltaalpha[j] * div(ilambdatwo,q-indexalpha[j])  );
	ndroite++;
      }
      if(ndroite!=0)(*para)(i,j+1)=(valtwo+valone)/ndroite;      
    }
  }
}


// interface subroutines 
void GetSimulData(FBSinogramme *p)
{
  int onemore; // logical 
  int i;
  char *response,*filename,car;
  FILE *pfile;
  int nbell,nbdisk,nbtri; // number of ellipsis, number of disks,
  //number of triangles (trinagles are not considered)
  response=new char[80];
  filename=new char[80];
  Scalaire er, etheta, epsi, ea, eb, edensity; // ellipse parameters
  float fr, ftheta, fpsi, fa, fb, fdensity; // ellipse parameters
  int retour ;

  printf(" Manual or File (*/f) ? ");
  retour=scanf("%s",response);
  printf(" votre reponse %s \n",response);
  if(response[0]=='f')
    {
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
      for(i=0;i<nbell;i++)
        {
          retour=fscanf(pfile,"%f %f %f %f %f %f ",
		 &fr,&ftheta,&fa,&fb,&fpsi,&fdensity);
	  er=fr; etheta=ftheta; epsi=fpsi; ea=fa; eb=fb; edensity=fdensity;
	  printf("r, theta, psi, a, b, density = :%e %e %e %e %e %e\n ", 
		 er, etheta, epsi, ea, eb,edensity);
	  p->FBAddEllipse(er, etheta, epsi, ea, eb, edensity,5);
        }
      for(i=0;i<nbdisk;i++)
        {
          retour=fscanf(pfile,"%f %f %f %f ",
		 &fr,&ftheta,&fa,&fdensity);
	  er=fr; etheta=ftheta; epsi=0.; ea=fa; eb=fa; edensity=fdensity;
	  printf("r, theta, psi, a, b, density = :%e %e %e %e %e %e\n ", 
		 er, etheta, epsi, ea, eb,edensity);
	  p->FBAddEllipse(er, etheta, epsi, ea, eb,edensity,5);
        }
      fclose(pfile);
    }
  else
    {
      onemore=1;
      while(onemore)
	{
	  GetEllParam(er,etheta, epsi, ea, eb, edensity);
	  p->FBAddEllipse(er, etheta, epsi, ea, eb, edensity,5);
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

