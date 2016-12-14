// Mathematical tools 
//
// LD Jan 2014 : creation

// (c) Copyright TIMC 2014

#include <mathtools.h>



//using namespace std;

void Moments(double *M,int n, double *f, double *x, int p, int q )
// Compute the Moments of order "n" of the functions f(i,.)
// in the form M(i)= sum_j f_i(j)x_i(j)^n (but use a trapezoidal rule) 
// f_i(j) is supposed to be stored at f(i*q+j) j=0,...,q-1 i=0,...,p-1
// x_i(j) is supposed to be stored at x(i*q+j) j=0,...,q-1 i=0,...,p-1
// (for the Data Consitency Conditions of order n)
{
  double momn;
  for(int i=0;i<p;i++){
    momn=0.;
    //formule des trapèzes si échantillonnage irrégulier
    for(int j=1;j<q-1;j++)
      momn+=f[i*q+j]*pow(x[i*q+j],(double)n)*(x[i*q+j+1]-x[i*q+j-1])/2;
    momn+=f[i*q]*pow(x[i*q],(double)n)*(x[i*q+1]-x[i*q])/2;
    momn+=f[i*q+q-1]*pow(x[i*q+q-1],(double)n)*(x[i*q+q-1]-x[i*q+q-2])/2;
    M[i]=momn;
  }
}


Scalaire DistancePoly(Scalaire *t,Scalaire *f,int n, int degree, int verbose)
// measure the norme 2 of the difference of the fonction "f" to the polynomial
// set of degree "degree"
// f is described by it sampled graph (t[i], f[i]), i=0,...,n-1
// the function is given only for degree =0 , 1 , 2 or 3
// the t[i] are supposed to be equidistant
// more important t[i] are supposed to sample an even intervalle [-b,b]
{
  if(degree<0){
    std::cout<< " ERR in DistancePoly: degree must be >= 0 ; degree = " 
	     << degree << std::endl;
    exit(0);
  }
  if(degree>3){
    std::cout<< " ERR in DistancePoly: degree must be <= 2 ; degree = " 
	     << degree << std::endl;
    exit(0);
  }
  Scalaire mean;
  mean=0.;
  for(int i=0;i<n;i++)
    mean+=f[i];
  mean-=(f[0]+f[n-1])/2; // trapezoidal integral discretization
  mean/=(n-1); // the mean is now computed 

  Scalaire normeinf=0;
  Scalaire polyi=0;
  Scalaire normeone=0;
  switch (degree) {
  case 0 :{
    Scalaire distancezero=0.;
    normeinf=0;normeone=0;
    for(int i=0;i<n;i++){
      distancezero+=(f[i]-mean)*(f[i]-mean);
      if(fabs(f[i]-mean) > normeinf)normeinf=fabs(f[i]-mean);
      normeone+=fabs(f[i]-mean);
    }
    distancezero-=((f[0]-mean)*(f[0]-mean) + (f[n-1]-mean)*(f[n-1]-mean))/2.;
    distancezero/=(n-1);
    normeone-=(fabs(f[0]-mean)+fabs(f[n-1]-mean))/2.;
    normeone/=(n-1);
    if(verbose==1){
      std::cout<<" Projection sur PO (mean) : "<<mean<<std::endl; 
      std::cout<<" Norm Inf (Max_i |f(i)-p(i)|) = "<<normeinf <<std::endl; 
      std::cout<<" Norme One (1/(n-1)(Sum_{i=1,n-2) |f(i)-p(i)|+(|f(0)-p(0)|+|f(n-1)-p(n-1)|)/2) = "<<normeone <<std::endl; 
    }
    return(distancezero);
    break;}
  case 1 : {
    Scalaire tf=0,tt=0;
    for(int i=0;i<n;i++){
      tf+=f[i]*t[i];
      tt+=t[i]*t[i];
    }
    tf-=(t[0]*f[0]+t[n-1]*f[n-1])/2; // trapezoidal integral discretization
    tt-=(t[0]*t[0]+t[n-1]*t[n-1])/2; 
    Scalaire coefft=tf/tt;
    Scalaire distanceun=0.;
    normeinf=0;normeone=0;
    for(int i=0;i<n;i++){
      polyi=coefft*t[i]+mean;
      distanceun+=(f[i]-polyi)*(f[i]-polyi);
      if(fabs(f[i]-polyi) > normeinf)normeinf=fabs(f[i]-polyi);
      normeone+=fabs(f[i]-polyi);
    }
    normeone-=(fabs(f[0]-(coefft*t[0]+mean))+fabs(f[n-1]-(coefft*t[n-1]+mean)))/2.;
    normeone/=(n-1);

    if(verbose==1){
      std::cout<<" Projection sur P1"<<std::endl; 
      for(int i=0;i<n;i++)
	std::cout<<coefft*t[i]+mean<<" "<<std::flush;
      std::cout<<std::endl;
      std::cout<<" Polynome "<< coefft << " x + "<< mean <<std::endl; 
      std::cout<<" Norm Inf (Max_i |f(i)-p(i)|) = "<<normeinf <<std::endl; 
      std::cout<<" Norme One (1/(n-1)(Sum_{i=1,n-2) |f(i)-p(i)|+(|f(0)-p(0)|+|f(n-1)-p(n-1)|)/2) = "<<normeone <<std::endl; 

      std::cout<< "#############################################"<<std::endl; 
    }
    distanceun-=((f[0]-(coefft*t[0]+mean))*(f[0]-(coefft*t[0]+mean))
		+(f[n-1]-(coefft*t[n-1]+mean))*(f[n-1]-(coefft*t[n-1]+mean)))/2;
    distanceun/=(n-1);

    return(distanceun);
    break;    }
  case 2 : {
    Scalaire tdeux;
    Scalaire tf=0,tt=0, ttf=0,tttt=0;
    for(int i=0;i<n;i++){
      tf+=f[i]*t[i];
      tdeux=t[i]*t[i];
      tt+=tdeux;
      ttf+=tdeux*f[i];
      tttt+=tdeux*tdeux;
    }
    tf-=(t[0]*f[0]+t[n-1]*f[n-1])/2; // trapezoidal integral discretization
    tt-=(t[0]*t[0]+t[n-1]*t[n-1])/2; 
    ttf-=(t[0]*t[0]*f[0]+t[n-1]*t[n-1]*f[n-1])/2; 
    tttt-=(t[0]*t[0]*t[0]*t[0]+t[n-1]*t[n-1]*t[n-1]*t[n-1])/2; 
    Scalaire coefft=tf/tt;
    Scalaire coefftt=(ttf-mean*tt)/(tttt-tt*tt/(n-1));
    Scalaire distancedeux=0.,disti=0.;
    normeinf=0;normeone=0;
    for(int i=0;i<n;i++){
      polyi=coefftt*(t[i]*t[i]-tt/(n-1))+coefft*t[i]+mean;
      disti=f[i]-polyi;
      distancedeux+=disti*disti;
      if(fabs(disti) > normeinf)normeinf=fabs(disti);
      normeone+=fabs(disti);     
    }
    normeone-=(fabs(f[0]-(coefftt*(t[0]*t[0]-tt/(n-1))+coefft*t[0]+mean))
	       +fabs(f[n-1]-(coefftt*(t[n-1]*t[n-1]-tt/(n-1))+coefft*t[n-1]+mean)))/2.;
    normeone/=(n-1);
    if(verbose==1){
      std::cout<<" Projection sur P2"<<std::endl; 
      for(int i=0;i<n;i++)
	std::cout<<coefftt*(t[i]*t[i]-tt/(n-1))+coefft*t[i]+mean<<" "<<std::flush;
      std::cout<<std::endl;
      std::cout<<" Polynome "<< coefftt << " (x^2 - "<<tt/(n-1)<<") + "<< coefft << " x + "<< mean <<std::endl; 
      std::cout<<" <=> Polynome "<< coefftt << " x^2 + "<< coefft << " x + "<< 
	mean-coefftt*tt/(n-1) <<std::endl; 
      std::cout<<" Norm Inf (Max_i |f(i)-p(i)|) = "<<normeinf <<std::endl; 
      std::cout<<" Norme One (1/(n-1)(Sum_{i=1,n-2) |f(i)-p(i)|+(|f(0)-p(0)|+|f(n-1)-p(n-1)|)/2) = "<<normeone <<std::endl; 

      std::cout<< "#############################################"<<std::endl; 
    }
    disti=(f[0]-(coefftt*(t[0]*t[0]-tt/(n-1))+coefft*t[0]+mean));
    distancedeux-= disti*disti/2;   
    disti=(f[n-1]-(coefftt*(t[n-1]*t[n-1]-tt/(n-1))+coefft*t[n-1]+mean));
    distancedeux-= disti*disti/2;   
    distancedeux/=(n-1);
    return(distancedeux);
    break;    }
  case 3 : {
    Scalaire tdeux,ttrois;
    Scalaire tf=0,tt=0, ttf=0,tttt=0, tttf=0,tttttt=0;
    for(int i=0;i<n;i++){
      tf+=f[i]*t[i];
      tdeux=t[i]*t[i];
      ttrois=tdeux*t[i];
      tt+=tdeux;
      ttf+=tdeux*f[i];
      tttt+=tdeux*tdeux;
      tttf+=ttrois*f[i];
      tttttt+=ttrois*ttrois;
    }
    tf-=(t[0]*f[0]+t[n-1]*f[n-1])/2; // trapezoidal integral discretization
    tt-=(t[0]*t[0]+t[n-1]*t[n-1])/2; 
    ttf-=(t[0]*t[0]*f[0]+t[n-1]*t[n-1]*f[n-1])/2; 
    tttt-=(t[0]*t[0]*t[0]*t[0]+t[n-1]*t[n-1]*t[n-1]*t[n-1])/2; 
    tttf-=(t[0]*t[0]*t[0]*f[0]+t[n-1]*t[n-1]*t[n-1]*f[n-1])/2; 
    tttttt-=(t[0]*t[0]*t[0]*t[0]*t[0]*t[0]
	     +t[n-1]*t[n-1]*t[n-1]*t[n-1]*t[n-1]*t[n-1])/2; 
    Scalaire coefft=tf/tt;
    Scalaire coefftt=(ttf-mean*tt)/(tttt-tt*tt/(n-1));
    Scalaire coeffttt=(tttf-tf*tttt/tt)/(tttttt-tttt*tttt/tt);
    Scalaire distancetrois=0.,disti=0.;
    normeinf=0;normeone=0;
    for(int i=0;i<n;i++){
      polyi=coeffttt*(t[i]*t[i]*t[i]-(tttt/tt)*t[i])
	+ coefftt*(t[i]*t[i]-tt/(n-1))+coefft*t[i]+mean;
      disti=f[i]-polyi;
      distancetrois+=disti*disti;
      if(fabs(disti) > normeinf)normeinf=fabs(disti);
      normeone+=fabs(disti);
    }
    normeone-=(fabs(f[0]-(coeffttt*(t[0]*t[0]*t[0]-(tttt/tt)*t[0])
			  + coefftt*(t[0]*t[0]-tt/(n-1))+coefft*t[0]+mean))
	       +fabs(f[n-1]-(coeffttt*(t[n-1]*t[n-1]*t[n-1]-(tttt/tt)*t[n-1])
			  + coefftt*(t[n-1]*t[n-1]-tt/(n-1))+coefft*t[n-1]+mean)))/2.;
    normeone/=(n-1);
    if(verbose==1){
      std::cout<<" Projection sur P3"<<std::endl; 
      for(int i=0;i<n;i++)
	std::cout<<coeffttt*(t[i]*t[i]*t[i]-(tttt/tt)*t[i])
		    + coefftt*(t[i]*t[i]-tt/(n-1))+coefft*t[i]+mean<<" "<<std::flush;
      std::cout<<std::endl;
      std::cout<<" Polynome "<< coeffttt << " (x^3 - "<<tttt/tt<<" x) + "<< coefftt << " (x^2 - "<<tt/(n-1)<<") + "<< coefft << " x + "<< mean <<std::endl;
      std::cout<<" <=> Polynome "<<coeffttt << " x^3 + " << coefftt << " x^ + "
	       << (coefft-coeffttt*tttt/tt)  << " x + "<< 
	mean-coefftt*tt/(n-1) <<std::endl; 
      std::cout<<" Norm Inf (Max_i |f(i)-p(i)|) = "<<normeinf <<std::endl; 
      std::cout<<" Norme One (1/(n-1)(Sum_{i=1,n-2) |f(i)-p(i)|+(|f(0)-p(0)|+|f(n-1)-p(n-1)|)/2) = "<<normeone <<std::endl; 
      std::cout<< "##########################"<<std::endl; 
    }
    disti=(f[0]-(coeffttt*(t[0]*t[0]*t[0]-(tttt/tt)*t[0])
		 + coefftt*(t[0]*t[0]-tt/(n-1))+coefft*t[0]+mean));
    distancetrois-=disti*disti/2;
    disti=(f[n-1]-(coeffttt*(t[n-1]*t[n-1]*t[n-1]-(tttt/tt)*t[n-1])
		 + coefftt*(t[n-1]*t[n-1]-tt/(n-1))+coefft*t[n-1]+mean));
    distancetrois-=disti*disti/2;
    distancetrois/=(n-1);
    return(distancetrois);
    break;    }
  default :{
    std::cout<< " ERR in DistancePoly: degree must be >0 and <=3, degree = " 
	     << degree << std::endl;    }
  }
}
