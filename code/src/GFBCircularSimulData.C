 /* 
test of the GFBSinogramme class data
LD mardi 17 décembre 2013, 15:25:09 (UTC+0100)
circular trajectory
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <iostream> 

// #include <sinogramme.h>
// #include <tomotypes.h>
// #include <realimage.h>
#include <mathtools.h>
#include <nr.h> // numerical receipes
#include <GFBsinogramme.h>


#define CHECK 0
#define PGMCHECK 1

int main()
{
  int i;
  int nbtrans, nbrot;
  Scalaire circleradius, FOVradius, circleangle, fanangle;
  char *filename;
  
  GFBSinogramme *p;
  Scalaire *a;

  filename=new char[80];
   
  std::cout<< " GFBSino in circular geometry :  "<<std::endl;
  std::cout<< " radius of source-detector circle ?:  "<<std::flush;
  std::cin>> circleradius;
  std::cout<< " Centered FOV radius  ?:  "<<std::flush;
  std::cin>> FOVradius;

  Scalaire circleanglebegin=0.;
  // for Rolf
  Scalaire yzero=5.;
  Scalaire muL=acos(yzero/circleradius),muR=-muL;
  // circleangle = M_PI/2+mu

  std::cout<<" angular begin of the source trajectory  : "<< M_PI/2+muR << " for Lambda_S Rolf test  "<<std::endl;
  std::cout<<" angular begin of the source trajectory  : "<< M_PI/2+muR+2*muL << " for Lambda_B Rolf test  "<<std::endl;
  std::cout<<" angular begin of the source trajectory ? : "<<std::flush;
  std::cin>>circleanglebegin;

//  fanangle=M_PI;
  fanangle=2*asin(FOVradius/circleradius);
  std::cout<<" angular interval of the source trajectory    : "<< M_PI+ fanangle << " for short scan "<<std::endl;
  std::cout<<" angular interval of the source trajectory    : "<< 2*M_PI << " for full scan "<<std::endl;
  std::cout<<" angular interval of the source trajectory    : "<< 2*muL << " for Lambda_S Rolf test  "<<std::endl;
  std::cout<<" angular interval of the source trajectory    : "<< 2*M_PI-2*muL << " for Lambda_B Rolf test  "<<std::endl;
  std::cout<<" angular interval of the source trajectory ? : "<<std::flush;
  std::cin>>circleangle;

  std::cout<< "   source position number ?:  "<<std::flush;
  std::cin>>nbrot;
  nbtrans=(int) (2.*(fanangle/M_PI)*(circleradius/FOVradius)*nbrot/4.); // nbtrans = 2q+1 cf Natterer
  if(2*(nbtrans/2)==nbtrans)
    nbtrans=nbtrans+1;
  std::cout<<" detector position number   ("<< nbtrans<< ") ???: " << std::flush;
  std::cin>>nbtrans;
  std::cout<< " ################################ "<<std::endl;
  std::cout<< " ########### geometry : ########### "<<std::endl;
  std::cout<< " radius of source-detector circle:  "<< circleradius <<std::endl;  
  std::cout<< " Centered FOV radius: "<< FOVradius <<std::flush;  
  std::cout<< " (fan angle=  "<< fanangle <<" ) "<<std::endl;
  std::cout<<" angular begin of the source trajectory: " <<  circleanglebegin <<std::endl;
  std::cout<<"  angular interval of the source trajectory : "<< circleangle <<std::endl;
  std::cout<< " source position number:  "<<nbrot <<std::endl;
  std::cout<<"  detector position number: "<< nbtrans << std::endl;
  std::cout<< " ################################ "<<std::endl;
  std::cout<< " ################################ "<<std::endl;

  //GFBSinogramme(int np,int nq, Scalaire radius, 
  //         Scalaire circleanglewidth=2*M_PI,Scalaire circleanglebegin=0.,
  //			     Scalaire FOVradius);
  p= new GFBSinogramme(nbrot,nbtrans,circleradius,circleangle,circleanglebegin,FOVradius);
 

  std::cout<< " ################################ "<<std::endl;
  std::cout<< "       Simulation of data: "<<std::endl;
  std::cout<< " ################################ "<<std::endl;
  GetSimulData(p);
  std::cout<< " ################################ "<<std::endl;
  std::cout<< " ################################ "<<std::endl;
  

  double noiselevel;
  std::cout<< " Noise level of the data (in %, 1 means 100% ; if <0 no noise) : "<<std::flush;
  std::cin>>noiselevel;
  if(noiselevel>0){
    int idum=(-13);
    for (int ip=0;ip<nbrot;ip++)
      for (int iq=0;iq<nbtrans;iq++)
	(*p)(ip,iq)+= noiselevel*(*p)(ip,iq)*NR::gasdev(idum);
  }


  //  p->Examine();
  std::cout<<" EcritSino sinogramme file name ???:  "<<std::flush;
  std::cin>>filename; 
  p->EcritSino(filename);
  std::cout<<" PGM sinogramme file name ???:  "<<std::flush;
  std::cin>>filename; 
  p->PGMWrite(filename);

  delete p;
  return(0);
}
  
