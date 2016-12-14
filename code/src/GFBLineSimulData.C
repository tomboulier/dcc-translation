 /* 
Line trajectory given by Rolf for FB DCC
Use of the GFBSinogramme class data
LD mardi 17 décembre 2013, 15:25:09 (UTC+0100)
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <iostream> 

#include <sinogramme.h>
// #include <tomotypes.h>
// #include <realimage.h>
#include <GFBsinogramme.h>

using namespace std;

int main()
{
  int i;
  int nbtrans, nbrot; // for parallel geometry
  Scalaire recradius; // Centered disk radius of the reconstruction region
  Sinogramme *p; // parallel projection
  char *filename;
  
  std::cout<< " geometry :  "<<std::endl;
  std::cout<< "   rotation number ?:  "<<std::flush;
  std::cin>>nbrot;
  nbtrans=(int)(nbrot*2/M_PI);
  std::cout<<" translation (detector position) number   ("<< nbtrans<< ") ???: " << std::flush;
  std::cin>>nbtrans;

  std::cout<< "   Centered disk radius of the reconstruction region ?:  "<<std::flush;
  std::cin>>recradius;

  p= new Sinogramme(nbrot,nbtrans, -M_PI/2, M_PI, recradius);

  std::cout<< " ################################ "<<std::endl;
  std::cout<< " Generation of simulated data: "<<std::endl;
  std::cout<< " ################################ "<<std::endl;
  GetSimulData(p);
  std::cout<< " ################################ "<<std::endl;
  std::cout<< " ################################ "<<std::endl;
  
  filename=new char[80];
  std::cout<<" PGM sinogramme file name ???:  "<<std::flush;
  std::cin>>filename; 
  p->PGMWrite(filename);



  Scalaire phimin=-M_PI/2+M_PI/nbrot,phimax=M_PI/2-M_PI/nbrot;
  Scalaire xima,yima=5;
  int nxima=101;
  Scalaire ximadeb=-8, xlength=16,dx=xlength/(nxima-1),ximafin=ximadeb+xlength;
  int n=1;
  Scalaire vecxima [nxima];
  Scalaire DCC [nxima];

  for (xima=ximadeb,i=0; i<nxima; xima+=dx, i++){
    DCC[i]=RolfRetroprojection(p, xima, yima,
			       phimin,phimax,
			       n);
    vecxima[i]=xima;
  }

  std::cout<< "vecxima: " << std::endl; 
  for (int i=0; i<nxima;  i++)
    std::cout<< vecxima[i] << "  "<<std::flush;
  std::cout<< std::endl; 
  std::cout<< "DCC n= " << n << std::endl; 
  for (int i=0; i<nxima;  i++)
    std::cout<< DCC[i] << "  "<<std::flush;
  std::cout<< std::endl; 
  std::cout<< std::endl; 

  delete p;
  delete filename;
   
}
  
