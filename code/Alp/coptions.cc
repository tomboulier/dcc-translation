// Linear algebra with C++
//
// $Header: /home/bainvil/Modules/alp/RCS/coptions.C,v 2.1 1995/02/18 14:10:10 bainvil Exp bainvil $

#include <alp.h>
#include <iostream>

// using namespace std;

std::ostream& AlpLibraryOptions(std::ostream& o)
{
  o<<"ALP library compilation options"<<std::endl;
  o<<"Revision $Revision: 2.1 $"<<std::endl;
  o<<"Architecture        : ";
#ifdef i486_dos_ARCH
  o<<"i486/dos"<<std::endl;
#endif
#ifdef i486_linux_ARCH
  o<<"i486/linux"<<std::endl;
#endif
#ifdef alpha_osf_ARCH
  o<<"DECalpha/osf"<<std::endl;
#endif
#ifdef mips_ultrix_ARCH
  o<<"DECmips/ultrix"<<std::endl;
#endif
  return o;
}
