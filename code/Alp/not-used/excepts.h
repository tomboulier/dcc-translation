// Linear algebra with C++
// exceptions
//
// This file is no more used in the librairy (June 1995)
//
// $Header: /home/bainvil/Modules/alp/RCS/excepts.h,v 1.2 1995/02/13 21:40:01 bainvil Exp bainvil $

#ifndef __EXCEPTS_H
#define __EXCEPTS_H

#ifdef USEEXCEPTIONS

#include <string.h>

// This class is thrown for every dimension
// error in the arguments of the functions of
// the library

class AlpDimensionError
    {
  protected:
    char *msg;			// Error message (function + diagnostic)
  public:
    AlpDimensionError(const char *error_function,
		      const char *error_diag="incompatible dimensions of the arguments")
	{
	msg=new char [strlen(error_function)+3+strlen(error_diag)+1];
	msg[0]=0;
	strcat(msg,error_function);
	strcat(msg," : ");
	strcat(msg,error_diag);
	}
    ~AlpDimensionError()
	{ delete [] msg; }
    const char *Msg()		// Get the message
	{ return msg; }
    };

// This class is thrown by the non-implemented functions
#define ALPNOTIMPLEMENTEDMSG "Function not implemented : "
class AlpNotImplemented
    {
  protected:
    char *msg;			// Error message (function)
  public:
    AlpNotImplemented(const char *error_function)
	{
	msg=new char [strlen(ALPNOTIMPLEMENTEDMSG)+strlen(error_function)+1];
	msg[0]=0;
	strcat(msg,ALPNOTIMPLEMENTEDMSG);
	strcat(msg,error_function);
	}
    ~AlpNotImplemented()
	{ delete [] msg; }
    const char *Msg()		// Get the message
	{ return msg; }
    };
#undef ALPNOTIMPLEMENTEDMSG

#endif // USEEXCEPTIONS
#endif
