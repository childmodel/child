/******************************************************************\
**  errors.cpp
**
**  Error-handling routines.
**
**  Created Dec. 97 from earlier routine embedded in child.cpp
**  $Id: errors.cpp,v 1.2 1999-01-12 21:05:51 gtucker Exp $
\******************************************************************/

#include "errors.h"


/*****************************************************************************\
**
**  ReportFatalError:  This is an error-handling routine that prints the
**                     message errMsg then halts the program.
**
**      Parameters:     errMsg -- error message
**      Called by:
**      Created: 12/96 GT
**
\*****************************************************************************/
void ReportFatalError( char *errMsg )
{
  cout << errMsg <<endl;
  cout << "That was a fatal error, my friend!" <<endl;
  exit(1);
}



