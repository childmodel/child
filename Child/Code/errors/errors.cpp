/******************************************************************\
**  errors.cpp
**
**  Error-handling routines.
**
**  Created Dec. 97 from earlier routine embedded in child.cpp
**  $Id: errors.cpp,v 1.3 2002-04-10 15:13:40 gtucker Exp $
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
void ReportFatalError( const char *errMsg )
{
  cout << errMsg <<endl;
  cout << "That was a fatal error, my friend!" <<endl;
  exit(1);
}



