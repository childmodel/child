/******************************************************************/
/**
**  @file errors.cpp
**  @brief Error-handling routines.
**
**  Created Dec. 97 from earlier routine embedded in child.cpp
**  $Id: errors.cpp,v 1.9 2004-06-16 13:37:29 childcvs Exp $
*/
/******************************************************************/

#include "errors.h"

#include <stdlib.h>
#include <iostream>

#define CHILD_ABORT_ON_ERROR "CHILD_ABORT_ON_ERROR"
#define CHILD_ABORT_ON_WARNING "CHILD_ABORT_ON_WARNING"

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
  std::cout << errMsg <<std::endl;
  std::cout << "That was a fatal error, my friend!" <<std::endl;
  if (getenv(CHILD_ABORT_ON_ERROR) != NULL)
    abort();
  std::cout << "(Set \"" CHILD_ABORT_ON_ERROR "\" to generate a crash.)"
	    <<std::endl;
  exit(1);
}



/*****************************************************************************\
**
**  ReportWarning:  This is an error-handling routine that prints the
**                     message errMsg but does not halt the program unless
**                     the environment variable CHILD_ABORT_ON_WARNING is set.
**
**      Parameters:     errMsg -- error message
**      Called by:
**      Created: 9/03 SL
**
\*****************************************************************************/
void ReportWarning( const char *errMsg )
{
  std::cout << "WARNING: " << errMsg <<std::endl;
  if (getenv(CHILD_ABORT_ON_WARNING) != NULL)
    abort();
  std::cout <<
    "(Set \"" CHILD_ABORT_ON_WARNING "\" to generate a crash."
    "\n.e.g. \"env " CHILD_ABORT_ON_WARNING "=1 child ...\")"
	    <<std::endl;
}

