/******************************************************************/
/**
**  @file errors.cpp
**  @brief Error-handling routines.
**
**  Created Dec. 97 from earlier routine embedded in child.cpp
**  $Id: errors.cpp,v 1.6 2003-02-10 12:21:02 childcvs Exp $
*/
/******************************************************************/

#include "errors.h"

#include <stdlib.h>
#if !defined(HAVE_NO_NAMESPACE)
# include <iostream>
using namespace std;
#else
# include <iostream.h>
#endif

#define CHILD_ABORT_ON_ERROR "CHILD_ABORT_ON_ERROR"

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
  if (getenv(CHILD_ABORT_ON_ERROR) != NULL)
    abort();
  cout << "(Set \"" CHILD_ABORT_ON_ERROR "\" to generate a crash.)" <<endl;
  exit(1);
}



