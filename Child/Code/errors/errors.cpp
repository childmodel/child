/******************************************************************/
/**
**  @file errors.cpp
**  @brief Error-handling routines.
**
**  Created Dec. 97 from earlier routine embedded in child.cpp
**  $Id: errors.cpp,v 1.5 2003-01-17 17:30:25 childcvs Exp $
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



