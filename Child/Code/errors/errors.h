/****************************************************************\
**  errors.h
**
**  Header file for CHILD error-handling routines.
**
**  Created Dec. 97
**  $Id: errors.h,v 1.5 2002-04-11 10:53:25 arnaud Exp $
\****************************************************************/

#ifndef ERRORS_H
#define ERRORS_H

#include <stdlib.h>
#include <iostream.h>

void ReportFatalError( const char *errStr )
#ifdef __GNUC__
 __attribute__ ((noreturn))
#endif
;

#endif
