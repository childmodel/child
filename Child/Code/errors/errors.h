//-*-c++-*- 

/****************************************************************\
**  errors.h
**
**  Header file for CHILD error-handling routines.
**
**  Created Dec. 97
**  $Id: errors.h,v 1.6 2002-07-08 17:21:50 arnaud Exp $
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
