//-*-c++-*- 

/****************************************************************/
/**
**  @file errors.h
**  @brief Header file for CHILD error-handling routines.
**
**  Created Dec. 97
**  $Id: errors.h,v 1.11 2004-01-07 10:53:25 childcvs Exp $
*/
/****************************************************************/

#ifndef ERRORS_H
#define ERRORS_H

#include "../compiler.h"

void ReportFatalError( const char *errStr ) ATTRIBUTE_NORETURN;

void ReportWarning( const char *errstr );

#endif
