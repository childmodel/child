//-*-c++-*- 

/****************************************************************/
/**
**  @file errors.h
**  @brief Header file for CHILD error-handling routines.
**
**  Created Dec. 97
**  $Id: errors.h,v 1.10 2003-08-06 16:12:22 childcvs Exp $
*/
/****************************************************************/

#ifndef ERRORS_H
#define ERRORS_H

#include "../compiler.h"

void ReportFatalError( const char *errStr ) ATTRIBUTE_NORETURN;

#endif
