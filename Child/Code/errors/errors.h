//-*-c++-*- 

/****************************************************************/
/**
**  @file errors.h
**  @brief Header file for CHILD error-handling routines.
**
**  Created Dec. 97
**  $Id: errors.h,v 1.9 2003-01-28 13:56:34 childcvs Exp $
*/
/****************************************************************/

#ifndef ERRORS_H
#define ERRORS_H

#ifndef ATTRIBUTE_NORETURN
# if defined(__GNUC__) || defined(__INTEL_COMPILER)
#  define ATTRIBUTE_NORETURN __attribute__ ((noreturn))
# else
#  define ATTRIBUTE_NORETURN
# endif
#endif

void ReportFatalError( const char *errStr ) ATTRIBUTE_NORETURN;

#endif
