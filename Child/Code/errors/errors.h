//-*-c++-*- 

/****************************************************************/
/**
**  @file errors.h
**  @brief Header file for CHILD error-handling routines.
**
**  Created Dec. 97
**  $Id: errors.h,v 1.8 2003-01-17 17:30:26 childcvs Exp $
*/
/****************************************************************/

#ifndef ERRORS_H
#define ERRORS_H

void ReportFatalError( const char *errStr )
#ifdef __GNUC__
 __attribute__ ((noreturn))
#endif
;

#endif
