//-*-c++-*- 

/****************************************************************\
**  errors.h
**
**  Header file for CHILD error-handling routines.
**
**  Created Dec. 97
**  $Id: errors.h,v 1.7 2002-09-23 12:11:48 arnaud Exp $
\****************************************************************/

#ifndef ERRORS_H
#define ERRORS_H

void ReportFatalError( const char *errStr )
#ifdef __GNUC__
 __attribute__ ((noreturn))
#endif
;

#endif
