//-*-c++-*- 

/*********************************************************************/
/**
**  @file mathutil.h
**  @brief Header file for special math utilities not in math.h.
**         All or most routines from Numerical Recipes in C by
**         Press et al.
**
**  $Id: mathutil.h,v 1.7 2003-07-15 17:24:58 childcvs Exp $
*/
/*********************************************************************/

#ifndef MATHUTIL_H
#define MATHUTIL_H

#define PI 3.14159265358979323846

double ran3( long *idum );

void init_genrand(unsigned long s);
void init_by_array(unsigned long init_key[], int key_length);
unsigned long genrand_int32(void);
long genrand_int31(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_res53(void);

#endif
