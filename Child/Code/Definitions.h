//-*-c++-*- 

/****************************************************************************/
/**
**   @file Definitions.h
**   @brief Header file for defining global constants
**
**   $Id: Definitions.h,v 1.39 2003-11-14 17:59:26 childcvs Exp $
*/
/****************************************************************************/

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

/** DEFINITIONS *************************************************************/
#define VERSION "2.2.0, July 2002"
#define TRUE 1
#define FALSE 0
typedef enum {      // method of grid construction
  kUniformMesh = 0,
  kPerturbedMesh = 1,
  kRandomMesh = 2
} tMeshType_t;
typedef enum { // type of open boundary
  kCornerOutlet = 0,
  kOpenSide = 1,
  kOppositeSidesOpen = 2,
  kAllSidesOpen = 3,
  kSpecifyOutlet = 4
} tOpenBoundary_t;
typedef enum { // type of boundary condition
  kClosedBoundary = 1,
  kOpenBoundary = 2,
  kNonBoundary = 0
} tBoundary_t;
typedef enum {
  kFlowAllowed = 1,
  kFlowNotAllowed = 0
} tEdgeBoundary_t;
#define kMeanderNode true
#define kNonMeanderNode false
#define RHO 1000.0      /* Density of water (kg/m^3) */
#define RHOSED 2650. /* density of sediment, [kg/m^3] */
#define GRAV 9.81       /* Gravitational acceleration, m/s^2 */
#define POROSITY 0.3        /* porosity of sediment on bed */
#define VISC .00000112      /* viscosity of water [m^2/s] */
#define SECPERYEAR 31557600.00  /*number of seconds in a year*/

// Macros
#define ROUND(x)    static_cast<int>((x)+0.5)
#define SIGN(x)     ( (x)>0 ? 1 : 0 )

#ifndef BOOL
# define BOOL(x) (x)
#endif
#if defined(__SUNPRO_CC)
# if __SUNPRO_CC==0x420 && !defined(ENUM_BOOL_DEFINED)
#  define ENUM_BOOL_DEFINED 1
typedef enum { false=0, true } bool;
#  undef BOOL
#  define BOOL(x) ((x)?true:false)
# endif
#endif

#define INT_TO_ENUM(type,x) static_cast<type>(x)
#if defined(__SUNPRO_CC)
# if __SUNPRO_CC==0x420
#  undef INT_TO_ENUM
#  define INT_TO_ENUM(type,x) (type)(x)
# endif
#endif

// file suffixes
#define SNODES ".nodes"
#define SEDGES ".edges"
#define STRI ".tri"
#define SZ ".z"
#define SRANDOM ".random"
#define SVAREA ".varea"

#define OPTREADINPUT_PREVIOUS 1

#endif
