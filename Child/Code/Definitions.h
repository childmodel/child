//-*-c++-*- 

/****************************************************************************/
/**
**   @file Definitions.h
**   @brief Header file for defining global constants
**
**   $Id: Definitions.h,v 1.46 2004-04-27 11:09:02 childcvs Exp $
*/
/****************************************************************************/

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

/** DEFINITIONS *************************************************************/
#define VERSION "2.3.0, April 2004"
typedef enum {      // method of grid construction
  kUniformMesh = 0,
  kPerturbedMesh = 1,
  kRandomMesh = 2
} tMeshType_t;
typedef enum { // type of open boundary (used in tMesh::MakeMeshFromScratch)
  kCornerOutlet = 0,       // corner outlet (lower left)
  kOpenSide = 1,           // one open side (lower)
  kOppositeSidesOpen = 2,  // two opposite sides (upper and lower)
  kAllSidesOpen = 3,       // all sides
  kSpecifyOutlet = 4 ,     // specify outlet coordinates
  kAllSideClosed = -1      // all sides closed
} tOpenBoundary_t;

typedef enum { // type of boundary condition
  kClosedBoundary = 1,
  kOpenBoundary = 2,
  kNonBoundary = 0
} tBoundary_t;
inline const char* BoundName( tBoundary_t b ){
  switch(b){
  case kClosedBoundary:
    return "1-Closed";
  case kOpenBoundary:
    return "2-Open";
  case kNonBoundary:
    return "0-Non";
  }
}
inline tBoundary_t IntToBound( int b_ ){
  return static_cast<tBoundary_t>(b_);
}
inline int BoundToInt( tBoundary_t b ){
  return static_cast<int>(b);
}

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
#define SVEG ".veg"

#define OPTREADINPUT_PREVIOUS 1

#endif
