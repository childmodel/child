/*****************************************************************************/
/**
**   @file predicates.cpp
**   @brief Routines for Arbitrary Precision Floating-point Arithmetic
**          and Fast Robust Geometric Predicates
*/
//  These functions are added as members of the class Predicates; should
//  ideally be a namespace, but some compilers, including Sun's, do not
//  support namespaces yet.
//  Presently, these functions are called from some member functions of
//  tGrid to check for line segment intersection because inexact arithmetic
//  can lead to the wrong answer.
//  Note that the group of functions in class Predicates is a subset of
//  the original "predicates" package, except that I've added two new
//  functions, DifferenceOfProductsOfDifferences(...) and
//  AdaptDiffOfProdsOfDiffs(...) to do segment intersection detection.
//  --Stephen Lancaster, 1/99
//  $Id: predicates.cpp,v 1.13 2004-06-16 13:32:06 childcvs Exp $
/*****************************************************************************/
/*                                                                           */ 
/*  Routines for Arbitrary Precision Floating-point Arithmetic               */ 
/*  and Fast Robust Geometric Predicates                                     */ 
/*  (predicates.c)                                                           */ 
/*                                                                           */ 
/*  May 18, 1996                                                             */ 
/*                                                                           */ 
/*  Placed in the public domain by                                           */ 
/*  Jonathan Richard Shewchuk                                                */ 
/*  School of Computer Science                                               */ 
/*  Carnegie Mellon University                                               */ 
/*  5000 Forbes Avenue                                                       */ 
/*  Pittsburgh, Pennsylvania  15213-3891                                     */ 
/*  jrs@cs.cmu.edu                                                           */ 
/*                                                                           */ 
/*  This file contains C implementation of algorithms for exact addition     */ 
/*    and multiplication of floating-point numbers, and predicates for       */ 
/*    robustly performing the orientation and incircle tests used in         */ 
/*    computational geometry.  The algorithms and underlying theory are      */ 
/*    described in Jonathan Richard Shewchuk.  "Adaptive Precision Floating- */ 
/*    Point Arithmetic and Fast Robust Geometric Predicates."  Technical     */ 
/*    Report CMU-CS-96-140, School of Computer Science, Carnegie Mellon      */ 
/*    University, Pittsburgh, Pennsylvania, May 1996.  (Submitted to         */ 
/*    Discrete & Computational Geometry.)                                    */ 
/*                                                                           */ 
/*  This file, the paper listed above, and other information are available   */ 
/*    from the Web page http://www.cs.cmu.edu/~quake/robust.html .           */ 
/*                                                                           */ 
/*****************************************************************************/ 
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  Using this code:                                                         */ 
/*                                                                           */ 
/*  First, read the short or long version of the paper (from the Web page    */ 
/*    above).                                                                */ 
/*                                                                           */ 
/*  Be sure to call exactinit() once, before calling any of the arithmetic   */ 
/*    functions or geometric predicates.  Also be sure to turn on the        */ 
/*    optimizer when compiling this file.                                    */ 
/*                                                                           */ 
/*                                                                           */ 
/*  Several geometric predicates are defined.  Their parameters are all      */ 
/*    points.  Each point is an array of two or three floating-point         */ 
/*    numbers.  The geometric predicates, described in the papers, are       */ 
/*                                                                           */ 
/*    orient2d(pa, pb, pc)                                                   */ 
/*    orient2dfast(pa, pb, pc)                                               */ 
/*    orient3d(pa, pb, pc, pd)                                               */ 
/*    orient3dfast(pa, pb, pc, pd)                                           */ 
/*    incircle(pa, pb, pc, pd)                                               */ 
/*    incirclefast(pa, pb, pc, pd)                                           */ 
/*    insphere(pa, pb, pc, pd, pe)                                           */ 
/*    inspherefast(pa, pb, pc, pd, pe)                                       */ 
/*                                                                           */ 
/*  Those with suffix "fast" are approximate, non-robust versions.  Those    */ 
/*    without the suffix are adaptive precision, robust versions.  There     */ 
/*    are also versions with the suffices "exact" and "slow", which are      */ 
/*    non-adaptive, exact arithmetic versions, which I use only for timings  */ 
/*    in my arithmetic papers.                                               */ 
/*                                                                           */ 
/*                                                                           */ 
/*  An expansion is represented by an array of floating-point numbers,       */ 
/*    sorted from smallest to largest magnitude (possibly with interspersed  */ 
/*    zeros).  The length of each expansion is stored as a separate integer, */ 
/*    and each arithmetic function returns an integer which is the length    */ 
/*    of the expansion it created.                                           */ 
/*                                                                           */ 
/*  Several arithmetic functions are defined.  Their parameters are          */ 
/*                                                                           */ 
/*    e, f           Input expansions                                        */ 
/*    elen, flen     Lengths of input expansions (must be >= 1)              */ 
/*    h              Output expansion                                        */ 
/*    b              Input scalar                                            */ 
/*                                                                           */ 
/*  The arithmetic functions are                                             */ 
/*                                                                           */ 
/*    grow_expansion(elen, e, b, h)                                          */ 
/*    grow_expansion_zeroelim(elen, e, b, h)                                 */ 
/*    expansion_sum(elen, e, flen, f, h)                                     */ 
/*    expansion_sum_zeroelim1(elen, e, flen, f, h)                           */ 
/*    expansion_sum_zeroelim2(elen, e, flen, f, h)                           */ 
/*    fast_expansion_sum(elen, e, flen, f, h)                                */ 
/*    fast_expansion_sum_zeroelim(elen, e, flen, f, h)                       */ 
/*    linear_expansion_sum(elen, e, flen, f, h)                              */ 
/*    linear_expansion_sum_zeroelim(elen, e, flen, f, h)                     */ 
/*    scale_expansion(elen, e, b, h)                                         */ 
/*    scale_expansion_zeroelim(elen, e, b, h)                                */ 
/*    compress(elen, e, h)                                                   */ 
/*                                                                           */ 
/*  All of these are described in the long version of the paper; some are    */ 
/*    described in the short version.  All return an integer that is the     */ 
/*    length of h.  Those with suffix _zeroelim perform zero elimination,    */ 
/*    and are recommended over their counterparts.  The procedure            */ 
/*    fast_expansion_sum_zeroelim() (or linear_expansion_sum_zeroelim() on   */ 
/*    processors that do not use the round-to-even tiebreaking rule) is      */ 
/*    recommended over expansion_sum_zeroelim().  Each procedure has a       */ 
/*    little note next to it (in the code below) that tells you whether or   */ 
/*    not the output expansion may be the same array as one of the input     */ 
/*    expansions.                                                            */ 
/*                                                                           */ 
/*                                                                           */ 
/*  If you look around below, you'll also find macros for a bunch of         */ 
/*    simple unrolled arithmetic operations, and procedures for printing     */ 
/*    expansions (commented out because they don't work with all C           */ 
/*    compilers) and for generating random floating-point numbers whose      */ 
/*    significand bits are all random.  Most of the macros have undocumented */ 
/*    requirements that certain of their parameters should not be the same   */ 
/*    variable; for safety, better to make sure all the parameters are       */ 
/*    distinct variables.  Feel free to send email to jrs@cs.cmu.edu if you  */ 
/*    have questions.                                                        */ 
/*                                                                           */ 
/*****************************************************************************/ 

#include "predicates.h"

// The algorithms below fail on processors using extended precision.
// Consequently, on Intel x86, we set the control word of the x87 device
// in double precision mode.
// The present implementation works on x86 for Linux, Cygwin, FreeBSD and
// Win32.
#if defined(i386) && (defined(linux) || defined(__CYGWIN__))
# if defined(linux)
#  include <fpu_control.h>
# elif defined(__CYGWIN__)
// _FPU_EXTENDED is as well the mask of the precision control
#  define _FPU_EXTENDED 0x300
#  define _FPU_DOUBLE   0x200
#  if defined(_lint)
#   define _FPU_GETCW(cw) (cw) = 0
#   define _FPU_SETCW(cw)
#  else
#   define _FPU_GETCW(cw) __asm__ __volatile__("fnstcw %0" : "=m" (*&cw))
#   define _FPU_SETCW(cw) __asm__ __volatile__("fldcw %0" : : "m" (*&cw))
#  endif
#  define fpu_control_t unsigned int
# else
#  error Platform not supported
# endif
class Set_local_fpu_precision {
  fpu_control_t oldcw_pc;
public:
  Set_local_fpu_precision() {
    fpu_control_t cw;
    _FPU_GETCW(cw);
    oldcw_pc = cw & _FPU_EXTENDED; // Extract precision mode
    cw = (cw & ~_FPU_EXTENDED) | _FPU_DOUBLE; // Set double precision mode
    _FPU_SETCW(cw);
  }
  ~Set_local_fpu_precision() {
    fpu_control_t cw;
    _FPU_GETCW(cw);
    cw = (cw & ~_FPU_EXTENDED) | oldcw_pc;
    _FPU_SETCW(cw); // Restore control word precision
  }
};
# define SET_DOUBLE_PRECISION_MODE Set_local_fpu_precision _local_precision

#elif defined(i386) && defined(__FreeBSD__)
# include <floatingpoint.h>
class Set_local_fpu_precision {
  fp_prec_t oldcw_pc;
public:
  Set_local_fpu_precision() {
    oldcw_pc = fpgetprec(); // Extract precision mode
    fpsetprec(FP_PD);  // Set double precision mode
  }
  ~Set_local_fpu_precision() {
    fpsetprec(oldcw_pc); // Restore control word precision
  }
};
# define SET_DOUBLE_PRECISION_MODE Set_local_fpu_precision _local_precision

#elif defined(_WIN32)
# include <float.h>
class Set_local_fpu_precision {
  unsigned int oldcw_pc;
public:
  Set_local_fpu_precision() {
    oldcw_pc = _control87(0,0) & MCW_PC; // Extract precision mode
    _control87(PC_53,MCW_PC); // Set double precision mode
  }
  ~Set_local_fpu_precision() {
    _control87(oldcw_pc,MCW_PC); // Restore control word precision
  }
};
# define SET_DOUBLE_PRECISION_MODE Set_local_fpu_precision _local_precision
#else
# define SET_DOUBLE_PRECISION_MODE
#endif

// constructor; just calls exactinit() (SL, 10/98):
Predicates::Predicates() 
{
   SET_DOUBLE_PRECISION_MODE;
   exactinit();
}
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  exactinit()   Initialize the variables used for exact arithmetic.        */ 
/*                                                                           */ 
/*  `epsilon' is the largest power of two such that 1.0 + epsilon = 1.0 in   */ 
/*  floating-point arithmetic.  `epsilon' bounds the relative roundoff       */ 
/*  error.  It is used for floating-point error analysis.                    */ 
/*                                                                           */ 
/*  `splitter' is used to split floating-point numbers into two half-        */ 
/*  length significands for exact multiplication.                            */ 
/*                                                                           */ 
/*  I imagine that a highly optimizing compiler might be too smart for its   */ 
/*  own good, and somehow cause this routine to fail, if it pretends that    */ 
/*  floating-point arithmetic is too much like real arithmetic.              */ 
/*                                                                           */ 
/*  Don't change this routine unless you fully understand it.                */ 
/*                                                                           */ 
/*****************************************************************************/ 
// with Predicates a class, this is basically the constructor;
// rather than explicitly make it so, call it from the constructor.
void Predicates::exactinit() 
{ 
  REAL half; 
  REAL check, lastcheck; 
  bool every_other; 
 
  every_other = true; 
  half = 0.5; 
  epsilon = 1.0; // Predicates private data member
  splitter = 1.0;  // Predicates private data member
  check = 1.0; 
  /* Repeatedly divide `epsilon' by two until it is too small to add to    */ 
  /*   one without causing roundoff.  (Also check if the sum is equal to   */ 
  /*   the previous sum, for machines that round up instead of using exact */ 
  /*   rounding.  Not that this library will work on such machines anyway. */ 
  do { 
    lastcheck = check; 
    epsilon *= half; 
    if (every_other) { 
      splitter *= 2.0; 
    } 
    every_other = !every_other; 
    check = 1.0 + epsilon; 
  } while ((check != 1.0) && (check != lastcheck)); 
  splitter += 1.0; 
 
  /* Error bounds for orientation and incircle tests. */
// Predicates private data members:
  resulterrbound = (3.0 + 8.0 * epsilon) * epsilon; 
  ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon; 
  ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon; 
  ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon; 
  o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon; 
  o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon; 
  o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon; 
  iccerrboundA = (10.0 + 96.0 * epsilon) * epsilon; 
  iccerrboundB = (4.0 + 48.0 * epsilon) * epsilon; 
  iccerrboundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon; 
  isperrboundA = (16.0 + 224.0 * epsilon) * epsilon; 
  isperrboundB = (5.0 + 72.0 * epsilon) * epsilon; 
  isperrboundC = (71.0 + 1408.0 * epsilon) * epsilon * epsilon; 
} 
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  grow_expansion()   Add a scalar to an expansion.                         */ 
/*                                                                           */ 
/*  Sets h = e + b.  See the long version of my paper for details.           */ 
/*                                                                           */ 
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */ 
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */ 
/*  properties as well.  (That is, if e has one of these properties, so      */ 
/*  will h.)                                                                 */ 
/*                                                                           */ 
/*****************************************************************************/ 
 
int Predicates::grow_expansion(int elen, const REAL* e, REAL b, REAL* h)
  /* e and h can be the same. */ 
{ 
  REAL Q; 
  INEXACT REAL Qnew; 
  int eindex; 
  REAL enow; 
  INEXACT REAL bvirt; 
  REAL avirt, bround, around; 
 
  Q = b; 
  for (eindex = 0; eindex < elen; eindex++) { 
    enow = e[eindex]; 
    Two_Sum(Q, enow, Qnew, h[eindex]); 
    Q = Qnew; 
  } 
  h[eindex] = Q; 
  return eindex + 1; 
} 
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  grow_expansion_zeroelim()   Add a scalar to an expansion, eliminating    */ 
/*                              zero components from the output expansion.   */ 
/*                                                                           */ 
/*  Sets h = e + b.  See the long version of my paper for details.           */ 
/*                                                                           */ 
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */ 
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */ 
/*  properties as well.  (That is, if e has one of these properties, so      */ 
/*  will h.)                                                                 */ 
/*                                                                           */ 
/*****************************************************************************/ 
 
int Predicates::grow_expansion_zeroelim(int elen, const REAL* e, REAL b, REAL* h)
  /* e and h can be the same. */ 
{ 
  REAL Q, hh; 
  INEXACT REAL Qnew; 
  int eindex, hindex; 
  REAL enow; 
  INEXACT REAL bvirt; 
  REAL avirt, bround, around; 
 
  hindex = 0; 
  Q = b; 
  for (eindex = 0; eindex < elen; eindex++) { 
    enow = e[eindex]; 
    Two_Sum(Q, enow, Qnew, hh); 
    Q = Qnew; 
    if (hh != 0.0) { 
      h[hindex++] = hh; 
    } 
  } 
  if ((Q != 0.0) || (hindex == 0)) { 
    h[hindex++] = Q; 
  } 
  return hindex; 
} 
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  expansion_sum()   Sum two expansions.                                    */ 
/*                                                                           */ 
/*  Sets h = e + f.  See the long version of my paper for details.           */ 
/*                                                                           */ 
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */ 
/*  with IEEE 754), maintains the nonadjacent property as well.  (That is,   */ 
/*  if e has one of these properties, so will h.)  Does NOT maintain the     */ 
/*  strongly nonoverlapping property.                                        */ 
/*                                                                           */ 
/*****************************************************************************/ 
 
int Predicates::expansion_sum(int elen, const REAL* e, int flen, const REAL* f, REAL* h)
  /* e and h can be the same, but f and h cannot. */ 
{ 
  REAL Q; 
  INEXACT REAL Qnew; 
  int findex, hindex, hlast; 
  REAL hnow; 
  INEXACT REAL bvirt; 
  REAL avirt, bround, around; 
 
  Q = f[0]; 
  for (hindex = 0; hindex < elen; hindex++) { 
    hnow = e[hindex]; 
    Two_Sum(Q, hnow, Qnew, h[hindex]); 
    Q = Qnew; 
  } 
  h[hindex] = Q; 
  hlast = hindex; 
  for (findex = 1; findex < flen; findex++) { 
    Q = f[findex]; 
    for (hindex = findex; hindex <= hlast; hindex++) { 
      hnow = h[hindex]; 
      Two_Sum(Q, hnow, Qnew, h[hindex]); 
      Q = Qnew; 
    } 
    h[++hlast] = Q; 
  } 
  return hlast + 1; 
} 
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  expansion_sum_zeroelim1()   Sum two expansions, eliminating zero         */ 
/*                              components from the output expansion.        */ 
/*                                                                           */ 
/*  Sets h = e + f.  See the long version of my paper for details.           */ 
/*                                                                           */ 
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */ 
/*  with IEEE 754), maintains the nonadjacent property as well.  (That is,   */ 
/*  if e has one of these properties, so will h.)  Does NOT maintain the     */ 
/*  strongly nonoverlapping property.                                        */ 
/*                                                                           */ 
/*****************************************************************************/ 
 
int Predicates::expansion_sum_zeroelim1(int elen, const REAL* e, int flen,
                                        const REAL* f, REAL* h) 
/* e and h can be the same, but f and h cannot. */ 
{ 
  REAL Q; 
  INEXACT REAL Qnew; 
  int index, findex, hindex, hlast; 
  REAL hnow; 
  INEXACT REAL bvirt; 
  REAL avirt, bround, around; 
 
  Q = f[0]; 
  for (hindex = 0; hindex < elen; hindex++) { 
    hnow = e[hindex]; 
    Two_Sum(Q, hnow, Qnew, h[hindex]); 
    Q = Qnew; 
  } 
  h[hindex] = Q; 
  hlast = hindex; 
  for (findex = 1; findex < flen; findex++) { 
    Q = f[findex]; 
    for (hindex = findex; hindex <= hlast; hindex++) { 
      hnow = h[hindex]; 
      Two_Sum(Q, hnow, Qnew, h[hindex]); 
      Q = Qnew; 
    } 
    h[++hlast] = Q; 
  } 
  hindex = -1; 
  for (index = 0; index <= hlast; index++) { 
    hnow = h[index]; 
    if (hnow != 0.0) { 
      h[++hindex] = hnow; 
    } 
  } 
  if (hindex == -1) { 
    return 1; 
  } else { 
    return hindex + 1; 
  } 
} 
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  expansion_sum_zeroelim2()   Sum two expansions, eliminating zero         */ 
/*                              components from the output expansion.        */ 
/*                                                                           */ 
/*  Sets h = e + f.  See the long version of my paper for details.           */ 
/*                                                                           */ 
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */ 
/*  with IEEE 754), maintains the nonadjacent property as well.  (That is,   */ 
/*  if e has one of these properties, so will h.)  Does NOT maintain the     */ 
/*  strongly nonoverlapping property.                                        */ 
/*                                                                           */ 
/*****************************************************************************/ 
 
int Predicates::expansion_sum_zeroelim2(int elen, const REAL* e,
                                        int flen, const REAL* f, REAL* h) 
/* e and h can be the same, but f and h cannot. */ 
{ 
  REAL Q, hh; 
  INEXACT REAL Qnew; 
  int eindex, findex, hindex, hlast; 
  REAL enow; 
  INEXACT REAL bvirt; 
  REAL avirt, bround, around; 
 
  hindex = 0; 
  Q = f[0]; 
  for (eindex = 0; eindex < elen; eindex++) { 
    enow = e[eindex]; 
    Two_Sum(Q, enow, Qnew, hh); 
    Q = Qnew; 
    if (hh != 0.0) { 
      h[hindex++] = hh; 
    } 
  } 
  h[hindex] = Q; 
  hlast = hindex; 
  for (findex = 1; findex < flen; findex++) { 
    hindex = 0; 
    Q = f[findex]; 
    for (eindex = 0; eindex <= hlast; eindex++) { 
      enow = h[eindex]; 
      Two_Sum(Q, enow, Qnew, hh); 
      Q = Qnew; 
      if (hh != 0) { 
        h[hindex++] = hh; 
      } 
    } 
    h[hindex] = Q; 
    hlast = hindex; 
  } 
  return hlast + 1; 
} 
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  fast_expansion_sum()   Sum two expansions.                               */ 
/*                                                                           */ 
/*  Sets h = e + f.  See the long version of my paper for details.           */ 
/*                                                                           */ 
/*  If round-to-even is used (as with IEEE 754), maintains the strongly      */ 
/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   */ 
/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      */ 
/*  properties.                                                              */ 
/*                                                                           */ 
/*****************************************************************************/ 
 
int Predicates::fast_expansion_sum(int elen, const REAL* e,
                                   int flen, const REAL* f, REAL* h)
/* h cannot be e or f. */ 
{ 
  REAL Q; 
  INEXACT REAL Qnew; 
  INEXACT REAL bvirt; 
  REAL avirt, bround, around; 
  int eindex, findex, hindex; 
  REAL enow, fnow; 
 
  enow = e[0]; 
  fnow = f[0]; 
  eindex = findex = 0; 
  if ((fnow > enow) == (fnow > -enow)) { 
    Q = enow; 
    enow = e[++eindex]; 
  } else { 
    Q = fnow; 
    fnow = f[++findex]; 
  } 
  hindex = 0; 
  if ((eindex < elen) && (findex < flen)) { 
    if ((fnow > enow) == (fnow > -enow)) { 
      Fast_Two_Sum(enow, Q, Qnew, h[0]); 
      enow = e[++eindex]; 
    } else { 
      Fast_Two_Sum(fnow, Q, Qnew, h[0]); 
      fnow = f[++findex]; 
    } 
    Q = Qnew; 
    hindex = 1; 
    while ((eindex < elen) && (findex < flen)) { 
      if ((fnow > enow) == (fnow > -enow)) { 
        Two_Sum(Q, enow, Qnew, h[hindex]); 
        enow = e[++eindex]; 
      } else { 
        Two_Sum(Q, fnow, Qnew, h[hindex]); 
        fnow = f[++findex]; 
      } 
      Q = Qnew; 
      hindex++; 
    } 
  } 
  while (eindex < elen) { 
    Two_Sum(Q, enow, Qnew, h[hindex]); 
    enow = e[++eindex]; 
    Q = Qnew; 
    hindex++; 
  } 
  while (findex < flen) { 
    Two_Sum(Q, fnow, Qnew, h[hindex]); 
    fnow = f[++findex]; 
    Q = Qnew; 
    hindex++; 
  } 
  h[hindex] = Q; 
  return hindex + 1; 
} 
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  fast_expansion_sum_zeroelim()   Sum two expansions, eliminating zero     */ 
/*                                  components from the output expansion.    */ 
/*                                                                           */ 
/*  Sets h = e + f.  See the long version of my paper for details.           */ 
/*                                                                           */ 
/*  If round-to-even is used (as with IEEE 754), maintains the strongly      */ 
/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   */ 
/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      */ 
/*  properties.                                                              */ 
/*                                                                           */ 
/*****************************************************************************/ 
 
int Predicates::fast_expansion_sum_zeroelim(int elen, const REAL* e,
                                int flen, const REAL* f, REAL* h)
/* h cannot be e or f. */ 
{ 
  REAL Q; 
  INEXACT REAL Qnew; 
  INEXACT REAL hh; 
  INEXACT REAL bvirt; 
  REAL avirt, bround, around; 
  int eindex, findex, hindex; 
  REAL enow, fnow; 
 
  enow = e[0]; 
  fnow = f[0]; 
  eindex = findex = 0; 
  if ((fnow > enow) == (fnow > -enow)) { 
    Q = enow; 
    enow = e[++eindex]; 
  } else { 
    Q = fnow; 
    fnow = f[++findex]; 
  } 
  hindex = 0; 
  if ((eindex < elen) && (findex < flen)) { 
    if ((fnow > enow) == (fnow > -enow)) { 
      Fast_Two_Sum(enow, Q, Qnew, hh); 
      enow = e[++eindex]; 
    } else { 
      Fast_Two_Sum(fnow, Q, Qnew, hh); 
      fnow = f[++findex]; 
    } 
    Q = Qnew; 
    if (hh != 0.0) { 
      h[hindex++] = hh; 
    } 
    while ((eindex < elen) && (findex < flen)) { 
      if ((fnow > enow) == (fnow > -enow)) { 
        Two_Sum(Q, enow, Qnew, hh); 
        enow = e[++eindex]; 
      } else { 
        Two_Sum(Q, fnow, Qnew, hh); 
        fnow = f[++findex]; 
      } 
      Q = Qnew; 
      if (hh != 0.0) { 
        h[hindex++] = hh; 
      } 
    } 
  } 
  while (eindex < elen) { 
    Two_Sum(Q, enow, Qnew, hh); 
    enow = e[++eindex]; 
    Q = Qnew; 
    if (hh != 0.0) { 
      h[hindex++] = hh; 
    } 
  } 
  while (findex < flen) { 
    Two_Sum(Q, fnow, Qnew, hh); 
    fnow = f[++findex]; 
    Q = Qnew; 
    if (hh != 0.0) { 
      h[hindex++] = hh; 
    } 
  } 
  if ((Q != 0.0) || (hindex == 0)) { 
    h[hindex++] = Q; 
  } 
  return hindex; 
} 
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  linear_expansion_sum()   Sum two expansions.                             */ 
/*                                                                           */ 
/*  Sets h = e + f.  See either version of my paper for details.             */ 
/*                                                                           */ 
/*  Maintains the nonoverlapping property.  (That is, if e is                */ 
/*  nonoverlapping, h will be also.)                                         */ 
/*                                                                           */ 
/*****************************************************************************/ 
 
int Predicates::linear_expansion_sum(int elen, const REAL* e,
                                     int flen, const REAL* f, REAL* h)
/* h cannot be e or f. */ 
{ 
  REAL Q, q; 
  INEXACT REAL Qnew; 
  INEXACT REAL R; 
  INEXACT REAL bvirt; 
  REAL avirt, bround, around; 
  int eindex, findex, hindex; 
  REAL enow, fnow; 
  REAL g0; 
 
  enow = e[0]; 
  fnow = f[0]; 
  eindex = findex = 0; 
  if ((fnow > enow) == (fnow > -enow)) { 
    g0 = enow; 
    enow = e[++eindex]; 
  } else { 
    g0 = fnow; 
    fnow = f[++findex]; 
  } 
  if ((eindex < elen) && ((findex >= flen) 
                          || ((fnow > enow) == (fnow > -enow)))) { 
    Fast_Two_Sum(enow, g0, Qnew, q); 
    enow = e[++eindex]; 
  } else { 
    Fast_Two_Sum(fnow, g0, Qnew, q); 
    fnow = f[++findex]; 
  } 
  Q = Qnew; 
  for (hindex = 0; hindex < elen + flen - 2; hindex++) { 
    if ((eindex < elen) && ((findex >= flen) 
                            || ((fnow > enow) == (fnow > -enow)))) { 
      Fast_Two_Sum(enow, q, R, h[hindex]); 
      enow = e[++eindex]; 
    } else { 
      Fast_Two_Sum(fnow, q, R, h[hindex]); 
      fnow = f[++findex]; 
    } 
    Two_Sum(Q, R, Qnew, q); 
    Q = Qnew; 
  } 
  h[hindex] = q; 
  h[hindex + 1] = Q; 
  return hindex + 2; 
} 
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  linear_expansion_sum_zeroelim()   Sum two expansions, eliminating zero   */ 
/*                                    components from the output expansion.  */ 
/*                                                                           */ 
/*  Sets h = e + f.  See either version of my paper for details.             */ 
/*                                                                           */ 
/*  Maintains the nonoverlapping property.  (That is, if e is                */ 
/*  nonoverlapping, h will be also.)                                         */ 
/*                                                                           */ 
/*****************************************************************************/ 
 
int Predicates::linear_expansion_sum_zeroelim(int elen, const REAL* e,
                                  int flen, const REAL* f, REAL* h)
/* h cannot be e or f. */ 
{ 
  REAL Q, q, hh; 
  INEXACT REAL Qnew; 
  INEXACT REAL R; 
  INEXACT REAL bvirt; 
  REAL avirt, bround, around; 
  int eindex, findex, hindex; 
  int count; 
  REAL enow, fnow; 
  REAL g0; 
 
  enow = e[0]; 
  fnow = f[0]; 
  eindex = findex = 0; 
  hindex = 0; 
  if ((fnow > enow) == (fnow > -enow)) { 
    g0 = enow; 
    enow = e[++eindex]; 
  } else { 
    g0 = fnow; 
    fnow = f[++findex]; 
  } 
  if ((eindex < elen) && ((findex >= flen) 
                          || ((fnow > enow) == (fnow > -enow)))) { 
    Fast_Two_Sum(enow, g0, Qnew, q); 
    enow = e[++eindex]; 
  } else { 
    Fast_Two_Sum(fnow, g0, Qnew, q); 
    fnow = f[++findex]; 
  } 
  Q = Qnew; 
  for (count = 2; count < elen + flen; count++) { 
    if ((eindex < elen) && ((findex >= flen) 
                            || ((fnow > enow) == (fnow > -enow)))) { 
      Fast_Two_Sum(enow, q, R, hh); 
      enow = e[++eindex]; 
    } else { 
      Fast_Two_Sum(fnow, q, R, hh); 
      fnow = f[++findex]; 
    } 
    Two_Sum(Q, R, Qnew, q); 
    Q = Qnew; 
    if (hh != 0) { 
      h[hindex++] = hh; 
    } 
  } 
  if (q != 0) { 
    h[hindex++] = q; 
  } 
  if ((Q != 0.0) || (hindex == 0)) { 
    h[hindex++] = Q; 
  } 
  return hindex; 
} 
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  scale_expansion()   Multiply an expansion by a scalar.                   */ 
/*                                                                           */ 
/*  Sets h = be.  See either version of my paper for details.                */ 
/*                                                                           */ 
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */ 
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */ 
/*  properties as well.  (That is, if e has one of these properties, so      */ 
/*  will h.)                                                                 */ 
/*                                                                           */ 
/*****************************************************************************/ 
 
int Predicates::scale_expansion(int elen, const REAL* e, REAL b, REAL* h)
/* e and h cannot be the same. */ 
{ 
  INEXACT REAL Q; 
  INEXACT REAL sum; 
  INEXACT REAL product1; 
  REAL product0; 
  int eindex, hindex; 
  REAL enow; 
  INEXACT REAL bvirt; 
  REAL avirt, bround, around; 
  INEXACT REAL c; 
  INEXACT REAL abig; 
  REAL ahi, alo, bhi, blo; 
  REAL err1, err2, err3; 
 
  Split(b, bhi, blo); 
  Two_Product_Presplit(e[0], b, bhi, blo, Q, h[0]); 
  hindex = 1; 
  for (eindex = 1; eindex < elen; eindex++) { 
    enow = e[eindex]; 
    Two_Product_Presplit(enow, b, bhi, blo, product1, product0); 
    Two_Sum(Q, product0, sum, h[hindex]); 
    hindex++; 
    Two_Sum(product1, sum, Q, h[hindex]); 
    hindex++; 
  } 
  h[hindex] = Q; 
  return elen + elen; 
} 
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  scale_expansion_zeroelim()   Multiply an expansion by a scalar,          */ 
/*                               eliminating zero components from the        */ 
/*                               output expansion.                           */ 
/*                                                                           */ 
/*  Sets h = be.  See either version of my paper for details.                */ 
/*                                                                           */ 
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */ 
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */ 
/*  properties as well.  (That is, if e has one of these properties, so      */ 
/*  will h.)                                                                 */ 
/*                                                                           */ 
/*****************************************************************************/ 
 
int Predicates::scale_expansion_zeroelim(int elen, const REAL* e, REAL b, REAL* h)
/* e and h cannot be the same. */ 
{ 
  INEXACT REAL Q, sum; 
  REAL hh; 
  INEXACT REAL product1; 
  REAL product0; 
  int eindex, hindex; 
  REAL enow; 
  INEXACT REAL bvirt; 
  REAL avirt, bround, around; 
  INEXACT REAL c; 
  INEXACT REAL abig; 
  REAL ahi, alo, bhi, blo; 
  REAL err1, err2, err3; 
 
  Split(b, bhi, blo); 
  Two_Product_Presplit(e[0], b, bhi, blo, Q, hh); 
  hindex = 0; 
  if (hh != 0) { 
    h[hindex++] = hh; 
  } 
  for (eindex = 1; eindex < elen; eindex++) { 
    enow = e[eindex]; 
    Two_Product_Presplit(enow, b, bhi, blo, product1, product0); 
    Two_Sum(Q, product0, sum, hh); 
    if (hh != 0) { 
      h[hindex++] = hh; 
    } 
    Fast_Two_Sum(product1, sum, Q, hh); 
    if (hh != 0) { 
      h[hindex++] = hh; 
    } 
  } 
  if ((Q != 0.0) || (hindex == 0)) { 
    h[hindex++] = Q; 
  } 
  return hindex; 
} 
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  compress()   Compress an expansion.                                      */ 
/*                                                                           */ 
/*  See the long version of my paper for details.                            */ 
/*                                                                           */ 
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */ 
/*  with IEEE 754), then any nonoverlapping expansion is converted to a      */ 
/*  nonadjacent expansion.                                                   */ 
/*                                                                           */ 
/*****************************************************************************/ 
 
int Predicates::compress(int elen, const REAL* e, REAL* h)
/* e and h may be the same. */ 
{ 
  REAL Q, q; 
  INEXACT REAL Qnew; 
  int eindex, hindex; 
  INEXACT REAL bvirt; 
  REAL enow, hnow; 
  int top, bottom; 
 
  bottom = elen - 1; 
  Q = e[bottom]; 
  for (eindex = elen - 2; eindex >= 0; eindex--) { 
    enow = e[eindex]; 
    Fast_Two_Sum(Q, enow, Qnew, q); 
    if (q != 0) { 
      h[bottom--] = Qnew; 
      Q = q; 
    } else { 
      Q = Qnew; 
    } 
  } 
  top = 0; 
  for (hindex = bottom + 1; hindex < elen; hindex++) { 
    hnow = h[hindex]; 
    Fast_Two_Sum(hnow, Q, Qnew, q); 
    if (q != 0) { 
      h[top++] = q; 
    } 
    Q = Qnew; 
  } 
  h[top] = Q; 
  return top + 1; 
} 
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  estimate()   Produce a one-word estimate of an expansion's value.        */ 
/*                                                                           */ 
/*  See either version of my paper for details.                              */ 
/*                                                                           */ 
/*****************************************************************************/ 
 
REAL Predicates::estimate( int elen, const REAL* e ) 
{ 
  REAL Q; 
  int eindex; 
 
  Q = e[0]; 
  for (eindex = 1; eindex < elen; eindex++) { 
    Q += e[eindex]; 
  } 
  return Q; 
} 
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  DifferenceOfProductsOfDifferences(...), AdaptDiffOfProdsOfDiffs(...)     */ 
/*                                                                           */ 
/*  Adaptations of orient2d and orient2dadapt, SL 10/98                      */ 
/*                                                                           */ 
/*****************************************************************************/ 
double Predicates::DifferenceOfProductsOfDifferences( double a, double b,
                                                      double c, double d,
                                                      double e, double f,
                                                      double g, double h )
{
   SET_DOUBLE_PRECISION_MODE;
   double left;
   double right;
   double diff, errbound;
   double sum;

   left = ( a - b ) * ( c - d );
   right = ( e - f ) * ( g - h );
   diff = left - right;
   if( left > 0.0 )
   {
      if( right <= 0.0 ) return diff;
      else sum = left + right;
   }
   else if( left < 0.0 )
   {
      if( right >= 0.0 ) return diff;
      else sum = -left - right;
   }
   else return diff;
   errbound = ccwerrboundA * sum;
   if( diff >= errbound || -diff >= errbound ) return diff;
   return AdaptDiffOfProdsOfDiffs( a, b, c, d, e, f, g, h, sum );
}

double Predicates::AdaptDiffOfProdsOfDiffs( double terma, double termb,
                                            double termc, double termd,
                                            double terme, double termf,
                                            double termg, double termh,
                                            double sum )
{
   INEXACT REAL diff1, diff2, diff3, diff4;
   REAL diff1tail, diff2tail, diff3tail, diff4tail;
   INEXACT REAL leftprod, rightprod;
   REAL leftprodtail, rightprodtail;
   REAL diff, errbound;
   REAL B[4], C1[8], C2[12], D[16];
   INEXACT REAL B3;
   int C1length, C2length, Dlength;
   REAL u[4];
   INEXACT REAL u3;
   INEXACT REAL s1, t1;
   REAL s0, t0;

   INEXACT REAL bvirt;
   REAL avirt, bround, around;
   INEXACT REAL c;
   INEXACT REAL abig;
   REAL ahi, alo, bhi, blo;
   REAL err1, err2, err3;
   INEXACT REAL _i, _j;
   REAL _0;

   diff1 = static_cast<REAL>( terma - termb );
   diff2 = static_cast<REAL>( termc - termd );
   diff3 = static_cast<REAL>( terme - termf );
   diff4 = static_cast<REAL>( termg - termh );

   Two_Product( diff1, diff2, leftprod, leftprodtail );
   Two_Product( diff3, diff4, rightprod, rightprodtail );

   Two_Two_Diff( leftprod, leftprodtail, rightprod, rightprodtail,
                 B3, B[2], B[1], B[0] );
   B[3] = B3;

   diff = estimate( 4, B);
   errbound = ccwerrboundB * sum;
   if( diff >= errbound || -diff >= errbound ) return diff;

   Two_Diff_Tail( terma, termb, diff1, diff1tail );
   Two_Diff_Tail( termc, termd, diff2, diff2tail );
   Two_Diff_Tail( terme, termf, diff3, diff3tail );
   Two_Diff_Tail( termg, termh, diff4, diff4tail );

   if( diff1tail == 0.0 && diff2tail == 0.0 &&
       diff3tail == 0.0 && diff4tail == 0.0 )
       return diff;

   errbound = ccwerrboundC * sum + resulterrbound * Absolute( diff );
   diff += ( diff1 * diff2tail + diff2 * diff1tail )
       - ( diff3 * diff4tail + diff4 * diff3tail );
   if( diff >= errbound || -diff >= errbound ) return diff;

   Two_Product( diff1tail, diff2, s1, s0 );
   Two_Product( diff3tail, diff4, t1, t0 );
   Two_Two_Diff( s1, s0, t1, t0, u3, u[2], u[1], u[0] );
   u[3] = u3;
   C1length = fast_expansion_sum_zeroelim( 4, B, 4, u, C1 );

   Two_Product( diff1, diff2tail, s1, s0 );
   Two_Product( diff3, diff4tail, t1, t0 );
   Two_Two_Diff( s1, s0, t1, t0, u3, u[2], u[1], u[0] );
   u[3] = u3;
   C2length = fast_expansion_sum_zeroelim( C1length, C1, 4, u, C2 );

   Two_Product( diff1tail, diff2tail, s1, s0 );
   Two_Product( diff3tail, diff4tail, t1, t0 );
   Two_Two_Diff( s1, s0, t1, t0, u3, u[2], u[1], u[0] );
   u[3] = u3;
   Dlength = fast_expansion_sum_zeroelim( C2length, C2, 4, u, D );

   return( D[Dlength - 1] );
}
/*****************************************************************************/ 
/*                                                                           */ 
/*  orient2dfast()   Approximate 2D orientation test.  Nonrobust.            */ 
/*  orient2dexact()   Exact 2D orientation test.  Robust.                    */ 
/*  orient2dslow()   Another exact 2D orientation test.  Robust.             */ 
/*  orient2d()   Adaptive exact 2D orientation test.  Robust.                */ 
/*                                                                           */ 
/*               Return a positive value if the points pa, pb, and pc occur  */ 
/*               in counterclockwise order; a negative value if they occur   */ 
/*               in clockwise order; and zero if they are collinear.  The    */ 
/*               result is also a rough approximation of twice the signed    */ 
/*               area of the triangle defined by the three points.           */ 
/*                                                                           */ 
/*  Only the first and last routine should be used; the middle two are for   */ 
/*  timings.                                                                 */ 
/*                                                                           */ 
/*  The last three use exact arithmetic to ensure a correct answer.  The     */ 
/*  result returned is the determinant of a matrix.  In orient2d() only,     */ 
/*  this determinant is computed adaptively, in the sense that exact         */ 
/*  arithmetic is used only to the degree it is needed to ensure that the    */ 
/*  returned value has the correct sign.  Hence, orient2d() is usually quite */ 
/*  fast, but will run more slowly when the input points are collinear or    */ 
/*  nearly so.                                                               */ 
/*                                                                           */ 
/*****************************************************************************/ 
 
REAL Predicates::orient2dfast(const REAL *pa, const REAL *pb, const REAL *pc) 
{ 
  REAL acx, bcx, acy, bcy; 
 
  acx = pa[0] - pc[0]; 
  bcx = pb[0] - pc[0]; 
  acy = pa[1] - pc[1]; 
  bcy = pb[1] - pc[1]; 
  return acx * bcy - acy * bcx; 
} 
 
 
REAL Predicates::orient2dadapt(const REAL *pa, const REAL *pb, const REAL *pc, REAL detsum) 
{ 
  SET_DOUBLE_PRECISION_MODE;
  INEXACT REAL acx, acy, bcx, bcy; 
  REAL acxtail, acytail, bcxtail, bcytail; 
  INEXACT REAL detleft, detright; 
  REAL detlefttail, detrighttail; 
  REAL det, errbound; 
  REAL B[4], C1[8], C2[12], D[16]; 
  INEXACT REAL B3; 
  int C1length, C2length, Dlength; 
  REAL u[4]; 
  INEXACT REAL u3; 
  INEXACT REAL s1, t1; 
  REAL s0, t0; 
 
  INEXACT REAL bvirt; 
  REAL avirt, bround, around; 
  INEXACT REAL c; 
  INEXACT REAL abig; 
  REAL ahi, alo, bhi, blo; 
  REAL err1, err2, err3; 
  INEXACT REAL _i, _j; 
  REAL _0; 
 
  acx = static_cast<REAL>(pa[0] - pc[0]); 
  bcx = static_cast<REAL>(pb[0] - pc[0]); 
  acy = static_cast<REAL>(pa[1] - pc[1]); 
  bcy = static_cast<REAL>(pb[1] - pc[1]); 
 
  Two_Product(acx, bcy, detleft, detlefttail); 
  Two_Product(acy, bcx, detright, detrighttail); 
 
  Two_Two_Diff(detleft, detlefttail, detright, detrighttail, 
               B3, B[2], B[1], B[0]); 
  B[3] = B3; 
 
  det = estimate(4, B); 
  errbound = ccwerrboundB * detsum; 
  if ((det >= errbound) || (-det >= errbound)) { 
    return det; 
  } 
 
  Two_Diff_Tail(pa[0], pc[0], acx, acxtail); 
  Two_Diff_Tail(pb[0], pc[0], bcx, bcxtail); 
  Two_Diff_Tail(pa[1], pc[1], acy, acytail); 
  Two_Diff_Tail(pb[1], pc[1], bcy, bcytail); 
 
  if ((acxtail == 0.0) && (acytail == 0.0) 
      && (bcxtail == 0.0) && (bcytail == 0.0)) { 
    return det; 
  } 
 
  errbound = ccwerrboundC * detsum + resulterrbound * Absolute(det); 
  det += (acx * bcytail + bcy * acxtail) 
       - (acy * bcxtail + bcx * acytail); 
  if ((det >= errbound) || (-det >= errbound)) { 
    return det; 
  } 
 
  Two_Product(acxtail, bcy, s1, s0); 
  Two_Product(acytail, bcx, t1, t0); 
  Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]); 
  u[3] = u3; 
  C1length = fast_expansion_sum_zeroelim(4, B, 4, u, C1); 
 
  Two_Product(acx, bcytail, s1, s0); 
  Two_Product(acy, bcxtail, t1, t0); 
  Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]); 
  u[3] = u3; 
  C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, C2); 
 
  Two_Product(acxtail, bcytail, s1, s0); 
  Two_Product(acytail, bcxtail, t1, t0); 
  Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]); 
  u[3] = u3; 
  Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, D); 
 
  return(D[Dlength - 1]); 
} 
 
REAL Predicates::orient2d(const REAL *pa, const REAL *pb, const REAL *pc)
{ 
  SET_DOUBLE_PRECISION_MODE;
  REAL detleft, detright, det; 
  REAL detsum, errbound; 
 
  detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]); 
  detright = (pa[1] - pc[1]) * (pb[0] - pc[0]); 
  det = detleft - detright; 
 
  if (detleft > 0.0)
  { 
    if (detright <= 0.0) return det; 
    else detsum = detleft + detright;
  }
  else if (detleft < 0.0)
  { 
    if (detright >= 0.0) return det; 
    else detsum = -detleft - detright; 
  }
  else return det;
  
  errbound = ccwerrboundA * detsum; 
  if ((det >= errbound) || (-det >= errbound))
  { 
    return det; 
  } 
 
  return orient2dadapt(pa, pb, pc, detsum); 
} 
 
 
/*****************************************************************************/ 
/*                                                                           */ 
/*  incirclefast()   Approximate 2D incircle test.  Nonrobust.               */ 
/*  incircleexact()   Exact 2D incircle test.  Robust.                       */ 
/*  incircleslow()   Another exact 2D incircle test.  Robust.                */ 
/*  incircle()   Adaptive exact 2D incircle test.  Robust.                   */ 
/*                                                                           */ 
/*               Return a positive value if the point pd lies inside the     */ 
/*               circle passing through pa, pb, and pc; a negative value if  */ 
/*               it lies outside; and zero if the four points are cocircular.*/ 
/*               The points pa, pb, and pc must be in counterclockwise       */ 
/*               order, or the sign of the result will be reversed.          */ 
/*                                                                           */ 
/*  Only the first and last routine should be used; the middle two are for   */ 
/*  timings.                                                                 */ 
/*                                                                           */ 
/*  The last three use exact arithmetic to ensure a correct answer.  The     */ 
/*  result returned is the determinant of a matrix.  In incircle() only,     */ 
/*  this determinant is computed adaptively, in the sense that exact         */ 
/*  arithmetic is used only to the degree it is needed to ensure that the    */ 
/*  returned value has the correct sign.  Hence, incircle() is usually quite */ 
/*  fast, but will run more slowly when the input points are cocircular or   */ 
/*  nearly so.                                                               */ 
/*                                                                           */ 
/*****************************************************************************/ 
 
REAL Predicates::incirclefast(const REAL *pa, const REAL *pb, const REAL *pc, const REAL *pd) 
{ 
  REAL adx, ady, bdx, bdy, cdx, cdy; 
  REAL abdet, bcdet, cadet; 
  REAL alift, blift, clift; 
 
  adx = pa[0] - pd[0]; 
  ady = pa[1] - pd[1]; 
  bdx = pb[0] - pd[0]; 
  bdy = pb[1] - pd[1]; 
  cdx = pc[0] - pd[0]; 
  cdy = pc[1] - pd[1]; 
 
  abdet = adx * bdy - bdx * ady; 
  bcdet = bdx * cdy - cdx * bdy; 
  cadet = cdx * ady - adx * cdy; 
  alift = adx * adx + ady * ady; 
  blift = bdx * bdx + bdy * bdy; 
  clift = cdx * cdx + cdy * cdy; 
 
  return alift * bcdet + blift * cadet + clift * abdet; 
} 
 
 
REAL Predicates::incircleadapt(const REAL *pa, const REAL *pb, const REAL *pc, const REAL *pd,
                   REAL permanent) 
{ 
  SET_DOUBLE_PRECISION_MODE;
  INEXACT REAL adx, bdx, cdx, ady, bdy, cdy; 
  REAL det, errbound; 
 
  INEXACT REAL bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1; 
  REAL bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0; 
  REAL bc[4], ca[4], ab[4]; 
  INEXACT REAL bc3, ca3, ab3; 
  REAL axbc[8], axxbc[16], aybc[8], ayybc[16], adet[32]; 
  int axbclen, axxbclen, aybclen, ayybclen, alen; 
  REAL bxca[8], bxxca[16], byca[8], byyca[16], bdet[32]; 
  int bxcalen, bxxcalen, bycalen, byycalen, blen; 
  REAL cxab[8], cxxab[16], cyab[8], cyyab[16], cdet[32]; 
  int cxablen, cxxablen, cyablen, cyyablen, clen; 
  REAL abdet[64]; 
  int ablen; 
  REAL fin1[1152], fin2[1152]; 
  REAL *finnow, *finother, *finswap; 
  int finlength; 
 
  REAL adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail; 
  INEXACT REAL adxadx1, adyady1, bdxbdx1, bdybdy1, cdxcdx1, cdycdy1; 
  REAL adxadx0, adyady0, bdxbdx0, bdybdy0, cdxcdx0, cdycdy0; 
  REAL aa[4], bb[4], cc[4]; 
  INEXACT REAL aa3, bb3, cc3; 
  INEXACT REAL ti1, tj1; 
  REAL ti0, tj0; 
  REAL u[4], v[4]; 
  INEXACT REAL u3, v3; 
  REAL temp8[8], temp16a[16], temp16b[16], temp16c[16]; 
  REAL temp32a[32], temp32b[32], temp48[48], temp64[64]; 
  int temp8len, temp16alen, temp16blen, temp16clen; 
  int temp32alen, temp32blen, temp48len, temp64len; 
  REAL axtbb[8], axtcc[8], aytbb[8], aytcc[8]; 
  int axtbblen, axtcclen, aytbblen, aytcclen; 
  REAL bxtaa[8], bxtcc[8], bytaa[8], bytcc[8]; 
  int bxtaalen, bxtcclen, bytaalen, bytcclen; 
  REAL cxtaa[8], cxtbb[8], cytaa[8], cytbb[8]; 
  int cxtaalen, cxtbblen, cytaalen, cytbblen; 
  REAL axtbc[8], aytbc[8], bxtca[8], bytca[8], cxtab[8], cytab[8]; 
  int axtbclen=0, aytbclen=0, bxtcalen=0, bytcalen=0, cxtablen=0, cytablen=0; 
  REAL axtbct[16], aytbct[16], bxtcat[16], bytcat[16], cxtabt[16], cytabt[16]; 
  int axtbctlen, aytbctlen, bxtcatlen, bytcatlen, cxtabtlen, cytabtlen; 
  REAL axtbctt[8], aytbctt[8], bxtcatt[8]; 
  REAL bytcatt[8], cxtabtt[8], cytabtt[8]; 
  int axtbcttlen, aytbcttlen, bxtcattlen, bytcattlen, cxtabttlen, cytabttlen; 
  REAL abt[8], bct[8], cat[8]; 
  int abtlen, bctlen, catlen; 
  REAL abtt[4], bctt[4], catt[4]; 
  int abttlen, bcttlen, cattlen; 
  INEXACT REAL abtt3, bctt3, catt3; 
  REAL negate; 
 
  INEXACT REAL bvirt; 
  REAL avirt, bround, around; 
  INEXACT REAL c; 
  INEXACT REAL abig; 
  REAL ahi, alo, bhi, blo; 
  REAL err1, err2, err3; 
  INEXACT REAL _i, _j; 
  REAL _0; 
 
  adx = static_cast<REAL>(pa[0] - pd[0]); 
  bdx = static_cast<REAL>(pb[0] - pd[0]); 
  cdx = static_cast<REAL>(pc[0] - pd[0]); 
  ady = static_cast<REAL>(pa[1] - pd[1]); 
  bdy = static_cast<REAL>(pb[1] - pd[1]); 
  cdy = static_cast<REAL>(pc[1] - pd[1]); 
 
  Two_Product(bdx, cdy, bdxcdy1, bdxcdy0); 
  Two_Product(cdx, bdy, cdxbdy1, cdxbdy0); 
  Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]); 
  bc[3] = bc3; 
  axbclen = scale_expansion_zeroelim(4, bc, adx, axbc); 
  axxbclen = scale_expansion_zeroelim(axbclen, axbc, adx, axxbc); 
  aybclen = scale_expansion_zeroelim(4, bc, ady, aybc); 
  ayybclen = scale_expansion_zeroelim(aybclen, aybc, ady, ayybc); 
  alen = fast_expansion_sum_zeroelim(axxbclen, axxbc, ayybclen, ayybc, adet); 
 
  Two_Product(cdx, ady, cdxady1, cdxady0); 
  Two_Product(adx, cdy, adxcdy1, adxcdy0); 
  Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]); 
  ca[3] = ca3; 
  bxcalen = scale_expansion_zeroelim(4, ca, bdx, bxca); 
  bxxcalen = scale_expansion_zeroelim(bxcalen, bxca, bdx, bxxca); 
  bycalen = scale_expansion_zeroelim(4, ca, bdy, byca); 
  byycalen = scale_expansion_zeroelim(bycalen, byca, bdy, byyca); 
  blen = fast_expansion_sum_zeroelim(bxxcalen, bxxca, byycalen, byyca, bdet); 
 
  Two_Product(adx, bdy, adxbdy1, adxbdy0); 
  Two_Product(bdx, ady, bdxady1, bdxady0); 
  Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]); 
  ab[3] = ab3; 
  cxablen = scale_expansion_zeroelim(4, ab, cdx, cxab); 
  cxxablen = scale_expansion_zeroelim(cxablen, cxab, cdx, cxxab); 
  cyablen = scale_expansion_zeroelim(4, ab, cdy, cyab); 
  cyyablen = scale_expansion_zeroelim(cyablen, cyab, cdy, cyyab); 
  clen = fast_expansion_sum_zeroelim(cxxablen, cxxab, cyyablen, cyyab, cdet); 
 
  ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet); 
  finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1); 
 
  det = estimate(finlength, fin1); 
  errbound = iccerrboundB * permanent; 
  if ((det >= errbound) || (-det >= errbound)) { 
    return det; 
  } 
 
  Two_Diff_Tail(pa[0], pd[0], adx, adxtail); 
  Two_Diff_Tail(pa[1], pd[1], ady, adytail); 
  Two_Diff_Tail(pb[0], pd[0], bdx, bdxtail); 
  Two_Diff_Tail(pb[1], pd[1], bdy, bdytail); 
  Two_Diff_Tail(pc[0], pd[0], cdx, cdxtail); 
  Two_Diff_Tail(pc[1], pd[1], cdy, cdytail); 
  if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0) 
      && (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)) { 
    return det; 
  } 
 
  errbound = iccerrboundC * permanent + resulterrbound * Absolute(det); 
  det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail) 
                                     - (bdy * cdxtail + cdx * bdytail)) 
          + 2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx)) 
       + ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail) 
                                     - (cdy * adxtail + adx * cdytail)) 
          + 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx)) 
       + ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail) 
                                     - (ady * bdxtail + bdx * adytail)) 
          + 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx)); 
  if ((det >= errbound) || (-det >= errbound)) { 
    return det; 
  } 
 
  finnow = fin1; 
  finother = fin2; 
 
  if ((bdxtail != 0.0) || (bdytail != 0.0) 
      || (cdxtail != 0.0) || (cdytail != 0.0)) { 
    Square(adx, adxadx1, adxadx0); 
    Square(ady, adyady1, adyady0); 
    Two_Two_Sum(adxadx1, adxadx0, adyady1, adyady0, aa3, aa[2], aa[1], aa[0]); 
    aa[3] = aa3; 
  } 
  if ((cdxtail != 0.0) || (cdytail != 0.0) 
      || (adxtail != 0.0) || (adytail != 0.0)) { 
    Square(bdx, bdxbdx1, bdxbdx0); 
    Square(bdy, bdybdy1, bdybdy0); 
    Two_Two_Sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0, bb3, bb[2], bb[1], bb[0]); 
    bb[3] = bb3; 
  } 
  if ((adxtail != 0.0) || (adytail != 0.0) 
      || (bdxtail != 0.0) || (bdytail != 0.0)) { 
    Square(cdx, cdxcdx1, cdxcdx0); 
    Square(cdy, cdycdy1, cdycdy0); 
    Two_Two_Sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0, cc3, cc[2], cc[1], cc[0]); 
    cc[3] = cc3; 
  } 
 
  if (adxtail != 0.0) { 
    axtbclen = scale_expansion_zeroelim(4, bc, adxtail, axtbc); 
    temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, 2.0 * adx, 
                                          temp16a); 
 
    axtcclen = scale_expansion_zeroelim(4, cc, adxtail, axtcc); 
    temp16blen = scale_expansion_zeroelim(axtcclen, axtcc, bdy, temp16b); 
 
    axtbblen = scale_expansion_zeroelim(4, bb, adxtail, axtbb); 
    temp16clen = scale_expansion_zeroelim(axtbblen, axtbb, -cdy, temp16c); 
 
    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                            temp16blen, temp16b, temp32a); 
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, 
                                            temp32alen, temp32a, temp48); 
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
                                            temp48, finother); 
    finswap = finnow; finnow = finother; finother = finswap; 
  } 
  if (adytail != 0.0) { 
    aytbclen = scale_expansion_zeroelim(4, bc, adytail, aytbc); 
    temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, 2.0 * ady, 
                                          temp16a); 
 
    aytbblen = scale_expansion_zeroelim(4, bb, adytail, aytbb); 
    temp16blen = scale_expansion_zeroelim(aytbblen, aytbb, cdx, temp16b); 
 
    aytcclen = scale_expansion_zeroelim(4, cc, adytail, aytcc); 
    temp16clen = scale_expansion_zeroelim(aytcclen, aytcc, -bdx, temp16c); 
 
    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                            temp16blen, temp16b, temp32a); 
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, 
                                            temp32alen, temp32a, temp48); 
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
                                            temp48, finother); 
    finswap = finnow; finnow = finother; finother = finswap; 
  } 
  if (bdxtail != 0.0) { 
    bxtcalen = scale_expansion_zeroelim(4, ca, bdxtail, bxtca); 
    temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, 2.0 * bdx, 
                                          temp16a); 
 
    bxtaalen = scale_expansion_zeroelim(4, aa, bdxtail, bxtaa); 
    temp16blen = scale_expansion_zeroelim(bxtaalen, bxtaa, cdy, temp16b); 
 
    bxtcclen = scale_expansion_zeroelim(4, cc, bdxtail, bxtcc); 
    temp16clen = scale_expansion_zeroelim(bxtcclen, bxtcc, -ady, temp16c); 
 
    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                            temp16blen, temp16b, temp32a); 
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, 
                                            temp32alen, temp32a, temp48); 
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
                                            temp48, finother); 
    finswap = finnow; finnow = finother; finother = finswap; 
  } 
  if (bdytail != 0.0) { 
    bytcalen = scale_expansion_zeroelim(4, ca, bdytail, bytca); 
    temp16alen = scale_expansion_zeroelim(bytcalen, bytca, 2.0 * bdy, 
                                          temp16a); 
 
    bytcclen = scale_expansion_zeroelim(4, cc, bdytail, bytcc); 
    temp16blen = scale_expansion_zeroelim(bytcclen, bytcc, adx, temp16b); 
 
    bytaalen = scale_expansion_zeroelim(4, aa, bdytail, bytaa); 
    temp16clen = scale_expansion_zeroelim(bytaalen, bytaa, -cdx, temp16c); 
 
    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                            temp16blen, temp16b, temp32a); 
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, 
                                            temp32alen, temp32a, temp48); 
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
                                            temp48, finother); 
    finswap = finnow; finnow = finother; finother = finswap; 
  } 
  if (cdxtail != 0.0) { 
    cxtablen = scale_expansion_zeroelim(4, ab, cdxtail, cxtab); 
    temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, 2.0 * cdx, 
                                          temp16a); 
 
    cxtbblen = scale_expansion_zeroelim(4, bb, cdxtail, cxtbb); 
    temp16blen = scale_expansion_zeroelim(cxtbblen, cxtbb, ady, temp16b); 
 
    cxtaalen = scale_expansion_zeroelim(4, aa, cdxtail, cxtaa); 
    temp16clen = scale_expansion_zeroelim(cxtaalen, cxtaa, -bdy, temp16c); 
 
    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                            temp16blen, temp16b, temp32a); 
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, 
                                            temp32alen, temp32a, temp48); 
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
                                            temp48, finother); 
    finswap = finnow; finnow = finother; finother = finswap; 
  } 
  if (cdytail != 0.0) { 
    cytablen = scale_expansion_zeroelim(4, ab, cdytail, cytab); 
    temp16alen = scale_expansion_zeroelim(cytablen, cytab, 2.0 * cdy, 
                                          temp16a); 
 
    cytaalen = scale_expansion_zeroelim(4, aa, cdytail, cytaa); 
    temp16blen = scale_expansion_zeroelim(cytaalen, cytaa, bdx, temp16b); 
 
    cytbblen = scale_expansion_zeroelim(4, bb, cdytail, cytbb); 
    temp16clen = scale_expansion_zeroelim(cytbblen, cytbb, -adx, temp16c); 
 
    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                            temp16blen, temp16b, temp32a); 
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, 
                                            temp32alen, temp32a, temp48); 
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
                                            temp48, finother); 
    finswap = finnow; finnow = finother; finother = finswap; 
  } 
 
  if ((adxtail != 0.0) || (adytail != 0.0)) { 
    if ((bdxtail != 0.0) || (bdytail != 0.0) 
        || (cdxtail != 0.0) || (cdytail != 0.0)) { 
      Two_Product(bdxtail, cdy, ti1, ti0); 
      Two_Product(bdx, cdytail, tj1, tj0); 
      Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]); 
      u[3] = u3; 
      negate = -bdy; 
      Two_Product(cdxtail, negate, ti1, ti0); 
      negate = -bdytail; 
      Two_Product(cdx, negate, tj1, tj0); 
      Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]); 
      v[3] = v3; 
      bctlen = fast_expansion_sum_zeroelim(4, u, 4, v, bct); 
 
      Two_Product(bdxtail, cdytail, ti1, ti0); 
      Two_Product(cdxtail, bdytail, tj1, tj0); 
      Two_Two_Diff(ti1, ti0, tj1, tj0, bctt3, bctt[2], bctt[1], bctt[0]); 
      bctt[3] = bctt3; 
      bcttlen = 4; 
    } else { 
      bct[0] = 0.0; 
      bctlen = 1; 
      bctt[0] = 0.0; 
      bcttlen = 1; 
    } 
 
    if (adxtail != 0.0) { 
      temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, adxtail, temp16a); 
      axtbctlen = scale_expansion_zeroelim(bctlen, bct, adxtail, axtbct); 
      temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, 2.0 * adx, 
                                            temp32a); 
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                              temp32alen, temp32a, temp48); 
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
                                              temp48, finother); 
      finswap = finnow; finnow = finother; finother = finswap; 
      if (bdytail != 0.0) { 
        temp8len = scale_expansion_zeroelim(4, cc, adxtail, temp8); 
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail, 
                                              temp16a); 
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, 
                                                temp16a, finother); 
        finswap = finnow; finnow = finother; finother = finswap; 
      } 
      if (cdytail != 0.0) { 
        temp8len = scale_expansion_zeroelim(4, bb, -adxtail, temp8); 
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail, 
                                              temp16a); 
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, 
                                                temp16a, finother); 
        finswap = finnow; finnow = finother; finother = finswap; 
      } 
 
      temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, adxtail, 
                                            temp32a); 
      axtbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adxtail, axtbctt); 
      temp16alen = scale_expansion_zeroelim(axtbcttlen, axtbctt, 2.0 * adx, 
                                            temp16a); 
      temp16blen = scale_expansion_zeroelim(axtbcttlen, axtbctt, adxtail, 
                                            temp16b); 
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                              temp16blen, temp16b, temp32b); 
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, 
                                              temp32blen, temp32b, temp64); 
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, 
                                              temp64, finother); 
      finswap = finnow; finnow = finother; finother = finswap; 
    } 
    if (adytail != 0.0) { 
      temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, adytail, temp16a); 
      aytbctlen = scale_expansion_zeroelim(bctlen, bct, adytail, aytbct); 
      temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, 2.0 * ady, 
                                            temp32a); 
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                              temp32alen, temp32a, temp48); 
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
                                              temp48, finother); 
      finswap = finnow; finnow = finother; finother = finswap; 
 
 
      temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, adytail, 
                                            temp32a); 
      aytbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adytail, aytbctt); 
      temp16alen = scale_expansion_zeroelim(aytbcttlen, aytbctt, 2.0 * ady, 
                                            temp16a); 
      temp16blen = scale_expansion_zeroelim(aytbcttlen, aytbctt, adytail, 
                                            temp16b); 
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                              temp16blen, temp16b, temp32b); 
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, 
                                              temp32blen, temp32b, temp64); 
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, 
                                              temp64, finother); 
      finswap = finnow; finnow = finother; finother = finswap; 
    } 
  } 
  if ((bdxtail != 0.0) || (bdytail != 0.0)) { 
    if ((cdxtail != 0.0) || (cdytail != 0.0) 
        || (adxtail != 0.0) || (adytail != 0.0)) { 
      Two_Product(cdxtail, ady, ti1, ti0); 
      Two_Product(cdx, adytail, tj1, tj0); 
      Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]); 
      u[3] = u3; 
      negate = -cdy; 
      Two_Product(adxtail, negate, ti1, ti0); 
      negate = -cdytail; 
      Two_Product(adx, negate, tj1, tj0); 
      Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]); 
      v[3] = v3; 
      catlen = fast_expansion_sum_zeroelim(4, u, 4, v, cat); 
 
      Two_Product(cdxtail, adytail, ti1, ti0); 
      Two_Product(adxtail, cdytail, tj1, tj0); 
      Two_Two_Diff(ti1, ti0, tj1, tj0, catt3, catt[2], catt[1], catt[0]); 
      catt[3] = catt3; 
      cattlen = 4; 
    } else { 
      cat[0] = 0.0; 
      catlen = 1; 
      catt[0] = 0.0; 
      cattlen = 1; 
    } 
 
    if (bdxtail != 0.0) { 
      temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, bdxtail, temp16a); 
      bxtcatlen = scale_expansion_zeroelim(catlen, cat, bdxtail, bxtcat); 
      temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, 2.0 * bdx, 
                                            temp32a); 
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                              temp32alen, temp32a, temp48); 
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
                                              temp48, finother); 
      finswap = finnow; finnow = finother; finother = finswap; 
      if (cdytail != 0.0) { 
        temp8len = scale_expansion_zeroelim(4, aa, bdxtail, temp8); 
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail, 
                                              temp16a); 
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, 
                                                temp16a, finother); 
        finswap = finnow; finnow = finother; finother = finswap; 
      } 
      if (adytail != 0.0) { 
        temp8len = scale_expansion_zeroelim(4, cc, -bdxtail, temp8); 
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail, 
                                              temp16a); 
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, 
                                                temp16a, finother); 
        finswap = finnow; finnow = finother; finother = finswap; 
      } 
 
      temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, bdxtail, 
                                            temp32a); 
      bxtcattlen = scale_expansion_zeroelim(cattlen, catt, bdxtail, bxtcatt); 
      temp16alen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, 2.0 * bdx, 
                                            temp16a); 
      temp16blen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, bdxtail, 
                                            temp16b); 
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                              temp16blen, temp16b, temp32b); 
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, 
                                              temp32blen, temp32b, temp64); 
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, 
                                              temp64, finother); 
      finswap = finnow; finnow = finother; finother = finswap; 
    } 
    if (bdytail != 0.0) { 
      temp16alen = scale_expansion_zeroelim(bytcalen, bytca, bdytail, temp16a); 
      bytcatlen = scale_expansion_zeroelim(catlen, cat, bdytail, bytcat); 
      temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, 2.0 * bdy, 
                                            temp32a); 
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                              temp32alen, temp32a, temp48); 
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
                                              temp48, finother); 
      finswap = finnow; finnow = finother; finother = finswap; 
 
 
      temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, bdytail, 
                                            temp32a); 
      bytcattlen = scale_expansion_zeroelim(cattlen, catt, bdytail, bytcatt); 
      temp16alen = scale_expansion_zeroelim(bytcattlen, bytcatt, 2.0 * bdy, 
                                            temp16a); 
      temp16blen = scale_expansion_zeroelim(bytcattlen, bytcatt, bdytail, 
                                            temp16b); 
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                              temp16blen, temp16b, temp32b); 
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, 
                                              temp32blen, temp32b, temp64); 
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, 
                                              temp64, finother); 
      finswap = finnow; finnow = finother; finother = finswap; 
    } 
  } 
  if ((cdxtail != 0.0) || (cdytail != 0.0)) { 
    if ((adxtail != 0.0) || (adytail != 0.0) 
        || (bdxtail != 0.0) || (bdytail != 0.0)) { 
      Two_Product(adxtail, bdy, ti1, ti0); 
      Two_Product(adx, bdytail, tj1, tj0); 
      Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]); 
      u[3] = u3; 
      negate = -ady; 
      Two_Product(bdxtail, negate, ti1, ti0); 
      negate = -adytail; 
      Two_Product(bdx, negate, tj1, tj0); 
      Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]); 
      v[3] = v3; 
      abtlen = fast_expansion_sum_zeroelim(4, u, 4, v, abt); 
 
      Two_Product(adxtail, bdytail, ti1, ti0); 
      Two_Product(bdxtail, adytail, tj1, tj0); 
      Two_Two_Diff(ti1, ti0, tj1, tj0, abtt3, abtt[2], abtt[1], abtt[0]); 
      abtt[3] = abtt3; 
      abttlen = 4; 
    } else { 
      abt[0] = 0.0; 
      abtlen = 1; 
      abtt[0] = 0.0; 
      abttlen = 1; 
    } 
 
    if (cdxtail != 0.0) { 
      temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, cdxtail, temp16a); 
      cxtabtlen = scale_expansion_zeroelim(abtlen, abt, cdxtail, cxtabt); 
      temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, 2.0 * cdx, 
                                            temp32a); 
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                              temp32alen, temp32a, temp48); 
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
                                              temp48, finother); 
      finswap = finnow; finnow = finother; finother = finswap; 
      if (adytail != 0.0) { 
        temp8len = scale_expansion_zeroelim(4, bb, cdxtail, temp8); 
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail, 
                                              temp16a); 
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, 
                                                temp16a, finother); 
        finswap = finnow; finnow = finother; finother = finswap; 
      } 
      if (bdytail != 0.0) { 
        temp8len = scale_expansion_zeroelim(4, aa, -cdxtail, temp8); 
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail, 
                                              temp16a); 
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, 
                                                temp16a, finother); 
        finswap = finnow; finnow = finother; finother = finswap; 
      } 
 
      temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, cdxtail, 
                                            temp32a); 
      cxtabttlen = scale_expansion_zeroelim(abttlen, abtt, cdxtail, cxtabtt); 
      temp16alen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, 2.0 * cdx, 
                                            temp16a); 
      temp16blen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, cdxtail, 
                                            temp16b); 
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                              temp16blen, temp16b, temp32b); 
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, 
                                              temp32blen, temp32b, temp64); 
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, 
                                              temp64, finother); 
      finswap = finnow; finnow = finother; finother = finswap; 
    } 
    if (cdytail != 0.0) { 
      temp16alen = scale_expansion_zeroelim(cytablen, cytab, cdytail, temp16a); 
      cytabtlen = scale_expansion_zeroelim(abtlen, abt, cdytail, cytabt); 
      temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, 2.0 * cdy, 
                                            temp32a); 
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                              temp32alen, temp32a, temp48); 
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, 
                                              temp48, finother); 
      finswap = finnow; finnow = finother; finother = finswap; 
 
 
      temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, cdytail, 
                                            temp32a); 
      cytabttlen = scale_expansion_zeroelim(abttlen, abtt, cdytail, cytabtt); 
      temp16alen = scale_expansion_zeroelim(cytabttlen, cytabtt, 2.0 * cdy, 
                                            temp16a); 
      temp16blen = scale_expansion_zeroelim(cytabttlen, cytabtt, cdytail, 
                                            temp16b); 
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, 
                                              temp16blen, temp16b, temp32b); 
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, 
                                              temp32blen, temp32b, temp64); 
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, 
                                              temp64, finother); 
      finswap = finnow; finnow = finother; finother = finswap; 
    } 
  } 
 
  return finnow[finlength - 1]; 
} 
 
REAL Predicates::incircle(const REAL *pa, const REAL *pb, const REAL *pc, const REAL *pd)
{ 
  SET_DOUBLE_PRECISION_MODE;
  REAL adx, bdx, cdx, ady, bdy, cdy; 
  REAL bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady; 
  REAL alift, blift, clift; 
  REAL det; 
  REAL permanent, errbound; 
 
  adx = pa[0] - pd[0]; 
  bdx = pb[0] - pd[0]; 
  cdx = pc[0] - pd[0]; 
  ady = pa[1] - pd[1]; 
  bdy = pb[1] - pd[1]; 
  cdy = pc[1] - pd[1]; 
 
  bdxcdy = bdx * cdy; 
  cdxbdy = cdx * bdy; 
  alift = adx * adx + ady * ady; 
 
  cdxady = cdx * ady; 
  adxcdy = adx * cdy; 
  blift = bdx * bdx + bdy * bdy; 
 
  adxbdy = adx * bdy; 
  bdxady = bdx * ady; 
  clift = cdx * cdx + cdy * cdy; 
 
  det = alift * (bdxcdy - cdxbdy) 
      + blift * (cdxady - adxcdy) 
      + clift * (adxbdy - bdxady); 
 
  permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * alift 
            + (Absolute(cdxady) + Absolute(adxcdy)) * blift 
            + (Absolute(adxbdy) + Absolute(bdxady)) * clift; 
  errbound = iccerrboundA * permanent; 
  if ((det > errbound) || (-det > errbound)) { 
    return det; 
  } 
 
  return incircleadapt(pa, pb, pc, pd, permanent); 
} 
 
