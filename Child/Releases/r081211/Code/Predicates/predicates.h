//-*-c++-*- 

/*****************************************************************************/
/**
**   @file predicates.h
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
//  $Id: predicates.h,v 1.9 2003/07/18 17:49:55 childcvs Exp $
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

#ifndef PREDICATES_H
#define PREDICATES_H

#include <stdio.h> 
#include <math.h> 



/* On some machines, the exact arithmetic routines might be defeated by the  */ 
/*   use of internal extended precision floating-point registers.  Sometimes */ 
/*   this problem can be fixed by defining certain values to be volatile,    */ 
/*   thus forcing them to be stored to memory and rounded off.  This isn't   */ 
/*   a great solution, though, as it slows the arithmetic down.              */ 
/*                                                                           */ 
/* To try this out, write "#define INEXACT volatile" below.  Normally,       */ 
/*   however, INEXACT should be defined to be nothing.  ("#define INEXACT".) */ 
 
#define INEXACT                          /* Nothing */ 
/* #define INEXACT volatile */ 
 
#define REAL double                      /* float or double */ 
 
/* Which of the following two methods of finding the absolute values is      */ 
/*   fastest is compiler-dependent.  A few compilers can inline and optimize */ 
/*   the fabs() call; but most will incur the overhead of a function call,   */ 
/*   which is disastrously slow.  A faster way on IEEE machines might be to  */ 
/*   mask the appropriate bit, but that's difficult to do in C.              */ 
 
#define Absolute(a)  ((a) >= 0.0 ? (a) : -(a)) 
/* #define Absolute(a)  fabs(a) */
//#define Absolute(a)  abs(a) // defined in <macros.h>
 
/* Many of the operations are broken up into two pieces, a main part that    */ 
/*   performs an approximate operation, and a "tail" that computes the       */ 
/*   roundoff error of that operation.                                       */ 
/*                                                                           */ 
/* The operations Fast_Two_Sum(), Fast_Two_Diff(), Two_Sum(), Two_Diff(),    */ 
/*   Split(), and Two_Product() are all implemented as described in the      */ 
/*   reference.  Each of these macros requires certain variables to be       */ 
/*   defined in the calling routine.  The variables `bvirt', `c', `abig',    */ 
/*   `_i', `_j', `_k', `_l', `_m', and `_n' are declared `INEXACT' because   */ 
/*   they store the result of an operation that may incur roundoff error.    */ 
/*   The input parameter `x' (or the highest numbered `x_' parameter) must   */ 
/*   also be declared `INEXACT'.                                             */ 
 
#define Fast_Two_Sum_Tail(a, b, x, y) bvirt = x - a; y = b - bvirt 
 
#define Fast_Two_Sum(a, b, x, y) x = static_cast<REAL>(a + b); Fast_Two_Sum_Tail(a, b, x, y) 

#define Two_Sum_Tail(a, b, x, y) bvirt = static_cast<REAL>(x - a); avirt = x - bvirt; bround = b - bvirt; around = a - avirt; y = around + bround 
 
#define Two_Sum(a, b, x, y) x = static_cast<REAL>(a + b); Two_Sum_Tail(a, b, x, y) 
 
#define Two_Diff_Tail(a, b, x, y) bvirt = static_cast<REAL>(a - x); avirt = x + bvirt; bround = bvirt - b; around = a - avirt; y = around + bround 
 
#define Two_Diff(a, b, x, y) x = static_cast<REAL>(a - b); Two_Diff_Tail(a, b, x, y) 
 
#define Split(a, ahi, alo) c = static_cast<REAL>(splitter * a); abig = static_cast<REAL>(c - a); ahi = c - abig; alo = a - ahi 
 
#define Two_Product_Tail(a, b, x, y) Split(a, ahi, alo); Split(b, bhi, blo); err1 = x - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); y = (alo * blo) - err3 
 
#define Two_Product(a, b, x, y) x = static_cast<REAL>(a * b); Two_Product_Tail(a, b, x, y) 
 
/* Two_Product_Presplit() is Two_Product() where one of the inputs has       */ 
/*   already been split.  Avoids redundant splitting.                        */ 
 
#define Two_Product_Presplit(a, b, bhi, blo, x, y) x = static_cast<REAL>(a * b); Split(a, ahi, alo); err1 = x - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); y = (alo * blo) - err3 
 
/* Square() can be done more quickly than Two_Product().                     */ 
 
#define Square_Tail(a, x, y) Split(a, ahi, alo); err1 = x - (ahi * ahi); err3 = err1 - ((ahi + ahi) * alo); y = (alo * alo) - err3 
 
#define Square(a, x, y) x = static_cast<REAL>(a * a); Square_Tail(a, x, y) 
 
/* Macros for summing expansions of various fixed lengths.  These are all    */ 
/*   unrolled versions of Expansion_Sum().                                   */ 
 
#define Two_One_Sum(a1, a0, b, x2, x1, x0) Two_Sum(a0, b , _i, x0); Two_Sum(a1, _i, x2, x1) 
 
#define Two_One_Diff(a1, a0, b, x2, x1, x0) Two_Diff(a0, b , _i, x0); Two_Sum( a1, _i, x2, x1) 
 
#define Two_Two_Sum(a1, a0, b1, b0, x3, x2, x1, x0) Two_One_Sum(a1, a0, b0, _j, _0, x0); Two_One_Sum(_j, _0, b1, x3, x2, x1) 
 
#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0) Two_One_Diff(a1, a0, b0, _j, _0, x0); Two_One_Diff(_j, _0, b1, x3, x2, x1) 
 

//namespace Predicates
// wanted to make this a namespace, but Sun's compiler was last updated
// in Nov., 1996 (4.2), and pre-dates the introduction of namespaces;
// probably waiting for the final version of the standard to come out
// with a new version compiler; what a pain.
class Predicates
{
public:
   Predicates(); // just calls exactinit()
   ~Predicates() {} // doesn't do anything
   
private:
   // basically the constructor:
   void exactinit();
   
   // basic "exact" arithmetic:
   int grow_expansion(int elen, const REAL* e, REAL b, REAL* h);
   int grow_expansion_zeroelim(int elen, const REAL* e, REAL b, REAL* h);
   int expansion_sum(int elen, const REAL* e, int flen, const REAL* f, REAL* h);
   int expansion_sum_zeroelim1(int elen, const REAL* e, int flen, const REAL* f,
                               REAL* h);
   int expansion_sum_zeroelim2(int elen, const REAL* e, int flen, const REAL* f,
                               REAL* h);
   int fast_expansion_sum(int elen, const REAL* e, int flen, const REAL* f, REAL* h);
   int fast_expansion_sum_zeroelim(int elen, const REAL* e, int flen, const REAL* f,
                                   REAL* h);
   int linear_expansion_sum(int elen, const REAL* e, int flen, const REAL* f, REAL* h);
   int linear_expansion_sum_zeroelim(int elen, const REAL* e, int flen, const REAL* f,
                                     REAL* h);
   int scale_expansion(int elen, const REAL* e, REAL b, REAL* h);
   int scale_expansion_zeroelim(int elen, const REAL* e, REAL b, REAL* h);
   int compress(int elen, const REAL* e, REAL* h);
   REAL estimate( int elen, const REAL* e );
   
public:
   // two functions added by SL, 10/98, to deal with line
   // segment intersection; modeled after orient2d() and
   // orient2dadapt(), which do cross-products and tell you
   // whether the result is positive, negative, or zero:
   double DifferenceOfProductsOfDifferences( double a, double b,
                                             double c, double d,
                                             double e, double f,
                                             double g, double h );
private:
   double AdaptDiffOfProdsOfDiffs( double a, double b,
                                   double c, double d,
                                   double e, double f,
                                   double g, double h,
                                   double sum );

public:
   // orient...() and in...() functions; the ...adapt()
   // functions should not be called directly--they are called
   // by ...()'s (no suffix), e.g., orient2d(); the latter are
   // generally the ones to use; for CHILD purposes, only need
   // orient2d(), orient2dadapt(), and, potentially, incircle()
   // and incircleadapt():
   REAL orient2dfast(const REAL *pa, const REAL *pb, const REAL *pc);
   REAL orient2dadapt(const REAL *pa, const REAL *pb, const REAL *pc, REAL detsum);
   REAL orient2d(const REAL *pa, const REAL *pb, const REAL *pc);
   REAL incirclefast(const REAL *pa, const REAL *pb, const REAL *pc, const REAL *pd);
   REAL incircleadapt(const REAL *pa, const REAL *pb, const REAL *pc, const REAL *pd,
                      REAL permanent);
   REAL incircle(const REAL *pa, const REAL *pb, const REAL *pc, const REAL *pd);
   
private:
   REAL splitter;     /* = 2^ceiling(p / 2) + 1.  Used to split floats in half. */ 
   REAL epsilon;                /* = 2^(-p).  Used to estimate roundoff errors. */ 
/* A set of coefficients used to calculate maximum roundoff errors.          */ 
   REAL resulterrbound; 
   REAL ccwerrboundA, ccwerrboundB, ccwerrboundC; 
   REAL o3derrboundA, o3derrboundB, o3derrboundC; 
   REAL iccerrboundA, iccerrboundB, iccerrboundC; 
   REAL isperrboundA, isperrboundB, isperrboundC; 
};

#endif
