/***************************************************************************/
/**
**  @file
**  @brief Enable flaoting point traps on Linux
**
** Enable traps on invalid, div/0 and overflow exception
** This code is specific to glibc 2.2 and later.
** Use before any other declaration otherwise _GNU_SOURCE
** will be overridden. Link with -lm (done by default by g++).
*/
/***************************************************************************/

#define _GNU_SOURCE 1
/* to import __GLIBC__ */
#include <stdlib.h>

#if defined (__GLIBC__)
#include <fenv.h>
static void __attribute__ ((constructor)) trapfpe(void)
{
  /* Enable some exceptions.  At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
#elif defined(i386) && defined(__CYGWIN__) && !defined(_lint)
static void __attribute__ ((constructor)) trapfpe(void)
{
  /* x86 specific */
  unsigned int cw = 0x037f & ~(0x01 | 0x04 | 0x08);
  __asm__ ("fldcw %0" : : "m" (*&cw));
}
#endif


#if 0
/* This code enables traps on the same exceptions. This  *
 * should work with older versions of glibc.             */

#include <fpu_control.h>
static void __attribute__ ((constructor)) trapfpe(void)
{
  fpu_control_t cw =
  _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM);
  _FPU_SETCW(cw);
  /* On x86, this expands to: */
  /* unsigned int cw = 0x037f & ~(0x01 | 0x04 | 0x08); */
  /* __asm__ ("fldcw %0" : : "m" (*&cw));              */
}

#endif
