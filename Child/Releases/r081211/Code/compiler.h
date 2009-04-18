#ifndef CHILD_COMPILER_H
#define CHILD_COMPILER_H

#ifdef __GNUC_
/* Somewhere in the middle of the GCC 2.96 development cycle, we implemented
   a mechanism by which the user can annotate likely branch directions and
   expect the blocks to be reordered appropriately.  Define __builtin_expect
   to nothing for earlier compilers.  */

# if __GNUC__ == 2 && __GNUC_MINOR__ < 96
# define __builtin_expect(x, expected_value) (x)
# endif

# define likely(x) (__builtin_expect(!!(x),1))
# define unlikely(x) (__builtin_expect(!!(x),0))
#else

# define likely(x) (x)
# define unlikely(x) (x)
#endif

#ifndef ATTRIBUTE_NORETURN
# if defined(__GNUC__) || defined(__INTEL_COMPILER)
#  define ATTRIBUTE_NORETURN __attribute__ ((noreturn))
# else
#  define ATTRIBUTE_NORETURN
# endif
#endif

#endif

