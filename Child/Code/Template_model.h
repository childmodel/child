#ifndef TEMPLATE_MODEL_H
#define TEMPLATE_MODEL_H

/**
   @file
   @brief Model used for template organisation
*/

// Either, we put all code in headers or we keep the code definition
// in .cpp files with the same name as the corresponding .h file.

#if defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(_lint) \
    || defined(__BORLANDC__) || defined(_MSC_VER)
# define CHILD_TEMPLATE_IN_HEADER
#endif

#endif
