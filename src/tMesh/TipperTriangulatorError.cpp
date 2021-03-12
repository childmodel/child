/***************************************************************************/
/**
**  @file
**  @brief Error handler for Tipper triangulator
*/
/***************************************************************************/

#include "TipperTriangulator.h"

#if defined(TIPPER_TEST)
# include <stdlib.h>
#else
# include "../errors/errors.h"
#endif

void tt_error_handler(void){
#if defined(TIPPER_TEST)
  exit(1);
#else
  ReportFatalError( "Fatal error in Tipper Triangulator" );
#endif
}

