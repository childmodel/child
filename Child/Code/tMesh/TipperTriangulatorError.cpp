#include "TipperTriangulator.h"
#include "../errors/errors.h"

#include <stdlib.h>

void tt_error_handler(void){
#if defined(TIPPER_TEST)
  exit(1);
#else
  ReportFatalError( "Fatal error in Tipper Triangulator" );
#endif
}

