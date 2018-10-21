/**************************************************************************/
/**
**  childDriver.cpp: This provides a test and example of the CHILD
**  interface.
**
**  July 2008
**
**  For information regarding this program, please contact Greg Tucker at:
**
**     Cooperative Institute for Research in Environmental Sciences (CIRES)
**     and Department of Geological Sciences
**     University of Colorado
**     2200 Colorado Avenue, Campus Box 399
**     Boulder, CO 80309-0399
**
*/
/**************************************************************************/

#include "bmi_model.h"


int
main(int argc, char *argv[])
{
  bmi::Model child;
  string argstr;

  for (int i=1; i<argc; i++ ) {
    argstr.append(argv[i]);
    if (i < (argc - 1))
      argstr.append( " " );
  }

  child.Initialize(argstr.c_str());
  child.Run(0);

  child.Finalize();

  return 0;
}
