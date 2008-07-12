/**************************************************************************/
/**
**  ChildDriver.cpp: This provides a test and example of the CHILD
**  interface.
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

#include "childInterface.h"

int main( int argc, char **argv )
{
	childInterface myChildInterface;
	
	myChildInterface.Initialize( argc, argv );
	myChildInterface.RunOneStorm();
	myChildInterface.CleanUp();

	return 0;
}
