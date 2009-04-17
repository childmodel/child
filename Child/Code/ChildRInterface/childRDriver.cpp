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

#include "childRInterface.h"

int main( int argc, char **argv )
{
	childRInterface myChildRInterface;
	
	myChildRInterface.Initialize( argc, argv );

	if(1) // make this zero to use "example 2" below
	{
		// Example 1: using "Run" method, and setting run duration to zero so model reads duration from input file
		myChildRInterface.Run( 0 );
	}
	else
	{
		// Example 2: using "RunOneStorm" 
		double mytime = 0;
		double myrunduration = 100000;
		
		while( mytime<myrunduration )
		{
			mytime = myChildRInterface.RunOneStorm();
		}
	}
	
	// Note that calling CleanUp() isn't strictly necessary, as the destructor will automatically clean it
	// up when myChildInterface is deleted ... but it's nice to be able to do this at will (and free up
	// memory)
	myChildRInterface.CleanUp();
	
	return 0;
}
