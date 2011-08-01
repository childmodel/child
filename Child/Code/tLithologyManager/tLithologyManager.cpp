//-*-c++-*-

/**************************************************************************/
/**
 **  @file tLithologyManager.cpp
 **
 **  @brief Implementation of the tLithologyManager class.
 **
 **  For information regarding this program, please contact Greg Tucker at:
 **
 **     Cooperative Institute for Research in Environmental Sciences (CIRES)
 **     and Department of Geological Sciences
 **     University of Colorado
 **     2200 Colorado Avenue, Campus Box 399
 **     Boulder, CO 80309-0399
 */
/**************************************************************************/

#include <vector>
#include <sstream>
#include "tLithologyManager.h"

using namespace std;

/**************************************************************************/
/**
 **  Basic constructor
 **
 **  (does nothing at the moment)
 */
/**************************************************************************/
tLithologyManager::
tLithologyManager()
{
  if(1) cout << "tLithologyManager default constructor" << endl;
}


/**************************************************************************/
/**
 **  Destructor
 **
 **  Closes and deletes any remaining ofstream objects.
 */
/**************************************************************************/
tLithologyManager::
~tLithologyManager()
{
  if(1) cout << "tLithologyManager destructor" << endl;
}


/**************************************************************************/
/**
 **  InitializeFromInputFile
 **
 */
/**************************************************************************/
void tLithologyManager::
InitializeFromInputFile( tInputFile &inputFile, tMesh<tLNode> *mesh )
{
  if(1) cout << "tWaterSedTracker::InitializeFromInputFile" << endl;
  
  meshPtr_ = mesh;
  
  // Read from the input file and do what the user requests
  
  /* Notes:
  
  Things we ought to handle include:
  - read from a restart run
  - read from another CHILD-format lith file with given name
  - set lithology according to a strat sequence of uniform thicknesses
    - with or without a uniform or varying regolith layer on top
  - set lithology according to a map pattern (uniform in depth)
    - with or without a uniform or varying regolith layer on top
  - set regolith thickness according to a map pattern
  - set to simple "bedrock and regolith" (current default)
  - ?? set according to a set of "contact surfaces" that represent the lower
    boundary of a unit
  
  */
  
}


/**************************************************************************/
/**
 */
/**************************************************************************/


/**************************************************************************/
/**
 */
/**************************************************************************/

/**************************************************************************/
/**
 */
/**************************************************************************/

/**************************************************************************/
/**
 */
/**************************************************************************/

