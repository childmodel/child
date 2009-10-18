//-*-c++-*-

/**************************************************************************/
/**
**  @file tWaterSedTracker.cpp
**
**  @brief Implementation of the tWaterSedTracker class.
**
**  A tWaterSedTracker is an object that implements tracking and writing
**  of time series water and sediment discharge at specified nodes.
**  It works with help from tLNode (which needs a data item to keep track
**  of sediment) and tErosion (which needs to know whether it should track
**  fluvial sediment flux at each time step, and if so where).
**
**  The class was first written in October 2009 by GT in order to
**  facilitate coupling of CHILD and SedFlux. The latter needs point
**  sources of water and sediment flux, representing rivers. The idea
**  behind allowing a caller to request particular nodes to track is that
**  a coupling code can find out from SedFlux where its shoreline lies,
**  use CHILD's IElement interface to figure out which CHILD nodes lie
**  along the shoreline, and tell CHILD to keep track of water and sediment
**  discharge at those nodes.
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
#include "tWaterSedTracker.h"

using namespace std;

/**************************************************************************/
/**
**  Basic constructor
**
**  (does nothing at the moment)
*/
/**************************************************************************/
tWaterSedTracker::
tWaterSedTracker()
{}
  

/**************************************************************************/
/**
**  InitializeFromInputFile
**
**  This method queries the input file to find out which nodes the user
**  wants to track. It assumes that the input file has already been
**  check to determine that the user wants to track some nodes.
*/
/**************************************************************************/
void tWaterSedTracker::
InitializeFromInputFile( tInputFile &inputFile, tMesh<tLNode> *mesh, tErosion *erosion )
{
  int number_of_nodes_to_track;
  
  number_of_nodes_to_track = inputFile.ReadInt( "NUMBER_OF_NODES_TO_TRACK_WATER_AND_SED" );
  
}



/**************************************************************************/
/**
*/
/**************************************************************************/
void tWaterSedTracker::
ResetListOfNodesToTrack( vector<int> ids_of_nodes_to_track )
{}


/**************************************************************************/
/**
*/
/**************************************************************************/
void tWaterSedTracker::
WriteAndResetWaterSedTimeseriesData()
{}

