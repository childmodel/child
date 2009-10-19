//-*-c++-*-

/**************************************************************************/
/**
**  @file tWaterSedTracker.h
**
**  @brief Header file for the tWaterSedTracker class.
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

#ifndef TWATERSEDTRACKER_H
#define TWATERSEDTRACKER_H

#include <vector>
#include <string>
#include "../tInputFile/tInputFile.h"
#include "../tMesh/tMesh.h"
#include "../tLNode/tLNode.h"



class tWaterSedTracker
{
public:

  // Constructor for basic initialization
  tWaterSedTracker();
  
  // Destructor
  ~tWaterSedTracker();
  
  // Initialize nodes to track using data from input file
  void InitializeFromInputFile( tInputFile &inputFile, tMesh<tLNode> *mesh );
  
  // Reset the list of nodes to track
  void ResetListOfNodesToTrack( std::vector<int> ids_of_nodes_to_track );
  
  // Write data to file, and zero out the cumulative flux records within tLNode objects
  void WriteAndResetWaterSedTimeseriesData();
  
  // Add to the total cumulative volume at each node
  void AddSedVolumesAtTrackingNodes( double flux_duration );

private:
  std::vector<tLNode *> tracking_node_list_;
  std::vector<ofstream *> output_file_list_;
  std::string output_file_base_name_;
};


#endif
