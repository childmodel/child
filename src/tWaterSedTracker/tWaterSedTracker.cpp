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
#include <sstream>
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
  : active(false)
{
  if(0) cout << "tWaterSedTracker constructor" << endl;
}

// copy constructor for copying to new mesh:
tWaterSedTracker::tWaterSedTracker( const tWaterSedTracker& orig,
				    tMesh<tLNode>* nPtr )
  : tracking_node_list_( orig.tracking_node_list_.size() )
{
  output_file_base_name_ = orig.output_file_base_name_;
  const int number_of_nodes_to_track = tracking_node_list_.size();
  tMesh< tLNode >::nodeListIter_t nNI( nPtr->getNodeList() );
  for( int i=0; i<number_of_nodes_to_track; ++i )
    tracking_node_list_[i] = nNI.GetP( orig.tracking_node_list_[i]->getID() );
  // Create and open the output files
  CreateAndOpenWaterAndSedimentOutputFiles( number_of_nodes_to_track, 0.0 );
}

/**************************************************************************/
/**
**  Destructor
**
**  Closes and deletes any remaining ofstream objects.
*/
/**************************************************************************/
tWaterSedTracker::
~tWaterSedTracker()
{
  if(0) cout << "tWaterSedTracker destructor" << endl;
  
  for( unsigned i=0; i<output_file_list_.size(); i++ )
  {
    output_file_list_[i]->close();
    delete output_file_list_[i];
  }
}
  

/**************************************************************************/
/**
**  InitializeFromInputFile
**
**  This method queries the input file to find out which nodes the user
**  wants to track. It assumes that the input file has already been
**  check to determine that the user wants to track some nodes.
*/
/**************************************************************************/
#define kBigNumber 1e10
void tWaterSedTracker::
InitializeFromInputFile( tInputFile &inputFile, tMesh<tLNode> *mesh /*, tErosion *erosion*/ )
{
  active = true;
  int number_of_nodes_to_track;
  string input_string;
  stringstream input_stringstream;
  
  if(1) cout << "tWaterSedTracker::InitializeFromInputFile" << endl;
  
  // Read in the number of nodes to track and their approximate (x,y) coordinates
  number_of_nodes_to_track = inputFile.ReadInt( "NUMBER_OF_NODES_TO_TRACK_WATER_AND_SED" );
  tracking_node_list_.resize( number_of_nodes_to_track, NULL );
  vector<double> x( number_of_nodes_to_track );   // List of tracking node x-coords
  vector<double> y( number_of_nodes_to_track );   // List of tracking node y-coords
  input_string = inputFile.ReadString( "COORDS_OF_NODES_TO_TRACK" );
  input_stringstream.str( input_string );   // use the stringstream to get #s from the input line
  for( int i=0; i<number_of_nodes_to_track; i++ )
  {
    if( input_stringstream.eof() )
      ReportFatalError( "Reached end of line in reading COORDS_OF_NODES_TO_TRACK:\n make sure there is at least one x and y value for each tracking node." );
    input_stringstream >> x[i];
    if( input_stringstream.eof() )
      ReportFatalError( "Reached end of line in reading COORDS_OF_NODES_TO_TRACK:\n make sure there is at least one x and y value for each tracking node." );    
    input_stringstream >> y[i];
  }
  
  // Find the closest node to each (x,y) pair and put it on the tracking list.
  // We do this by sweeping through all interior nodes, checking the square of
  // the distance, and remembering which node has the smallest for each (x,y)
  // pair.
  vector<double> min_distance_squared( number_of_nodes_to_track, kBigNumber );   // Min dist^2 found so far
  tLNode *current_node;
  tMesh< tLNode >::nodeListIter_t ni( mesh->getNodeList() );
  for( current_node=ni.FirstP(); ni.IsActive(); current_node=ni.NextP() )
  {
    for( int i=0; i<number_of_nodes_to_track; i++ )
    {
      double dist_squared = pow(x[i]-current_node->getX(),2)+pow(y[i]-current_node->getY(),2);
      if( dist_squared < min_distance_squared[i] )
      {
        min_distance_squared[i] = dist_squared;  // if the distance is closer, remember it ...
        tracking_node_list_[i] = current_node;   // ... and set this as the closest node.
      }
    }
  }
  
  // Remember the base name for the sed/water flux output files
  output_file_base_name_ = inputFile.ReadString( "OUTFILENAME" );
  
  // Create and open the output files
  CreateAndOpenWaterAndSedimentOutputFiles( number_of_nodes_to_track, 0.0 );
  
  // Debugging error trap!
  cout << "Tracking node list:" << endl;
  for( int i=0; i<number_of_nodes_to_track; i++ )
  {
    assert( tracking_node_list_[i] != NULL );
    cout << i << " " << x[i] << " " << y[i] << " " << tracking_node_list_[i]->getPermID();
    cout << " " << tracking_node_list_[i]->getX() << " " << tracking_node_list_[i]->getY() << endl;
  }
}
#undef kBigNumber


/**************************************************************************/
/**
*/
/**************************************************************************/
void tWaterSedTracker::
ResetListOfNodesToTrack( vector<tLNode *> list_of_nodes_to_track,
                         double current_time )
{
  if(1) cout << "tWaterSedTracker::ResetListOfNodesToTrack" << endl;
  
  // Store a copy of the tracking list (uses vector's overloaded assignment operator)
  tracking_node_list_ = list_of_nodes_to_track;
  
  // Release memory used to store ofstreams
  for( unsigned i=0; i<output_file_list_.size(); i++ )
  {
    output_file_list_[i]->close();
    delete output_file_list_[i];
  }
    
  // Reset output file list
  output_file_base_name_.resize( 0 );
  
  // Create and open output files
  CreateAndOpenWaterAndSedimentOutputFiles( list_of_nodes_to_track.size(),
                                            current_time );

}


/**************************************************************************/
/**
*/
/**************************************************************************/
void tWaterSedTracker::
WriteAndResetWaterSedTimeseriesData( double period_starting_time, 
                                     double period_duration )
{
  if(1) cout << "tWaterSedTracker::WriteAndResetWaterSedTimeseriesData" << endl;
  
  for( unsigned i=0; i<tracking_node_list_.size(); i++ )
  {
    assert( tracking_node_list_[i] != NULL );
    assert( output_file_list_[i]->good() );
    tLNode * cn = tracking_node_list_[i];
    *output_file_list_[i] << period_starting_time << " " << period_duration
        << " " << cn->getQ()
        << " " << cn->getCumulativeSedXportVolume()/period_duration << endl;
    cn->ResetCumulativeSedXportVolume();
  }

}

/**************************************************************************/
/**
*/
/**************************************************************************/
void tWaterSedTracker::
AddSedVolumesAtTrackingNodes( double flux_duration )
{
  if(1) cout << "tWaterSedTracker::AddSedVolumesAtTrackingNodes" << endl;

  for( unsigned i=0; i<tracking_node_list_.size(); i++ )
  {
    tracking_node_list_[i]->AddInfluxToCumulativeSedXportVolume( flux_duration );
  }
}

/**************************************************************************/
/**
*/
/**************************************************************************/
void tWaterSedTracker::
CreateAndOpenWaterAndSedimentOutputFiles( int number_of_nodes_to_track,
                                          double current_time )
{
  for( int i=0; i<number_of_nodes_to_track; i++ )
  {
    stringstream ss;
    ss << output_file_base_name_;
    ss << "_node" << tracking_node_list_[i]->getPermID();
    ss << "_t" << current_time << ".water_sed";
    output_file_list_.push_back( new ofstream( ss.str().c_str() ) );
    /*if( output_file_list_[i]->good() ) cout << "little froglet" << endl;
    cout << "the size of ofl is " << output_file_list_.size() << endl;*/
    if( !output_file_list_[i]->good() )
    {
      cout << "When trying to create output file '" << ss.str() << endl;
      ReportFatalError( "Unable to create file." );
    }
    *output_file_list_[i] << "NODE " << tracking_node_list_[i]->getPermID() << endl
        << "X " << tracking_node_list_[i]->getX() << endl
        << "Y " << tracking_node_list_[i]->getY() << endl;
    *output_file_list_[i] << "Time_start Duration Discharge(m3/yr) Sedflux(m3/yr)" << endl;
  }
}

