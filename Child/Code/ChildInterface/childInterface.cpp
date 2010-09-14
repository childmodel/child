/**************************************************************************/
/**
**  ChildInterface.cpp: This is the C++ source file for ChildInterface, which
**    provides a generalized, OpenMI-style interface to the CHILD
**    landscape evolution model.
**
**  The first version is really just a sketch.
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

using namespace std;

Predicates predicate;

/**************************************************************************/
/**
**  Default constructor for childInterface
**
**  Sets pointers to NULL and initialized flag to false.
*/
/**************************************************************************/
childInterface::
childInterface() : element_set_id("CHILD_node_element_set"), 
                   element_set_description( "Element set interface for CHILD's voronoi nodes" ),
				   version(0)
{
	initialized = false;
	
	rand = NULL;
	mesh  = NULL;
	output = NULL;
	storm = NULL;
	strmNet = NULL;
	erosion = NULL;
	uplift = NULL;
	time = NULL;
	vegetation = NULL;
	floodplain = NULL;
	stratGrid = NULL;
	loess = NULL;
	strmMeander = NULL;
}


/**************************************************************************/
/**
**  Initialize
**
**  This method accesses the input file given on the command line (and
**  passed via the parameters argc and argv), and creates the necessary
**  objects.
*/
/**************************************************************************/
void childInterface::
Initialize( int argc, char **argv )
{
		
   /****************** INITIALIZATION *************************************\
   **  ALGORITHM
   **    Get command-line arguments (name of input file + any other opts)
   **    Set silent_mode flag
   **    Open main input file
   **    Create and initialize objects for...
   **      Mesh
   **      Output files
   **      Storm
   **      Stream network
   **      Erosion
   **      Uplift (or baselevel change)
   **      Run timer
   **    Write output for initial state
   **    Get options for erosion type, meandering, etc.
   \**********************************************************************/

   // Check command-line arguments
   tOption option( argc, argv );

   // Say hello
   option.version();

   // Open main input file
   tInputFile inputFile( option.inputFile );

   // Create a random number generator for the simulation itself
   rand = new tRand( inputFile );

   // Create and initialize objects:
   std::cout << "Creating mesh...\n";
   mesh = new tMesh<tLNode>( inputFile, option.checkMeshConsistency );

   std::cout << "Creating output files...\n";
   output = new tLOutput<tLNode>( mesh, inputFile, rand );

   storm = new tStorm( inputFile, rand );
   
   std::cout << "Creating stream network...\n";
   strmNet = new tStreamNet( *mesh, *storm, inputFile );
   erosion = new tErosion( mesh, inputFile );
   uplift = new tUplift( inputFile );

   // Get various options
   optDetachLim = inputFile.ReadBool( "OPTDETACHLIM" );
   optDiffuseDepo = inputFile.ReadBool( "OPTDIFFDEP" );
   optVegetation = inputFile.ReadBool( "OPTVEG" );
   optForest = inputFile.ReadBool( "OPTFOREST" );
   optFire = inputFile.ReadBool( "OPTFIRE" );
   optFloodplainDep = inputFile.ReadBool( "OPTFLOODPLAIN" );
   optLoessDep = inputFile.ReadBool( "OPTLOESSDEP" );
   optMeander = inputFile.ReadBool( "OPTMEANDER" );
   optStratGrid = inputFile.ReadBool( "OPTSTRATGRID" ,false);
   optNonlinearDiffusion = inputFile.ReadBool( "OPT_NONLINEAR_DIFFUSION", false );
   optDepthDependentDiffusion = 
     inputFile.ReadBool( "OPT_DEPTH_DEPENDENT_DIFFUSION", false );
   optTrackWaterSedTimeSeries = inputFile.ReadBool( "OPT_TRACK_WATER_SED_TIMESERIES", false );
   
   // create run timer object:
   time = new tRunTimer( inputFile, !option.silent_mode );

   // If applicable, create Vegetation object
   if( optVegetation )
     {
       if( optFire )
	 {
	   if( optForest )
	     vegetation = new tVegetation( mesh, inputFile, time, storm );
	   else
	     vegetation = new tVegetation( mesh, inputFile, time );
	 }
       else
	 vegetation = new tVegetation( mesh, inputFile );
     }
   
   // If applicable, create floodplain object
   if( optFloodplainDep )
       floodplain = new tFloodplain( inputFile, mesh );

   // If applicable, create eolian deposition object
   if( optLoessDep )
       loess = new tEolian( inputFile );

   // If applicable, create Stream Meander object
   if( optMeander )
     strmMeander = new tStreamMeander( *strmNet, *mesh, inputFile, rand );

   // If applicable, create Stratigraphy Grid object
   // and pass it to output
   if( optStratGrid ) {
     if (!optFloodplainDep)
       ReportFatalError("OPTFLOODPLAIN must be enabled.");
     stratGrid = new tStratGrid (inputFile, mesh);
     output->SetStratGrid( stratGrid, strmNet );
   }
   
   // If applicable, set up tracking of water and sediment flux
   if( optTrackWaterSedTimeSeries )
   {
     water_sed_tracker_.InitializeFromInputFile( inputFile, mesh );
     erosion->ActivateSedVolumeTracking( &water_sed_tracker_ );
   }

   std::cout << "Writing data for time zero...\n";
   output->WriteOutput( 0. );
   initialized = true;
   std::cout << "Initialization done.\n";

}


/**************************************************************************/
/**
**  RunOneStorm
**
**  This method executes the model for one storm, calling each of several
**  subroutines in turn. It returns the updated simulation time.
*/
/**************************************************************************/

double childInterface::
RunOneStorm()
{
  double stormDuration, stormPlusDryDuration;
	
  /**************** Run 1 storm ****************************************\
   **  ALGORITHM
   **    Generate storm
   **    Do storm...
   **      Update network (flow directions, drainage area, runoff)
   **      Water erosion/deposition (vertical)
   **      Meandering (if applicable)
   **      Floodplain deposition (if applicable)
   **    Do interstorm...
   **      Hillslope transport
   **      Vegetation (if applicable)
   **      Exposure history
   **      Mesh densification (if applicable)
   **      Eolian (loess) deposition (if applicable)
   **      Uplift (or baselevel change)
   **********************************************************************/
  if(0) //debug
    std::cout << "         " << std::endl;
  time->ReportTimeStatus();
	
  // Do storm...
  storm->GenerateStorm( time->getCurrentTime(),
			strmNet->getInfilt(), strmNet->getSoilStore() );
  stormDuration = min( storm->getStormDuration(), time->RemainingTime() );
  stormPlusDryDuration = min( storm->getStormDuration() + storm->interstormDur(),
			      time->RemainingTime() );
                                
  if(0) // DEBUG
    std::cout << "Remaining time: " << time->RemainingTime() << std::endl;
	
  if(0) //DEBUG
    std::cout<< "Storm: "<< storm->getRainrate() << " " << storm->getStormDuration() << " "
	     << stormDuration << " " << storm->interstormDur() << " " << stormPlusDryDuration 
	     << std::endl;
			
  strmNet->UpdateNet( time->getCurrentTime(), *storm );
  if(0) //DEBUG
    std::cout << "UpdateNet::Done.." << std::endl;
	
  if(0) //DEBUG
    {
      tMesh< tLNode >::nodeListIter_t mli( mesh->getNodeList() );  // gets nodes from the list
      tLNode * cn;
      for( cn=mli.FirstP(); mli.IsActive(); cn=mli.NextP() )
	{
	  if( cn->getX()>64 && cn->getX()<65 )
	    std::cout << "In main loop 1, node at " << cn->getX() << "," << cn->getY() << " has plain id " << cn->getID() << " and perm ID " << cn->getPermID() << std::endl;
	  if( cn->getID()==8121 || cn->getID()==8122 ) 
	    {
	      tEdge * debugedg = cn->getFlowEdg();
	      tLNode * nbr = static_cast<tLNode *>(debugedg->getDestinationPtrNC());
	      std::cout<<"Childmain 1: node "<<cn->getID()<<" edge "<<debugedg->getID()<<" downstream nbr "<<nbr->getID()<<std::endl;
	      std::cout<<"z "<<cn->getZ()<<" dsn z "<<nbr->getZ()<<std::endl;
	    }
	}
    }
	
  // Link tLNodes to StratNodes, adjust elevation StratNode to surrounding tLNodes
  if( optStratGrid )
    stratGrid->UpdateStratGrid(tStratGrid::k0, time->getCurrentTime());
	
  if(0) //DEBUG
    {
      tMesh< tLNode >::nodeListIter_t mli( mesh->getNodeList() );  // gets nodes from the list
      tLNode * cn;
      for( cn=mli.FirstP(); mli.IsActive(); cn=mli.NextP() )
	{
	  if( cn->getX()>64 && cn->getX()<65 )
	    std::cout << "In main loop 2, node at " << cn->getX() << "," << cn->getY() << " has plain id " << cn->getID() << " and perm ID " << cn->getPermID() << std::endl;
	  if( cn->getID()==8121 || cn->getID()==8122 ) 
	    {
	      tEdge * debugedg = cn->getFlowEdg();
	      tLNode * nbr = static_cast<tLNode *>(debugedg->getDestinationPtrNC());
	      std::cout<<"Childmain 2: node "<<cn->getID()<<" edge "<<debugedg->getID()<<" downstream nbr "<<nbr->getID()<<std::endl;
	      std::cout<<"z "<<cn->getZ()<<" dsn z "<<nbr->getZ()<<std::endl;
	    }
	}
    }

  // Do chemical weathering before physical weathering, which may be dependent on
  // the degree of chemical weathering 
  // what this does, whether it does anything,
  // is determined by the chemical weathering option, one of which is "None."
  // this option is read in the tErosion constructor:
  erosion->WeatherBedrock( stormPlusDryDuration );

  // Do physical weathering before diffusion, which may be dependent on thickness
  // and availability
  // what this does, whether it does anything,
  // is determined by the physical weathering option, one of which is "None."
  // this option is read in the tErosion constructor:
  erosion->ProduceRegolith( stormPlusDryDuration, time->getCurrentTime() );

  //Diffusion is now before fluvial erosion in case the tools
  //detachment laws are being used.
  if( optNonlinearDiffusion )
    {
      if( optDepthDependentDiffusion )
	// note: currently only nonlinear version of depth-dependent diffusion,
	// and this function does not allow non-deposition (optDiffuseDepo)
	erosion->DiffuseNonlinearDepthDep( stormPlusDryDuration,
					   time->getCurrentTime() );
      else
	erosion->DiffuseNonlinear( stormPlusDryDuration, optDiffuseDepo );
    }
  else
    erosion->Diffuse( stormPlusDryDuration, optDiffuseDepo );
	
	
  if( optDetachLim )
    erosion->ErodeDetachLim( stormDuration, strmNet,
			     vegetation );
  else
    {
      erosion->DetachErode( stormDuration, strmNet,
			    time->getCurrentTime(), vegetation );
      // To use tools rules, you must use DetachErode2 NMG 07/05
      //erosion.DetachErode2( stormDuration, &strmNet,
      //                      time.getCurrentTime(), vegetation );
    }
	
	
  if(0) //DEBUG
    std::cout << "Erosion::Done.." << std::endl;
	
  // Link tLNodes to StratNodes, adjust elevation StratNode to surrounding tLNodes
  if( optStratGrid )
    stratGrid->UpdateStratGrid(tStratGrid::k1,time->getCurrentTime() );
	
  if(0) //DEBUG
    {
      tMesh< tLNode >::nodeListIter_t mli( mesh->getNodeList() );  // gets nodes from the list
      tLNode * cn;
      for( cn=mli.FirstP(); mli.IsActive(); cn=mli.NextP() )
	{
	  if( cn->getX()>64 && cn->getX()<65 )
	    std::cout << "In main loop 3, node at " << cn->getX() << "," << cn->getY() << " has plain id " << cn->getID() << " and perm ID " << cn->getPermID() << std::endl;
	  if( cn->getID()==8121 || cn->getID()==8122 ) 
	    {
	      tEdge * debugedg = cn->getFlowEdg();
	      tLNode * nbr = static_cast<tLNode *>(debugedg->getDestinationPtrNC());
	      std::cout<<"Childmain 3: node "<<cn->getID()<<" edge "<<debugedg->getID()<<" downstream nbr "<<nbr->getID()<<std::endl;
	      std::cout<<"z "<<cn->getZ()<<" dsn z "<<nbr->getZ()<<std::endl;
	    }
	}	
    }
	
  // TODO: how does Migrate know how long to migrate??
  if( optMeander )
    strmMeander->Migrate( time->getCurrentTime() );
	
  if(0) //DEBUG
    std::cout << "Meander-Migrate::Done..\n";
	
  // Link tLNodes to StratNodes, adjust elevation StratNode to surrounding tLNodes
  if( optStratGrid )
    stratGrid->UpdateStratGrid(tStratGrid::k2,time->getCurrentTime());
	
  //----------------FLOODPLAIN---------------------------------
  if( optFloodplainDep )
    {
      if( floodplain->OptControlMainChan() )
	floodplain->UpdateMainChannelHeight( time->getCurrentTime(), strmNet->getInletNodePtrNC() );
      std::cout << "UpdateChannelHeight::Done..\n";
		
      if( optStratGrid ){
	stratGrid->UpdateStratGrid(tStratGrid::k3,time->getCurrentTime());
      }
		
      floodplain->DepositOverbank( storm->getRainrate(),
				   storm->getStormDuration(),
				   time->getCurrentTime() );
      std::cout << "tFloodplain::Done..\n";
		
      if( optStratGrid ){
	stratGrid->UpdateStratGrid(tStratGrid::k4,time->getCurrentTime());
      }
		
    } // end of floodplain stuff
	
	
#define NEWVEG 1
  if( optVegetation ) {
    if( NEWVEG )
       vegetation->GrowVegetation( mesh, stormPlusDryDuration );
    // previously only grew during interstorm:
//       vegetation->GrowVegetation( mesh, stormPlusDryDuration - stormDuration );
    else
      vegetation->UpdateVegetation( mesh, stormDuration,
				    stormPlusDryDuration - stormDuration );
  }
#undef NEWVEG
	
  // Do interstorm...
  //Diffusion has now been moved to before fluvial erosion NMG 07/05
  // erosion.Diffuse( storm.getStormDuration() + storm.interstormDur(),
  // 		       optDiffuseDepo );
	
  erosion->UpdateExposureTime( stormPlusDryDuration );
	
  if( optLoessDep )
    loess->DepositLoess( mesh,
			 stormPlusDryDuration,
			 time->getCurrentTime() );
	
  if( time->getCurrentTime() < uplift->getDuration() )
    uplift->DoUplift( mesh,
		      stormPlusDryDuration, 
		      time->getCurrentTime() );
	
  if( optTrackWaterSedTimeSeries )
    water_sed_tracker_.WriteAndResetWaterSedTimeseriesData( time->getCurrentTime(),
							    stormPlusDryDuration );
		
  time->Advance( stormPlusDryDuration );
	
  if( time->CheckOutputTime() )
    output->WriteOutput( time->getCurrentTime() );
	
  if( output->OptTSOutput() ) output->WriteTSOutput();
  
  return( time->getCurrentTime() );
}



/**************************************************************************/
/**
**  Run
**
**  This method runs CHILD for a duration specified on the command line or,
**  if that duration is <=0, for the duration specified in the previously
**  read input file.
*/
/**************************************************************************/

void childInterface::
Run( double run_duration )
{
	if( initialized==false )
		ReportFatalError( "childInterface must be initialized (with Initialize() method) before Run() method is called." );

	if( run_duration>0.0 )
		time->Start( time->getCurrentTime(), time->getCurrentTime()+run_duration );

   while( !time->IsFinished() )
		RunOneStorm();

}


/**************************************************************************/
/**
**  CleanUp
**
**  This method deletes the various objects created. Note that this 
**  function is called by the destructor, so if a user function calls it
**  explicitly, it will end up being called twice. 
*/
/**************************************************************************/

void childInterface::
CleanUp()
{
	if( rand ) {
		delete rand;
		rand = NULL;
	}
	if( mesh ) {
		delete mesh;
		mesh = NULL;
	}
	if( output ) {
		delete output;
		output = NULL;
	}
	if( storm ) {
		delete storm;
		storm = NULL;
	}
	if( strmNet ) {
		delete strmNet;
		strmNet = NULL;
	}
	if( erosion ) {
		delete erosion;
		erosion = NULL;
	}
	if( uplift ) {
		delete uplift;
		uplift = NULL;
	}
	if( time ) {
		delete time;
		time = NULL;
	}
	if( optVegetation && vegetation ) {
		delete vegetation;
		vegetation = NULL;
	}
	if( optFloodplainDep && floodplain ) {
		delete floodplain;
		floodplain = NULL;
	}
	if( optLoessDep && loess ) {
		delete loess;
		loess = NULL;
	}
	if( optMeander && strmMeander ) {
		delete strmMeander;
		strmMeander = NULL;
	}
	if( optStratGrid && stratGrid ) {
		delete stratGrid;
		stratGrid = NULL;
	}
	initialized = false;
}


childInterface::
~childInterface()
{
	CleanUp();
}


/**************************************************************************/
/**
**  childInterface::ExternalErosionAndDeposition
**
**  This function allows a calling program to tell CHILD to perform a
**  given depth of erosion and/or sedimentation at each node.
**
**  Presently, it is set up without any grain-size information. A future
**  version could also pass grain-size information.
**
**  GT, June 09 
*/
/**************************************************************************/

void childInterface::ExternalErosionAndDeposition( vector<double> dz )
{
  tMesh< tLNode >::nodeListIter_t mli( mesh->getNodeList() );  // gets nodes from the list
  tLNode * cn;
  for( cn=mli.FirstP(); mli.IsActive(); cn=mli.NextP() )
  {
    int current_id = cn->getPermID();
    cn->EroDep( dz[current_id] );
  }
}


/**************************************************************************/
/**
**  childInterface::ExternalErodeAndDeposit
**
**  This function allows a calling program to tell CHILD to perform
**  erosion and/or deposition sufficient to bring make the new elevation
**  field equal to z. The length of z must be equal to the total
**  number of nodes, and must be listed in ID order 0, 1, 2 ... N-1.
**  Note that although the vector is as long as the total number of nodes,
**  only the interior (non-boundary) nodes are affected.
**
**  Presently, it is set up without any grain-size information. A future
**  version could also pass grain-size information.
**
**  GT, Sep 2010
*/
/**************************************************************************/

void childInterface::ExternalErodeAndDepositToElevation( vector<double> z )
{
  tMesh< tLNode >::nodeListIter_t mli( mesh->getNodeList() );  // gets nodes from the list
  tLNode * cn;
  for( cn=mli.FirstP(); mli.IsActive(); cn=mli.NextP() )
  {
    int current_id = cn->getPermID();
    cn->EroDep( z[current_id] - cn->getZ() );
  }
}



/**************************************************************************/
/**
**  Functions that implement the OpenMI IElement interface
*/
/**************************************************************************/

string childInterface::
getID()
{
    return element_set_id;
}

string childInterface::
getDescription()
{
    return element_set_description;
}

ElementType childInterface::
getElementType()
{
    return XYPolygon;
}

int childInterface::
getElementCount()
{
    if( mesh!=NULL )
	    return mesh->getNodeList()->getSize();
    else
	    return(0);
}

int childInterface::
getVersion()
{
    return version;
}

int childInterface::
GetElementIndex( string elementID )
{
    int my_index;
	std::istringstream my_buffer( elementID );
    my_buffer >> my_index;
    return my_index;
}

string childInterface::
GetElementID( int element_index )
{
    std::string my_ID_string;
	std::stringstream my_buffer;
    my_buffer << element_index;
	my_ID_string = my_buffer.str();
    return my_ID_string;
}

int childInterface::
GetVertexCount( int element_index )
{
   tLNode *my_node;
   tMesh< tLNode >::nodeListIter_t ni( mesh->getNodeList() );
   tEdge * firstedg(0);   // ptr to first edg
   tEdge * curedg;     // pointer to current edge
   int num_vertices;   // Number of vertices
	
   my_node = ni.GetPByPermID( element_index );
   firstedg = my_node->getEdg();
   if( !firstedg ) return( 0 );
   curedg = firstedg->getCCWEdg();
   num_vertices = 1;
   while( curedg!=firstedg ) {
      curedg = curedg->getCCWEdg();
	  num_vertices++;
   }
   
   return( num_vertices );
   
}

// Note: with Voronoi polygons, # faces always = # of vertices
int childInterface::
GetFaceCount( int element_index )
{
   return( GetVertexCount( element_index ) );
}


// Numbering works like this:
//
// - each spoke (directed edge) corresponds to one face of the voronoi polygon
// - spokes and corresponding faces are considered to be numbered counter-clockwise,
//   starting from 0 at the first one (whichever getEdg() returns)
// - the vertex to the RIGHT of each spoke/face has the same index as that spoke/face
// - so, we return the index of the right-hand vertex in face_vertex_indices[0]
// - we return the index of the left-hand vertex, which is associated with the next
//   spoke and face counter-clockwise, in face_vertex_indices[1]
// - because numbers increase counter-clockwise, the left-hand is always equal to
//   either one more than the right (if the right isn't already at the highest index
//   number, which is num_faces-1), or zero.
//
std::vector<int> childInterface::
GetFaceVertexIndices( int element_index, int face_index )
{
   std::vector<int> face_vertex_indices(2);
   int num_faces = GetVertexCount( element_index );
   
   if( face_index>=num_faces ) 
   {
      face_vertex_indices[0] = -1;
	  face_vertex_indices[1] = -1;
   }
   else
   {
      face_vertex_indices[0] = face_index;   // This is the right-hand one
      face_vertex_indices[1] = ( face_index==(num_faces-1) ) ? 0 : face_index+1;
   }
   return( face_vertex_indices );

}


// Returns the right-hand voronoi vertex x-coordinate
double childInterface::
GetXCoordinate( int element_index, int vertex_index )
{
   int i;
   tLNode *my_node;
   tEdge *current_edge;
   tMesh< tLNode >::nodeListIter_t ni( mesh->getNodeList() );
   tArray2<double> right_hand_voronoi_vertex_coords;
   
   my_node = ni.GetPByPermID( element_index );
   if( my_node<=0 ) return( -9999 );   // Temporary hacked NODATA code
   
   current_edge = my_node->getEdg();
   if( !current_edge ) return( -9999 );
   
   for( i=0; i<vertex_index; i++ )
      current_edge = current_edge->getCCWEdg();
	  
   current_edge->getRVtx( right_hand_voronoi_vertex_coords );
   
   return( right_hand_voronoi_vertex_coords.at(0) );

}


// Returns the right-hand voronoi vertex y-coordinate
double childInterface::
GetYCoordinate( int element_index, int vertex_index )
{
   int i;
   tLNode *my_node;
   tEdge *current_edge;
   tMesh< tLNode >::nodeListIter_t ni( mesh->getNodeList() );
   tArray2<double> right_hand_voronoi_vertex_coords;
   
   my_node = ni.GetPByPermID( element_index );
   if( my_node<=0 ) return( -9999 );   // Temporary hacked NODATA code
   
   current_edge = my_node->getEdg();
   if( !current_edge ) return( -9999 );
   
   for( i=0; i<vertex_index; i++ )
      current_edge = current_edge->getCCWEdg();
	  
   current_edge->getRVtx( right_hand_voronoi_vertex_coords );
   
   return( right_hand_voronoi_vertex_coords.at(1) );

}


double childInterface::
GetZCoordinate( int element_index, int vertex_index )
{
   int i;
   tLNode *my_node;
   tEdge *current_edge;
   tMesh< tLNode >::nodeListIter_t ni( mesh->getNodeList() );
   tArray2<double> right_hand_voronoi_vertex_coords;
   double z_coord;
   tTriangle *my_triangle;
   tArray<double> p0, p1, p2;
   tNode * triangle_vertex_node;
   tArray<double> zs;
   
   // Get the node
   my_node = ni.GetPByPermID( element_index );
   if( my_node<=0 ) return( -9999 );   // Temporary hacked NODATA code
   
   // Get the right edge
   current_edge = my_node->getEdg();
   if( !current_edge ) return( -9999 );
   for( i=0; i<vertex_index; i++ )
      current_edge = current_edge->getCCWEdg();

   // Get the vertex coordinates for the right-hand voronoi vertex
   current_edge->getRVtx( right_hand_voronoi_vertex_coords );
	  
   // Get the triangle for which the Voronoi vertex is the circumcenter
   my_triangle = current_edge->TriWithEdgePtr();
   if( !my_triangle ) return( -9999 );
   
   // Find the z-coordinate by interpolation around the 3 points of the triangle
   zs.setSize(3);
   triangle_vertex_node = my_triangle->pPtr(0);
   triangle_vertex_node->get2DCoords( p0 );
   zs[0] = triangle_vertex_node->getZ();
   triangle_vertex_node = my_triangle->pPtr(1);
   triangle_vertex_node->get2DCoords( p1 );
   zs[1] = triangle_vertex_node->getZ();
   triangle_vertex_node = my_triangle->pPtr(2);
   triangle_vertex_node->get2DCoords( p2 );
   zs[2] = triangle_vertex_node->getZ();
   z_coord = PlaneFit( right_hand_voronoi_vertex_coords.at(0), right_hand_voronoi_vertex_coords.at(1),
                       p0, p1, p2, zs);

   return( z_coord );

}

// Returns TRUE if the node is an interior node (non-boundary), FALSE otherwise
bool childInterface::
IsInteriorNode( int element_index )
{
   tMesh< tLNode >::nodeListIter_t ni( mesh->getNodeList() );
   tLNode *my_node;
   
   my_node = ni.GetPByPermID( element_index );
   return( my_node->getBoundaryFlag()==kNonBoundary );
}


/**************************************************************************/
/**
**  childInterface::GetNodeCount
**
**  Returns the total number of nodes (including boundaries).
**
**  GT, Apr 2010
*/
/**************************************************************************/
long childInterface::GetNodeCount()
{
   return mesh->getNodeList()->getSize();
}

/**************************************************************************/
/**
**  childInterface::GetNodeCoords
**
**  Returns the (x,y,z) coordinates of each node in the form
**  x0,y0,z0,x1,y1,z1, ... xN-1,yN-1,zN-1
**
**  GT, Apr 2010
*/
/**************************************************************************/
std::vector<double> childInterface::GetNodeCoords()
{
   tLNode *current_node;
   tMesh<tLNode>::nodeListIter_t ni( mesh->getNodeList() );
   std::vector<double> coords( 3*mesh->getNodeList()->getSize() );
   
   for( current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP() )
   {
      coords[3*current_node->getPermID()+0] = current_node->getX();
      coords[3*current_node->getPermID()+1] = current_node->getY();
      coords[3*current_node->getPermID()+2] = current_node->getZ();
      if(1) std::cout << "Node " << current_node->getPermID()
                      << " x=" << current_node->getX()
                      << " y=" << current_node->getY()
                      << " z=" << current_node->getZ() << std::endl;
   }
   
   return coords;

}

/**************************************************************************/
/**
**  childInterface::GetTriangleCount
**
**  Returns the total number of triangles.
**
**  GT, Apr 2010
*/
/**************************************************************************/
long childInterface::GetTriangleCount()
{
   return mesh->getTriList()->getSize();
}

/**************************************************************************/
/**
**  childInterface::GetTriangleVertexIDs
**
**  Returns the ID numbers of the 3 triangle vertex nodes, in
**  counter-clockwise order.
**
**  Format of the list is p0_0,p0_1,p0_2,p1_0,p1_1,... etc., where the
**  first number is triangle position in list (not necessarily ID), and
**  the second is vertex number.
**
**  GT, Apr 2010
*/
/**************************************************************************/
std::vector<long> childInterface::GetTriangleVertexIDs()
{
   tTriangle *current_tri;
   tMesh<tLNode>::triListIter_t ti( mesh->getTriList() );
   std::vector<long> vertex_ids( 3*mesh->getTriList()->getSize() );
   
   long k=0;
   for( current_tri=ti.FirstP(); !ti.AtEnd(); current_tri=ti.NextP() )
   {
      vertex_ids[k+0] = current_tri->pPtr(0)->getPermID();
      vertex_ids[k+1] = current_tri->pPtr(1)->getPermID();
      vertex_ids[k+2] = current_tri->pPtr(2)->getPermID();
      k += 3;
   }
   
   return vertex_ids;

}

/**************************************************************************/
/**
**  childInterface::GetCurrentTime
**
**  Returns the current time in the simulation.
**
**  GT, Apr 2010
*/
/**************************************************************************/
double childInterface::
GetCurrentTime()
{
   assert( time != NULL );
   return time->getCurrentTime();
}

/**************************************************************************/
/**
**  childInterface::GetRemainingRunTime()
**
**  Returns the remaining simulation time before the end of the run.
**
**  GT, Apr 2010
*/
/**************************************************************************/
double childInterface::
GetRemainingRunTime()
{
   assert( time != NULL );
   return time->RemainingTime();
}

/**************************************************************************/
/**
**  childInterface::GetValueSet
**
**  This method finds and returns a set of values on nodes. It does this
**  by parsing its string argument to discover which quantity the caller
**  wants, and calling the appropriate helper function to get and return
**  the result.
**
**  GT, Sep/Oct 2009
*/
/**************************************************************************/
std::vector<double> childInterface::
GetValueSet( string var_name )
{
  if(1) std::cout << "childInterface::GetValueSet() here with request '"
                  << var_name << "'\n";
  if( var_name.compare( 0,4,"elev" )==0 )
  {
    if(1) std::cout << "request for elevs\n";
    return GetNodeElevationVector();
  }
  else if( var_name.compare( 0,2,"dz" )==0 || var_name.compare( 0,3,"ero" )==0 )
  {
    if(1) std::cout << "request for dz\n";
    return GetNodeErosionVector();
  }
  else if( var_name.compare( 0,5,"disch" )==0 || var_name.compare( 0,5,"water" )==0 )
  {
    if(1) std::cout << "request for Q\n";
    return GetNodeDischargeVector();
  }
  else if( var_name.compare( 0,3,"sed" )==0 )
  {
    if(1) std::cout << "request for Qs\n";
    return GetNodeSedimentFluxVector();
  }
  else
  {
    if(1) std::cout << "request for NOTHING!\n";
    std::vector<double> empty_vector;
    return empty_vector;
  };
}



/**************************************************************************/
/**
**  childInterface::GetNodeElevationVector
**
**  This private method is used to create and return a vector of elevation
**  values for the nodes. The vector is in order by permanent node ID.
**
**  GT, Sep/Oct 2009
*/
/**************************************************************************/
std::vector<double> childInterface::
GetNodeElevationVector()
{
   tLNode *current_node;
   tMesh<tLNode>::nodeListIter_t ni( mesh->getNodeList() );
   std::vector<double> elevations( mesh->getNodeList()->getSize() );
   
   for( current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP() )
   {
      if(1) std::cout << "node " << current_node->getPermID()
                      << " z=" << current_node->getZ() << std::endl;
      elevations[current_node->getPermID()] = current_node->getZ();
   }
   
   return elevations;

}

/**************************************************************************/
/**
**  childInterface::GetNodeErosionVector
**
**  This private method is used to create and return a vector of 
**  cumulative erosion/deposition depths since the last time the function
**  was invoked.
**
**  GT, Oct 2009
*/
/**************************************************************************/
std::vector<double> childInterface::
GetNodeErosionVector()
{
   tLNode *current_node;
   tMesh<tLNode>::nodeListIter_t ni( mesh->getNodeList() );
   std::vector<double> dz( mesh->getNodeList()->getSize() );
      
   for( current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP() )
   {
      dz[current_node->getPermID()] = current_node->getCumulativeEroDep();
      current_node->ResetCumulativeEroDep();
   }

   return dz;

}


/**************************************************************************/
/**
**  childInterface::GetNodeDischargeVector
**
**  This private method is used to create and return a vector of discharge
**  values for the nodes. The vector is in order by permanent node ID.
**  Warning: it is the discharge of the most recent storm, not necessarily
**  an average or steady value.
**
**  GT, Apr 2010
*/
/**************************************************************************/
std::vector<double> childInterface::
GetNodeDischargeVector()
{
   tLNode *current_node;
   tMesh<tLNode>::nodeListIter_t ni( mesh->getNodeList() );
   std::vector<double> discharge( mesh->getNodeList()->getSize() );
   
   for( current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP() )
   {
      if(1) std::cout << "NOde " << current_node->getPermID()
                      << " Q=" << current_node->getQ() << std::endl;
      discharge[current_node->getPermID()] = current_node->getQ();
   }
   
   return discharge;

}

/**************************************************************************/
/**
**  childInterface::GetNodeSedimentFluxVector
**
**  This private method is used to create and return a vector of sed flux
**  values for the nodes. The vector is in order by permanent node ID.
**  Warning: it is the most recent instantaneous value, not necessarily
**  an average or steady value. It is also just the water-borne flux,
**  not the hillslope mass flux.
**
**  GT, Apr 2010
*/
/**************************************************************************/
std::vector<double> childInterface::
GetNodeSedimentFluxVector()
{
   tLNode *current_node;
   tMesh<tLNode>::nodeListIter_t ni( mesh->getNodeList() );
   std::vector<double> sedflux( mesh->getNodeList()->getSize() );
   
   for( current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP() )
   {
      if(1) std::cout << "noDE " << current_node->getPermID()
                      << " Qs=" << current_node->getQs() << std::endl;
      sedflux[current_node->getPermID()] = current_node->getQs();
   }
   
   return sedflux;

}

/**************************************************************************/
/**
**  childInterface::SetValueSet
**
**  This method sets values of a CHILD variable at each node.
**  The variable to be changed is passed as a string, while the new values
**  are passed as a vector of doubles. Use with caution!
**  The nodes are assumed to be in order by permanent ID, starting from 0.
**  The caller MUST specify an elevation value for every node.
**
**  GT, May 2010
*/
/**************************************************************************/
void childInterface::
SetValueSet( string var_name, std::vector<double> value_set )
{
  if(1) std::cout << "childInterface::SetValueSet() here with request '"
                  << var_name << "'\n";
  if( var_name.compare( 0,4,"elev" )==0 )
  {
    if(1) std::cout << "request for elevs\n";
    return SetNodeElevations( value_set );
  }
  else
  {
    std::cout << "Warning: unrecognized value set '" << var_name << "'\n";
    std::cout << "Request to set values ignored\n";
  }
}


/**************************************************************************/
/**
**  childInterface::SetNodeElevations
**
**  Sets the elevations of nodes to the values specified in "elevations".
**  The nodes are assumed to be in order by permanent ID, starting from 0.
**  The caller MUST specify an elevation value for every node.
**
**  GT, May 2010
*/
/**************************************************************************/
void childInterface::
SetNodeElevations( std::vector<double> elevations )
{
   tLNode *current_node;
   tMesh<tLNode>::nodeListIter_t ni( mesh->getNodeList() );
   
   for( current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP() )
   {
      if(1) std::cout << "node " << current_node->getPermID()
                      << " changing z from " << current_node->getZ() 
                      << " to " << elevations[current_node->getPermID()] << std::endl;
      current_node->setZ( elevations[current_node->getPermID()] );
   }
   
}


/**************************************************************************/
/**
**  childInterface::AdjustElevations
**
**  This public function adjusts the elevations of the nodes, without
**  changing layers/stratigraphy, by amounts specified in "dz".
**  It is intended to be used, for example, when a tectonic model 
**  calculates a pattern of uplift and/or subsidence over a particular
**  time interval.
**  The nodes are assumed to be in order by permanent ID, starting from 0.
**  The caller MUST specify an elevation value for every node.
**
**  GT, May 2010
*/
/**************************************************************************/
void childInterface::
AdjustElevations( std::vector<double> dz )
{
   tLNode *current_node;
   tMesh<tLNode>::nodeListIter_t ni( mesh->getNodeList() );
   
   for( current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP() )
   {
      if(1) std::cout << "Adjust elevations: node " << current_node->getPermID()
                      << " changing z from " << current_node->getZ() 
                      << " to " 
                      << current_node->getZ() + dz[current_node->getPermID()] 
                      << std::endl;
      current_node->ChangeZ( dz[current_node->getPermID()] );
   }
   
}

/**************************************************************************/
/**
**  childInterface::AdjustInteriorElevations
**
**  This public function adjusts the elevations of the nodes, without
**  changing layers/stratigraphy, by amounts specified in "dz".
**  It is intended to be used, for example, when a tectonic model 
**  calculates a pattern of uplift and/or subsidence over a particular
**  time interval.
**  The nodes are assumed to be in order by permanent ID, starting from 0.
**  The caller MUST specify an elevation value for every node, but only
**  the interior (non-boundary) nodes are actually modified.
**
**  GT, May 2010
*/
/**************************************************************************/
void childInterface::
AdjustInteriorElevations( std::vector<double> dz )
{
   tLNode *current_node;
   tMesh<tLNode>::nodeListIter_t ni( mesh->getNodeList() );
   
   for( current_node=ni.FirstP(); ni.IsActive(); current_node=ni.NextP() )
   {
      if(1) std::cout << "AdjustInteriorElevations: node " 
                      << current_node->getPermID()
                      << " changing z from " << current_node->getZ() 
                      << " to " 
                      << current_node->getZ() + dz[current_node->getPermID()] 
                      << std::endl;
      current_node->ChangeZ( dz[current_node->getPermID()] );
   }
   
}


/**************************************************************************/
/**
**  childInterface::TrackWaterAndSedFluxAtNodes
**
**  This method tells CHILD to track and record water and sediment flux
**  time series at specified nodes. The method switches on tracking if it
**  is not already on. It uses the tWaterSedTracker interface functions to
**  reset the list of nodes to track.
**
**  GT, Oct 2009
*/
/**************************************************************************/
void childInterface::
TrackWaterAndSedFluxAtNodes( vector<int> ids_of_nodes_to_track,
                             double current_time )
{
  // Set tracking option on, and make sure Erosion knows about it too
  if( !optTrackWaterSedTimeSeries )
    optTrackWaterSedTimeSeries = true;
  erosion->ActivateSedVolumeTracking( &water_sed_tracker_ );
  
  // Compose a list of pointers to tracking nodes
  tMesh< tLNode >::nodeListIter_t ni( mesh->getNodeList() );
  vector<tLNode *> list_of_nodes_to_track;
  for( unsigned i=0; i<ids_of_nodes_to_track.size(); i++ )
  {
    list_of_nodes_to_track.push_back( ni.GetPByPermID( ids_of_nodes_to_track[i] ) );
  }

  // Tell the WaterSedTracker to reset its list
  water_sed_tracker_.ResetListOfNodesToTrack( list_of_nodes_to_track,
                                              current_time );
}


