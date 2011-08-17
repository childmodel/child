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
**  Initialize_Copy
**
**  Like a copy constructor or assignment operator, except that the members
**  of childInterface that refer to important Child objects are
**  pointers. A copy constructor or assignment operator would simply copy
**  those pointers. For separate copies of those objects, we need something
**  different, which is this. Before calling this function, construct a 
**  childInterface object with the default constructor, and then call this
**  function from that uninitialized object and pass the state you want to 
**  copy:
**    childInterface myChildInterface;
**    myChildInterface.Initialize( argc, argv );
**    childInterface childCopy;
**    childCopy.Initialize_Copy( myChildInterface );
**
**  I've made a strong attempt to make sure that all necessary copy 
**  constructors and assignment operators exist, and that pointer members
**  of the various objects are assigned the appropriate pointers that
**  point to members of the copy rather than simply copying the pointers
**  to members of the original. That said, be careful...
**
**  This function can be used to (a) store an initial condition, e.g.,
**  for comparison to a final result, and (b) create a "daughter"
**  simulation with all the parameter values of the "parent," where those
**  parameter values may then be changed, e.g., in a Monte Carlo 
**  optimization.
**
**  SL, 10/2010
*/
/**************************************************************************/
void childInterface::Initialize_Copy( const childInterface& orig )
{
  initialized = orig.initialized; // if you copy an unitialized object...
  optNoDiffusion = orig.optNoDiffusion;
  optNoFluvial = orig.optNoFluvial;
  optNoUplift = orig.optNoUplift;
  optDetachLim = orig.optDetachLim;
  optFloodplainDep = orig.optFloodplainDep;
  optLoessDep = orig.optLoessDep;
  optVegetation = orig.optVegetation;
  optFire = orig.optFire;
  optForest = orig.optForest;
  optMeander = orig.optMeander;
  optDiffuseDepo = orig.optDiffuseDepo;
  optStratGrid = orig.optStratGrid;
  optTrackWaterSedTimeSeries = orig.optTrackWaterSedTimeSeries;
  optNonlinearDiffusion = orig.optNonlinearDiffusion;
  optDepthDependentDiffusion = orig.optDepthDependentDiffusion;
  optLandslides = orig.optLandslides;
  opt3DLandslides = orig.opt3DLandslides;
  optChemicalWeathering = orig.optChemicalWeathering;
  optPhysicalWeathering = orig.optPhysicalWeathering;
  optStreamLineBoundary = orig.optStreamLineBoundary;
  
  if( orig.rand )
    rand = new tRand( *orig.rand );
  if( orig.mesh )
    mesh  = new tMesh<tLNode>( orig.mesh );
  // copying output objects too problematic:
  output = 0;
  //   if( orig.output )
  //     output = new tLOutput<tLNode>( *orig.output );
  if( orig.storm )
  {
    storm = new tStorm( *orig.storm );
    storm->setRand( rand );
  }
  if( orig.strmNet )
  {
    //       const int nA = orig.mesh->getNodeList()->getActiveSize();
    //       vector<int> flowDest( nA );
    //       tMesh<tLNode>::nodeListIter_t oNI( orig.mesh->getNodeList() );
    //       int i=0;
    //       for( tLNode* cn = oNI.FirstP(); oNI.IsActive(); cn = oNI.NextP() )
    // 	flowDest[i] = cn->getDownstrmNbr()->getID();
    //       strmNet = new tStreamNet( *orig.strmNet, storm, mesh, flowDest );
    strmNet = new tStreamNet( *orig.strmNet, storm, mesh );
  }
  if( orig.erosion )
    erosion = new tErosion( *orig.erosion, mesh );
  if( orig.uplift )
    uplift = new tUplift( *orig.uplift );
  if( orig.time )
    time = new tRunTimer( *orig.time );
  if( orig.vegetation )
  {
    vegetation = new tVegetation( *orig.vegetation );
    if( vegetation->FirePtr() )
      vegetation->FirePtr()->setTimePtr( time );
    if( vegetation->ForestPtr() )
    {
      vegetation->ForestPtr()->setMeshPtr( mesh ); // also sets pointers in tTrees
      vegetation->ForestPtr()->setStormPtr( storm );
    }
  }
  if( orig.floodplain )
  {
    floodplain = new tFloodplain( *orig.floodplain );
    floodplain->setMeshPtr( mesh );
  }
  if( orig.stratGrid )
  {
    stratGrid = new tStratGrid( *orig.stratGrid );
    stratGrid->setMesh( mesh );
    stratGrid->updateConnect();
  }
  if( orig.loess )
    loess = new tEolian( *orig.loess );
  if( orig.strmMeander )
  {
    strmMeander = new tStreamMeander( *orig.strmMeander );
    strmMeander->setMeshPtr( mesh );
    strmMeander->setNetPtr( strmNet );
    strmMeander->setRandPtr( rand );
  }
  
  element_set_id = orig.element_set_id;
  element_set_description = orig.element_set_description;
  version = orig.version;
  
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
Initialize( string argument_string )
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
  tOption option( argument_string );
  
  // Say hello
  option.version();
  
  // Open main input file
  tInputFile inputFile( option.inputFile );
  
  // Get various options
  optNoDiffusion = inputFile.ReadBool( "OPTNODIFFUSION", false );
  optNoFluvial = inputFile.ReadBool( "OPTNOFLUVIAL", false );
  optNoUplift = inputFile.ReadBool( "OPTNOUPLIFT", false );
  optDetachLim = inputFile.ReadBool( "OPTDETACHLIM" );
  optDiffuseDepo = inputFile.ReadBool( "OPTDIFFDEP" );
  optVegetation = inputFile.ReadBool( "OPTVEG" );
  optForest = inputFile.ReadBool( "OPTFOREST", false );
  optFire = inputFile.ReadBool( "OPTFIRE", false );
  optFloodplainDep = inputFile.ReadBool( "OPTFLOODPLAIN" );
  optLoessDep = inputFile.ReadBool( "OPTLOESSDEP" );
  optMeander = inputFile.ReadBool( "OPTMEANDER" );
  optStratGrid = inputFile.ReadBool( "OPTSTRATGRID" ,false);
  optNonlinearDiffusion = inputFile.ReadBool( "OPT_NONLINEAR_DIFFUSION", false );
  optDepthDependentDiffusion = 
  inputFile.ReadBool( "OPT_DEPTH_DEPENDENT_DIFFUSION", false );
  optLandslides = inputFile.ReadBool( "OPT_LANDSLIDES", false );
  opt3DLandslides = inputFile.ReadBool( "OPT_3D_LANDSLIDES", false );
  optChemicalWeathering = inputFile.ReadBool( "CHEM_WEATHERING_LAW", false );
  optPhysicalWeathering = inputFile.ReadBool( "PRODUCTION_LAW", false );
  optTrackWaterSedTimeSeries = 
  inputFile.ReadBool( "OPT_TRACK_WATER_SED_TIMESERIES", false );
  
  
  // Create a random number generator for the simulation itself
  rand = new tRand( inputFile );
  
  // CREATE AND INITIALIZE OBJECTS
  //
  // Create (or read) model mesh
  if( !option.silent_mode )
    std::cout << "Creating mesh...\n";
  mesh = new tMesh<tLNode>( inputFile, option.checkMeshConsistency );
  
  // Initialize the lithology manager
  lithology_manager_.InitializeFromInputFile( inputFile, mesh );
  
  // Create and initialize output object
  if( !option.no_write_mode )
  {
    if( !option.silent_mode )
      std::cout << "Creating output handler...\n";
    output = new tLOutput<tLNode>( mesh, inputFile, rand );
  }
  
  // Create and initialize storm object
  storm = new tStorm( inputFile, rand, option.no_write_mode );
  
  // Create and initialize stream network object
  if( !option.silent_mode )
    std::cout << "Initializing stream network...\n";
  strmNet = new tStreamNet( *mesh, *storm, inputFile );
  optStreamLineBoundary = inputFile.ReadBool( "OPTSTREAMLINEBNDY", false );
  // option to convert streamlines from specified points along streamlines
  // into open boundary nodes; requires flow edges be set:
  if( optStreamLineBoundary )
  {
    tPtrList< tLNode > streamList;
    strmNet->FindStreamLines( inputFile, streamList );
    mesh->MakeStreamLineBoundaries( streamList );
    // do stuff in tStreamNet::UpdateNet(), but also need to do 
    // tStreamNet::InitFlowDirs(), so call functions piecemeal:
    strmNet->CalcSlopes();
    strmNet->InitFlowDirs(); // TODO: should all be done in call to updatenet
    strmNet->FlowDirs();
    strmNet->CheckNetConsistency();
    strmNet->MakeFlow( 0.0 );
  }

  // Create and initialize erosion/sedimentation object
  if( !option.silent_mode )
    std::cout << "Initializing erosion/transport/sedimentation...\n";
  erosion = new tErosion( mesh, inputFile, option.no_write_mode );
  
  // Create and initialize tectonics/baselevel object
  if( !optNoUplift )
  {
    if( !option.silent_mode )
      std::cout << "Initializing tectonics/baselevel...\n";
    uplift = new tUplift( inputFile );
  }
  
  // Create and initialize run timer object:
  time = new tRunTimer( inputFile, !option.silent_mode );
  
  // If applicable, create and initialize Vegetation object
  if( optVegetation )
  {
    if( !option.silent_mode )
      std::cout << "Initializing Vegetation module...\n";
    if( optFire )
    {
      if( optForest )
        vegetation = new tVegetation( mesh, inputFile, option.no_write_mode, 
                                     time, storm );
      else
        vegetation = new tVegetation( mesh, inputFile, option.no_write_mode, 
                                     time );
    }
    else
      vegetation = new tVegetation( mesh, inputFile );
  }
  
  // If applicable, create and initialize floodplain object
  if( optFloodplainDep )
  {
    if( !option.silent_mode )
      std::cout << "Initializing Floodplain module...\n";
    floodplain = new tFloodplain( inputFile, mesh );
  }
  
  // If applicable, create and initialize eolian deposition object
  if( optLoessDep )
  {
    if( !option.silent_mode )
      std::cout << "Initializing Eolian deposition module...\n";
    loess = new tEolian( inputFile );
  }
  
  // If applicable, create and initialize Stream Meander object
  if( optMeander )
  {
    if( !option.silent_mode )
      std::cout << "Initializing Stream Meander module...\n";
    strmMeander = new tStreamMeander( *strmNet, *mesh, inputFile, rand );
  }
  
  // If applicable, create and initialize Stratigraphy Grid object
  // and pass it to output
  if( optStratGrid ) {
    if (!optFloodplainDep)
      ReportFatalError("OPTFLOODPLAIN must be enabled.");
    if( !option.silent_mode )
      std::cout << "Initializing StratGrid module...\n";
    stratGrid = new tStratGrid (inputFile, mesh);
    if( output )
      output->SetStratGrid( stratGrid, strmNet );
  }
  
  // If applicable, set up tracking of water and sediment flux
  if( optTrackWaterSedTimeSeries && !option.no_write_mode )
  {
    water_sed_tracker_.InitializeFromInputFile( inputFile, mesh );
    erosion->ActivateSedVolumeTracking( &water_sed_tracker_ );
  }
  
  // Write output for time zero
  if( !option.silent_mode )
    std::cout << "Writing data for time zero...\n";
  if( output )
    output->WriteOutput( 0. );
  
  // Finish up initialization
  initialized = true;
  if( !option.silent_mode )
    std::cout << "******* Initialization done *******\n";
  
}


/**************************************************************************/
/**
 **  Initialize (argc, argv version)
 **
 **  This version of Initialize turns (argc, argv) into a string, and
 **  passes it to Initialize( string ).
 */
/**************************************************************************/
void childInterface::
Initialize( int argc, char **argv )
{
  string argstr;
  for( int i=1; i<argc; i++ )
  {
    argstr.append( argv[i] );
    if( i<(argc-1) )
      argstr.append( " " );
  }
  Initialize( argstr );
}
  
  
/**************************************************************************/
/**
**  VaryParameters
**
**  This method varies the parameters by amounts given be the values
**  in the input file (typically the standard deviation of measured values)
**  and the step size (for multiplication by the values in the input file).
**
**  NOTE: Need to set variable parameters "by hand" in this function, i.e.,
**  if you want to vary a parameter, you need to include a 
**  tInputFile::ReadDouble (or whatever) for that variable (with keyword),
**  find what object that variable is a member of, and use the appropriate 
**  "get" and "set" functions to change the value. You may need to write
**  new "get" and "set" functions also, as many parameters are not meant
**  to be changed and therefore don't have associated "get" and "set"
**  functions. And, for some parameters, you'll need to go through every
**  active node in the mesh to change the associated quantity. 
**  Have fun!
**
**  -SL, 10/2010
*/
/**************************************************************************/
vector<double> childInterface::
VaryParameters( const tInputFile &optFile, const double &delta, 
		tRand &rand, bool yesVary /*=true*/ )
{
  tList<double> paramL;
  // precipitation rate:
  double tmpVal = optFile.ReadDouble( "ST_PMEAN", false );
  if( tmpVal > 0.0 )
    {
      storm->GenerateStorm( 0.0 );
      double newRainrate = storm->getRainrate();
      if(yesVary)
	{
	  do
	    newRainrate += tmpVal * delta * ( rand.ran3() - 0.5 ) * 2.0;
	  while( newRainrate <= 0.0 );
	  storm->setRainrate( newRainrate );
	}
      paramL.insertAtBack( newRainrate );
    }
  // infiltration rate, also saturated hydraulic conductivity:
  tmpVal = optFile.ReadDouble( "INFILTRATION", false );
  if( tmpVal > 0.0 )
    {
      double newInfilt = strmNet->getInfilt();
      if(yesVary)
	{
	  do
	    newInfilt += tmpVal * delta * ( rand.ran3() - 0.5 ) * 2.0;
	  while( newInfilt <= 0.0 );
	  strmNet->setInfilt( newInfilt );
	}
      paramL.insertAtBack( newInfilt );
    }
  // depth scale for depth-dependent diffusion, also used with 
  // root cohesion:
  tmpVal = optFile.ReadDouble( "DIFFDEPTHSCALE", false );
  if( tmpVal > 0.0 )
    {
      double newDiffDepth = erosion->getDiffusionH();
      if(yesVary)
	{
	  do
	    newDiffDepth += tmpVal * delta * ( rand.ran3() - 0.5 ) * 2.0;
	  while( newDiffDepth <= 0.0 );
	  erosion->setDiffusionH( newDiffDepth );
	}
      paramL.insertAtBack( newDiffDepth );
    }
  // soil bulk density; also resets wet bulk density (done in "set"):
  tmpVal = optFile.ReadDouble( "SOILBULKDENSITY", false );
  if( tmpVal > 0.0 )
    {
      double newSoilDens = erosion->getSoilBulkDensity();
      if(yesVary)
	{
	  do
	    newSoilDens += tmpVal * delta * ( rand.ran3() - 0.5 ) * 2.0;
	  while( newSoilDens <= 0.0 );
	  erosion->setSoilBulkDensity( newSoilDens );
 	}
      paramL.insertAtBack( newSoilDens );
   }
  // soil depth:
  tmpVal = optFile.ReadDouble( "REGINIT", false );
  if( tmpVal > 0.0 )
    {
      tMesh<tLNode>::nodeListIter_t nI( mesh->getNodeList() );
      double newRegDepth = nI.FirstP()->getRegolithDepth();
      if(yesVary)
	{
	  do
	    newRegDepth += tmpVal * delta * ( rand.ran3() - 0.5 ) * 2.0;
	  while( newRegDepth <= 0.0 );
	  for( tLNode* cn = nI.FirstP(); nI.IsActive(); cn = nI.NextP() )
	    cn->setLayerDepth( 0, newRegDepth );
  	}
      paramL.insertAtBack( newRegDepth );
   }
  // maximum (initial) root cohesion; assumes partitioning between
  // lateral and vertical is constant, so resets lateral cohesion
  // from new value of vertical:
  tmpVal = optFile.ReadDouble( "MAXVERTROOTCOHESION", false );
  if( tmpVal > 0.0 )
    {
      double newVertCohesion = vegetation->ForestPtr()->getMVRC();
      if(yesVary)
	{
	  do
	    newVertCohesion += tmpVal * delta * ( rand.ran3() - 0.5 ) * 2.0;
	  while( newVertCohesion <= 0.0 );
	  double newLatCohesion = 
	    newVertCohesion * vegetation->ForestPtr()->getRSPartition();
	  vegetation->ForestPtr()->setMVRC( newVertCohesion );
	  vegetation->ForestPtr()->setMLRC( newLatCohesion );
	  tMesh<tLNode>::nodeListIter_t nI( mesh->getNodeList() );
	  tTrees *trees=0;
	  for( tLNode* cn = nI.FirstP(); nI.IsActive(); cn = nI.NextP() )
	    {
	      trees = cn->getVegCover().getTrees();
	      trees->setMaxRootStrength( -1.0 );
	      trees->TreesInitialize( cn, vegetation->ForestPtr(), 0.0 );
	      // 	  vegetation->ForestPtr()->RootStrengthInit( trees );
	    }
  	}
      paramL.insertAtBack( newVertCohesion );
    }
  // friction slope (or tan(phi)) for landsliding:
  tmpVal = optFile.ReadDouble( "FRICSLOPE", false );
  if( tmpVal > 0.0 )
    {
      double newFricSlope = erosion->getFricSlope();
      if(yesVary)
	{
	  do
	    newFricSlope += tmpVal * delta * ( rand.ran3() - 0.5 ) * 2.0;
	  while( newFricSlope <= 0.0 );
	  erosion->setFricSlope( newFricSlope );
  	}
      paramL.insertAtBack( newFricSlope );
    }
  const int sz = paramL.getSize();
  vector<double> paramV(sz);
  double param(0.0);
  for( int i=0;i<sz;++i )
    {
      paramL.removeFromFront( param );
      paramV[i] = param;
    }
  return paramV;
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
			
  //-------------HYDROLOGY--------------------------------------------
  // This is where flow directions are updated and discharges are
  // calculated:
  if( !optNoFluvial || optLandslides)
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

  //-------------CHEMICAL WEATHERING------------------------------
  // Do chemical weathering before physical weathering, which may be dependent on
  // the degree of chemical weathering 
  // what this does, whether it does anything,
  // is determined by the chemical weathering option, one of which is "None."
  // this option is read in the tErosion constructor:
  if( optChemicalWeathering )
    erosion->WeatherBedrock( stormPlusDryDuration );

  //-------------PHYSICAL WEATHERING------------------------------
  // Do physical weathering before diffusion, which may be dependent on thickness
  // and availability
  // what this does, whether it does anything,
  // is determined by the physical weathering option, one of which is "None."
  // this option is read in the tErosion constructor:
  if( optPhysicalWeathering )
    erosion->ProduceRegolith( stormPlusDryDuration, time->getCurrentTime() );

  //-------------DIFFUSION----------------------------------------
  //Diffusion is now before fluvial erosion in case the tools
  //detachment laws are being used.
  if( !optNoDiffusion )
    {
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
    }
  
  //-------------VEGETATION---------------------------------------
  // Vegetation moved up before fluvial erosion and (new) landsliding; 
  // latter depends on vegetation
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

  //-------------FLUVIAL------------------------------------------
  if( !optNoFluvial )
    {
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
    }
  
  //-------------LANDSLIDES---------------------------------------
  if( optLandslides )
    {
      if( opt3DLandslides )
	erosion->LandslideClusters3D( storm->getRainrate(),
				      time->getCurrentTime() );
      else
	erosion->LandslideClusters( storm->getRainrate(),
				    time->getCurrentTime() );
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
	
  //-------------MEANDERING-------------------------------------
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
	
	
  // Do interstorm...
  //Diffusion has now been moved to before fluvial erosion NMG 07/05
  // erosion.Diffuse( storm.getStormDuration() + storm.interstormDur(),
  // 		       optDiffuseDepo );
	
  erosion->UpdateExposureTime( stormPlusDryDuration );
	
  //----------------EOLIAN------------------------------------
  if( optLoessDep )
    loess->DepositLoess( mesh,
			 stormPlusDryDuration,
			 time->getCurrentTime() );
	
  //----------------TECTONICS---------------------------------
  if( !optNoUplift )
    {
      if( time->getCurrentTime() < uplift->getDuration() )
	uplift->DoUplift( mesh,
			  stormPlusDryDuration, 
			  time->getCurrentTime() );
    }
  
  if( optTrackWaterSedTimeSeries && water_sed_tracker_.IsActive() )
    water_sed_tracker_.WriteAndResetWaterSedTimeseriesData( time->getCurrentTime(),
							    stormPlusDryDuration );
		
  time->Advance( stormPlusDryDuration );
	
  if( output > 0 && time->CheckOutputTime() )
    output->WriteOutput( time->getCurrentTime() );
	
  if( output > 0 && output->OptTSOutput() ) output->WriteTSOutput();
  
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
	// why check option flags? if the pointers are
	// non-null, then we want to delete those objects
	// regardless of whether the option flag is set 
	// (SL, 12/10)
// 	if( optVegetation && vegetation ) {
	if( vegetation ) {
		delete vegetation;
		vegetation = NULL;
	}
// 	if( optFloodplainDep && floodplain ) {
	if( floodplain ) {
		delete floodplain;
		floodplain = NULL;
	}
// 	if( optLoessDep && loess ) {
	if( loess ) {
		delete loess;
		loess = NULL;
	}
// 	if( optMeander && strmMeander ) {
	if( strmMeander ) {
		delete strmMeander;
		strmMeander = NULL;
	}
// 	if( optStratGrid && stratGrid ) {
	if( stratGrid ) {
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
**  WARNING: THIS ONLY CHANGES ELEVATIONS; NEED TO CALL ALTERNATE ERODEP
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
 **  WARNING: THIS ONLY CHANGES ELEVATIONS; NEED TO CALL ALTERNATE ERODEP
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
  if( var_name.compare( 0,1,"x" )==0 || var_name.compare( 0,5,"nodex" )==0)
  {
    if(1) std::cout << "request for node x coordinates\n";
    return GetNodeXCoords();
  }
  if( var_name.compare( 0,1,"y" )==0 || var_name.compare( 0,5,"nodey" )==0)
  {
    if(1) std::cout << "request for node y coordinates\n";
    return GetNodeYCoords();
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
  else if( var_name.compare( 0,4,"land" )==0 )
	{
		if(1) std::cout << "request for landslides\n";
		return GetLandslideAreasVector();
	}
  else if( var_name.compare( 0,4,"load" )==0 )
	{
		if(1) std::cout << "request for loads\n";
		return GetLoads();
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
**  childInterface::GetLandslideAreasVector
**
**  This private method is used to create and return a vector of landslide
**  area values from tErosion. The vector is in order of occurrence, i.e.,
**  it is a time series. Note that the other private methods query 
**  instantaneous values at a node.
**
**  SL, Oct 2010
*/
/**************************************************************************/
std::vector<double> childInterface::GetLandslideAreasVector()
{
  std::vector<double> slideAreas( erosion->landslideAreas.getSize() );
  tListIter<double> aI( erosion->landslideAreas );
  int i=0;
  for( double* Ptr = aI.FirstP(); !aI.AtEnd(); Ptr = aI.NextP() )
    {
      slideAreas[i] = *Ptr;
      ++i;
    }
  return slideAreas;
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
    if(1) std::cout << "request to set elevs\n";
    SetNodeElevations( value_set );
  }
  else if( var_name.compare( 0,2,"kr")==0 || var_name.compare( 0,2,"KR" )==0 )
  {
    if(1) std::cout << "request to set KR\n";
    lithology_manager_.SetRockErodibilityValuesAtAllDepths( value_set );
  }
  else
  {
    std::cout << "Warning: unrecognized value set '" << var_name << "'\n";
    std::cout << "Request to set values ignored\n";
  }
}

/**************************************************************************/
/**
**  childInterface::setWriteOption( bool )
**
**  This method sets the value of the write option. 
**
**  If changing from "no write" to "write" then creates a tLOutput object 
**  and turns on timeseries writing in objects such as erosion. 
**
**  If changing from "write" to "no write" then deletes tLOutput object
**  and turns off timeseries writing in objects such as erosion.
**
**  Takes boolean: if true then write, if false then don't write.
**
**  SL, Oct. 2010
*/
/**************************************************************************/
void childInterface::setWriteOption( bool write_mode, tInputFile& inputFile )
{
  if( write_mode && output ) return; // already have output -> do nothing
  if( !write_mode && !output ) return; // already no output -> do nothing
  if( write_mode )
    { // want output and have none -> create output
      output = new tLOutput<tLNode>( mesh, inputFile, rand );
      storm->TurnOnOutput( inputFile );
      erosion->TurnOnOutput( inputFile );
      if( optTrackWaterSedTimeSeries )  
	{
	  water_sed_tracker_.InitializeFromInputFile( inputFile, mesh );
	  erosion->ActivateSedVolumeTracking( &water_sed_tracker_ );
	}
      if( vegetation && vegetation->FirePtr() )
	vegetation->FirePtr()->TurnOnOutput( inputFile );
    }
  else
    { // don't want output and have it -> destroy output
      delete output;
      output = NULL;
      storm->TurnOffOutput();
      erosion->TurnOffOutput();
      if( optTrackWaterSedTimeSeries )
	{
	  water_sed_tracker_.IsActive() = false;
	  erosion->DeactivateSedVolumeTracking();
	}
      if( vegetation && vegetation->FirePtr() )
	vegetation->FirePtr()->TurnOffOutput();
    }
}

/**************************************************************************/
/**
**  childInterface::WriteChildStyleOutput
**
**  Simply writes output in usual Child formats for the current time. 
**  Provides a way to call a "write" from outside the childInterface scope.
**  Note that the output object has to exist; this function won't create it.
**
**  SL, Oct. 2010
*/
/**************************************************************************/

void childInterface::WriteChildStyleOutput()
{
  if( output )
    output->WriteOutput( time->getCurrentTime() );
}

/**************************************************************************/
/**
**  childInterface::ChangeOption
**
**  Method to change options/switches on the fly.
**
**  SL, Nov. 2010
*/
/**************************************************************************/
void childInterface::ChangeOption( string option, int val )
{
  if( option.compare( 0,7,"no diff" )==0 )
    optNoDiffusion = ( val > 0 );
  if( option.compare( 0,7,"no fluv" )==0 )
    optNoFluvial = ( val > 0 );
  if( option.compare( 0,5,"no up" )==0 )
    optNoUplift = ( val > 0 );
  if( option.compare( 0,6,"detach" )==0 )
    optDetachLim = ( val > 0 );
  if( option.compare( 0,5,"flood" )==0 )
    optFloodplainDep = ( val > 0 );
  if( option.compare( 0,5,"loess" )==0 )
    optLoessDep = ( val > 0 );
  if( option.compare( 0,3,"veg" )==0 )
    optVegetation = ( val > 0 );
  if( option.compare( 0,4,"fire" )==0 )
    optFire = ( val > 0 );
  if( option.compare( 0,6,"forest" )==0 )
    optForest = ( val > 0 );
  if( option.compare( 0,7,"meander" )==0 )
    optMeander = ( val > 0 );
  if( option.compare( 0,7,"no depo" )==0 )
    optDiffuseDepo = ( val > 0 );
  if( option.compare( 0,5,"strat" )==0 )
    optStratGrid = ( val > 0 );
  if( option.compare( 0,5,"track" )==0 )
    optTrackWaterSedTimeSeries = ( val > 0 );
  if( option.compare( 0,6,"nonlin" )==0 )
    optNonlinearDiffusion = ( val > 0 );
  if( option.compare( 0,9,"depth dep" )==0 )
    optDepthDependentDiffusion = ( val > 0 );
  if( option.compare( 0,6,"landsl" )==0 )
    optLandslides = ( val > 0 );
  if( option.compare( 0,6,"3D lan" )==0 )
    opt3DLandslides = ( val > 0 );
  if( option.compare( 0,4,"chem" )==0 )
    optChemicalWeathering = ( val > 0 );
  if( option.compare( 0,4,"phys" )==0 )
    optPhysicalWeathering = ( val > 0 );
  if( option.compare( 0,6,"stream" )==0 )
    optStreamLineBoundary = ( val > 0 );

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
      if(0) std::cout << "Adjust elevations: node " << current_node->getPermID()
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


/**************************************************************************/
/**
 **  childInterface::GetLoads
 **
 **  Calculates and returns the weight of each rock and sediment column
 **  at each interior node. The output can be used by a lithosphere 
 **  flexure model to compute isostatic deflection. Returns a vector of
 **  loads in Newtons. Assumes gravitational acceleration given by GRAV,
 **  so change this if you want to apply it to another planet!
 **
 **  GT, Aug 2011
 */
/**************************************************************************/
vector<double> childInterface::
GetLoads()
{
	vector<double> the_loads( mesh->getNodeList()->getSize() );

	tLNode *current_node;
	tMesh<tLNode>::nodeListIter_t ni( mesh->getNodeList() );
	
	for( current_node=ni.FirstP(); ni.IsActive(); current_node=ni.NextP() )
	{
		int node_id = current_node->getPermID();
		if(0) std::cout << "Getting load for node " << node_id << std::endl;
	  double load = 0.0;
		double varea = current_node->getVArea();
		tListIter<tLayer> li( current_node->getLayersRefNC() );
		//int j=0;
		for( tLayer * layp = li.FirstP(); !li.AtEnd(); layp=li.NextP() )
		{
			//std::cout << " Layer bulk density=" << layp->getBulkDensity();
      //std::cout << " Thickness=" << layp->getDepth();
			load += layp->getDepth() * layp->getBulkDensity() * varea * GRAV;
			//std::cout << " sum so far=" << load << std::endl;
		}
		the_loads[node_id] = load;
		//std::cout << "Total for node " << node_id << "=" << load << std::endl;
	}
	
	return the_loads;

}


/**************************************************************************/
/**
 **  childInterface::GetNodeXCoords
 **
 **  GT, Aug 2011
 */
/**************************************************************************/
void childInterface::
GetNodeXCoords( vector<double> & x )
{
	if( x.size() != mesh->getNodeList()->getSize() )
		x.resize( mesh->getNodeList()->getSize() );
	
	tLNode *current_node;
	tMesh<tLNode>::nodeListIter_t ni( mesh->getNodeList() );
	
	for( current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP() )
	{
		int node_id = current_node->getPermID();
		x[node_id] = current_node->getX();
	}
	
}

/**************************************************************************/
/**
 **  childInterface::GetNodeYCoords
 **
 **  GT, Aug 2011
 */
/**************************************************************************/
void childInterface::
GetNodeYCoords( vector<double> & y )
{
	if( y.size() != mesh->getNodeList()->getSize() )
		y.resize( mesh->getNodeList()->getSize() );
	
	tLNode *current_node;
	tMesh<tLNode>::nodeListIter_t ni( mesh->getNodeList() );
	
	for( current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP() )
	{
		int node_id = current_node->getPermID();
		y[node_id] = current_node->getY();
	}
	
}

/**************************************************************************/
/**
 **  childInterface::GetNodeXCoords
 **
 **  GT, Aug 2011
 */
/**************************************************************************/
std::vector<double> childInterface::
GetNodeXCoords()
{
	tLNode *current_node;
	tMesh<tLNode>::nodeListIter_t ni( mesh->getNodeList() );
	std::vector<double> x( mesh->getNodeList()->getSize() );
	
	for( current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP() )
	{
		int node_id = current_node->getPermID();
		if(1) std::cout << "  In GetNodeXCoords, node " << node_id << " has X coord " << current_node->getX() << std::endl;
		x[node_id] = current_node->getX();
	}
	return x;
	
}

/**************************************************************************/
/**
 **  childInterface::GetNodeYCoords
 **
 **  GT, Aug 2011
 */
/**************************************************************************/
std::vector<double> childInterface::
GetNodeYCoords()
{
	tLNode *current_node;
	tMesh<tLNode>::nodeListIter_t ni( mesh->getNodeList() );
	std::vector<double> y( mesh->getNodeList()->getSize() );
	
	for( current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP() )
	{
		int node_id = current_node->getPermID();
		if(1) std::cout << "  In GetNodeYCoords, node " << node_id << " has Y coord " << current_node->getY() << std::endl;
		y[node_id] = current_node->getY();
	}
	return y;
	
}

