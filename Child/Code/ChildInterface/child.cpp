#include "child.h"

#define VERBOSE (false)

void Child::Initialize (std::string argument_string) {
  
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
  else
    vegetation = NULL;
  
  // If applicable, create and initialize floodplain object
  if( optFloodplainDep )
  {
    if( !option.silent_mode )
      std::cout << "Initializing Floodplain module...\n";
    floodplain = new tFloodplain( inputFile, mesh );
  }
  else
    floodplain = NULL;
  
  // If applicable, create and initialize eolian deposition object
  if( optLoessDep )
  {
    if( !option.silent_mode )
      std::cout << "Initializing Eolian deposition module...\n";
    loess = new tEolian( inputFile );
  }
  else
    loess = NULL;
  
  // If applicable, create and initialize Stream Meander object
  if( optMeander )
  {
    if( !option.silent_mode )
      std::cout << "Initializing Stream Meander module...\n";
    strmMeander = new tStreamMeander( *strmNet, *mesh, inputFile, rand );
  }
  else
    strmMeander = NULL;
  
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
  else
    stratGrid = NULL;
  
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
 **  Run
 **
 **  This method runs CHILD for a duration specified on the command line or,
 **  if that duration is <=0, for the duration specified in the previously
 **  read input file.
 */
/**************************************************************************/

void Child::Run(double run_duration) {
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

void Child::CleanUp() {
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

/**************************************************************************/
/**
 **  RunOneStorm
 **
 **  This method executes the model for one storm, calling each of several
 **  subroutines in turn. It returns the updated simulation time.
 */
/**************************************************************************/

double Child::RunOneStorm() {
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
        erosion->DiffuseNonlinear( stormPlusDryDuration, optDiffuseDepo, time->getCurrentTime() );
    }
    else
      erosion->Diffuse( stormPlusDryDuration, optDiffuseDepo, time->getCurrentTime() );
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
  if(0) std::cout << "Calculating fluvial erosion and transport ...\n" << std::flush;
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
    std::cout << "Erosion::Done.." << std::flush;
	
  
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
**  childInterface::MaskNodesBelowElevation
**
**  This public function sets on or off the "mask" status of individual
**  CHILD nodes, according to whether they are above or below elevation
**  "elev". Its meant to be used to mask out nodes below sea level, so
**  submarine nodes can be handled by another model (like SedFlux).
**
**  GT, Sep 2010
*/
/**************************************************************************/
void Child::MaskNodesBelowElevation(double elev) {
   tLNode *current_node;
   tMesh<tLNode>::nodeListIter_t ni( mesh->getNodeList() );
   
   for( current_node=ni.FirstP(); ni.IsActive(); current_node=ni.NextP() )
   {
      if( current_node->getZ() <= elev )
         current_node->setMask( true );
      else
         current_node->setMask( false );
   }
   
}

void Child::CopyNodeElevations (double * const dest)
{
  tLNode *current_node;
  tMesh<tLNode>::nodeListIter_t ni (mesh->getNodeList());
   
  for (current_node = ni.FirstP(); !ni.AtEnd(); current_node = ni.NextP())
  {
    if (VERBOSE)
      std::cout << "node " << current_node->getPermID()
        << " z=" << current_node->getZ() << std::endl;

    dest[current_node->getPermID()] = current_node->getZ();
  }

  if (VERBOSE)
    fprintf (stderr, "Number of elevation nodes is %d\n",
        mesh->getNodeList ()->getSize ());
   
  return;
}

void Child::CopyNodeErosion (double * const dest)
{
   tLNode *current_node;
   tMesh<tLNode>::nodeListIter_t ni (mesh->getNodeList());
      
   for (current_node = ni.FirstP(); !ni.AtEnd(); current_node = ni.NextP())
   {
      dest[current_node->getPermID()] = current_node->getCumulativeEroDep();
      current_node->ResetCumulativeEroDep();
   }

   return;
}

void Child::CopyNodeDischarge (double * const dest)
{
  tLNode *current_node;
  tMesh<tLNode>::nodeListIter_t ni (mesh->getNodeList());
   
  for (current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP())
  {
    if (VERBOSE)
      std::cout << "NOde " << current_node->getPermID()
        << " Q=" << current_node->getQ() << std::endl;

     dest[current_node->getPermID()] = current_node->getQ ();
  }
   
  return;
}

void Child::CopyNodeSedimentFlux (double * const dest)
{
  tLNode *current_node;
  tMesh<tLNode>::nodeListIter_t ni (mesh->getNodeList());
   
  for (current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP())
  {
    if (VERBOSE)
      std::cout << "noDE " << current_node->getPermID()
        << " Qs=" << current_node->getQsin() << std::endl;

    dest[current_node->getPermID()] = current_node->getQsin ();
  }
   
  return;
}

void Child::SetNodeElevations (const double * elevations) {
   tLNode *current_node;
   tMesh<tLNode>::nodeListIter_t ni (mesh->getNodeList());
   
   for (current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP())
   {
      if (VERBOSE)
        std::cout << "node " << current_node->getPermID()
                  << " changing z from " << current_node->getZ() 
                  << " to " << elevations[current_node->getPermID()]
                  << std::endl;
      current_node->setZ (elevations[current_node->getPermID ()]);
   }
}

void Child::SetNodeUplift (const double * uplift) {
   tLNode *current_node;
   tMesh<tLNode>::nodeListIter_t ni (mesh->getNodeList());
   
   for (current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP())
   {
      if(VERBOSE) std::cout << "Adjust elevations: node " << current_node->getPermID()
                      << " changing z from " << current_node->getZ() 
                      << " to " 
                      << current_node->getZ() + uplift[current_node->getPermID()] 
                      << std::endl;
      current_node->ChangeZ (uplift[current_node->getPermID()]);
   }
}

