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

Predicates predicate;

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
   optFloodplainDep = inputFile.ReadBool( "OPTFLOODPLAIN" );
   optLoessDep = inputFile.ReadBool( "OPTLOESSDEP" );
   optMeander = inputFile.ReadBool( "OPTMEANDER" );
   optStratGrid = inputFile.ReadBool( "OPTSTRATGRID" ,false);
   optNonlinearDiffusion = inputFile.ReadBool( "OPT_NONLINEAR_DIFFUSION", false );
   
   // If applicable, create Vegetation object
   if( optVegetation )
       vegetation = new tVegetation( mesh, inputFile );

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

   std::cout << "Writing data for time zero...\n";
   time = new tRunTimer( inputFile, !option.silent_mode );
   output->WriteOutput( 0. );
   std::cout << "Initialization done.\n";

}



void childInterface::
RunOneStorm()
{
	
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
	if(0) //DEBUG
		std::cout<< "Storm: "<< storm->getRainrate() << " " << storm->getStormDuration() << " "
			<< storm->interstormDur() << std::endl;
	
	strmNet->UpdateNet( time->getCurrentTime(), *storm );
	if(0) //DEBUG
		std::cout << "UpdateNet::Done.." << std::endl;
	
	if(0) //DEBUG
	{
		tMesh< tLNode >::nodeListIter_t mli( mesh->getNodeList() );  // gets nodes from the list
		tLNode * cn;
		for( cn=mli.FirstP(); mli.IsActive(); cn=mli.NextP() )
		{
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
		    if( cn->getID()==8121 || cn->getID()==8122 ) 
			{
				tEdge * debugedg = cn->getFlowEdg();
				tLNode * nbr = static_cast<tLNode *>(debugedg->getDestinationPtrNC());
				std::cout<<"Childmain 2: node "<<cn->getID()<<" edge "<<debugedg->getID()<<" downstream nbr "<<nbr->getID()<<std::endl;
				std::cout<<"z "<<cn->getZ()<<" dsn z "<<nbr->getZ()<<std::endl;
			}
		}
		
	}
	
	//Diffusion is now before fluvial erosion in case the tools
	//detachment laws are being used.
	if( optNonlinearDiffusion )
		erosion->DiffuseNonlinear( storm->getStormDuration() + storm->interstormDur(),
								   optDiffuseDepo );
	else
		erosion->Diffuse( storm->getStormDuration() + storm->interstormDur(),
						  optDiffuseDepo );
	
	
	if( optDetachLim )
		erosion->ErodeDetachLim( storm->getStormDuration(), strmNet,
								 vegetation );
	else{
		erosion->DetachErode( storm->getStormDuration(), strmNet,
							  time->getCurrentTime(), vegetation );
		// To use tools rules, you must use DetachErode2 NMG 07/05
		//erosion.DetachErode2( storm.getStormDuration(), &strmNet,
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
		    if( cn->getID()==8121 || cn->getID()==8122 ) 
			{
				tEdge * debugedg = cn->getFlowEdg();
				tLNode * nbr = static_cast<tLNode *>(debugedg->getDestinationPtrNC());
				std::cout<<"Childmain 3: node "<<cn->getID()<<" edge "<<debugedg->getID()<<" downstream nbr "<<nbr->getID()<<std::endl;
				std::cout<<"z "<<cn->getZ()<<" dsn z "<<nbr->getZ()<<std::endl;
			}
		}
		
	}
	
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
	        floodplain->UpdateMainChannelHeight( time->getCurrentTime(),
												 strmNet->getInletNodePtrNC() );
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
			vegetation->GrowVegetation( mesh, storm->interstormDur() );
		else
			vegetation->UpdateVegetation( mesh, storm->getStormDuration(),
										  storm->interstormDur() );
	}
#undef NEWVEG
	
	// Do interstorm...
	//Diffusion has now been moved to before fluvial erosion NMG 07/05
	// erosion.Diffuse( storm.getStormDuration() + storm.interstormDur(),
	// 		       optDiffuseDepo );
	
	erosion->UpdateExposureTime( storm->getStormDuration() +
								 storm->interstormDur() );
	
	if( optLoessDep )
		loess->DepositLoess( mesh,
							 storm->getStormDuration()+storm->interstormDur(),
							 time->getCurrentTime() );
	
	if( time->getCurrentTime() < uplift->getDuration() )
		uplift->DoUplift( mesh,
						  storm->getStormDuration() + storm->interstormDur(), 
						  time->getCurrentTime() );
	
	time->Advance( storm->getStormDuration() + storm->interstormDur() );
	
	if( time->CheckOutputTime() )
		output->WriteOutput( time->getCurrentTime() );
	
	if( output->OptTSOutput() ) output->WriteTSOutput();
		
}

void childInterface::
CleanUp()
{
	delete rand;
	delete mesh;
	delete output;
	delete storm;
	delete strmNet;
	delete erosion;
	delete uplift;
	delete time;
	delete vegetation;
	delete floodplain;
	delete loess;
	delete strmMeander;
	delete stratGrid;
	
}


childInterface::
~childInterface()
{
	CleanUp();
}
