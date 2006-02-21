/**************************************************************************/
/**
**                            C H I L D
**
**       CHANNEL-HILLSLOPE INTEGRATED LANDSCAPE DEVELOPMENT MODEL
**
**				CU EXECUTABLE RELEASE VERSION 0601 alpha
**
**  Designed and created by Gregory E. Tucker, Stephen T. Lancaster,
**      Nicole M. Gasparini, and Rafael L. Bras
**
**
**  @file   childrmain.cpp
**  @brief  This file contains the main() routine that handles
**          top-level initialization and implements the main
**          time-loop. It is equivalent to childmain.cpp but lacks
**          vegetation, stratgrid, floodplain, and eolian.
**
**  NOTE: This source code is copyrighted material. It is distributed
**        solely for noncommercial research and educational purposes
**        only. Use in whole or in part for commercial purposes without
**        a written license from the copyright holder(s) is expressly
**        prohibited. Copies of this source code or any of its components
**        may not be transferred to any other individuals or organizations
**        without written consent.
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

/* set traps for some floating point exceptions on Linux */
#include "trapfpe.h"
#include "Inclusions.h"
#include "tOption/tOption.h"

#include "tMeshList/tMeshList.h"

Predicates predicate;

int main( int argc, char **argv )
{
   bool optDetachLim,      // Option for detachment-limited erosion only
        optDiffuseDepo;    // Option for deposition / no deposition by diff'n
   tVegetation *dummy_veg(0);  // -> vegetation object ... in this case, just a dummy

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
   tRand rand( inputFile );

   // Create and initialize objects:
   std::cout << "Creating mesh...\n";
   tMesh<tLNode> mesh( inputFile, option.checkMeshConsistency );

   std::cout << "Creating output files...\n";
   tLOutput<tLNode> output( &mesh, inputFile, &rand );

   tStorm storm( inputFile, &rand );
   std::cout << "Creating stream network...\n";
   tStreamNet strmNet( mesh, storm, inputFile );
   tErosion erosion( &mesh, inputFile );
   tUplift uplift( inputFile );

   // Get various options
   optDetachLim = inputFile.ReadBool( "OPTDETACHLIM" );
   optDiffuseDepo = inputFile.ReadBool( "OPTDIFFDEP" );

   std::cout << "Writing data for time zero...\n";
   tRunTimer time( inputFile, !option.silent_mode );
   output.WriteOutput( 0. );
   std::cout << "Initialization done.\n";


   /**************** MAIN LOOP ******************************************\
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
   while( !time.IsFinished() )
   {
      if(0) //debug
      	std::cout << "         " << std::endl;
      time.ReportTimeStatus();

      // Do storm...
      storm.GenerateStorm( time.getCurrentTime(),
                           strmNet.getInfilt(), strmNet.getSoilStore() );

      strmNet.UpdateNet( time.getCurrentTime(), storm );

      //Diffusion is now before fluvial erosion in case the tools
      //detachment laws are being used.
      erosion.Diffuse( storm.getStormDuration() + storm.interstormDur(),
                       optDiffuseDepo );

      if( optDetachLim )
          erosion.ErodeDetachLim( storm.getStormDuration(), &strmNet,
				  dummy_veg );
      else{
         erosion.DetachErode( storm.getStormDuration(), &strmNet,
                             time.getCurrentTime(), dummy_veg );
		 // To use tools rules, you must use DetachErode2 NMG 07/05
         //erosion.DetachErode2( storm.getStormDuration(), &strmNet,
         //                      time.getCurrentTime(), dummy_veg );
      }
      
      if( time.getCurrentTime() < uplift.getDuration() )
          uplift.DoUplift( &mesh,
                           storm.getStormDuration() + storm.interstormDur(), 
						   time.getCurrentTime() );

      time.Advance( storm.getStormDuration() + storm.interstormDur() );

      if( time.CheckOutputTime() )
          output.WriteOutput( time.getCurrentTime() );

      if( output.OptTSOutput() ) output.WriteTSOutput();

   } // end of main loop

   delete dummy_veg;

   return 0;
}

