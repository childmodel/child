/**************************************************************************\
**
**       TODDLER: Slimmed-down version of CHILD without Vegetation,
**                Stream Meandering, or mesh densification.
**
**                            C H I L D
**
**       CHANNEL-HILLSLOPE INTEGRATED LANDSCAPE DEVELOPMENT MODEL
**
**                      VERSION 2.1.2, JUNE 2000
**
**  Designed and created by Gregory E. Tucker, Stephen T. Lancaster,
**      Nicole M. Gasparini, and Rafael L. Bras, Department of Civil
**      and Environmental Engineering, Massachusetts Institute of
**      Technology, Cambridge, MA, USA.
**
**  MAIN.CPP -- This file contains the main() routine that handles
**              top-level initialization and implements the main
**              time-loop.
**
**  NOTE: This source code is copyrighted material. It is distributed
**        solely for noncommercial research and educational purposes
**        only. Use in whole or in part for commercial purposes without 
**        a written license from the copyright holder(s) is expressly
**        prohibited. Copies of this source code or any of its components 
**        may not be transferred to any other individuals or organizations
**        without written consent. Copyright (C) Massachusetts Institute
**        of Technology, 1997-2000. All rights reserved.
**
**  For information regarding this program, please contact Greg Tucker at:
**
**     Before July 1, 2000:               After July 1, 2000:
**       Dept. of Civil and                 School of Geography
**         Environmental Engr'ing           University of Oxford
**       MIT, Bldg 48, Room 429             Mansfield Road
**       Cambridge, MA 02139 USA            Oxford OX1 3TB United Kingdom
**
**  $Id: toddlermain.cpp,v 1.2 2000-06-19 17:39:35 gtucker Exp $
\**************************************************************************/


#include "Inclusions.h"
#include "tFloodplain/tFloodplain.h"
#include "tEolian/tEolian.h"

Predicates predicate;


main( int argc, char **argv )
{
   int silent_mode,       // Option for silent mode (no time output to stdout)
       optDetachLim,      // Option for detachment-limited erosion only
       optFloodplainDep,  // Option for floodplain (overbank) deposition
       optLoessDep,       // Option for eolian deposition
       optDiffuseDepo;    // Option for deposition / no deposition by diff'n
   tFloodplain *floodplain;  // -> floodplain object
   tEolian *loess;           // -> eolian deposition object
   
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
   if( argc<2 )
   {
      cerr << "Usage: " << argv[0] << " <input file>" << endl;
      ReportFatalError( "You need to give the name of an input file." );
   }

   // Check whether we're in silent mode
   silent_mode = ( argc>2 && argv[2][1]=='s' );
   
   // Say hello
   cout << "\nThis is TODDLER, version " << VERSION << endl << endl;
   
   // Open main input file
   tInputFile inputFile( argv[1] );

   // Create and initialize objects:
   cout << "Creating mesh...\n";
   tMesh<tLNode> mesh( inputFile );
   cout << "Creating output files...\n";
   tLOutput<tLNode> output( &mesh, inputFile );
   tStorm storm( inputFile );
   cout << "Creating stream network...\n";
   tStreamNet strmNet( mesh, storm, inputFile );
   tErosion erosion( &mesh, inputFile );
   tUplift uplift( inputFile );
   cout << "Writing data for time zero...\n";
   tRunTimer time( inputFile, !silent_mode );
   output.WriteOutput( 0 );
   cout << "Initialization done.\n";

   // Get various options
   optDetachLim = inputFile.ReadItem( optDetachLim, "OPTDETACHLIM" );
   optDiffuseDepo = inputFile.ReadItem( optDiffuseDepo, "OPTDIFFDEP" );
   optFloodplainDep = inputFile.ReadItem( optFloodplainDep, "OPTFLOODPLAIN" );
   optLoessDep = inputFile.ReadItem( optLoessDep, "OPTLOESSDEP" );

   // If applicable, create floodplain object
   if( optFloodplainDep )
       floodplain = new tFloodplain( inputFile, &mesh );

   // If applicable, create eolian deposition object
   if( optLoessDep )
       loess = new tEolian( inputFile );

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
      time.ReportTimeStatus();

      // Do storm...
      storm.GenerateStorm( time.getCurrentTime(),
                           strmNet.getInfilt(), strmNet.getSoilStore() );
      cout << storm.getRainrate() << " " << storm.getStormDuration() << " "
           << storm.interstormDur() << endl;

      strmNet.UpdateNet( time.getCurrentTime(), storm );
      
      if( optDetachLim )
          erosion.ErodeDetachLim( storm.getStormDuration() );
      else
          erosion.DetachErode( storm.getStormDuration(), &strmNet,
                               time.getCurrentTime() );

      if( optFloodplainDep )
          floodplain->DepositOverbank( storm.getRainrate(),
                                       storm.getStormDuration(),
                                       time.getCurrentTime() );

      // Do interstorm...
      erosion.Diffuse( storm.getStormDuration() + storm.interstormDur(),
      optDiffuseDepo );

      erosion.UpdateExposureTime( storm.getStormDuration() + 
                                      storm.interstormDur() );

      if( optLoessDep )
          loess->DepositLoess( &mesh, 
                               storm.getStormDuration()+storm.interstormDur(),
                               time.getCurrentTime() );

      if( time.getCurrentTime() < uplift.getDuration() )
          uplift.DoUplift( &mesh,
                           storm.getStormDuration() + storm.interstormDur() );

      time.Advance( storm.getStormDuration() + storm.interstormDur() );

      if( time.CheckOutputTime() )
          output.WriteOutput( time.getCurrentTime() );

   } // end of main loop
   
}

