/**************************************************************************/
/**
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
**  @file toddlermain.cpp
**  @brief This file contains the main() routine that handles
**         top-level initialization and implements the main
**         time-loop.
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
**       School of Geography and the Environment
**       University of Oxford
**       Mansfield Road
**       Oxford OX1 3TB United Kingdom
**
**  $Id: toddlermain.cpp,v 1.18 2004-03-24 14:54:32 childcvs Exp $
*/
/**************************************************************************/

/* set traps for some floating point exceptions on Linux */
#include "trapfpe.h"
#include "Inclusions.h"
#include "tFloodplain/tFloodplain.h"
#include "tEolian/tEolian.h"

Predicates predicate;


int main( int argc, char **argv )
{
   bool silent_mode;      // Option for silent mode (no time output to stdout)
   int optDetachLim,      // Option for detachment-limited erosion only
       optFloodplainDep,  // Option for floodplain (overbank) deposition
       optLoessDep,       // Option for eolian deposition
       optVegetation,     // Option for dynamic vegetation cover
       optDiffuseDepo;    // Option for deposition / no deposition by diff'n
   tVegetation *vegetation(0);  // -> vegetation object
   tFloodplain *floodplain(0);  // -> floodplain object
   tEolian *loess(0);           // -> eolian deposition object
   
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
   silent_mode = BOOL( argc>2 && argv[2][1]=='s' );
   
   // Say hello
   cout << "\nThis is TODDLER, version " << CHILD_VERSION 
	<< " (compiled " __DATE__ " " __TIME__ ")"
	<< endl << endl;
   
   // Open main input file
   tInputFile inputFile( argv[1] );

   // Create a random number generator for the simulation itself
   tRand rand( inputFile );
   // Create and initialize objects:
   cout << "Creating mesh...\n";
   tMesh<tLNode> mesh( inputFile );
   cout << "Creating output files...\n";
   tLOutput<tLNode> output( &mesh, inputFile, &rand );
   tStorm storm( inputFile, &rand );
   cout << "Creating stream network...\n";
   tStreamNet strmNet( mesh, storm, inputFile );
   tErosion erosion( &mesh, inputFile );
   tUplift uplift( inputFile );

   // Get various options
   optDetachLim = inputFile.ReadItem( optDetachLim, "OPTDETACHLIM" );
   optDiffuseDepo = inputFile.ReadItem( optDiffuseDepo, "OPTDIFFDEP" );
   optVegetation = inputFile.ReadItem( optVegetation, "OPTVEG" );
   optFloodplainDep = inputFile.ReadItem( optFloodplainDep, "OPTFLOODPLAIN" );
   optLoessDep = inputFile.ReadItem( optLoessDep, "OPTLOESSDEP" );

   // If applicable, create Vegetation object
   if( optVegetation )
       vegetation = new tVegetation( &mesh, inputFile );

   // If applicable, create floodplain object
   if( optFloodplainDep )
       floodplain = new tFloodplain( inputFile, &mesh );

   // If applicable, create eolian deposition object
   if( optLoessDep )
       loess = new tEolian( inputFile );

   cout << "Writing data for time zero...\n";
   tRunTimer time( inputFile, BOOL(!silent_mode) );
   output.WriteOutput( 0. );
   cout << "Initialization done.\n";

   // Option for time series output (IN PROGRESS)
   /*   switch( optTSOutput ){
   case 1:   // Volume output each N years.
     if( time.CheckTSOutputTime() )
       output.WriteVolOutput();
     break;
   case 2:   // Volume and vegetation cover output each N years.
     cout << "here" << endl;
     if( time.CheckTSOutputTime() ){
       cout << "there" << endl;
       output.WriteVolVegOutput();}
     break;
   case 3:   // All data at each storm.
     output.WriteTSOutput();
     break;      
   case 0:   // No additional timeseries output.
     break;
   default:  // Invalid option.
     ReportFatalError( "The input file contains an invalid value for 
OptTSOutput." );
   }   */

#if 0
   // Temporary -- for canyon range stuff
   cout << "Densifying around fold strike ...\n";
   int optFoldDens = inputFile.ReadItem( optFoldDens, "OPTFOLDDENS" );
   if( optFoldDens )
     {
       double foldDensYmin = inputFile.ReadItem( foldDensYmin, "FOLDDENSYMIN" );
       double foldDensYmax = inputFile.ReadItem( foldDensYmax, "FOLDDENSYMAX" );
       double fault = inputFile.ReadItem( fault, "FAULTPOS" );
       for( int i=1; i<=optFoldDens; i++ )
	 {
	   tPtrList<tLNode> foldPoints;
	   tMesh<tLNode>::nodeListIter_t ni( mesh.getNodeList() );
	   tLNode *cn;
	   for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
	     {
	     if( cn->getY()>=foldDensYmin && cn->getY()<=foldDensYmax &&
		 cn->getX()>=7500 && cn->getX()<=92500 )
	       foldPoints.insertAtBack( cn );
	     if( cn->getY()<fault )
	       cn->setAlluvThickness( 100000 );
	     }
	   tPtrListIter<tLNode> fpi( foldPoints );
	   for( cn=fpi.FirstP(); !(fpi.AtEnd()); cn=fpi.NextP() )
	     mesh.AddNodesAround( cn, 0.0 );
	   foldPoints.Flush();
	 }
#endif

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
      cout
	<< "Storm: "
	<< storm.getRainrate() << " " << storm.getStormDuration() << " "
	<< storm.interstormDur() << endl;

      strmNet.UpdateNet( time.getCurrentTime(), storm );
      
      if( optDetachLim )
          erosion.ErodeDetachLim( storm.getStormDuration(), &strmNet,
				  vegetation );
      else
          erosion.DetachErode( storm.getStormDuration(), &strmNet,
                               time.getCurrentTime(), vegetation );

      if( optFloodplainDep )
          floodplain->DepositOverbank( storm.getRainrate(),
                                       storm.getStormDuration(),
                                       time.getCurrentTime() );

      if( optVegetation )
	  vegetation->UpdateVegetation( &mesh, storm.getStormDuration(),
					storm.interstormDur() );


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

      if( output.OptTSOutput() ) output.WriteTSOutput();

      /* IN PROGRESS
      switch( optTSOutput ){
      case 1:   // Volume output each N years.
        if( time.CheckTSOutputTime() )
  output.WriteVolOutput();
break;
      case 2:   // Volume and vegetation cover output each N years.
if( time.CheckTSOutputTime() )
  output.WriteVolVegOutput();
break;
      case 3:   // All data at each storm.
output.WriteTSOutput();
break;      
      case 0:   // No additional timeseries output.
break;
      default:  // Invalid option.
ReportFatalError( "The input file contains an invalid value for OptTSOutput." 
*/ 

     /*nodeListIter_t ni( mesh.getNodeList() );
      tLNode *cn;
      for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
	{
	  if( cn->getY()<25 && cn->getX()>250 && cn->getDrArea()>1000 )
	    cn->TellAll();
	}*/

   } // end of main loop

   delete vegetation;
   delete floodplain;
   delete loess;
   
   return 0;
}

