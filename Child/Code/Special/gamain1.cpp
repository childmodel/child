
/******************************************************************

  gamain.cpp: specialized main file for Geoarchaeology simulations

  Greg Tucker, Sept. 1999

*****************************************************************/

/*#include "tList/tList.h"
#include "tGridList/tGridList.h"
#include "tArray/tArray.h"
#include "tPtrList/tPtrList.h"
#include "GridElements/gridElements.h"
#include "tListInputData/tListInputData.h"
#include "tGrid/tGrid.h"
#include "tLNode/tLNode.h"
#include "Erosion/erosion.h"
#include "tUplift/tUplift.h"
#include "tOutput/tOutput.h"
#include "tStorm/tStorm.h"
#include "tStreamNet/tStreamNet.h"
#include "tRunTimer/tRunTimer.h"
#include <iostream.h>*/
#include "Inclusions.h"
#include "tFloodplain/tFloodplain.h"
#include "tEolian/tEolian.h"



Predicates predicate;


main( int argc, char **argv )
{
   int silent_mode,       // Option for silent mode (no time output to stdout)
       optDetachLim,      // Option for detachment-limited erosion only
       optMeander,        // Option for stream meandering
       optFloodplainDep,  // Option for floodplain (overbank) deposition
       optLoessDep,       // Option for eolian deposition
       optDiffuseDepo;    // Option for deposition / no deposition by diff'n
   tStreamMeander *meander;  // -> meander object
   tFloodplain *floodplain;  // -> floodplain object
   tEolian *loess;           // -> eolian deposition object

   ofstream oefile;
   
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
   cout << "\nThis is CHILD, version " << VERSION << endl << 
       "Geoarchaeology special version 1.0" << endl << endl;
   
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
   optMeander = inputFile.ReadItem( optMeander, "OPTMNDR" );
   optDiffuseDepo = inputFile.ReadItem( optDiffuseDepo, "OPTDIFFDEP" );
   optFloodplainDep = inputFile.ReadItem( optFloodplainDep, "OPTFLOODPLAIN" );
   optLoessDep = inputFile.ReadItem( optLoessDep, "OPTLOESSDEP" );

   // If applicable, create stream meander object
   if( optMeander )
       meander = new tStreamMeander( strmNet, mesh, inputFile );

   // If applicable, create floodplain object
   if( optFloodplainDep )
       floodplain = new tFloodplain( inputFile, &mesh );

   // If applicable, create eolian deposition object
   if( optLoessDep )
       loess = new tEolian( inputFile );

   // For Geoarchaeology special application
   double kr = inputFile.ReadItem( kr, "KR" );
   double drop = inputFile.ReadItem( drop, "GA_VALDROP" );
   double inletElev = inputFile.ReadItem( inletElev, "GA_INLETELEV" );
   double meanInletElev = inletElev;
   double period = inputFile.ReadItem( period, "GA_PERIOD" );
   int optwave = inputFile.ReadItem( optwave, "GA_OPTWAVE" );
   double amplitude = inputFile.ReadItem( amplitude, "GA_AMPLITUDE" );
   double tpeak, mr, mf, ttime, oldar1, noise0, noise1;
   int numpts; //number of points in the floodplain curve
   int fpindex=0;//where are you in the floodplain data
   double fpslope;//slope of floodplain curve
   double chanslp;//slope of channel, read in if optwave==2
   tArray<double> fpht;
   tArray<double> fptime;
   if( optwave==0 ) period = 2.0 * PI / period;
   else if( optwave==1 )
   {
      tpeak = inputFile.ReadItem( tpeak, "GA_TPEAK" );
      if( tpeak<=0.0 || tpeak>=1.0 )
          ReportFatalError("GA_TPEAK must be between 0 and 1 (not inclusive");
      tpeak = tpeak*period;
      mr = amplitude/tpeak;
      mf = amplitude/(period-tpeak);
      oldar1=0;
      noise0=0;
      noise1=0;
      oefile.open("Geoarch/outletelev");
      
   } 
   else if( optwave==2){
      numpts=inputFile.ReadItem( numpts, "NUMFLDPLNPTS" );
      fpht.setSize( numpts );
      fptime.setSize( numpts );
      int i=0;
      char add='1';
      char add2='0';
      char name[30];
      double help;
      double inittime;
      chanslp=inputFile.ReadItem(chanslp, "CHANSLOPE" );
      cout<<"channel slope is "<<chanslp<<endl;
      inittime=inputFile.ReadItem(inittime, "INPUTTIME" );
      while (i<numpts){
         if(i<9){
            strcpy(name, "FLDPLNTIME" );
            strcat(name, &add );
            help=inputFile.ReadItem(help,name);
            fptime[i]=help+inittime;
            cout<<"index "<<i<<" fldplntime "<<fptime[i];
            strcpy(name, "FLDPLNHT" );
            strcat(name, &add );
            help=inputFile.ReadItem(help,name);
            fpht[i]=help;
            cout<<" fldplnht "<<fpht[i]<<endl;
            i++;
            add++;
         }
         if(i>=9){
            add='1';
            strcpy(name, "FLDPLNTIME" );
            strcat(name, &add );
            strcat(name, &add2 );
            help=inputFile.ReadItem(help,name);
            fptime[i]=help;
            cout<<"index "<<i<<" fldplntime "<<fptime[i];
            strcpy(name, "FLDPLNHT" );
            strcat(name, &add );
            strcat(name, &add2 );
            help=inputFile.ReadItem(help,name);
            fpht[i]=help;
            cout<<" fldplnht "<<fpht[i]<<endl;
            i++;
            add2++;
         }
            
            
      }
      fpslope=(fpht[fpindex+1]-fpht[fpindex])/(fptime[fpindex+1]-fptime[fpindex]);
      oefile.open("Terraces/outletelev");
      
   }
   
      
   int numg = inputFile.ReadItem( numg, "NUMGRNSIZE" );
   //if( numg<2 ) ReportFatalError("Must use at least 2 sizes with GA." );
   tArray<double> deparr( numg );
   int i;
   for( i=0; i<numg; i++ ) deparr[i] = 0.0;
   assert( strmNet.getInletNodePtr() != 0 );

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
   **      Eolian (loess) deposition (if applicable)
   **      Uplift (or baselevel change)
   **********************************************************************/
   while( !time.IsFinished() )
   {
      time.ReportTimeStatus();

      // Do storm...
      storm.GenerateStorm( time.getCurrentTime(),
                           strmNet.getInfilt(), strmNet.getSoilStore() );
      //cout << storm.getRainrate() << " " << storm.getStormDuration() << " "
      //   << storm.interstormDur() << endl;
      //cin >> dbg;
      strmNet.UpdateNet( time.getCurrentTime(), storm );
      
      // Addition for Geoarchaeology model: set erodibility of main
      // stream to zero and set its profile elevations as a boundary
      // condition
      tMeshListIter<tLNode> nodeIter( mesh.getNodeList() );
      tLNode *cn;
      double elev, totlen;
      // Start by resetting erodibility for all nodes; will be overridden
      // to zero for main channel nodes
      for( cn=nodeIter.FirstP(); nodeIter.IsActive(); cn=nodeIter.NextP() )
          cn->setLayerErody( 0, kr );
      // Set the new drop elevation
      if( optwave==0 ) {
          inletElev = drop + amplitude * sin( period*time.getCurrentTime() );
          //cout << "Inlet " << inletElev << " at " << time.getCurrentTime() << endl;
      }
      else if( optwave==1 )
      {
         noise1=(0.99*noise0+drand48()-0.5)*0.9;
         ttime = fmod( time.getCurrentTime(), period );
         if( ttime<=tpeak )
             inletElev = drop + mr*ttime + noise1;
         else
             inletElev = drop + amplitude - mf*(ttime-tpeak) + noise1;
         noise0=noise1;
         oefile<<noise1<<endl;
         
      }
      else if( optwave==2 )
      {
         if(time.getCurrentTime()>=fptime[fpindex] && time.getCurrentTime()<fptime[fpindex+1]){
            //slope and index don't need to be changed, just calculate elev
            inletElev = fpht[fpindex] + fpslope*(time.getCurrentTime()-fptime[fpindex]);
         }
         else{
            //calculate new slope and update index
            fpindex++;
            fpslope=(fpht[fpindex+1]-fpht[fpindex])/(fptime[fpindex+1]-fptime[fpindex]);
            inletElev = fpht[fpindex] + fpslope*(time.getCurrentTime()-fptime[fpindex]);
         }
         cout<<"fhpt = "<<fpht[fpindex]<<" fpslope "<<fpslope<<" current time "<<time.getCurrentTime()<<" fptime "<<fptime[fpindex]<<endl;
         oefile<<inletElev<<endl;
      }
      

         
      // Find the total length of the main channel and compute slope
      if( optwave<2){
         cn = strmNet.getInletNodePtr();
         totlen = 0.0;
         do
         {
            cn->setLayerErody( 0, 0.0 );  // Main channel elev is a B.C., thus unerodible
            totlen += cn->getFlowEdg()->getLength();
            cn = cn->getDownstrmNbr();
         }
         while( cn->getBoundaryFlag()==kNonBoundary );
         
         chanslp = drop/totlen;
      }
      
      
      // Now set elevations along main channel
      elev = inletElev;  // starting at elevation at inlet
      cn = strmNet.getInletNodePtr(); // begin at inlet
      cn->setZ( inletElev );  // set inlet's elev
      do // work downstream along main channel, setting elevations
      {
         double delz;
         elev = elev - chanslp * cn->getFlowEdg()->getLength();
         cn = cn->getDownstrmNbr();
         delz = elev - cn->getZ();
         while( delz < -0.1 ) {  // Test: erode one active layer thick at a tm
            deparr[0] = -0.1;
            cn->EroDep( 0, deparr, time.getCurrentTime() );
            delz += 0.1;
         }
         deparr[0] = delz;
         cn->EroDep( 0, deparr, time.getCurrentTime() );
      }
      while( cn->getBoundaryFlag()==kNonBoundary );

      //cout << "eroding...\n";
      if( optDetachLim )
          erosion.ErodeDetachLim( storm.getStormDuration() );
      else
          erosion.DetachErode( storm.getStormDuration(), &strmNet,
                               time.getCurrentTime() );
      //cout << "meandering...\n";
      if( optMeander )
          meander->Migrate( time.getCurrentTime() );

      //cout << "overbanking...\n";
      if( optFloodplainDep )
          floodplain->DepositOverbank( storm.getRainrate(),
                                       storm.getStormDuration(),
                                       time.getCurrentTime() );

      // Do interstorm...
      //cout << "Doing diffusion\n";
      erosion.Diffuse( storm.getStormDuration() + storm.interstormDur(),
      optDiffuseDepo );

      //cout << "exposure time...\n";
      erosion.UpdateExposureTime( storm.getStormDuration() + 
                                      storm.interstormDur() );

      if( optLoessDep )
          loess->DepositLoess( &mesh, 
                               storm.getStormDuration()+storm.interstormDur(),
                               time.getCurrentTime() );
      //cout << "Uplift\n";
      if( time.getCurrentTime() < uplift.getDuration() )
          uplift.DoUplift( &mesh,
                           storm.getStormDuration() + storm.interstormDur() );
      time.Advance( storm.getStormDuration() + storm.interstormDur() );
      //cout << "Output\n";
      if( time.CheckOutputTime() )
          output.WriteOutput( time.getCurrentTime() );
      
   }
   
}




