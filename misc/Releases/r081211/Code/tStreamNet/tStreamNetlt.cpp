/**************************************************************************/
/**
**  @file tStreamNet.cpp
**  @brief Functions for class tStreamNet and related class tInlet.
**
**    Modifications:
**     - 2/00 added DensifyMeshDrArea plus modifications to tStreamNet
**       constructor to implement dynamic node addition in regions of
**       high drainage area (ie, main channels; see below) GT
**     - 6/01 added functions to implement Parker-Paola self-formed
**       channel model GT
**     - 2/02 changes to tParkerChannels, tInlet GT
**
**  $Id: tStreamNet.cpp,v 1.83 2005/02/20 20:24:07 childcvs Exp $
*/
/**************************************************************************/

#include <assert.h>
//#include <string>
#include "../errors/errors.h"
#include "tStreamNet.h"

tStreamNet::kChannelType_t tStreamNet::IntToChannelType( int c ){
  switch(c){
  case 1: return kRegimeChannels;
  case 2: return kParkerChannels;
  default:
    std::cout << "You asked for channel geometry model number " << c
	 << " but there is no such thing.\n"
      "Available models are:\n"
      " 1. Regime theory (empirical power-law scaling)\n"
      " 2. Parker-Paola self-formed channel theory\n";
    ReportFatalError( "Unrecognized channel geometry model code.\n" );
  }
}

/*****************************************************************************\
**
**      DistanceToLine: given x,y coords, finds distance to the line
**              defined by given a, b, and c (ax + by + c = 0)
**      Global function.
**
**      Data members updated:
**      Called by:
**      Calls:
**
\*****************************************************************************/
double DistanceToLine( double x2, double y2, double a, double b, double c )
{
   double f, g, h, x, y, d;

   f = -b;
   g = a;
   h = b * x2 - a * y2;
   if( fabs(b) > 0 && fabs(a) > 0 )
   {
      x = (g * c / b - h) / (f - g * a / b);
      y = (-c - a * x) / b;
   }
   else
   {
      if( fabs(a) == 0.0 )
      {
         y = -c / b;
         x = -h / f;
      }
      else
      {
         y = -h / g;
         x = -c / a;
      }
   }
   d = sqrt( (x2 - x) * (x2 - x) + (y2 - y) * (y2 - y) );
   return d;
}


/*****************************************************************************\
**
**      DistanceToLine: given x,y coords, finds distance to the line
**              formed by points  p0->(x, y) and p1->(x, y)
**      Global function.
**
**      Data members updated:
**      Called by:
**      Calls:
**
\*****************************************************************************/
double DistanceToLine( double x2, double y2, tNode const *p0, tNode const *p1 )
{
   double a, b, c, x0, y0, x1, y1, d;

   x0 = p0->getX();
   y0 = p0->getY();
   x1 = p1->getX();
   y1 = p1->getY();
   a = y1 - y0;
   b = x0 - x1;
   c = -( a * x0 + b * y0 );

   d = DistanceToLine( x2, y2, a, b, c );
   return d;
}


/**************************************************************************\
**  Functions for class tStreamNet here:
\**************************************************************************/

/**************************************************************************\
**
**  Constructor
**
**     The constructor asks for references to a mesh (i.e., _the_
**     mesh), a storm, and an input file). It sets the _meshPtr_ and
**     _stormPtr_, initializes the rainfall rate, and reads in various
**     hydrologic options and parameters from the input file.
**     It then initializes flow directions, drainage areas, etc.
**
**     Inputs:  meshRef -- reference to tMesh object. This ref is copied
**                         and used to communicate w/ the mesh.
**              storm   -- ref to storm object, copied and used to obtain
**                         rainfall input.
**              infile  -- ref to input file used to read parameters
**
**     Modifications:
**       - GT added input of parameters for sinusoidal variation in
**         infiltration capacity, 7/20/98.
**       - GT commented out mndrchngprob, which appears to be unused, 6/99
\**************************************************************************/

tStreamNet::tStreamNet( tMesh< tLNode > &meshRef, tStorm &storm,
                        const tInputFile &infile )
  :
  meshPtr(&meshRef),
  stormPtr(&storm),
  trans(0), infilt(0),
  inlet( &meshRef, infile ),
  optSinVarInfilt(false),
  mpParkerChannels(0)
{
   if (0) //DEBUG
     std::cout << "tStreamNet(...)...";
   assert( meshPtr != 0 );
   assert( stormPtr != 0 );

   // Read option for runoff generation/routing and get relevant parameters
   knds = infile.ReadItem( knds, "HYDR_ROUGH_COEFF_DS" );

   {
     int tmp_;
     tmp_ = infile.ReadItem( tmp_, "FLOWGEN" );
     miOptFlowgen = static_cast<kFlowGen_t>(tmp_);
   }
   filllakes = infile.ReadBool( "LAKEFILL" );
   if( miOptFlowgen == kSaturatedFlow1 || miOptFlowgen==kSaturatedFlow2 )
      trans = infile.ReadItem( trans, "TRANSMISSIVITY" );
   if( miOptFlowgen==kSaturatedFlow2 || miOptFlowgen==kConstSoilStore
       || miOptFlowgen==kHortonian )
   {
      infilt = infile.ReadItem( infilt, "INFILTRATION" );
      optSinVarInfilt = infile.ReadBool( "OPTSINVARINFILT" );
      if( optSinVarInfilt )
      {
         twoPiLam = (2.0*PI)/(infile.ReadItem( twoPiLam, "PERIOD_INFILT" ));
         infilt_dev = infile.ReadItem( infilt_dev, "MAXICMEAN" ) - infilt;
         infilt0 = infilt;
      }
   }
   else infilt = 0.0;
   if( miOptFlowgen == kConstSoilStore )
       soilStore = infile.ReadItem( soilStore, "SOILSTORE" );
   else soilStore = 0.0;
   if( miOptFlowgen == k2DKinematicWave )
   {
       mdKinWaveExp = infile.ReadItem( mdKinWaveExp, "KINWAVE_HQEXP" );
       mdKinWaveRough = knds / SECPERYEAR;
       std::cout << "mdKinWaveRough " << mdKinWaveRough << std::endl;
   }
   else mdKinWaveExp = mdKinWaveRough = 0.0;
   if( miOptFlowgen == kHydrographPeakMethod )
     {
       mdFlowVelocity = infile.ReadItem( mdFlowVelocity, "FLOWVELOCITY" );
       mdHydrgrphShapeFac = infile.ReadItem( mdHydrgrphShapeFac, "HYDROSHAPEFAC" );
     }


   // Get the initial rainfall rate from the storm object, and read in option
   // for stochastic variation in rainfall
   rainrate = stormPtr->getRainrate();
   optrainvar = infile.ReadBool( "OPTVAR" );
   // Arnaud: useless
   // Options related to stream meandering
   // int itMeanders = infile.ReadItem( itMeanders, "OPTMEANDER" );

   // Read hydraulic geometry parameters: first, those used (potentially)
   // by both "regime" and "parker" hydraulic geometry models
   bankfullevent = infile.ReadItem( bankfullevent, "BANKFULLEVENT" );
   bankfullevent /= SECPERYEAR;
   if( bankfullevent<=0.0 )
       ReportFatalError(
           "Input error: BANKFULLEVENT must be greater than zero" );

   ewstn = infile.ReadItem( ewstn, "HYDR_WID_EXP_STN" );
   //std::cout << "knds: " << knds << std::endl;
   //assert( knds > 0 );
   ends = infile.ReadItem( ends, "HYDR_ROUGH_EXP_DS" );
   //std::cout << "ends: " << ends << std::endl;
   enstn = infile.ReadItem( enstn, "HYDR_ROUGH_EXP_STN" );
   //std::cout << "enstn: " << enstn << std::endl;
   klambda = infile.ReadItem( klambda, "BANK_ROUGH_COEFF" );
   elambda = infile.ReadItem( elambda, "BANK_ROUGH_EXP" );

   // Read remaining hydraulic geometry parameters according to which
   // (of currently 2) model is chosen
   {
     int cread;
     cread = infile.ReadItem( cread, "CHAN_GEOM_MODEL" );
     miChannelType = IntToChannelType(cread);
   }
   switch(miChannelType){
   case kRegimeChannels:
     {
       kwds = infile.ReadItem( kwds, "HYDR_WID_COEFF_DS" );
       //std::cout << "kwds: " << kwds << std::endl;
       assert( kwds > 0 );
       kdds = infile.ReadItem( kdds, "HYDR_DEP_COEFF_DS" );
       ewds = infile.ReadItem( ewds, "HYDR_WID_EXP_DS" );
       //std::cout << "ewds: " << ewds << std::endl;
       edds = infile.ReadItem( edds, "HYDR_DEP_EXP_DS" );
       edstn = infile.ReadItem( edstn, "HYDR_DEP_EXP_STN" );
     }
     break;
   case kParkerChannels:
     mpParkerChannels = new tParkerChannels( infile );
     break;
   }

   // Option for adaptive meshing related to drainage area
   bool optMeshAdapt = infile.ReadBool( "OPTMESHADAPTAREA" );
   if( optMeshAdapt )
   {
      mdMeshAdaptMinArea = infile.ReadItem( mdMeshAdaptMinArea,
                                            "MESHADAPTAREA_MINAREA" );
      mdMeshAdaptMaxVArea = infile.ReadItem( mdMeshAdaptMaxVArea,
                                             "MESHADAPTAREA_MAXVAREA" );
   }

   // Initialize the network by calculating slopes, flow directions,
   // drainage areas, and discharge
   CalcSlopes();  // TODO: should be in tMesh
   InitFlowDirs(); // TODO: should all be done in call to updatenet
   FlowDirs();
   CheckNetConsistency();
   MakeFlow( 0.0 );
   if (0) //DEBUG
     std::cout << "finished" << std::endl;
}

//necessary?
tStreamNet::~tStreamNet()
{
   meshPtr = 0;
   stormPtr = 0;
   delete mpParkerChannels;
   if (0) //DEBUG
     std::cout << "~tStreamNet()" << std::endl;
}


/**************************************************************************\
**
**  tStreamNet get/set functions
**
\**************************************************************************/
void tStreamNet::ResetMesh( tMesh< tLNode > &meshRef )
{
   assert( &meshRef != 0 );
   meshPtr = &meshRef;
   assert( meshPtr != 0 );
}

const tMesh< tLNode > *
tStreamNet::getMeshPtr() const {return meshPtr;}

tMesh< tLNode > *
tStreamNet::getMeshPtrNC() {return meshPtr;}

const tStorm *tStreamNet::getStormPtr() const {return stormPtr;}

tStorm *tStreamNet::getStormPtrNC() {return stormPtr;}

inline tStreamNet::kFlowGen_t tStreamNet::getFlowGenOpt() const {return miOptFlowgen;}

bool tStreamNet::getFillLakesOpt() const {return filllakes;}

double tStreamNet::getRainRate() const {return rainrate;}

double tStreamNet::getTransmissivity() const {return trans;}

double tStreamNet::getInfilt() const {return infilt;}

double tStreamNet::getSoilStore() const {return soilStore;}

double tStreamNet::getInDrArea() const {return inlet.inDrArea;}

double tStreamNet::getInSedLoad() const {return inlet.inSedLoad;}

tArray< double > tStreamNet::getInSedLoadm( ) const {return inlet.inSedLoadm;}

tLNode const *tStreamNet::getInletNodePtr() const {return inlet.innode;}
tLNode *tStreamNet::getInletNodePtrNC() {   return inlet.innode;}

double tStreamNet::getInletSlope() const { return inlet.inletSlope; }

tArray< double > tStreamNet::getInletSedSizeFraction() const {return inlet.inletSedSizeFraction; }

//double tStreamNet::getMndrDirChngProb() const {return mndrDirChngProb;}

// TODO: the value checks are nice, but will hurt performance. Should
// probably be removed.
void tStreamNet::setFlowGenOpt( kFlowGen_t val )
{miOptFlowgen=val;
/*miOptFlowgen = ( val == 0 || val == 1 ) ? val : 0;*/}

void tStreamNet::setFillLakesOpt( bool val )
{filllakes = val;}

void tStreamNet::setRainRate( double val ) {rainrate = val;
/*( val >= 0 ) ? val : 0;*/}

void tStreamNet::setTransmissivity( double val )
{trans = ( val >= 0 ) ? val : 0;}

void tStreamNet::setInfilt( double val )
{infilt = ( val >= 0 ) ? val : 0;}

void tStreamNet::setInDrArea( double val )
{inlet.inDrArea = ( val >= 0 ) ? val : 0;}

void tStreamNet::setInSedLoad( double val )
{inlet.inSedLoad = ( val >= 0 ) ? val : 0;}

void tStreamNet::setInSedLoadm( size_t i, double val )
{
   //shouldn't this be an assert?
   if( i >= inlet.inSedLoadm.getSize() )
       ReportFatalError("Tried to index a size that doesn't exist in tstreamnet");
   inlet.inSedLoadm[i]=val;
}

void tStreamNet::setInletNodePtr( tLNode *Ptr )
{ inlet.innode = Ptr;}

//void tStreamNet::setMndrDirChngProb( double val )
//{mndrDirChngProb = ( val >= 0.0 && val <= 1.0 ) ? val : 1.0;}



/**************************************************************************\
**
**  tStreamNet::UpdateNet
**
**  This function updates the network and flow field by calling various
**  helper functions. There are two versions; the second takes a reference
**  to a storm object as a parameter and uses it to update the current
**  rainfall rate.
**
**  Inputs:  time -- current time in simulation; passed to MakeFlow,
**                   where it is optionally used to implement
**                   sinusoidal time-variations in hydrologic parameters
**                   like infiltration capacity.
**           storm (2) -- reference to storm object, used to update
**                        current rainfall rate.
**
**  Modifications:
**    - 7/20/98: now takes current time as a param, to use for updating
**      sine-varying infilt cap if applicable. Also, "storm" version now
**      simply calls "regular" version after doing its own thing. GT
**
**  TODO: move mesh-related routines -- slopes, voronoi areas, etc --
**         to tMesh
**
\**************************************************************************/
void tStreamNet::UpdateNet( double time )
{
   if (0) //DEBUG
     std::cout << "UpdateNet()...";
   CalcSlopes();          // TODO: should be in tMesh
   FlowDirs();
   
  if(0) //DEBUG
  {
         tMesh< tLNode >::nodeListIter_t mli( meshPtr->getNodeList() );  // gets nodes from the list
		 tLNode * cn;
		 for( cn=mli.FirstP(); mli.IsActive(); cn=mli.NextP() )
		 {
		    if( cn->getID()==8121 ) 
			{
		       tEdge * debugedg = cn->getFlowEdg();
		       tLNode * nbr = static_cast<tLNode *>(debugedg->getDestinationPtrNC());
		       std::cout<<"UpdateNet 1: node "<<cn->getID()<<" edge "<<debugedg->getID()<<" downstream nbr "<<nbr->getID()<<std::endl;
		       std::cout<<"z "<<cn->getZ()<<" dsn z "<<nbr->getZ() << " mdr "<<cn->Meanders()<<" nbr mdr "<<nbr->Meanders()<<std::endl;
            }
		}
   }

   MakeFlow( time );
   
  if(0) //DEBUG
  {
         tMesh< tLNode >::nodeListIter_t mli( meshPtr->getNodeList() );  // gets nodes from the list
		 tLNode * cn;
		 for( cn=mli.FirstP(); mli.IsActive(); cn=mli.NextP() )
		 {
		    if( cn->getID()==8121 ) 
			{
		       tEdge * debugedg = cn->getFlowEdg();
		       tLNode * nbr = static_cast<tLNode *>(debugedg->getDestinationPtrNC());
		       std::cout<<"UpdateNet 2: node "<<cn->getID()<<" edge "<<debugedg->getID()<<" downstream nbr "<<nbr->getID()<<std::endl;
		       std::cout<<"z "<<cn->getZ()<<" dsn z "<<nbr->getZ() << " mdr "<<cn->Meanders()<<" nbr mdr "<<nbr->Meanders()<<std::endl;
            }
		}

   }

   CheckNetConsistency();
   
  if(0) //DEBUG
  {
         tMesh< tLNode >::nodeListIter_t mli( meshPtr->getNodeList() );  // gets nodes from the list
		 tLNode * cn;
		 for( cn=mli.FirstP(); mli.IsActive(); cn=mli.NextP() )
		 {
		    if( cn->getID()==8121 ) 
			{
		       tEdge * debugedg = cn->getFlowEdg();
		       tLNode * nbr = static_cast<tLNode *>(debugedg->getDestinationPtrNC());
		       std::cout<<"UpdateNet 3: node "<<cn->getID()<<" edge "<<debugedg->getID()<<" downstream nbr "<<nbr->getID()<<std::endl;
		       std::cout<<"z "<<cn->getZ()<<" dsn z "<<nbr->getZ() << " mdr "<<cn->Meanders()<<" nbr mdr "<<nbr->Meanders()<<std::endl;
            }
		}

   }

   if (0) //DEBUG
     std::cout << "UpdateNet() finished" << std::endl;
}

void tStreamNet::UpdateNet( double time, tStorm &storm )
{
   if (0) //DEBUG
     std::cout << "UpdateNet(...)...";
   stormPtr = &storm;
   assert( stormPtr != 0 );
   rainrate = stormPtr->getRainrate();
   UpdateNet( time );

   //debug test
   //SortNodesByNetOrder( 1 );
   //SortNodesByNetOrder();
}


/**************************************************************************\
**
**  tStreamNet::CheckNetConsistency
**
**  This is a debugging function that tests the integrity of the stream
**  network. Tests include making sure all active nodes have a valid
**  downstream neighbor; all active nodes drain to an outlet (or a sink);
**  no node drains to a closed boundary; etc. No data are modified.
**
**  Modifications:
**    - disabled O(nnodes^2) loop. 8/03 AD
**    - some code moved to CheckNetConsistencyFlowPath 8/03 AD
**
\**************************************************************************/
void tStreamNet::CheckNetConsistency()
{
   tLNode *cn, *dn;
   tMesh< tLNode >::nodeListIter_t nI( meshPtr->getNodeList() ),
       tI( meshPtr->getNodeList() );
   enum{
     expensive_test = false
   };

   for( cn = nI.FirstP(); nI.IsActive(); cn = nI.NextP() )
   {
      //make sure each node has a viable downstrm nbr:
      dn = cn->getDownstrmNbr();
      if( dn == 0 )
      {
         std::cerr << "NODE #" << cn->getID() << " has no downstrm nbr\n";
         goto error;
      }
      if( dn->getID() == cn->getID() ){
         std::cerr<< "NODE #" << cn->getID() << " flows to itself!\n";
         goto error;
      }
      if (expensive_test) { // This is a O(nnodes^2) loop.
	tLNode *ln;
	if( ( ln = tI.GetP( dn->getID() ) ) != NULL )
	  {
	    if( ln != dn )
	      {
		std::cerr << "NODE #" << cn->getID()
		     << " downstrm nbr not id. to node in nodeList with same ID\n";
		goto error;
	      }
	  }
	else
	  {
	    std::cerr << "NODE #" << cn->getID()
		 << " downstrm nbr is not in nodeList\n";
	    goto error;
	  }
      }
      if( cn->Meanders() )
	{
	  if( cn->calcSlope() < 0.0 )
	    {
	      std::cerr << "error in tStreamNet::CheckNetConsistency at NODE #" << cn->getID()
	      << ",meanders and returns negative getSlope of "<<cn->calcSlope()<<'\n';
	      goto error;
	    }
	}
   }

   if (CheckNetConsistencyFlowPath(&cn))
     goto error;

   if (0) //DEBUG
     std::cout << "NETWORK PASSED\n";

   return;

  error:
   std::cerr << "Problem detected at the following node:\n";
   cn->TellAll();
   ReportFatalError( "Error in network consistency." );
}


/**************************************************************************\
**
**  tStreamNet::CheckNetConsistencyFlowPath
**
**  This is a debugging function that tests the integrity of the stream
**  network. all active nodes drain to an outlet (or a sink);
**  no node drains to a closed boundary; etc. No data are modified
**  (Except the tLNode->public1)
**
**  Modifications:
**    - added "sink" case to outlet-path test, 8/98 GT
**    - fixed bug wherein use of a fixed "kLargeNumber" to test for
**      infinite loop in flow routing resulted in detection of false
**      loops for very large meshes! 3/02 GT
**    - speed-up using tLNode::public1 8/03 AD
**
\**************************************************************************/
int tStreamNet::CheckNetConsistencyFlowPath( tLNode **pcn )
{
   tLNode *cn;
   tMesh< tLNode >::nodeListIter_t nI( meshPtr->getNodeList() );
   const long nodesInMesh = meshPtr->getNodeList()->getSize();
   enum{ notdone, valid };

   // mark all active nodes as not done.
   for( cn = nI.FirstP(); nI.IsActive(); cn = nI.NextP() )
     cn->public1 = notdone;
   for( cn = nI.FirstP(); nI.IsActive(); cn = nI.NextP() )
   {
      if (cn->public1 == valid) continue;

      // Make sure each node has path to outlet (or to a sink):
      bool metExistingFlowPath = false;
      int ctr = 1;
      tLNode *dn = cn->getDownstrmNbr();
      while( dn->getBoundaryFlag() == kNonBoundary
             && dn->getFloodStatus() != tLNode::kSink )
      {
         dn = dn->getDownstrmNbr();
         ctr++;
	 if (dn->public1 == valid ) {
	   metExistingFlowPath = true;
	   break;
	 }
         if( ctr > nodesInMesh ) {
             dn->TellAll();
	     if( ctr > nodesInMesh+4 )
	       {
		 std::cerr << "NODE #" << cn->getID()
		      << " has infinite loop in path downstream\n";
		 goto error;
	       }
	 }
      }
      if( !metExistingFlowPath &&
	  dn->getBoundaryFlag() != kOpenBoundary
          && dn->getFloodStatus() != tLNode::kSink )
      {
         std::cerr << "NODE #" << cn->getID()
              << " does not flow to outlet\n"
	      << "Flow ends up at the following node:\n";
         dn->TellAll();
         goto error;
      }
      // no error, mark the path as valid up to the first valid node
      cn->public1 = valid;
      dn = cn->getDownstrmNbr();
      if (metExistingFlowPath) {
	while( dn->public1 != valid )
	  {
	    dn->public1 = valid;
	    dn = dn->getDownstrmNbr();
	  }
      } else {
	while( dn->getBoundaryFlag() == kNonBoundary
	       && dn->getFloodStatus() != tLNode::kSink )
	  {
	    dn->public1 = valid;
	    dn = dn->getDownstrmNbr();
	  }
      }
   }
   *pcn = 0;
   return 0;
 error:
   *pcn = cn;
   return 1;
}


void tStreamNet::CheckMeander() const
{

   tLNode const *cn, *secondnode, *thirdnode;
   tEdge const * fe;
   double totlen;


   std::cout<<"                              \n";
   std::cout<<" --tStreamNet::CheckMeander-- \n";

   totlen=0.;
   cn= getInletNodePtr();
   int counter=0;
   while (counter < 150 && cn != NULL)		// check the first to hundred meander cells
    {
      fe = cn->getFlowEdg();
      if(fe==0) break;
      assert( fe );
      totlen += fe->getLength();


      secondnode = cn->getDownstrmNbr();
      thirdnode  = secondnode->getDownstrmNbr();

      //----------------------------AVULSION---------------------------------------------------
      // Check for avulsion
      // If an avulsion happened in the the experiments with the Thames DEM as template
      // the y-coordinate of a channel cells should be > 20250 or < 19500
      if( (cn->getY() > 22000 || cn->getY() < 18000) && cn->getX() > 7100){
      	std::cout<<" --AVULSION ERROR--"<<std::endl;
      	std::cout<<"Printing the whole meander list in function CheckMeander:"<<std::endl;
      	std::cout<<"Locations "<<cn->getID()<<" "<<cn->getX()<<" "<<cn->getY()<<" are too far away from the channel/ basin axis "<<std::endl;
      	std::cout<<"Find the avulsion site...."<<std::endl;

      	std::cout<<"Showing the existing meander list:  "<<std::endl;
        cn= getInletNodePtr();
        counter=0;
        while (counter < 150 && cn->getX() < 11900){
   	      if(!cn->Meanders()) std::cout<<" *|*| Avulsion here?"<<std::endl;
   	      std::cout <<cn->getID()<<" "
		   <<cn->getX()<<" "<<cn->getY()<<" "<<cn->getZ()
		   <<" M: "<<cn->Meanders()
		   <<" F: "<< tLNode::FloodName(cn->getFloodStatus())<<std::endl;
   	      cn=cn->getDownstrmNbr();
   	      counter++;
         }


        //exit(1);

      } // end avulsion..

      //-------------------------RECURSION-----------------------------------------------------

      //Check For Recursion
      if(thirdnode==cn && cn != NULL){
      	std::cout<<"   "<<std::endl;
        std::cout<<"ERROR: in Checkmeander, recursion in channel"<<std::endl;
        std::cout<<"Node "<<cn->getID()<<' '<<cn->getX()<<' '<<cn->getY()<<' '<<cn->getZ()<<std::endl;
        std::cout<<"points to "<<std::endl;
        std::cout<<"Node "<<secondnode->getID()<<' '<<secondnode->getX()<<' '<<secondnode->getY()<<' '<<secondnode->getZ()<<std::endl;
        std::cout<<"and back, ID of dwnstrnbt is: "<<std::endl;
        std::cout<<thirdnode->getID()<<std::endl;
        std::cout<<"   "<<std::endl;
        std::cout<<"The nodes surrounding "<< cn->getID()<<" are: "<<std::endl;
        tEdge *ce;
        tSpkIter spokIter( const_cast<tLNode*>(cn) );
        for( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() )
        {
          tLNode* surnode = static_cast<tLNode *>( ce->getDestinationPtrNC() );
          std::cout <<surnode->getID()<<" "
	       <<surnode->getX()<<" "<<surnode->getY()<<" "<<surnode->getZ()
	       <<" M: "<<surnode->Meanders()
	       <<" F: "<< tLNode::FloodName(surnode->getFloodStatus()) <<std::endl;
        }

        std::cout<<"   "<<std::endl;
        std::cout<<"The nodes surrounding its dwnstr nbr "<< secondnode->getID()<<" are: "<<std::endl;
        tSpkIter spokIter2( const_cast<tLNode*>(secondnode) );
        for( ce = spokIter2.FirstP(); !( spokIter2.AtEnd() ); ce = spokIter2.NextP() )
        {
          tLNode* surnode = static_cast<tLNode *>( ce->getDestinationPtrNC() );
          std::cout<<surnode->getID()<<" "
	      <<surnode->getX()<<" "<<surnode->getY()<<" "<<surnode->getZ()
	      <<" M: "<< surnode->Meanders()
	      <<" F: "<< tLNode::FloodName(surnode->getFloodStatus())<<std::endl;
        }

        std::cout<<"     "<<std::endl;
        std::cout<<"Showing the existing meander list:  "<<std::endl;
        cn= getInletNodePtr();
        counter=0;
        while (counter < 150){

   	      std::cout<<cn->getID()<<" "
		  <<cn->getX()<<" "<<cn->getY()<<" "<<cn->getZ()
		  <<" M: "<<cn->Meanders()
		  <<" F: "<<tLNode::FloodName(cn->getFloodStatus())<<std::endl;
   	      cn=cn->getDownstrmNbr();
   	      counter++;
         }


        //exit(1);

       }




      //---------------------------CHANNEL BUMPS-------------------------------------------------

      // Remove bumps
      if(secondnode->getZ() > cn->getZ()){

#if 0
	 double newelev;
	 if(thirdnode !=NULL) {
           newelev = ( thirdnode->getZ() + cn->getZ() )/2.0;
         } else {
	   newelev = cn->getZ();
         }
#endif
         std::cout<<"    "<<std::endl;
         std::cout<<"ERROR: in CheckMeander, bump in channel"<<std::endl;
         std::cout<<"Node "<<cn->getID()<<' '<<cn->getX()<<' '<<cn->getY()<<' '<<cn->getZ()<<std::endl;
         std::cout<<"is lower than its downstreamNeighbour"<<std::endl;
         std::cout<<"Node "<<secondnode->getID()<<' '<<secondnode->getX()<<' '<<secondnode->getY()<<' '<<secondnode->getZ()<<std::endl;


         tEdge *ce;
         tSpkIter spokIter( const_cast<tLNode*>(cn) );
         for( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() )
         {
          tLNode* surnode = static_cast<tLNode *>( ce->getDestinationPtrNC() );
          if(surnode !=NULL){
            std::cout<<surnode->getID()<<" "
		<<surnode->getX()<<" "<<surnode->getY()<<" "<<surnode->getZ()
		<<" M: "<<surnode->Meanders()
		<<" F: "<<tLNode::FloodName(surnode->getFloodStatus())<<std::endl;
          }
         }

        std::cout<<"    "<<std::endl;
        std::cout<<"The nodes surrounding its dwnstr nbr "<< secondnode->getID()<<" are: "<<std::endl;
        tSpkIter spokIter3( const_cast<tLNode*>(secondnode) );
        for( ce = spokIter3.FirstP(); !( spokIter3.AtEnd() ); ce = spokIter3.NextP() )
        {
          tLNode* surnode = static_cast<tLNode *>( ce->getDestinationPtrNC() );
          std::cout<<surnode->getID()<<" "
	      <<surnode->getX()<<" "<<surnode->getY()<<" "<<surnode->getZ()
	      <<" M: "<<surnode->Meanders()
	      <<" F: "<<tLNode::FloodName(surnode->getFloodStatus())<<std::endl;
        }

	std::cout<<"     "<<std::endl;
	std::cout<<"Showing the existing meander list:  "<<std::endl;
	cn= getInletNodePtr();
	counter=0;
	while (counter < 150 && cn != NULL){
	  std::cout<<cn->getID()<<" "
	      <<cn->getX()<<" "<<cn->getY()<<" "<<cn->getZ()
	      <<" M: "<<cn->Meanders()
	      <<" F: "<<tLNode::FloodName(cn->getFloodStatus())<<std::endl;
	  cn=cn->getDownstrmNbr();

	  counter++;
	}

	std::cout<<"Stop Run due to bump "<<std::endl;
	//exit(1);
      }

      if(cn != NULL){
        if(cn->getFloodStatus() != tLNode::kNotFlooded){
      	  std::cout<<"Warning: node "<<cn->getID()<<' '
	      <<cn->getX()<<' '<<cn->getY()<<' '<<cn->getZ()<<std::endl;
      	  std::cout<<"has flood status: "<<tLNode::FloodName(cn->getFloodStatus())<<std::endl;
      	}
      }

      counter++;

      if(cn == NULL) break;

      cn = cn->getDownstrmNbr();

      if(cn == NULL) break;

      if(cn->getBoundaryFlag() != kNonBoundary) break;

      }

    //---------------end new correction stuff------------------


   //std::cout<< "M-Nodes " <<cn->getID()<<' '<<cn->getX()<< ' ' << cn->getY()<<' '<<cn->getZ()<<" A= "<<cn->getDrArea()<<" Q= "<<cn->getQ()<<" W= "<<cn->getChanWidth()<< " Mndr="<<cn->Meanders()<<" Fld= "<<cn->getFloodStatus()<<std::endl;

  std::cout<<" --tStreamNet::CheckMeander, Finished -- \n";
  std::cout<<"   \n";

} // CheckMeander

//---------------------------------------------------------
void tStreamNet::ShowMeanderNeighbours(int debugID) const
{

         tLNode const *cn = getInletNodePtr();
         int counter = 0;
         while (counter < 100){

   	      //std::cout<<cn->getID()<<" "<<cn->getX()<<" "<<cn->getY()<<" "<<cn->getZ()<<" M: "<<cn->Meanders()<<" F: "<<cn->getFloodStatus()<<std::endl;
   	      if( cn->getID() ==debugID){

   	      	tEdge *ce;
                tSpkIter spokIter( const_cast<tLNode*>(cn) );

                for( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() )
                {
                 tLNode* surnode = static_cast<tLNode *>( ce->getDestinationPtrNC() );
                 std::cout<<surnode->getID()<<" "
		     <<surnode->getX()<<" "<<surnode->getY()<<" "<<surnode->getZ()
		     <<" M: "<<surnode->Meanders()
		     <<" F: "<<tLNode::FloodName(surnode->getFloodStatus())<<std::endl;
                }
   	      }

   	      cn=cn->getDownstrmNbr();
   	      counter++;
         }



}

/****************************************************************************\
**
**  CalcSlopes
**
**  Computes the slope of each edge in the mesh. This is done by computing
**  the slope for the first of each pair of complementary edges, then
**  assigning the negative of that value to the second of the pair.
**
**  Memory allocated: none
**  Notes:
**   - if edge length is zero, it is calculated
**   - complementary edges on the list are assumed to be organized pairwise;
**     that is, edges AB and BA are always together, for example.
**
**  TODO: should be a member of tMesh!
**
\****************************************************************************/
void tStreamNet::CalcSlopes()
{
  assert( meshPtr != 0 );
  tEdge *curedg;
  tMesh< tLNode >::edgeListIter_t i( meshPtr->getEdgeList() );

  if (0) //DEBUG
    std::cout << "CalcSlopes()...";

  // Loop through each pair of edges on the list
  for( curedg = i.FirstP(); !( i.AtEnd() ); curedg = i.NextP() )
    {
      // Make sure edge is valid, and length is nonzero
      assert( curedg != 0 );
      assert( curedg->getLength() > 0 );

      // Compute the slope and assign it to the current edge
      const double slp =
	( curedg->getOrgZ() - curedg->getDestZ() )
	/ curedg->getLength();
      curedg->setSlope( slp );

      // Advance to the edge's complement, and assign it -slp
      curedg = i.NextP();
      assert( !( i.AtEnd() ) );
      curedg->setSlope( -slp );
      //curedg->setLength( length );
      assert( curedg->getLength() > 0 );
    }
  if (0) //DEBUG
    std::cout << "CalcSlopes() finished" << std::endl;
}


/****************************************************************************\
**
**  tStreamNet::InitFlowDirs
**
**  Initialize flow directions such that each active (non-boundary) node
**  flows to another active node (or open boundary node). This initialization
**  process allows the FlowDirs function to assume that the previous flowedg
**  is always valid in the sense that it doesn't point to a closed boundary
**  (which otherwise could happen when the mesh is first read in).
**
**  Modifies: node flow directions (flowedg)
**  Written 12/1/97 gt.
**
\****************************************************************************/
#define kMaxSpokes 100
void tStreamNet::InitFlowDirs()
{
   tMesh< tLNode >::nodeListIter_t i( meshPtr->getNodeList() );
   tLNode * curnode;
   tEdge * flowedg;
   int ctr;

   if (0) //DEBUG
     std::cout << "InitFlowDirs()...\n";

   // For every active (non-boundary) node, initialize it to flow to a
   // non-boundary node (ie, along a "flowAllowed" edge)
   curnode = i.FirstP();
   while( i.IsActive() )
   {
      // Start with the node's default edge
      assert( curnode!=0 );
      flowedg = curnode->getEdg();
      assert( flowedg!=0 );

      // As long as the current edge is a no-flow edge, advance to the next one
      // counter-clockwise
      ctr = 0;
      while( !flowedg->FlowAllowed() )
      {
         flowedg = flowedg->getCCWEdg();
         assert( flowedg!=0 );
         ctr++;
         if( ctr>kMaxSpokes ) // Make sure to prevent std::endless loops
         {
            std::cerr << "Mesh error: node " << curnode->getID()
                 << " appears to be surrounded by closed boundary nodes"
                 << std::endl;
            ReportFatalError( "Bailing out of InitFlowDirs()" );
         }
      }
      curnode->setFlowEdg( flowedg );
      assert( curnode->getFlowEdg() != 0 );
      curnode = i.NextP();
   }

   if (0) //DEBUG
     std::cout << "finished\n";

}
#undef kMaxSpokes


/****************************************************************************\
**
**  tStreamNet::ReInitFlowDirs
**
**  Initialize flow directions such that each active (non-boundary) node
**  flows to another active node (or open boundary node). This initialization
**  process allows the FlowDirs function to assume that the previous flowedg
**  is always valid in the sense that it doesn't point to a closed boundary
**  (which otherwise could happen when the mesh is first read in).
**
**  Modifies: node flow directions (flowedg)
**  Written 12/1/97 gt.
**
\****************************************************************************/
#define kMaxSpokes 100
void tStreamNet::ReInitFlowDirs()
{
   if (0) //DEBUG
     std::cout << "ReInitFlowDirs()...\n";
   // For every active (non-boundary) node, initialize it to flow to a
   // non-boundary node (ie, along a "flowAllowed" edge)
   tMesh< tLNode >::nodeListIter_t i( meshPtr->getNodeList() );
   tLNode * curnode = i.FirstP();
   while( i.IsActive() )
   {
      // reinitialize flowedge only if flowedge currently NULL or is no-flow edge:
      if( curnode->getFlowEdg() == NULL || !curnode->getFlowEdg()->FlowAllowed() )
      {
         // Start with the node's default edge
         assert( curnode!=0 );
         tEdge * flowedg = curnode->getEdg();
         assert( flowedg!=0 );

         // As long as the current edge is a no-flow edge, advance to the next one
         // counter-clockwise
         int ctr = 0;
         while( !flowedg->FlowAllowed() )
         {
            flowedg = flowedg->getCCWEdg();
            assert( flowedg!=0 );
            ctr++;
            if( ctr>kMaxSpokes ) // Make sure to prevent std::endless loops
            {
               std::cerr << "Mesh error: node " << curnode->getID()
                    << " appears to be surrounded by closed boundary nodes"
                    << std::endl;
               ReportFatalError( "Bailing out of InitFlowDirs()" );
            }
         }
         curnode->setFlowEdg( flowedg );
         assert( curnode->getFlowEdg() != 0 );
      }
      curnode = i.NextP();
   }

   if (0) //DEBUG
     std::cout << "finished\n";

}
#undef kMaxSpokes

/****************************************************************************\
**
**  TryDamBypass: See whether we can avoid flowing uphill by deleting nearby
**     points.
**
**      Parameters:     none
**      Called by:
**      Assumes:
**         node points to a valid flow edge which connects the lowest nbr
**      Created: SL
**
\****************************************************************************/
/*Xint tStreamNet::DamBypass( tLNode *snknod )
{
   std::cout << "DamBypass" << std::endl;
   int nv, nvopp, nvother, maxcnt = 10, cntr = 0;
   double cz, nz, nnz, dis0, dis1, slp;
   tTriangle *mytri, *opptri;
   int uphill = 1;
   tEdge *fedg, *ce;
   tLNode *nbr;
   tLNode *nxtnbr;
   tLNode *othernbr;
   tLNode *pointtodelete;
   tArray< double > cxy, nxy, nnxy;
   tSpkIter sI( snknod );
   cxy = snknod->get2DCoords();

   cz = snknod->getZ();
   //while( uphill && cntr <= maxcnt )
   //{
      fedg = snknod->getFlowEdg();
      nbr = snknod->getDownstrmNbr();
      nxtnbr = nbr->getDownstrmNbr();
      if( nxtnbr == 0 ) return 0;
      if ( snknod == nxtnbr ) return 0;
      nnz = nxtnbr->getZ();
      if( nnz > cz ) return 0; //proceed no further if next nbr is not downhill
      nxy = nbr->get2DCoords();
      nnxy = nxtnbr->get2DCoords();
      if( PointsCCW( cxy, nxy, nnxy ) )
      {
         mytri = meshPtr->TriWithEdgePtr( fedg );
         nv = mytri->nVtx( snknod );
         othernbr = (tLNode *) mytri->pPtr( (nv+2)%3 );
      }
      else
      {
         mytri = meshPtr->TriWithEdgePtr( meshPtr->getEdgeComplement( fedg ) );
         nv = mytri->nVtx( snknod );
         othernbr = (tLNode *) mytri->pPtr( (nv+1)%3 );
      }
      opptri = mytri->tPtr(nv);
      dis0 = nbr->Dist( snknod, nxtnbr );
      dis1 = othernbr->Dist( snknod, nxtnbr );
      if( dis0 < dis1 ) pointtodelete = nbr;
      else pointtodelete = othernbr;
      if( pointtodelete->getBoundaryFlag() == kNonBoundary )
          meshPtr->DeleteNode( pointtodelete );
      else return 0;
      cntr++;
      for( ce = sI.FirstP(); !(sI.AtEnd()); ce = sI.NextP() )
      {
         nz = ce->getDestinationPtrNC()->getZ();
         slp = ( cz - nz ) / ce->getLength();
         ce->setSlope( slp );
         meshPtr->getEdgeComplement( ce )->setSlope( -slp );
         //if( ce->getDestinationPtrNC()->getZ() < cz ) uphill = 0;
      }
      //}
      std::cout << "DamBypass: deleted " << cntr << " nodes" << std::endl;
   return 1;
}*/


/****************************************************************************\
**
**  tStreamNet::FlowDirs
**
**  Computes flow directions for each node using a
**  steepest-direction routing algorithm. Flow proceeds along the
**  steepest of the directed edges (spokes) around a given node.
**  This edge is assigned as the node's flowedg. Nodes at which
**  the steepest slope is negative (uphill) are flagged as sinks.
**  (If the lake-fill option is active, sink nodes are processed
**  in order to find an outlet; otherwise, all runoff reaching a
**  sink node is assumed to evaporate and/or infiltrate).
**
**      Parameters:     none
**      Called by: UpdateNet
**      Modifies: node flow directions and flood status
**      Assumes:
**       - each node points to a valid flow edge (returned by getFlowEdg())
**       - each edge has a valid counter-clockwise edge
**       - edge slopes are up to date
**      Updated: 12/19/97 SL; 12/30/97 GT
**
\****************************************************************************/
#define kMaxSpokes 100
void tStreamNet::FlowDirs()
{
   tMesh< tLNode >::nodeListIter_t i( meshPtr->getNodeList() );  // gets nodes from the list
   double slp=0;                          // steepest slope found so far
   double meanderslp = 0;		// steepest meander slope found so far
   double selectslope;			// value of the selected slope
   tLNode *curnode;                     // ptr to the current node
   //XtLNode *newnode;  // ptr to new downstream node
   //XtLNode * nbr;
   tEdge * firstedg(0);   // ptr to first edg
   tEdge * curedg;     // pointer to current edge
   tEdge * nbredg(0);     // steepest neighbouring edge so far
   tEdge * meanderedg = NULL;  // steepest meander edge so far

   int ctr;

   // Find the connected edge with the steepest slope
   curnode = i.FirstP();
   while( i.IsActive() )  // DO for each non-boundary (active) node   = LOOP OVER NODES !
   {
      selectslope = 0.0;
      curnode->setFloodStatus( tLNode::kNotFlooded );  // Init flood status flag
      firstedg =  curnode->getFlowEdg();
      if( unlikely(firstedg == 0) ) {
          curnode->TellAll();
	  assert( 0 );
      }
      slp = firstedg->getSlope();
      nbredg = firstedg;
	  if(0) //DEBUG
	  {
        if(curnode->getID()==8121 /*|| curnode->getID()==213*/) {
           tLNode * nbr = static_cast<tLNode *>(firstedg->getDestinationPtrNC());
           std::cout<<"FlowDirs 1: node "<<curnode->getID()<<" edge "<<nbredg->getID()<<" slp "<<slp<<" downstream nbr "<<nbr->getID()<<std::endl;
           std::cout<<"z "<<curnode->getZ()<<" dsn z "<<nbr->getZ();
           std::cout<<" meander "<<curnode->Meanders()<<" nbr mndr "<<nbr->Meanders()<<std::endl;
        }
	  }
      curedg = firstedg->getCCWEdg();			// Go to the next counter clockwise edge
      ctr = 0;

      /*******************************************************************\
      ** MEANDER - SPECIFIC
      ** If the node meanders, please check whether it is still connected
      ** to another downstream meander node. If yes, get the  connecting
      ** spoke and its downstream slope
      \*******************************************************************/
#define FIXINLETMEANDERBUG 1
#if FIXINLETMEANDERBUG
      if( curnode->Meanders() || curnode==inlet.innode ){
#else
      if( curnode->Meanders() ){
#endif
      	// if the current node meeanders, check whether its current downstream neighbor also meanders
      	// (should be..)
      	tLNode *NodeAlongEdge =
		   static_cast<tLNode *>(firstedg->getDestinationPtrNC());
	    meanderslp = 0.0;
	    meanderedg = NULL;
      	if( NodeAlongEdge->Meanders()) {
	       meanderslp = firstedg->getSlope();
	       meanderedg = firstedg;
		   if(0) //DEBUG
	       {
              if(curnode->getID()==8121 /*|| curnode->getID()==213*/)
                 std::cout<<"FlowDirs: just set meanderslp+edg = "
			       <<" meanderslp "<<meanderslp<<" meanderedg "<<meanderedg->getID()<<std::endl;
           }
      	}
      }
	  if(0) //DEBUG
	  {
        if(curnode->getID()==8121 /*|| curnode==inlet.innode*/ ) {
           tLNode * nbr = static_cast<tLNode *>(firstedg->getDestinationPtrNC());
           std::cout<<"FlowDirs 2: node "<<curnode->getID()<<" edge "<<nbredg->getID()<<" slp "<<slp<<" downstream nbr "<<nbr->getID()<<std::endl;
           std::cout<<"z "<<curnode->getZ()<<" dsn z "<<nbr->getZ();
           std::cout<<" meander "<<curnode->Meanders()<<" nbr mndr "<<nbr->Meanders()
			<<" meanderslp "<<meanderslp<<" meanderedg ";
		   if( meanderedg!=NULL ) std::cout<<meanderedg->getID()<<std::endl;
		   else std::cout<<"NULL\n";
        }
	  }

     /*************************************************************\
     ** Standard: Check all existing spokes for the steepest
     ** downstream direction
     \**************************************************************/

      while( curedg!=firstedg )
      {
         assert( curedg != 0 );
         if ( curedg->getSlope() > slp && curedg->FlowAllowed())

         {
            slp = curedg->getSlope();
            nbredg = curedg;

         }
         curedg = curedg->getCCWEdg();
         ctr++;
         if( unlikely(ctr>kMaxSpokes) ) // Make sure to prevent endless loops
         {
            std::cerr << "Mesh error: node " << curnode->getID()
                 << " going round and round"
                 << std::endl;
            ReportFatalError( "Bailing out of FlowDirs()" );
         }
      }

      /***************************************************************************************\
      ** MEANDER - SPECIFIC
      ** Now make a choice. Compare the steepest descent spoke with the existing meander spoke
      ** if both go to meandering nodes, select the one with the steepest spoke
      ** if not, give preference to the meandering one, also if the normal spoke is steeper
      \***************************************************************************************/
#if FIXINLETMEANDERBUG
	  if( ( curnode->Meanders() || curnode==inlet.innode ) && meanderedg != NULL ){
#else
      if(curnode->Meanders() && meanderedg != NULL ){
#endif
      	 tLNode *SteepestDescentNode =
	        static_cast<tLNode *>(nbredg->getDestinationPtrNC());

	     // the steepest descent one is meandering
	     if( SteepestDescentNode->Meanders() ){
	        if(slp > meanderslp){
	           if(0) //DEBUG
		          if( curnode->getID()==8121 || curnode->getID()==8122 ) std::cout << "FlowDirs: steepest desc mnds, change dir\n";
	           curnode->setFlowEdg( nbredg);
	           selectslope = slp;
			 }
			 else if(slp <= meanderslp){
				 if(0) //DEBUG
					 if( curnode->getID()==8121 || curnode->getID()==8122 ) std::cout << "FlowDirs: cur mndr IS steepest\n";
				 curnode->setFlowEdg( meanderedg);
				 selectslope = meanderslp;
			 }
		 } // end i
		 // the steepest descent one is not meandering
		 else if ( !SteepestDescentNode->Meanders() ){
			 // pick the meander edge if it is positive and the steeper choice
			 // does not lead to an open boundary
#define TESTFIX 0
			 if( TESTFIX )
			 {
				 if( nbredg->getDestinationPtr()->getBoundaryFlag() != kOpenBoundary )
				 {
					 curnode->setFlowEdg( meanderedg);
					 selectslope = meanderslp;
				 }
				 else{
					 curnode->setFlowEdg( nbredg);
					 selectslope = slp;
					 if(0) //debug
					 {
						 std::cout << "Case meand->nonmeand invoked at node " << curnode->getX() << " " << curnode->getY() << std::endl;
						 std::cout << "meanderslp = " << meanderslp << std::endl;
					 }
				 }
			 }
			 else // NOT TESTFIX
			 {
				 if(meanderslp > 0.0 && 
					nbredg->getDestinationPtr()->getBoundaryFlag() != kOpenBoundary)
				 {
					 if(0) //DEBUG
						 if( curnode->getID()==8121 || curnode->getID()==8122 ) std::cout << "FlowDirs: steepest doesn't mdr, staying w/ current dir\n";
					 curnode->setFlowEdg( meanderedg);
					 selectslope = meanderslp;
					 if(0) //DEBUG
					 {
						 if(curnode->getID()==8121 || curnode->getID()==8122 ) {
							 tLNode * nbr = static_cast<tLNode *>(meanderedg->getDestinationPtrNC());
							 std::cout<<"FlowDirs 2A: node "<<curnode->getID()<<" edge "<<meanderedg->getID()<<" slp "<<meanderslp<<" downstream nbr "<<nbr->getID()<<std::endl;
							 std::cout<<"z "<<curnode->getZ()<<" dsn z "<<nbr->getZ();
							 std::cout<<" meander "<<curnode->Meanders()<<" nbr mndr "<<nbr->Meanders()<<std::endl;
						 }
					 }
				 }
				 else
				 {
					 curnode->setFlowEdg( nbredg);
					 selectslope = slp;
					 if(0) //debug
					 {
						 std::cout << "FlowDirs: Case meand->nonmeand invoked at node " << curnode->getX() << " " << curnode->getY() << " ";
						 std::cout << "meanderslp = " << meanderslp << std::endl;
					 } // end if
				 } // end else
			 } // end else
		 } // end else if
		 
	  } // end if
      else{ // all other cases, no menadering nodes involved
		  curnode->setFlowEdg( nbredg );
		  selectslope = slp;
      }
	  
	  if(0) //DEBUG
	  {
		  if(curnode->getID()==8121 || curnode->getID()==8122 ) {
		      tEdge * debugedg = curnode->getFlowEdg();
			  tLNode * nbr = static_cast<tLNode *>(debugedg->getDestinationPtrNC());
			  std::cout<<"FlowDirs 3: node "<<curnode->getID()<<" edge "<<debugedg->getID()<<" slp "<<selectslope<<" downstream nbr "<<nbr->getID()<<std::endl;
			  std::cout<<"z "<<curnode->getZ()<<" dsn z "<<nbr->getZ();
			  std::cout<<" meander "<<curnode->Meanders()<<" nbr mndr "<<nbr->Meanders()<<std::endl;
		  }
	  }
	  
	  
#if 1
      // ocasionally there are bumps in the meandering channel.
      // Even when all normal erosion and deposition functions
      // are swithched off (k's are 0)! Modify the downstream elevation
      // by flattening the occasional bump.

      //tLNode *secondnode = curnode->getDownstrmNbr();
      //tLNode *thirdnode  = secondnode->getDownstrmNbr();

      //if(  curnode->Meanders() && curnode != thirdnode ){

      	//if( secondnode->Meanders() && secondnode->getZ() < curnode->getZ() ){
      	  //curnode->setFlowEdg( firstedg);
      	//}
      	//else if ( secondnode->Meanders() && secondnode->getZ() >= curnode->getZ() ){
      	  //double newelev = (curnode->getZ() + thirdnode->getZ() )/2.0;
      	  //secondnode->setZ(newelev);
      	  //curnode->setFlowEdg( firstedg);
      	  //curedg->getSlope();
      	  //std::cout<<"in FlowDirs: Flattening a bump in the channel: "<<std::endl;
      	  //std::cout<<curnode->getX()<<' '<<curnode->getY()<<' '<<curnode->getZ()<<std::endl;
      	  //std::cout<<secondnode->getX()<<' '<<secondnode->getY()<<' '<<secondnode->getZ()<<std::endl;
      	  //std::cout<<thirdnode->getX()<<' '<<thirdnode->getY()<<' '<<thirdnode->getZ()<<std::endl;
      	  //exit(1);
      	//}
      	//else{
      	  //curnode->setFlowEdg( nbredg );
      	//}
      //}
      //else{
       // curnode->setFlowEdg( nbredg );
      //}
#endif

      //add a wrinkle: if node is a meander node and presently flows
      //to another meander node and the new 'nbredg' does not lead to a
      //meander node, then choose a random number and
      //compare it to the probability that a meander node will change
      //flow direction to a non-meander node
      /*if( mndrDirChngProb != 1.0 )
        {
        newnode = (tLNode *) nbredg->getDestinationPtrNC();
        if( curnode->getDownstrmNbr()->Meanders() &&
        curnode->getDownstrmNbr()->getZ() < curnode->getZ() &&
        !(newnode->Meanders()) )
        {
        chngnum = ran3( &seed );
        if( chngnum <= mndrDirChngProb ) curnode->setFlowEdg( nbredg );
        }
        else curnode->setFlowEdg( nbredg );
        }
        else curnode->setFlowEdg( nbredg );*/




      if(selectslope <= 0.0 && curnode->Meanders()){
      	std::cout<<"WARNING-Type 1, from tStreamNet::CalcSlopes....detected a meander node without positive drainage"<<std::endl;
        std::cout<<"ID= "<<curnode->getID()<<", X= "<<curnode->getX()<<", Y= "<<curnode->getY()<<", Z= "<<curnode->getZ()<<std::endl;

        //DebugShowNbrs( curnode );
        //exit(1);
      }

       // If the selected node has a positve slope
      if( (selectslope>0) && (curnode->getBoundaryFlag() != kClosedBoundary) ){
          curnode->setFloodStatus( tLNode::kNotFlooded );

      }
      else{
         curnode->setFloodStatus( tLNode::kSink );
         if( curnode->Meanders() ){
          std::cout<<"WARNING-Type 2, from tStreamNet::CalcSlopes....detected a meander node without positive drainage"<<std::endl;
          std::cout<<"ID= "<<curnode->getID()<<", X= "<<curnode->getX()<<", Y= "<<curnode->getY()<<", Z= "<<curnode->getZ()<<std::endl;

          //DebugShowNbrs( curnode );
          //exit(1);
         }

      }
	  
      if(0) //DEBUG
      {
	     if(curnode->getID()==8121 || curnode->getID()==8122 ) {
		   tEdge * debugedg = curnode->getFlowEdg();
		   tLNode * nbr = static_cast<tLNode *>(debugedg->getDestinationPtrNC());
		   std::cout<<"FlowDirs 4: node "<<curnode->getID()<<" edge "<<debugedg->getID()<<" slp "<<selectslope<<" downstream nbr "<<nbr->getID()<<std::endl;
		   std::cout<<"z "<<curnode->getZ()<<" dsn z "<<nbr->getZ();
		   std::cout<<" meander "<<curnode->Meanders()<<" nbr mndr "<<nbr->Meanders()<<std::endl;
	     }
      }
	  
      curnode = i.NextP();                // Continue with the next node on the active list
	  
   } // end of node loop

   if (0) //DEBUG
     std::cout << "FlowDirs() finished" << std::endl;
}
#undef kMaxSpokes
#undef FIXINLETMEANDERBUG


/*****************************************************************************\
**
**  tStreamNet::DrainAreaVoronoi
**
**  Computes drainage area for each node by summing the Voronoi areas of all
**  nodes that drain to it, using the following algorithm:
**
**    Reset drainage area for all active nodes to zero
**    FOR each active node
**      Cascade downstream, adding starting node's Voronoi area to each
**         downstream node's drainage area, until an outlet or sink is reached
**    IF there is an inlet, add its associated drainage area to all nodes
**         downstream
**
**  Note that each node's drainage area includes its own Voronoi area.
**
**    Calls: RouteFlowArea, tLNode::setDrArea, tInlet::FindNewInlet
**    Modifies:  node drainage area
**
\*****************************************************************************/
void tStreamNet::DrainAreaVoronoi()
{
   if (0) //DEBUG
     std::cout << "DrainAreaVoronoi()..." << std::endl;

   tLNode * curnode;
   tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );

   // Reset drainage areas to zero
   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
       curnode->setDrArea( 0. );

   // send voronoi area for each node to the node at the other end of the
   // flowedge and downstream
   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
   {

      // Debug Quintijn
      if(curnode->getVArea() < 0.0){
        std::cout<< "Voronoi area <  0.0 at \n";
        std::cout<< curnode->getX() <<' '<<curnode->getY()<<std::endl;
        std::cout<< "Area= "<<curnode->getVArea()<<std::endl;
      	exit(1);
      }

      RouteFlowArea( curnode, curnode->getVArea() );
   }
   if( inlet.innode != 0 )
   {
#if 1
      //inlet.FindNewInlet();            QC commented this out -> avoid channel hopping around
#endif
      RouteFlowArea( inlet.innode, inlet.inDrArea );
   }
   if (0) //DEBUG
     std::cout << "DrainAreaVoronoi() finished" << std::endl;
}


/*****************************************************************************\
**
**  tStreamNet::FlowPathLength
**
**  Computes the longest flow path length from divide to a node, for each
**  node on the mesh. This is used to approximate peak discharge.
**
\*****************************************************************************/
void tStreamNet::FlowPathLength()
{
  // Local variables
  tLNode * curnode,    // Pointer to the current node
    * downstreamNode;  // Pointer to current node's downstream neighbor
  double localPathLength;  // Potential flow path length to downstream nbr

  // Get list of nodes and node iterator
  tMesh< tLNode >::nodeListIter_t nodeIter( meshPtr->getNodeList() );

  // Sort nodes in upstream-to-downstream order
  SortNodesByNetOrder( false );

  // Reset all flow path lengths to zero
  for( curnode = nodeIter.FirstP(); nodeIter.IsActive();
       curnode = nodeIter.NextP() )
    curnode->setFlowPathLength( 0.0 );

  // Work through all active nodes, from upstream to downstream, setting
  // the flow path length of each node's downstream neighbor to the
  // maximum of (a) the current node's flow path length plus the length of
  // the flow edge, or (b) the downstream node's existing flow path length.
  for( curnode = nodeIter.FirstP(); nodeIter.IsActive();
       curnode = nodeIter.NextP() )
    {
      // Compute "local" flow path length to the downstream neighbor --
      // equal to flow path length at the current node plus the length of
      // the edge connecting current node to its downstream neighbor
      localPathLength = curnode->getFlowPathLength()
	+ ( curnode->getFlowEdg() )->getLength();

      // Get a pointer to the downstream neighbor
      downstreamNode = curnode->getDownstrmNbr();

      // Compare the local flow path length to the downstream neighbor's
      // current value of flow path length, which might have been set via
      // another of its upstream nodes -- if the "local" route is longer,
      // set downstream's flow path length to localPathLength.
      if( localPathLength > downstreamNode->getFlowPathLength() )
	downstreamNode->setFlowPathLength( localPathLength );
    }

}

/*****************************************************************************\
**
**  tStreamNet::RouteFlowHydrographPeak
**
**
**
\*****************************************************************************/
void tStreamNet::RouteFlowHydrographPeak()
{
  // Local variables
  tLNode * curnode;
  double runoff = rainrate - infilt;
  double stormdur = stormPtr->getStormDuration();
  double travelTime,  // Travel time from divide to this point
    Qp;               // Peak discharge

  assert( runoff >= 0.0 );

  // Get iterator for list of nodes
  tMesh< tLNode >::nodeListIter_t nodeIter( meshPtr->getNodeList() );

  // Set peak discharge for each node
  for( curnode=nodeIter.FirstP(); nodeIter.IsActive(); curnode=nodeIter.NextP() )
    {
      travelTime = curnode->getFlowPathLength() / mdFlowVelocity;
      Qp = ( runoff * curnode->getDrArea() * stormdur ) /
	( mdHydrgrphShapeFac * ( stormdur + travelTime ) );
      curnode->setDischarge( Qp );
    }

}



/*****************************************************************************\
**
**  tStreamNet::RouteFlowArea
**
**  Starting with the current node, this routine increments
**  the drainage area of the node and each node downstream by _addedArea_.
**
\*****************************************************************************/
inline
void tStreamNet::RouteFlowArea( tLNode *curnode, double addedArea )
{
   if (0) //DEBUG
     std::cout << "RouteFlowArea()..." << std::endl;
//#if DEBUG
   int niterations=0;  // Safety feature: prevents std::endless loops
//#endif

   // As long as the current node is neither a boundary nor a sink, add
   // _addedArea_ to its total drainage area and advance to the next node
   // downstream
   while( (curnode->getBoundaryFlag() == kNonBoundary) &&
	  (curnode->getFloodStatus()!=tLNode::kSink) )
   {
      curnode->AddDrArea( addedArea );
      curnode = curnode->getDownstrmNbr();
//#if DEBUG
      niterations++;
      if( unlikely(niterations>9990) )
	RouteError( curnode );
//#endif
   }
   if (0) //DEBUG
     std::cout << "RouteFlowArea() finished" << std::endl;
}


/*****************************************************************************\
**
**  tStreamNet::RouteRunoff
**
**  This is a variant of RouteFlowArea that adds addedArea and addedRunoff
**  to the drainage area and discharge, respectively, of each downstream node.
**
**  Parameters: curnode -- node to start with
**              addedArea -- drainage area to add to downstream nodes
**              addedRunoff -- runoff rate to add to downstream nodes
**  Created: 6/98 GT
**
\*****************************************************************************/
inline
void tStreamNet::RouteRunoff( tLNode *curnode, double addedArea,
                              double addedRunoff )
{
   if (0) //DEBUG
     std::cout << "RouteRunoff()..." << std::endl;
//#if DEBUG
   int niterations=0;  // Safety feature: prevents std::endless loops
//#endif

   // As long as the current node is neither a boundary nor a sink, add
   // _addedArea_ to its total drainage area and _addedRunoff_ to its total
   // discharge and then advance to the next node downstream.
   while( (curnode->getBoundaryFlag() == kNonBoundary)
	  && (curnode->getFloodStatus()!=tLNode::kSink) )
   {
      curnode->AddDrArea( addedArea );
      curnode->AddDischarge( addedRunoff );
      curnode = curnode->getDownstrmNbr();
//#if DEBUG
      niterations++;
      if( unlikely(niterations>9990) )
	RouteError( curnode );
//#endif
   }
   if (0) //DEBUG
     std::cout << "RouteRunoff() finished" << std::endl;
}

/*****************************************************************************\
**
**  tStreamNet::RouteError
**
**  Error handler
**  Created: 8/2003 AD
**
\*****************************************************************************/
void tStreamNet::RouteError( tLNode *curnode )
{
  while( (curnode->getBoundaryFlag() == kNonBoundary)
	 && (curnode->getFloodStatus()!=tLNode::kSink) )
    {
      curnode->TellAll();
      std::cout << std::flush;
      curnode = curnode->getDownstrmNbr();
    }
  ReportFatalError("Maximum number of iterations exceeded in RouteError.");
}


/*****************************************************************************\
**
**      MakeFlow : flow routing functions
**
**      Data members updated:
**      Called by:
**      Calls: FillLakes, DrainAreaVoronoi, FlowSaturated, FlowUniform
**
**      Created:  YC
**      Added:   YC
**      Modified:  GT made it a mbr func of tMesh 7/97
**         updated 12/19/97 SL
**         - GT added sinusoidal variation in infiltration capacity, and tm
**             parameter to go along w/ it
**
\*****************************************************************************/
void tStreamNet::MakeFlow( double tm )
{
   if (0) //DEBUG
     std::cout << "MakeFlow()..."<<std::flush;

   if( filllakes ) FillLakes();
   
   if(0) //DEBUG
   {
         tMesh< tLNode >::nodeListIter_t mli( meshPtr->getNodeList() );  // gets nodes from the list
		 tLNode * cn;
		 for( cn=mli.FirstP(); mli.IsActive(); cn=mli.NextP() )
		 {
		    if( cn->getID()==8121 ) 
			{
		       tEdge * debugedg = cn->getFlowEdg();
		       tLNode * nbr = static_cast<tLNode *>(debugedg->getDestinationPtrNC());
		       std::cout<<"MakeFlow 1: node "<<cn->getID()<<" edge "<<debugedg->getID()<<" downstream nbr "<<nbr->getID()<<std::endl;
		       std::cout<<"z "<<cn->getZ()<<" dsn z "<<nbr->getZ() << " mdr "<<cn->Meanders()<<" nbr mdr "<<nbr->Meanders()<<std::endl;
            }
		}
   }

   DrainAreaVoronoi();

   // If a hydrologic parameter varies through time, update it here
   // (currently, only infiltration capacity varies)
   if( optSinVarInfilt && infilt>0 )
       infilt = infilt0 + infilt_dev * sin( tm*twoPiLam );

   // Call appropriate function for runoff generation option
   switch( miOptFlowgen )
   {
      case kSaturatedFlow1:   // Topmodel-like runoff model with steady-state
          FlowSaturated1();   //   surface flow and possible return flow.
          break;
      case kSaturatedFlow2:   // Topmodel-type runoff model without return
          FlowSaturated2();   //   flow (infiltrated water does not contribute
          break;              //   to storm flow)
      case kConstSoilStore:   // "Bucket" model: spatially uniform soil storage
          FlowBucket();       //   capacity; any excess contributes to runoff.
          break;
      case kHydrographPeakMethod:
	  FlowPathLength();
	  RouteFlowHydrographPeak();
	  break;
      case k2DKinematicWave:
	  RouteFlowKinWave( rainrate );
	  break;
      case kHortonian:
          FlowUniform();      // Spatially uniform infiltration-excess runoff
	  break;
   }

   if (0) //DEBUG
     std::cout << "MakeFlow() finished" << std::endl;
}


/*****************************************************************************\
**
**  tStreamNet::FlowUniform
**
**  Computes discharge as the product of precip and drainage area. 3/97 GT
**
**  Called by:  main
**  Modifications:
**    - added infiltration parameter (default zero) (8/97 GT)
**    - updated 12/19/97 SL
**
\*****************************************************************************/
void tStreamNet::FlowUniform()
{
   if (0) //DEBUG
     std::cout << "FlowUniform..." << std::endl;
   tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );
   tLNode *curnode;
   double runoff = rainrate - infilt;
   double discharge;

   if( runoff<0 ) runoff = 0;  // Make sure runoff isn't negative
   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
   {
      discharge = curnode->getDrArea() * runoff;
      curnode->setDischarge( discharge );
   }
   if (0) //DEBUG
     std::cout << "FlowUniform finished" << std::endl;
}


/*****************************************************************************\
**
**  tStreamNet::FlowSaturated1
**
**  Computes surface runoff using Topmodel concept. Steady-state
**  subsurface flow capacity is transmissivity times slope. Total
**  steady-state runoff is rainfall rate times drainage area. Surface
**  runoff is total minus subsurface. The assumption of steady-state
**  subsurface flow implies that subsurface storm flow can contribute
**  significantly to erosion (in other words, both saturation-excess
**  and return flow are modeled). To model cases in which return flow
**  is insigificant, use FlowSaturated2.
**
**  Transmissivity has units of L^2/T; it is hydraulic conductivity
**  times depth. To get a volume discharge, we must multiply the
**  transmissivity by the Voronoi edge length for the flow edge.
**
**  Called by:  main
**  Modifies: node discharge
**  Modifications:
**    - updated 12/19/97 SL
**    - 5/2/98 SL
**
\*****************************************************************************/
void tStreamNet::FlowSaturated1()
{
   if (0) //DEBUG
     std::cout << "FlowSaturated1...";
   tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );
   tLNode *curnode;
   tEdge *fedg;
   double total_discharge, surface_discharge, subsurf_discharge;

   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
   {
      fedg = curnode->getFlowEdg();
      total_discharge = curnode->getDrArea() * rainrate;
      subsurf_discharge = fedg->getSlope() * fedg->getVEdgLen() * trans;
      surface_discharge = total_discharge - subsurf_discharge;
      if( surface_discharge < 0. ) {
	surface_discharge = 0.;
	subsurf_discharge = total_discharge;
      }
      curnode->setDischarge( surface_discharge );
      curnode->setSubSurfaceDischarge( subsurf_discharge );
   }
   if (0) //DEBUG
     std::cout << "finished" << std::endl;
}


/*****************************************************************************\
**
**  tStreamNet::FlowSaturated2
**
**  Computes surface runoff using a more precise formulation of the
**  Topmodel concept, based on Ijjasz-Vasquez et al (1992). Total
**  runoff from a node is the sum of infiltration-excess (if any)
**  and saturation-excess (if any). Saturation-excess is determined
**  by computing a saturation deficit that depends on slope, area,
**  and a parameter (trans), and subtracting the deficit from the storm
**  depth to obtain the saturation-excess runoff depth, which is
**  converted to a rate by multiplying by the storm duration.
**
**  Parameters:  none
**  Called by:  main
**  Modifies:  node discharges
**
\*****************************************************************************/
void tStreamNet::FlowSaturated2()
{
   tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );
   tLNode *curnode;
   //XtEdge *fedg;
   //Xdouble discharge;
   double infiltExRunoff,  // Infiltration-excess runoff rate (L/T)
       sd,                 // Storm depth minus any infilt-excess runoff (L)
       asRatio,            // Slope-area-per-width ratio (ie, A/Sb) (L),
                           //   where b is Voronoi edge length
       satDeficit,         // Saturation deficit (L)
       rsat,               // Saturation-excess runoff depth (L)
       runoff,             // Total runoff from node (L/T)
       stormDur = stormPtr->getStormDuration();  // Storm duration
   int nsat=0,nsr=0,nhort=0,nflat=0; // 4dbg

   if (0) //DEBUG
     std::cout << "FlowSaturated1" << std::endl;

  // Reset drainage areas and discharges to zero
  for ( curnode=nodIter.FirstP(); nodIter.IsActive(); curnode=nodIter.NextP() )
  {
     //if( curnode->drarea>1.15e7 )
     //    std::cout << "Q(" << curnode->id << ") " << curnode->q << std::endl;
     curnode->setDischarge( 0. );
  }

  for ( curnode=nodIter.FirstP(); nodIter.IsActive(); curnode=nodIter.NextP() )
  {
     //std::cout<<"Node " <<curnode->getID();
     assert( curnode->getFlowEdg() != 0 );
     satDeficit = 0;
     infiltExRunoff = rainrate - infilt;
     if( infiltExRunoff<0 ) infiltExRunoff = 0;
     else nhort++;
     sd = (rainrate - infiltExRunoff)*stormDur;
     if( curnode->calcSlope()>0 )
     {
         asRatio = curnode->getDrArea() /
             ( curnode->calcSlope() * curnode->getFlowEdg()->getVEdgLen() );
         assert( asRatio>0 );
         if( asRatio < trans  ) // if true, then NOT saturated => deficit
         {
            assert( trans>0 && asRatio>0 );
            satDeficit = log( trans / asRatio );
         }
         else nsat++;
         //std::cout << " satdef " << satDeficit;
     }
     else {
        //std::cout<<" zero slp or ve";
        nflat++;
        nsat++;
     }
     rsat = sd - satDeficit;
     if( rsat<0 ) rsat = 0;
     else nsr++;
     runoff = infiltExRunoff + rsat/stormDur;
     //std::cout<<" sat excess " << rsat << " total " << runoff << std::endl;
     RouteRunoff( curnode, curnode->getVArea(), runoff*curnode->getVArea() );
  }

  if (0) //DEBUG
    std::cout << nhort << " generate Horton runoff, " << nsat
	 << " pre-saturated, (" << nflat << "flat) "
	 << nsr << " saturate during storm.\n";

}


/*****************************************************************************\
**
**  tStreamNet::FlowBucket
**
**  Computes runoff rate and discharge assuming a spatially and temporally
**  uniform "bucket" soil storage capacity. Runoff rate is the sum of
**  infiltration excess runoff ( = rainrate - infilt ), if any, and
**  saturation excess runoff, if any. Saturation excess runoff is
**  calculated as ( total water depth of infiltrated storm water
**  minus soil storage capacity ) divided by storm duration, where the
**  infiltration depth is just the rainfall depth minus any infiltration-
**  excess runoff. Discharge is runoff rate times contributing area.
**    We save a multiplication operation by noting that:
**      satEx = ( infiltDepth - soilStore ) / stormDuration
**            = (rainrate - infiltEx) - (soilStore/stormDuration)
**
**  Parameters: none
**  Called by:  MakeFlow
**  Assumes:  stormPtr is valid -- ie, to use this function, StreamNet
**                                  must be init'd with a storm object,
**                                  or UpdateNet must be called w/ a
**                                  storm object.
**  Created:  5/98 GT
**  Modifications:
**
\*****************************************************************************/
void tStreamNet::FlowBucket()
{
   assert( stormPtr!=0 );
   if (0) //DEBUG
     std::cout << "FlowBucket..." << std::endl;
   tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );
   tLNode *curnode;
   double infiltEx=0.0,  /* infiltration excess runoff (L/T) */
       satEx=0.0,        /* saturation excess runoff (L/T) */
       infiltRate,       /* soil store / storm duration */
       runoff,           /* total runoff */
       discharge;

   // Compute total runoff as the sum of infiltration excess (if any)
   // and saturation excess (if any)
   rainrate = stormPtr->getRainrate();
   if( rainrate > infilt )
       infiltEx = rainrate - infilt;
   if( (rainrate-infiltEx) >
       (infiltRate=(soilStore/(stormPtr->getStormDuration()) ) ) )
       satEx = (rainrate - infiltEx) - infiltRate;
   runoff = infiltEx + satEx;
   //std::cout << "  R " << runoff << " = IEx " << infiltEx << " + SEx " << satEx << std::endl;

   // Compute and assign discharge for each node
   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
   {
      discharge = curnode->getDrArea() * runoff;
      curnode->setDischarge( discharge );
   }
   if (0) //DEBUG
     std::cout << "FlowBucket finished" << std::endl;
}


/*****************************************************************************\
**
**  tStreamNet::FillLakes
**
**  Finds drainage for closed depressions. The algorithm assumes
**  that sinks (nodes that are lower than any of their neighbors)
**  have already been identified during the flow directions
**  procedure. For each sink, the algorithm creates a list of
**  nodes in the current lake, which initially is just the sink
**  itself. It then iteratively looks for the lowest node on the
**  perimeter of the current lake. That node is checked to see
**  whether it can be an outlet, meaning that one of its
**  neighbors is both lower than itself and is not already
**  flooded (or is an open boundary). If the low node is not an
**  outlet, it is added to the current lake and the process
**  is repeated. If it is an outlet, then all of the nodes on the
**  current-lake list are identified draining it. The list is then
**  cleared, and the next sink is processed. If during the search
**  for the lowest node on the perimeter a flooded node is found
**  that isn't already part of the current lake (i.e., it was
**  flagged as a lake node when a previous sink was processed),
**  then it is simply added to the current-lake list --- in other
**  words, the "new" lake absorbs any "old" ones that are encountered.
**
**  Once an outlet has been found, flow directions for nodes in the
**  lake are resolved in order to create a contiguous path through
**  the lake.
**
**    Calls: FindLakeNodeOutlet
**    Called by: MakeFlow
**    Modifies:  flow direction and flood status flag of affected nodes
**    Created: 6/97 GT
**    Modifications:
**     - fixed memory leak on deletion of lakenodes 8/5/97 GT
**     - updated: 12/19/97 SL
**     - exploded and optimized: 08/2003 AD
**
\*****************************************************************************/
void tStreamNet::FillLakes()
{
   if (0) //DEBUG
   {
     std::cout << "FillLakes()..." << std::endl;
   }
   int debugcount=0; //DEBUG

   tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() ); // node iterator

   // Check each active node to see whether it is a sink
   for( tLNode *cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP() )
   {
      if( cn->getFloodStatus() == tLNode::kSink )
      {
	     // Create a new lake-list, initially containing just the sink node.
	     tPtrList< tLNode > lakeList;

	     // Build lakeList and iteratively search for an outlet along the
	     // perimeter of the lake
	     tLNode *lowestNode = BuildLakeList( lakeList, cn );

   
         if(0) //DEBUG
         {
		    debugcount++;
            tMesh< tLNode >::nodeListIter_t mli( meshPtr->getNodeList() );  // gets nodes from the list
		    tLNode * cn;
		    for( cn=mli.FirstP(); mli.IsActive(); cn=mli.NextP() )
		    {
		       if( cn->getID()==8121 ) break;
		    }
		    tEdge * debugedg = cn->getFlowEdg();
		    tLNode * nbr = static_cast<tLNode *>(debugedg->getDestinationPtrNC());
		    std::cout<<"FillLakes 1 iter " << debugcount << ": node "<<cn->getID()<<" edge "<<debugedg->getID()<<" downstream nbr "<<nbr->getID()<<std::endl;
			std::cout<<"z "<<cn->getZ()<<" dsn z "<<nbr->getZ() << " mdr "<<cn->Meanders()<<" nbr mdr "<<nbr->Meanders()<<std::endl;
         }

	     // Assign a flow path
	     tPtrListIter< tLNode > lakeIter( lakeList ); // Iterator for lake list
	     FillLakesFlowDirs(lakeIter, lowestNode);

	     // Finally, flag all of the nodes in it as "kFlooded"
		 if(0) //DEBUG
		    std::cout<<"FillLakes: " << lakeList.getSize() << " flooded nodes:" <<std::endl;
	     for( tLNode *cln = lakeIter.FirstP(); !( lakeIter.AtEnd() ); cln = lakeIter.NextP() )
		 {
	        cln->setFloodStatus( tLNode::kFlooded );
			if(0) //DEBUG
			{
			   std::cout<<cln->getID();
			   if( cln->Meanders() ) std::cout<<"*";
			   std::cout<<std::endl;
			}
		 }
      } /* END if Sink */
   } /* END Active Nodes */

   if (0) //DEBUG
     std::cout << "FillLakes() finished" << std::endl;

} // end of tStreamNet::FillLakes


/*****************************************************************************\
**
**  tStreamNet::BuildLakeList
**
**  Build LakeList and iteratively search for an outlet along the perimeter of
**  the lake
**
**  Moved from FillLakes (AD 09/2003)
**
\*****************************************************************************/
// insert a new node in lakelist
inline
void insertInLakeList(tPtrList< tLNode > &lakeList, tLNode *cn )
{
  lakeList.insertAtBack( cn );
  cn->setFloodStatus( tLNode::kCurrentLake );
}

tLNode *tStreamNet::BuildLakeList( tPtrList< tLNode > &lakeList, tLNode *cn )
{
  tLNode *lowestNode;      // Lowest node on perimeter found so far
  tPtrListIter< tLNode > lakeIter( lakeList ); // Iterator for lake list

  // insert the first node
  insertInLakeList( lakeList, cn );

  bool done = false;       // Flag indicating whether outlet has been found
  do
    {
      lowestNode = lakeIter.FirstP();
      double lowestElev = kVeryHigh; // Initialize lowest elev to very high val.

      // Check the neighbors of every node on the lake-list
      for( tLNode *cln = lakeIter.FirstP(); !( lakeIter.AtEnd() );
	   cln = lakeIter.NextP() )
	{
	  // Check all the neighbors of the node
	  tEdge *ce = cln->getEdg();
	  do
	    {
	      tLNode *thenode =
		static_cast<tLNode *>(ce->getDestinationPtrNC());
	      switch( thenode->getFloodStatus() ){
		// Is it a potential outlet (ie, not flooded and not
		// a boundary)?
	      case tLNode::kNotFlooded:
		if( ce->FlowAllowed() != tEdge::kFlowNotAllowed)
		  {
		    // Is it lower than the lowest found so far?
		    if( thenode->getZ() < lowestElev )
		      {
			lowestNode = thenode;
			lowestElev = thenode->getZ();
		      }
		  }
		break;
	      // If it's a previous lake node or a sink, add it to the list
	      case tLNode::kFlooded:
	      case tLNode::kSink:
		insertInLakeList( lakeList, thenode );
		break;
	      default:
		break;
	      }
	    } while( ( ce=ce->getCCWEdg() ) != cln->getEdg() );// END spokes


	} /* END lakeList */

      // Now we've found the lowest point on the perimeter. Now test
      // to see whether it's an outlet. If it's an open boundary, it's
      // an outlet...
      if( lowestNode->getBoundaryFlag() == kOpenBoundary ) done = true;
      else // ...it's also an outlet if it can drain to a "dry" location.
	{
	  // Can lowestNode drain to a non-flooded location?
	  if( FindLakeNodeOutlet( lowestNode ) ) done = true;
	  // no, it can't, so add it to the list and continue:
	  else
	    insertInLakeList( lakeList, lowestNode );
	}
      if( unlikely(lakeList.getSize() >
		   meshPtr->getNodeList()->getActiveSize()) )
	{
	  std::cout << "LAKE LIST SIZE=" << lakeList.getSize() << "\n"
	    "active node size=" << meshPtr->getNodeList()->getActiveSize()
	       << std::endl;
	  std::cerr <<
	    "Error in Lake Filling algorithm: "
	    "Unable to find a drainage outlet.\n"
	    "This error can occur when open boundary node(s) "
	    "are isolated from the interior of the mesh.\n"
	    "This is especially common when a single outlet point "
	    "(open boundary) is used.\n"
	    "Re-check mesh configuration or try changing SEED.\n";
	  ReportFatalError( "No drainage outlet found for one or more interior nodes." );
	}
    } while( !done );
  return lowestNode;
}


/*****************************************************************************\
**
**  tStreamNet::FillLakesFlowDirs
**
**  Moved from FillLakes (AD 08/2003)
**  - AD 08/2003: Speed-up: break early from inner loop and we allow to flow
**    to a kOutletPreFlag node.
**
\*****************************************************************************/
void tStreamNet::FillLakesFlowDirs(tPtrListIter< tLNode > &lakeIter,
				   tLNode *lowestNode) const
{
  // Now we've found an outlet for the current lake.
  // This next bit of code assigns a flowsTo for each node so there's
  // a complete flow path toward the lake's outlet. This isn't strictly
  // necessary --- the nodes could all point directly to the outlet,
  // skipping anything in between --- but it prevents potential problems
  // in ordering the list by network order. This also works by pointing
  // each node toward the first neighboring node they happen to find
  // that has been flagged as having its flow direction resolved.
  // Initially, the low node is thus flagged, and the algorithm repeats
  // until all the lake nodes are flagged as having a flow direction.
  // The algorithm isn't unique---there are many paths that could be
  // taken; this simply finds the most convenient one.
  lowestNode->setFloodStatus( tLNode::kOutletFlag );

  // Test for error in mesh: if the lowestNode is a closed boundary, it
  // means no outlet can be found.
  bool done;  // Flag indicating whether outlet has been found
  do
    {
      bool newOutLetFound = false;
      done = true;  // assume done until proven otherwise
      tLNode *cln;  // current lake node
      for( cln = lakeIter.FirstP(); !( lakeIter.AtEnd() );
	   cln = lakeIter.NextP() )
	{
	  const tLNode::tFlood_t fs = cln->getFloodStatus();
	  assert( fs != tLNode::kOutletPreFlag );
	  if( fs != tLNode::kOutletFlag )
	    {
	      done = false;

	      // Check each neighbor
	      // take the first node with tLNode::kOutletFlag, if not found
	      // take the first node with kOutletPreFlag
	      tEdge * const ce1 = cln->getEdg(); // first edge
	      tEdge *ce = ce1;                   // current edge
	      tEdge *cePreFlag = 0;              // edge with kOutletPreFlag
	      do
		{
		  const tLNode::tFlood_t fsnode =
		    static_cast<tLNode const *>(ce->getDestinationPtr())
		    ->getFloodStatus();
		  if( fsnode == tLNode::kOutletFlag )
		    {     // found one!
		      newOutLetFound = true;
		      cln->setFloodStatus( tLNode::kOutletPreFlag );
		      cln->setFlowEdg( ce );
		      goto foundOutlet;
		    }
		  if( cePreFlag == 0 && fsnode == tLNode::kOutletPreFlag  )
		    cePreFlag = ce;
		} while( ( ce=ce->getCCWEdg() ) != ce1 );
	      if ( cePreFlag != 0 )
		{
		  newOutLetFound = true;
		  cln->setFloodStatus( tLNode::kOutletPreFlag );
		  cln->setFlowEdg( cePreFlag );
		}
	    foundOutlet: ;	      ;
	    } // END if node not flagged as outlet
	} // END for each lake node

      // Now flag all the "preflagged" lake nodes as outlets
      if (newOutLetFound)
	for( cln = lakeIter.FirstP(); !( lakeIter.AtEnd() );
	     cln = lakeIter.NextP() )
	  if( cln->getFloodStatus() == tLNode::kOutletPreFlag )
	    cln->setFloodStatus( tLNode::kOutletFlag );
    } while( !done );
  lowestNode->setFloodStatus( tLNode::kNotFlooded );
}


/*****************************************************************************\
**
**  FindLakeNodeOutlet
**
**  This function is part of the lake-filling algorithm. It checks to see
**  whether there is a valid outlet for the current node, and if so it
**  assigns that outlet. An "outlet" essentially means a downhill neighbor
**  that isn't already flooded to the level of the current node. The function
**  performs basically the same operation as FlowDirs, but with stricter
**  criteria. The criteria for a valid outlet are:
**
**  (1) It must be lower than the current node (slope > 0)
**  (2) It must not be part of the current lake (a lake can't outlet to itself)
**  (3) It must not be a closed boundary (_flowAllowed_ must be TRUE)
**  (4) If the outlet is itself part of a different lake, the water surface
**      elevation of that lake must be lower than the current node.
**
**  Returns: TRUE if a valid outlet is found, FALSE otherwise
**  Calls: (none)
**  Called by: FillLakes
**  Created: 6/97 GT
**  Updated: 12/19/97 SL; 1/15/98 gt bug fix (open boundary condition)
**
\*****************************************************************************/
bool tStreamNet::FindLakeNodeOutlet( tLNode *node ) const
{
   double maxslp = 0;  // Maximum slope found so far
   tEdge * ce;        // Current edge
   tLNode *dn,        // Potential outlet
       *an;           // Node ptr used to find outlet of a previously
                      // identified lake

   // Check all the neighbors
   ce = node->getEdg();
   do
   {
      // If it passes this test, it's a valid outlet
      dn = static_cast<tLNode *>(ce->getDestinationPtrNC());
      assert( dn!=0 );
/*X      if( ce->getSlope() > maxslp &&
          dn->getFloodStatus() != kCurrentLake &&
          ce->FlowAllowed() &&
          ( dn->getBoundaryFlag()==kOpenBoundary ||
            dn->getDownstrmNbr()->getZ() < node->getZ() ) )*/
#define NEWFIX 1
#if NEWFIX
      if( ce->getSlope() > maxslp &&
          dn->getFloodStatus() != tLNode::kCurrentLake &&
          ce->FlowAllowed() && !FlowDirBreaksMeanderChannel( node, ce ) )
#else
      if( ce->getSlope() > maxslp &&
          dn->getFloodStatus() != tLNode::kCurrentLake &&
          ce->FlowAllowed() )
#endif
#undef NEWFIX
      {
         // Handle a very special and rare case: if the "target" node dn is
         // part of a previous lake, it's still a valid exit as long as its
         // water surface elevation is lower than the current lake (whose
         // wse, assuming an outlet is found, would be equal to _node_'s
         // elevation). It can sometimes happen that the target lake's wse is
         // exactly equal in elevation to _node_, in which case
         // the point is not considered an outlet---if it were, infinite loops
         // could result. (This fix added 4/98)
         if( dn->getFloodStatus() == tLNode::kFlooded )
         {
            // Iterate "downstream" through the "old" lake until reaching the
            // outlet, then test its elevation. If the elevation is exactly
            // equal to _node_, skip the rest and go on to the next iteration.
            an = dn;
            while( an->getFloodStatus() != tLNode::kNotFlooded )
                an = an->getDownstrmNbr();
            if( an->getZ()==node->getZ() ) continue;
         }

         // Assign the new max slope and set the flow edge accordingly
         maxslp = ce->getSlope();
         node->setFlowEdg( ce );
         // std::cout << "Node " << node->getID() << " flows to "
         //      << node->getDownstrmNbr()->getID() << std::endl;
      }
   } while( ( ce=ce->getCCWEdg() ) != node->getEdg() );

   return BOOL( maxslp > 0 );
}


/*****************************************************************************\
**
**  FlowDirBreaksMeander()
**
**  Tests whether a potential change in flow direction would break an existing
**  path from one meander node to another. (Note: it returns FALSE even if
**  the potential new path would also connect to a meandering node)
**
**  Returns: TRUE if the new direction breaks an existing meander path,
**           FALSE otherwise.
**  Calls: (none)
**  Called by: FindLakeNodeOutlet
**  Created: 12/04
**  Updated: 
**
\*****************************************************************************/
bool tStreamNet::FlowDirBreaksMeanderChannel( tLNode * theNode, tEdge * potentialFlowEdg ) const
{
   if( !theNode->Meanders() && theNode!=inlet.innode ) return false;   // no problem unless it's a meander node or the inlet
   if( potentialFlowEdg==theNode->getFlowEdg() ) return false; // no prob if no change in dir
   tLNode * nbr = theNode->getDownstrmNbr();
   if( !nbr->Meanders() ) return false;   // no prob if neighbor doesn't meander
   else return true;
}



/*****************************************************************************\
**
**  SortNodesByNetOrder:
**
**  This function sorts the list of nodes according to their order in the
**  network (upstream to downstream), in preparation for computing erosion &
**  deposition. (Note that this is only necessary when the sediment output
**  from a given node depends on the input, e.g. for mixed bedrock-alluvial
**  mode in which case the channel type [br or alluvial] depends on the
**  difference between sediment influx and carrying capacity). The sorting
**  algorithm is based on the "cascade" algorithm of Braun and Sambridge
**  (Basin Research, 1997, vol. 9, p. 27).
**    The single-direction sorting algorithm works by initially assigning a
**  tracer (like a packet of water) to each node. At each iteration, a tracer
**  from each node issent downstream. Any nodes that have zero tracers left
**  are moved to the bottom of the list (a FIFO stack), so that for example
**  the very first node moved will be the first node on the list when the
**  sorting is completed. The process continues until no unsorted nodes
**  remain.
**    The multi-flow option was added to allow for multiple flow directions
**  and kinematic-wave routing. The algorithm is slightly different. At
**  each pass, the unsorted nodes are "de-flagged" by setting their
**  tracer variables to zero. Then each unsorted node flags _all_ of the
**  adjacent nodes that are downhill by setting their tracer to 1 (this
**  is accomplished through a call to ActivateSortTracer). Finally, any
**  unflagged nodes are moved to the back of the list, and the process
**  is repeated until all nodes have been sorted.
**
**  Modifications:
**   - adapted from previous CHILD code by GT, 12/97
**   - multiflow sort capability added 1/2000, GT
**
**  TODO: it is possible that the "flagging" method, as opposed to the
**  "tracer movement" method, is more efficient even for single-flow
**  routing. The added cost lies in unflagging the unsorted nodes at
**  each pass. However, the gain comes from the fact that you never have
**  multiple "tracers" at a node that need to be removed one by one.
**  The two methods should be tested and compared.
**
\*****************************************************************************/
void tStreamNet::SortNodesByNetOrder( bool optMultiFlow )
{
   int nThisPass;                      // Number moved in current iteration
   int i;
   bool done;
   tLNode * cn;
   tMesh< tLNode >::nodeList_t *nodeList = meshPtr->getNodeList();
   int nUnsortedNodes = nodeList->getActiveSize();  // Number not yet sorted
   tMesh< tLNode >::nodeListIter_t listIter( nodeList );

   //std::cout << "SortNodesByNetOrder, optMultiFlow=" << optMultiFlow << std::endl;

   //test
   /*std::cout << "BEFORE: " << std::endl;
   for( cn=listIter.FirstP(); listIter.IsActive(); cn=listIter.NextP() )
   std::cout << cn->getID() << " " << cn->getZ() << " " << cn->getDrArea() << std::endl;*/

   // Iterate: move tracers downstream and sort until no nodes with tracers
   // are left. The only difference with the multiflow option is that we
   // call MoveTracersDownstreamMulti to send tracers down to all lower
   // elevation nodes, not just the steepest
   if( !optMultiFlow )
   {
      // Assign initial tracers: use "qs" field, which contains garbage at
      // this stage.
      for( cn=listIter.FirstP(); listIter.IsActive(); cn=listIter.NextP() )
          cn->ActivateSortTracer();
      do
      {
         // Send tracers downstream
         cn = listIter.FirstP();
         for( i=1; i<=nUnsortedNodes; i++ )
         {
            assert( cn!=0 );
            cn->MoveSortTracerDownstream();
            cn = listIter.NextP();
         }

         // Scan for any nodes that have no tracers, and move them to the
         // bottom of the list.
         tMesh<tLNode>::nodeListNode_t  *nodeToMove;
         nThisPass = 0;
         done = true;
         cn = listIter.FirstP();
         for( i=1; i<=nUnsortedNodes; i++ )
         {
            if( cn->NoMoreTracers() ) // If no tracers, move node to bottom
            {
               nodeToMove = listIter.NodePtr();
               cn = listIter.NextP();
               nodeList->moveToActiveBack( nodeToMove );
               nThisPass++;
            }
            else
            {
               cn = listIter.NextP();
               done = false;
            }
         }

         nUnsortedNodes -= nThisPass;

          /*std::cout << "NO. UNSORTED: " << nUnsortedNodes << std::endl;
            for( cn=listIter.FirstP(); listIter.IsActive(); cn=listIter.NextP() )
            std::cout << cn->getID() << " " << cn->getQ() << " " << cn->getQs()
            << std::endl;*/

      } while( !done );
   }

   else  // For multiple flow directions (e.g., kinematic wave)
       do
       {
          // Reinitialize by unflagging the nodes
          cn = listIter.FirstP();
          for( i=1; i<=nUnsortedNodes; i++ )
          {
             assert( cn!=0 );
             cn->DeactivateSortTracer();
             cn = listIter.NextP();
          }

          // Send tracers downstream
          cn = listIter.FirstP();
          for( i=1; i<=nUnsortedNodes; i++ )
          {
             assert( cn!=0 );
             cn->FlagDownhillNodes();
             cn = listIter.NextP();
          }

          // Scan for any nodes that are unflagged, and move them to the
          // bottom of the list.
          tMesh<tLNode>::nodeListNode_t *nodeToMove;
          nThisPass = 0;
          done = true;
          cn = listIter.FirstP();
          for( i=1; i<=nUnsortedNodes; i++ )
          {
             if( cn->NoMoreTracers() )  // If no tracers, move node to bottom
             {
                nodeToMove = listIter.NodePtr();
                cn = listIter.NextP();
                nodeList->moveToActiveBack( nodeToMove );
                nThisPass++;
             }
             else
             {
                cn = listIter.NextP();
                done = false;
             }
          }

          nUnsortedNodes -= nThisPass;

          /*std::cout << "NO. UNSORTED: " << nUnsortedNodes << std::endl;
            for( cn=listIter.FirstP(); listIter.IsActive(); cn=listIter.NextP() )
            std::cout << cn->getID() << " " << cn->getZ() << " " << cn->getQs()
            << " " << cn->getDrArea() << " " << std::endl;*/

       } while( !done );

   /*std::cout << "AFTER: " << std::endl;
   cn = listIter.FirstP();
   std::cout << "First node:\n";
   cn->TellAll();
   for( cn=listIter.FirstP(); listIter.IsActive(); cn=listIter.NextP() )
       std::cout << cn->getID() << " " << cn->getZ() << " " << cn->getDrArea() << std::endl;
       std::cout << "Leaving Sort\n" << flush;*/


}


/*****************************************************************************\
**
**       FindHydrGeom: goes through reach nodes and calculates/assigns
**                     values for hydraulic geometry.
**
**
**		  Parameters: kwds, ewds, ewstn, knds, ends, enstn, optrainvar
**      Data members updated: tLNode->chan.hydrwidth, hydrnrough, hydrdepth
**      Called by: FindReaches (needs to know how long to make reach "tails"
**      Calls: no "major" functions
**
**      Created: SL 1/98
**      Modified:
**       - calculation of power terms modified to handle case of zeros
**         GT 6/99
**       - GT 6/01 -- added statement to call alternative ParkerChannels
**         model instead of "regime" hydraulic geometry.
**       - 7/03: Parker geometry function now called from FindChanGeom().
**               This allows us to retain the power-law regime functions for
**               at-a-station values while using Parker for bankfull
**               conditions. Also moved the if-else out of the node loop. (GT)
**
\*****************************************************************************/
void tStreamNet::FindHydrGeom()
{

   double kwdspow, kndspow, kddspow,
       widpow, deppow, npow, qpsec;
   double width, depth, rough=0.0, slope;
   tLNode *cn;
   tMesh< tLNode >::nodeListIter_t nIter( meshPtr->getNodeList() );

   // If rainfall and hence discharge varies in time, set flow width, depth
   // and roughness using power law functions of their bankfull values
   if( optrainvar )
     {
       // Set up exponents and coefficients
       // NB: Calculate these just once and store--TODO
       if( edds>0.0 )
	 {
	   kddspow = pow(kdds, edstn / edds);
	   deppow = 1.0 - edstn / edds;
	 }
       else
	 {
	   kddspow = kdds;
	   deppow = 0.0;
	 }
       if( ewds>0.0 )
	 {
	   kwdspow = pow(kwds, ewstn / ewds);
	   widpow = 1.0 - ewstn / ewds;
	 }
       else
	 {
	   kwdspow = kwds;
	   widpow = 0.0;
	 }
       if( ends!=0.0 )
	 {
	   kndspow = pow(knds, enstn / ends);
	   npow = 1.0 - enstn / ends;
	 }
       else
	 {
	   kndspow = knds;
	   npow = 0.0;
	 }

       // Now loop over nodes, using at-a-station power law to set
       // width, depth & roughness
       for( cn = nIter.FirstP(); nIter.IsActive(); cn = nIter.NextP() )
	 {
	   //removed an if cn->Meanders(), so stuff calculated everywhere
	   //if rainfall varies, find hydraulic width "at-a-station"
	   //based on the channel width "downstream":
	   if( cn->getQ()>0.0 )
	     {
	       qpsec = cn->getQ()/SECPERYEAR;
	       width = pow(cn->getChanWidth(), widpow) * kwdspow
		 * pow(qpsec, ewstn);
	       cn->setHydrWidth( width );
	       depth = pow(cn->getChanDepth(), deppow) * kddspow
		 * pow(qpsec, edstn);
	       cn->setHydrDepth( depth );
	       rough = pow(cn->getChanRough(), npow) * kndspow
		 * pow(qpsec, enstn);
	       cn->setHydrRough( rough );
	       slope = cn->getChanSlope();
	       assert( slope >= 0. ); // slope can be 0 -- changed assert 3/99
	       //Depth now calculated as above - done to be consistent
	       //with changes made in FindChanGeom
	       //radfactor = qpsec * rough / width / sqrt(slope);
	       //hradius = pow(radfactor, 0.6);
	       //depth = hradius;
	       //depth = width / ( width / hradius - 2.0 );
	       //cn->setHydrDepth( depth );
	       cn->setHydrSlope( slope );
	     }
	 }
     }
   //if rainfall does not vary, set hydraulic geom. = channel geom.
   else
     {
       for( cn = nIter.FirstP(); nIter.IsActive(); cn = nIter.NextP() )
	 {
	   width = cn->getChanWidth();
	   rough = cn->getChanRough();
	   depth = cn->getChanDepth();
	   slope = cn->getChanSlope();
	   cn->setHydrWidth( width );
	   cn->setHydrRough( rough );
	   cn->setHydrSlope( slope );
	   cn->setHydrDepth( depth );
	 }
     }

   if (0) //DEBUG
     std::cout << "done tStreamNet::FindHydrGeom" << std::endl;
}


/*****************************************************************************\
**
**  tStreamNet::FindChanGeom
**
**  Calculates/assigns values for channel hydraulic geometry (width, depth,
**  roughness). Originally done only for meandering nodes; now done for
**  all nodes but only when invoked by the meander routines.
**
**		  Parameters: kwds, ewds, ewstn, knds, ends, enstn
**      Data members updated: tLNode->chan.hydrwidth, hydrnrough, hydrdepth
**      Called by: FindReaches (needs to know how long to make reach "tails"
**      Calls: no "major" functions
**
**      Created: SL 1/98
**      Modifications:
**       - 3/99: changed from calculation of bankfull event based on
**               recurrence interval to simply using an input parameter
**               that defines bankfull precip event magnitude. This was
**               done to prevent errors when the interstorm duration
**               exceeds 1.5 years. The variable _bankfullevent_ is the
**               magnitude of runoff which, when applied uniformly over
**               catchment and allowed to come to steady state, will
**               produce a bankfull event. GT
**       - 6/01 GT: added statement at beginning to ignore this fn if
**               alternative Parker-Paola models is used.
**       - 7/03: Parker-Paola model now called directly from here to set
**               bankfull geometry (GT)
**
\*****************************************************************************/
void tStreamNet::FindChanGeom()
{
   if( miChannelType==kParkerChannels )
     {
       assert( mpParkerChannels != 0 );
       mpParkerChannels->CalcChanGeom( meshPtr );
       return;
     }

   double qbf,      // Bankfull discharge in m3/s
       width,       // Channel width, m
       depth,       // Channel depth, m
       rough;       // Roughness
   double slope;
   double lambda;
   tLNode *cn;
   tMesh< tLNode >::nodeListIter_t nIter( meshPtr->getNodeList() );

   // the following modification made by gt, 3/99, to avoid hydraulic geom
   // errors during runs w/ long storms:
   //gt3/99 if (isdmn > 0 )  qbffactor = pmn * log(1.5 / isdmn);

   for( cn = nIter.FirstP(); nIter.IsActive(); cn = nIter.NextP() )
   {
      //took out an if cn->Meanders() so stuff will be calculated at all nodes
      //gt3/99 qbf = cn->getDrArea() * qbffactor;
      // Here we compute bankfull discharge and use it to compute width, depth,
      // etc. Note that if the user enters 0 for BANKFULLEVENT, the actual
      // current discharge will be used instead.
      qbf = cn->getDrArea()*bankfullevent;
      if( !qbf ) qbf = cn->getQ()/SECPERYEAR;  // q is in m^3/s
      width = kwds * pow(qbf, ewds);
      depth = kdds * pow(qbf, edds);
      rough = knds * pow(qbf, ends);
      lambda = klambda * pow(qbf, elambda);
      cn->setChanWidth( width );
      cn->setChanDepth( depth );
      cn->setChanRough( rough );
      cn->setBankRough( lambda );
      slope = cn->calcSlope();
      cn->setChanSlope( slope );
	  
	  if(0) //DEBUG
	  {
	     if( cn->getID()==8121 ) std::cout << "FindChanGeom 1: new slope at " << cn->getID() << " is " << slope << std::endl;
	  }
 
      //std::cout<<"FCG node "<<cn->getID()<<" new dep "<<depth;

      //Nic changed below, thinks it was causing problems
      //just using discharge relation to calculate depth instead
      //make sure slope will produce a positive depth:
      //double critS = qbf * qbf * rough * rough * 8.0 * pow( 2.0, 0.333 )/
      //( width * width * width * width * width * pow( width, 0.333 ) );
      //if( slope > critS ) //should also catch negative slope flag
      //{
      // std::cout << "in FindChanGeom, slope = " << slope << std::endl;
//            cn->setChanSlope( slope );
      //radfactor = qbf * rough / width / sqrt(slope);
      // hradius = pow(radfactor, 0.6);
         //depth = hradius;
      // depth = width / (width / hradius - 2.0);
//            cn->setChanDepth( depth );
      // std::cout<<" old dep "<<depth<<std::endl;
      //}
//         else cn->setMeanderStatus( kNonMeanderNode );
      //nmg
#define TESTFIX 0
      if( cn->Meanders() && slope <= 0.00000001 ){ // added "meanders" clause Gt 12/04
		 if(0) //debug
			std::cout << "Meander node de-activation in FindChanGeom, at node " << cn->getID() << " " << cn->getX() << " " << cn->getY() << std::endl;
         if( !TESTFIX ) cn->setMeanderStatus( kNonMeanderNode );
      }
      if( unlikely(slope < 0.0) )
      {
#if 1
         std::cout << "Negative slope = " << slope<<" at x,y "<< cn->getX() << ' ' << cn->getY()<< '\n';
         std::cout << "MeanderStatus= "<< cn->Meanders() << '\n';
         std::cout << " probably from infinite loop in tLNode::GetSlope()" << std::endl;
#else
         std::cout << "negative slope,"
              << " probably from infinite loop in tLNode::GetSlope()" << std::endl;
#endif
         ReportFatalError("negative slope in tStreamNet::FindChanGeom");
      }
   }
   if (0) //DEBUG
     std::cout << "done tStreamNet::FindChanGeom" << std::endl;
}


/**************************************************************************\
**
**  tStreamNet::RouteFlowKinWave
**
**  This routine was created to solve for 2D flow across a surface using
**  the kinematic wave approximation, rather than assuming channelized,
**  unidirectional (1D) flow. The algorithm obtains a steady state
**  solution under temporally constant runoff and discharge.
**    The algorithm begins by calling a routine to sort all active nodes
**  in upstream-to-downstream order, accounting for multiple downslope
**  flow pathways. This way, we always know the incoming flow from
**  upstream when we go to calculate depth and flow routing for each
**  node.
**    The velocity of downslope flow from a node i to one of its
**  neighbors j is computed from the general Chezy/Manning equation,
**  using the surface slope Sij between the two nodes and the water depth,
**  Yi, at node i:
**      Uij = (1/Kr) Yi^m Sij^(0.5)
**  To write this in terms of total discharge, Qij, instead of velocity,
**  we can multiply both sides by the depth, Yi, and by the width of
**  the Voronoi polygon face across which water is flowing, Wij:
**      Qij = (1/Kr) Yi^(m+1) Sij^(0.5) Wij
**  Rearranging, we can solve for Yi:
**      Yi = [ Kr Qij / (Sij^0.5 Wij) ] ^ (1/(m+1))
**  It remains now to find the discharge in a given direction, Qij.
**    To do so, we start by noting
**  that the total outgoing discharge (equal to discharge coming in
**  plus local runoff, both of which we know) is equal to the sum of
**  Qij over all Voronoi faces:
**      Qi = SUM[ (1/Kr) Yi^(m+1) Sij^(0.5) Wij ] for j=1 to M faces, or
**      Qi = (1/Kr) Yi^(m+1) SUM[ Sij^(0.5) Wij ]
**  From this, we see that the ratio Qij / Qi is given by
**      Qij/Qi = Sij^(0.5) Wij / SUM[ Sij^(0.5) Wij ], or
**      Qij = Qi * Sij^(0.5) Wij / SUM[ Sij^(0.5) Wij ]
**  Substituting into the expression for Yi above, we solve for
**      Yi = [ Kr Qi / SUM[ Sij^(0.5) Wij ] ] ^ (1/(m+1))
**  and we've nailed it!
**    In the algorithm below, we start by computing the sum term, then
**  find the depth, and then apportion discharge to each lower neighbor
**  using the expression above for Qij (Qij is added to the neighbor's
**  total discharge; since we're working upstream-to-downstream, we
**  know that the neighbor will have all its flow accounted for by the
**  time we come to calculate its flow depth and routing).
**    The method uses two parameters: the exponent 1/(m+1) (mdKinWaveExp),
**  and the roughness factor Kr (mdKinWaveRough), which is calculated
**  from the roughness parameter HYDR_ROUGH_COEF_DS (knds), but with
**  a conversion factor included to convert from seconds to years.
**
**  Created by GT, Jan 2000
**
\**************************************************************************/
void tStreamNet::RouteFlowKinWave( double rainrate_ )
{
   tLNode * cn;
   tEdge * ce;
   tMesh< tLNode >::nodeListIter_t niter( meshPtr->getNodeList() );
   double sum;                         // Sum used in to apportion flow
   double runoff = rainrate_ - infilt;  // Local runoff rate at node

   if (0) //DEBUG
     std::cout << "RouteFlowKinWave\n";

   if( runoff <= 0.0 ) return;

   // Reset discharge to zero everywhere
   for( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
       cn->setDischarge( 0.0 );

   // Sort nodes uphill-to-downhill
   SortNodesByNetOrder( true );

   // Route flow and compute water depths
   for( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
   {
      // Add local runoff to total incoming discharge
      cn->AddDischarge( runoff * cn->getVArea() );

      if( cn->getFloodStatus() == tLNode::kNotFlooded )
      {
         // Flow is apportioned among downhill neighbors according to
         // slope and Voronoi edge width, so first we perform the summation
         // of the product of Voronoi edge width and sqrt of slope in each dir
         //if( cn->getQ()>50000.0 ) cn->TellAll();
         sum = 0.0;
         ce = cn->getEdg();
         do
         {
            if( cn->getZ() > ce->getDestinationPtr()->getZ() &&
                ce->FlowAllowed() )
                sum += ce->getVEdgLen() * sqrt( ce->getSlope() );
            ce = ce->getCCWEdg();
         }
         while( ce!=cn->getEdg() );

         // Compute the flow depth
         //std::cout << "Q: " << cn->getQ() << " sum: " << sum << " DEPTH: " << cn->getHydrDepth() << std::endl;
         assert( cn->getQ()>0.0 );
         if( sum>0.0 )
             cn->setHydrDepth( pow( cn->getQ() * mdKinWaveRough / sum,
                                    mdKinWaveExp ) );
         else
             cn->setHydrDepth( 0.0 );

         // Route flow downhill
         if( sum>0.0 )
         {
            ce = cn->getEdg();
            do
            {
	       tLNode * dn = static_cast<tLNode *>(ce->getDestinationPtrNC());
               if( cn->getZ() > dn->getZ() && ce->FlowAllowed() )
                   dn->AddDischarge( cn->getQ()
                                     * (sqrt(ce->getSlope())*ce->getVEdgLen())/sum );
               ce = ce->getCCWEdg();
            }
            while( ce!=cn->getEdg() );
         }
      }
      else
          cn->setHydrDepth( 0.0 );

   } // end of for loop


} // End of tStreamNet::RouteFlowKinWave


/*****************************************************************************\
**
**  tStreamNet::DensifyMeshDrArea
**
**  Densifies the locally along channels/valleys by adding new nodes
**  wherever (a) drainage area exceeds a user-specified threshold
**  (stored in mdMeshAdaptMinArea) and (b) Voronoi area is larger than
**  a user specified maximum (stored in mdMeshAdaptMaxVArea). When a
**  node is found that meets these criteria, the tMesh function
**  AddNodesAround is called to add new nodes at the Voronoi vertices of
**  the node in question.
**
**		  Inputs: time -- current simulation time (if applicable; defaults
**                      to zero)
**      Other parameters: mdMeshAdaptMinArea, mdMeshAdaptMaxVArea
**      Called by: functions outside tStreamNet (eg, main)
**      Calls: tMesh::AddNodesAround
**      Created: GT 2/2000
**      Modifications:
**
\*****************************************************************************/
void tStreamNet::DensifyMeshDrArea( double time )
{
   tMesh< tLNode >::nodeListIter_t niter( meshPtr->getNodeList() );  // node list iter.
   tLNode *cn;              // Current node being checked

   // Check all active nodes
   for( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
   {
      // If drainage area and Voronoi area exceed thresholds, add new nodes
      if( cn->getDrArea()>mdMeshAdaptMinArea
          && cn->getVArea() > mdMeshAdaptMaxVArea )
      {
         meshPtr->AddNodesAround( cn, time );
      }
   }

}



/**************************************************************************\
**  FUNCTIONS FOR CLASS tInlet.
\**************************************************************************/

/**************************************************************************\
**
**  tInlet Constructors
**
**  (1) default sets all values to zero
**  (2) reads location and other variables from input file and
**      initializes (if applicable; the user may not want an inlet)
**
**  Modifications:
**    - Bug: if a new inlet node is added, using AddNodeAt won't assign
**  correct variables for layers, regolith, etc. Fix: new node is assigned
**  properties of the nearest neighbor (not elevation; that's done by
**  interpolation). GT 7/98
**    - 2/02 fixed bug: for single-size, inSedLoadm[0] wasn't being
**      initialized, even though it's used in DetachErode; also,
**      changed reading of inputs to use stl string class instead
**      of C-style character manipulation (GT)
**    - 7/03 the elevation of the inlet node is now found by interpolating
**      only the non-boundary nodes of the triangle the inlet is placed into.
**      ,and not interpolating using the nodes of the surrounding triangles
**      because this may lead to a resultant elevation of the inlet which
**      is lower than its direct surroundings, making it a pond. (QC)
**    - 5/06 added code to tInlet to handle option of having sediment
**      influx computed for each storm based on prescribed slope and
**      bed grain-size distribution at the inlet (GT)
**
\**************************************************************************/

tInlet::tInlet() :
  innode(0), inDrArea(0.), inSedLoad(0.), inSedLoadm(), meshPtr(0),
  optCalcSedFeed(false), inletSlope(0.), inletSedSizeFraction(0)
{}

#define LARGE_DISTANCE 1e9
tInlet::tInlet( tMesh< tLNode > *gPtr, const tInputFile &infile )
  :
  innode(0), inDrArea(0.),
  inSedLoad(0.),
  meshPtr(gPtr)
{
   bool inletbc = infile.ReadBool( "OPTINLET" );
   int i;
   int numg = infile.ReadItem( numg, "NUMGRNSIZE" );
   int add = 1;       // Use AddNode at for placing the inlet
   //char end, name[20];
   double xin, yin,   // Coords of inlet node
       mindist,       // Minimum distance above which new node will be added
       dist,          // Distance btwn inlet and a nearby node
       x, y,          // Location of a nearby node
       zin = 0,       // Elevation of inlet
       suminvdist = 0,// Sum of 1/dist for all nearby non-boundary nodes
       minDistFound;  // Smallest distance to nearby node found so far
   tTriangle *intri, *ntri;
   tLNode *cn,
       *closestNode(0);  // -> to closest nearby node found so far
   tPtrList< tLNode > nPL;            // List of nearby non-boundary nodes
   tPtrListIter< tLNode > itr( nPL ); // Iterator for the above list

   assert( meshPtr != 0 );
   if( inletbc )
   {
      // Read drainage area and sediment load at inlet. If more than one
      // grain size is simulated, read in a sediment load for each size
      // individually
     //std::string taglinebase = "INSEDLOAD";
     //std::string digits = "123456789";
     //std::string tagline;
     char tagline[12], lastdigit;
     strcpy( tagline, "INSEDLOAD0" );
     lastdigit = '0';
      inDrArea = infile.ReadItem( inDrArea, "INDRAREA" );
      
      // Stuff related to option for dynamic calculation of sediment influx at inlet
      optCalcSedFeed = infile.ReadItem( optCalcSedFeed, "INLET_OPTCALCSEDFEED" );
      if( optCalcSedFeed ) 
      {
        inletSlope = infile.ReadItem( inletSlope, "INLET_SLOPE" );
        inletSedSizeFraction.setSize(numg);
      }
      else
      {
        inletSlope = 0;
        inletSedSizeFraction.setSize(0);
      }
      
      // Read parameters related to sediment influx (either option)
      if(numg <= 1)
	{
          inSedLoad = infile.ReadItem( inSedLoad, "INSEDLOAD" );
	  inSedLoadm[0] = inSedLoad;
          if( optCalcSedFeed ) inletSedSizeFraction[0] = 1.0;
	}
      else{
         inSedLoadm.setSize(numg);
         inSedLoad=0.0;
         //i=1;
         //end='1';
         for( i=0; i<numg; i++ ) {
	   double help;
	   //strcpy( name, "INSEDLOAD");
           // strcat( name, &end );
           // help = infile.ReadItem( help, name);
	   //tagline = taglinebase + digits.substr( i, i );
	   // std::cout << tagline << std::endl;
	   // help = infile.ReadItem( help, tagline.c_str() );
            lastdigit++;
            tagline[9] = lastdigit;
            help = infile.ReadItem( help, tagline );
            inSedLoadm[i] = help;
            inSedLoad += help;
	    if (0) //DEBUG
	      std::cout<<"insedload of "<<i-1<<" is "<<inSedLoadm[i]<<std::endl;
            //i++;
            //end++;
            if( optCalcSedFeed )
                inletSedSizeFraction[i] = inSedLoadm[i];
         }
      }

      // Read in the location of the inlet node. If the specified coordinates
      // are "close" to an existing non-boundary node, assign that node as
      // the inlet; otherwise, create a new node. The elevation for the new
      // node is found by interpolation.
      xin = infile.ReadItem( xin, "INLET_X" );
      yin = infile.ReadItem( yin, "INLET_Y" );
      intri = meshPtr->LocateTriangle( xin, yin );
      assert( intri != 0 );  //TODO: should be error-check not assert
      for( i=0; i<3; i++ )
      {
       	 // Surrounding nodes of the Triangle
        cn = static_cast<tLNode *>(intri->pPtr(i));
        if( cn->getBoundaryFlag() == kNonBoundary ) nPL.insertAtBack( cn );

        // Nodes of the Triangles surrounding the one selected
        if(add==0){
	  ntri = intri->tPtr(i);
	  if( ntri != 0 )
	    {
	      cn = static_cast<tLNode *>(ntri->pPtr( ntri->nVOp( intri ) ));
	      if( cn->getBoundaryFlag() == kNonBoundary ) nPL.insertAtBack( cn );
	    }
	}
      }
      minDistFound = LARGE_DISTANCE;
      mindist = 0.000001;
      std::cout<< "x,y,z of nodes surrounding inlet are: \n";
      for( cn = itr.FirstP(); !(itr.AtEnd()); cn = itr.NextP() )
      {
         x = cn->getX();
         y = cn->getY();
         dist = sqrt( (xin - x) * (xin - x) + (yin - y) * (yin - y) );
         if( dist < minDistFound )
         {
            minDistFound = dist;
            closestNode = cn;
         }
         if( dist > mindist )          // find elev by interpolation
         {
            zin += cn->getZ() / dist;
            suminvdist += 1 / dist;
         }
         else                         // we're on a existing node, re-use
         {
            innode = cn;
            add = 0;
         }


         // DEBUG
         std::cout <<cn->getX()<<' '<<cn->getY()<<' '<<cn->getZ()
	      << ' ' << BoundName(cn->getBoundaryFlag())
	      <<" dist= "<<dist<<std::endl;
        /*
         if( dist < mindist )
         {
            dist = mindist;
            innode = cn;
         }*/
      }

      // Debug
      if( add == 0){
       std::cout<< "ADDING INLET by resetting an existing node:" << std::endl;
       std::cout <<"Innode boundary flag ="
	    << BoundName(innode->getBoundaryFlag()) << '\n';
      }

      if( add ) // fix here:
      {
	assert( closestNode != 0 );
 	 std::cout <<' '<<std::endl;
         std::cout << "ADDING INLET by calling AddNode():" << std::endl;
         std::cout <<' '<<std::endl;
         //closestNode->TellAll();
         //BE AWARE
         //The following commented line caused many problems for the code:
         //tLNode newnode( *closestNode );
         //This did not call the copy constructor which was created for
         //tLNode.  Most likely it called some default copy constructor
         //which could not handle copying of user created classes.
         //This created errors in the closestNode's layerlist.
         //The following calls the correct copy constructor.
	 // The comment above is not correct anymore. AD
         tLNode newnode( *closestNode );
         newnode.setZ( zin / suminvdist );
         newnode.setX( xin );
         newnode.setY( yin );
         //zin = zin / suminvdist;
	 //tArray< double > xyz(3);
         //xyz[0] = xin;
         //xyz[1] = yin;
         //xyz[2] = zin;
         //innode = meshPtr->AddNodeAt( xyz );
         innode = meshPtr->AddNode( newnode, kUpdateMesh ); // true means update mesh
         std::cout <<"Innode boundary flag ="
	      << BoundName(innode->getBoundaryFlag())<< '\n';
         std::cout <<"Innode x,y,z= "<<innode->getX()<<' '<<innode->getY()<<' '<<innode->getZ()<<std::endl;
         std::cout <<' '<<std::endl;
         //std::cout << "INLET NODE IS:\n";
         //innode->TellAll();
      }
   }
}

tInlet::~tInlet()
{
   innode = 0;
   meshPtr = 0;
}


/**************************************************************************\
**
**  tInlet::FindNewInlet
**
**  Search for points 'up-valley' of present inlet;
**  of those points and the present inlet, set new inlet to one with the
**  lowest elevation.
**
**  5/18/98 SL: Try something less arbitrary. Find active nodes
**
\**************************************************************************/
void tInlet::FindNewInlet()
{
   double yin, zmin, dmn, dmnn, dmin;
   tLNode *cn, *newinnode, *mn;
   tNode *bn0(0), *bn1(0), *mnn;
   tEdge *ce, *me;
   tMesh< tLNode >::nodeListIter_t nI( meshPtr->getNodeList() );
   int n;
     //tPtrList< tLNode > bList;
     //tPtrListIter< tLNode > bI( bList );
   const tArray< double > xyz( innode->get3DCoords() );
   yin = xyz[1];
   zmin = xyz[2];
   newinnode = innode;
     //for( cn = nI.FirstP(); nI.IsActive(); cn = nI.NextP() )
     //go through boundary nodes:
   /*cn =*/ nI.LastActiveP();
   for( cn = nI.NextP(); !(nI.AtEnd()); cn = nI.NextP() )
   {
        //select for 'northern' bndy nodes:
      if( cn->getY() > yin ) //(cn was originally any active node)
      {
           //go through bndy node's nbrs...
    	 tSpkIter sI( cn );
         for( ce = sI.FirstP(); !(sI.AtEnd()); ce = sI.NextP() )
         {
            mn = static_cast<tLNode *>(ce->getDestinationPtrNC());
              //easier to check node's elevation at this point to make
              //sure it's worth going further with all this logic;
              //so, find an active neighbor with elevation lower than
              //present inlet:
            // DEBUG: try to catch floating excep in if below:
            /*assert( mn->getZ()>0.0 );
            assert( mn->getZ()<10000.0 );
            assert( zmin>0.0 );
            assert( zmin<10000.0 );*/
            if( mn->getBoundaryFlag() == kNonBoundary && mn->getZ() < zmin )
            {
	       tSpkIter msI( mn );
               n = 0;
                 //go through the active node's nbrs to find and count
                 //'northern' bndy nodes:
               for(  me = msI.FirstP(); !(msI.AtEnd()); me = msI.NextP() )
               {
                  mnn = me->getDestinationPtrNC();
                  if( mnn->getBoundaryFlag() != kNonBoundary && mnn->getY() > yin )
                  {
                     if( n == 0 ) bn0 = mnn;
                     else if( n == 1 ) bn1 = mnn;
                     n++;
                  }
               }
                 //if active node has >1 northern bndy nbr...
               if( n > 1 )
               {
                    //find node's distance to boundary:
		  assert( bn0 != 0 );
		  assert( bn1 != 0 );
                  dmn = DistanceToLine( mn->getX(), mn->getY(), bn0, bn1 );
                  dmin = dmn;
                    //find it's active nbrs' distances:
                  for(  me = msI.FirstP(); !(msI.AtEnd()); me = msI.NextP() )
                  {
                     mnn = static_cast<tLNode *>(me->getDestinationPtrNC());
                     if( mnn->getBoundaryFlag() == kNonBoundary )
                     {
                        dmnn = DistanceToLine( mnn->getX(), mnn->getY(),
                                                    bn0, bn1 );
                        if( dmnn < dmin ) dmin = dmnn;
                     }
                  }
                    //if none of active nbrs' distances smaller,
                    //node is the new inlet if no better candidates are found:
                  if( dmin == dmn )
                  {
                     zmin = mn->getZ();
                     newinnode = mn;
                  }
               }
            }
         }
           /*if( cn->getZ() < zmin )
         {
            zmin = cn->getZ();
            newinnode = cn;
         }*/
      }
   }
     //finally reset inlet node; if no new inlet was found, we're just
     //setting it equal to itself:
   if( innode!=newinnode )
       std::cout << "*** MOVING INLET from " << innode->getX()
            << "," << innode->getY() << " to "
            << newinnode->getX() << ","
            << newinnode->getY() << std::endl;
   innode = newinnode;
}


/**************************************************************************\
**
**  tInlet "get" and "set" functions
**
\**************************************************************************/
double tInlet::getInSedLoad() const {return inSedLoad;}

double tInlet::getInSedLoad( size_t i )
{
   if(i>=inSedLoadm.getSize())
        ReportFatalError( "Trying to set size in sediment load that doesn't exist");
   return inSedLoadm[i];
}

tArray< double >
tInlet::getInSedLoadm( ) const
{
   return inSedLoadm;
}

void tInlet::setInSedLoad( double val )
{inSedLoad = ( val > 0.0 ) ? val : 0.0;}

void tInlet::setInSedLoad( size_t i, double val )
{
   if(i>=inSedLoadm.getSize())
        ReportFatalError( "Trying to set size in sediment load that doesn't exist");
   inSedLoadm[i]=val;
}

double tInlet::getInDrArea() const {return inDrArea;}
void tInlet::setInDrArea( double val ) {inDrArea = ( val > 0.0 ) ? val : 0.0;}
tLNode *tInlet::getInNodePtr() {return innode;}
void tInlet::setInNodePtr( tLNode *ptr ) {innode = ptr;}


/**************************************************************************\
**  FUNCTIONS FOR CLASS tParkerChannels.
\**************************************************************************/

/**************************************************************************\
**
**  tParkerChannel Constructor
**
**  The constructor reads in necessary parameters to implement the
**  Parker-Paola self-formed channel geometry model. A brief derivation
**  is below:
**
**    hypothesis: taub / taucrit = P = 1.2 to 1.4
**
**  where taub is bankfull bed shear stress and taucrit is critical
**  shear stress for motion of d50 bed material.
**
**    steady-uniform, wide-channel flow:  taub = kt (Qb/W)^alpha S^beta
**
**  Combining gives
**
**    W = Qb [ (kt S^beta) / (taucrit P) ] ^ alpha
**
**  From Shields' relation,
**
**    taucrit = thetac ( sigma - rho ) g d50   (see below for defs)
**    thetac = constant for fully turbulent flow ~ 0.045
**      (see Slingerland, Harbaugh, and Furlong, 1994)
**
**  Here, width is recalculated for any given (not necessarily bankfull)
**  discharge. Multiplicative constant and exponent are then:
**
**    W = Q * uconv * mdPPfac * S^(beta/alpha)
**
**  "uconv" is unit conversion from dishcarge in m3/yr to m3/s, consistent
**  with the use of SI units in kt and taucrit.
**
**      Created 6/01, GT
**
**    Modifications:
**
\**************************************************************************/
tParkerChannels::tParkerChannels( const tInputFile &infile )
{
  double kt,              // Shear-stress coefficient (SI units)
    taucrit,              // Critical shear stress for d50 (Pa)
    thetac,               // Critical shields stress
    shearRatio,           // Parker-Paola constant (tau/taucrit ratio)
    d50,                  // Median grain diameter (m)
    alpha,                // Specific-discharge exponent in shear stress eqn
    beta;                 // Slope exponent in shear stress eqn
  int i;                    // Counter
  const double
    sigma = RHOSED,       // Sediment density
    rho = RHO,            // Water density
    grav = GRAV,          // Gravitational acceleration
    secPerYear=SECPERYEAR;  // # of seconds in one year
  char astring[12];   // string var used in reading grain-size classes

   if (0) //DEBUG
     std::cout << "tParkerChannels::tParkerChannels\n";

   // Calculate coefficient and slope exponent for width equation (see above)
   // If more than one grain size is used, then grain sizes can vary
   // dynamically and we need to read and store the median diameter of each
   // size-class.
   kt = infile.ReadItem( kt, "KT" );
   alpha = infile.ReadItem( alpha, "MF" );
   beta = infile.ReadItem( beta, "NF" );
   thetac = infile.ReadItem( thetac, "THETAC" );
   shearRatio = infile.ReadItem( shearRatio, "SHEAR_RATIO" );
   miNumGrainSizeClasses = infile.ReadItem( miNumGrainSizeClasses,
					    "NUMGRNSIZE" );
   if(  miNumGrainSizeClasses < 1 ||  miNumGrainSizeClasses>=10 )
     ReportFatalError( "Number of grain-size classes must be between 1 and 9 inclusive." );

   mdPPexp1 = beta / alpha;
   if( miNumGrainSizeClasses>1 ) // In this case we will calc D50 on the fly
     {
       mdPPfac = ( 1.0 / secPerYear )
	 * pow( kt / (  thetac*(sigma-rho)*grav*shearRatio ), 1.0 / alpha );
       mdPPexp2 = -1.0 / alpha;
       mD50BySizeClass.setSize( miNumGrainSizeClasses );
       strcpy( astring, "GRAINDIAM0" );
       for( i=0; i<miNumGrainSizeClasses; i++ )
	 {
	   astring[9]++;  // Next number (1, 2, etc.)
	   //strcpy( astring, "GRAINDIAM" );
	   //strcat( astring, nstr );
	   mD50BySizeClass[i] = infile.ReadItem( d50, astring );
	 }
     }
   else // In this case we incorporate constant d50 into the parameters
     {
       d50 = infile.ReadItem( d50, "GRAINDIAM0" );
       taucrit = thetac*(sigma-rho)*grav*d50;
       if (0) //DEBUG
	 std::cout << "Tau crit = " << taucrit << std::endl;
       mdPPfac = ( 1.0 / secPerYear )
	 * pow( kt / ( taucrit*shearRatio ), 1.0 / alpha );
       mdPPexp2 = 0.0; // Not used in this case
       if (0) //DEBUG
	 std::cout << "mdPPfac=" << mdPPfac << "  mdPPexp1=" << mdPPexp1 << std::endl;
     }

  // Depth could be calculated from width plus Manning equation
   // (however, at the moment this is disabled)
  mdRough = infile.ReadItem( mdRough, "HYDR_ROUGH_COEFF_DS" ); // Manning's n
  mdRough = mdRough / secPerYear;
  mdDepthexp = 0.6;  // From Manning eqn; for Chezy / Darcy would be 2/3

}

/**************************************************************************\
**
**  tParkerChannels::CalcChanGeom
**
**  This function updates channel width and depth at each node according
**  to the Parker-Paola self-formed channel hypothesis. See the constructor
**  function for derivation.
**
**  Note: channels "self-form" for every discharge value in this
**  implementation; this is not an ideal approximation.
**
**  Note: chanwidth and chandepth are not set. Side effects unknown.
**
**      Created 6/01, GT
**
**      Modifications:
**         - correct potential error in depth-setting via division by
**           zero; set chan and hydr depth; 02/02 GT
**
\**************************************************************************/
void tParkerChannels::CalcChanGeom( tMesh<tLNode> *meshPtr )
{
  tMesh< tLNode >::nodeListIter_t ni( meshPtr->getNodeList() );
  tLNode *cn;

  if (0) //DEBUG
    std::cout << "tParkerChannels::CalcChanGeom\n";

  if( miNumGrainSizeClasses==1 )
    for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
      {
	cn->setChanWidth( mdPPfac * cn->getQ() * pow(cn->calcSlope(),mdPPexp1 ) );
	/* double denom;
	   if( ( denom = cn->getHydrWidth() * sqrt( cn->calcSlope() ) ) > 0.0 )
	   cn->setHydrDepth( pow( ( cn->getQ() * mdRough ) / denom,
	   mdDepthexp ) );
	   else
	   cn->setHydrDepth( 0.0 );*/
	cn->setHydrDepth( 1. );
	cn->setChanDepth( cn->getHydrDepth() );
      }
  else
    for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
      {
	// Calc mean (NOT median) grain size
	double d50 = 0.0;
	for( int i=0; i<miNumGrainSizeClasses; i++ )
	  {
	    d50 += cn->getLayerDgrade(0,i) * mD50BySizeClass[i];
	  }
	assert( cn->getLayerDepth(0)>0. );
	d50 = d50 / cn->getLayerDepth(0);
	assert( d50>0. );
	cn->setChanWidth( mdPPfac * cn->getQ() * pow(cn->calcSlope(),mdPPexp1 )
			  * pow( d50, mdPPexp2 ) );
	if(0) { // debug
	  std::cout << mdPPfac << " " << cn->getQ() << " " << cn->calcSlope()
	     << " " << mdPPexp1 << " " << d50 << " " << mdPPexp2 << std::endl;
	}
	cn->setHydrDepth( 1. );
	cn->setChanDepth( cn->getHydrDepth() );
	if( cn->getChanWidth()==0. && cn->getFloodStatus() == tLNode::kNotFlooded )
	  cn->TellAll();
      }

}




void tStreamNet::DebugShowNbrs( tLNode *theNode ) const
{
  if(theNode->getX() < 12000 && theNode->getX() > 7000){	// stay within boundary
    std::cout<<"Look around node number: "
	<< theNode->getID() << " "
	<< theNode->getX() << " " << theNode->getY() <<" "<<theNode->getZ()
	<<" M: "<< theNode->Meanders()
	<<" F: "<< tLNode::FloodName(theNode->getFloodStatus()) <<std::endl;
    tLNode *nextNode = theNode->getDownstrmNbr();
    if(nextNode != NULL){
      std::cout<<"with dstr nbr:		 "
	  << nextNode->getID() << " "
	  << nextNode->getX() << " " << nextNode->getY() <<" "<<nextNode->getZ()
	  <<" M: "<< nextNode->Meanders()
	  <<" F: "<< tLNode::FloodName(nextNode->getFloodStatus()) <<std::endl;
    }
    std::cout << "flowedge flippable ? " << (theNode->getFlowEdg()->isFlippable()
					? "yes": "no") << std::endl;

    std::cout<<"other nbrs are:          "<<std::endl;

    tEdge *ce;
    tSpkIter spokIter( theNode );

    int counter=0;
    for( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() )
      {
	counter++;
	tLNode* surnode = static_cast<tLNode *>( ce->getDestinationPtrNC() );
	std::cout<<"    "<<std::endl;
	std::cout<<"neighbour "<<counter<< ":"<<std::endl;
	std::cout<<surnode->getID()<<" "
	    <<surnode->getX()<<" "<<surnode->getY()<<" "<<surnode->getZ()
	    <<" M: " << surnode->Meanders()
	    <<" F: " << tLNode::FloodName(surnode->getFloodStatus())
	    <<" ce-slope: "<<ce->getSlope()<<std::endl;
	std::cout<<"has nbrs itself: "<<std::endl;
	tEdge *cee;
	tSpkIter spokIter2(surnode);
	if(0){
	  for( cee = spokIter2.FirstP(); !( spokIter2.AtEnd() ); cee = spokIter2.NextP() ){
	    tLNode* morenode = static_cast<tLNode *>( cee->getDestinationPtrNC() );
	    std::cout<<morenode->getID()<<" "
		<<morenode->getX()<<" "<<morenode->getY()<<" "<<morenode->getZ()
		<<" M: "<<morenode->Meanders()
		<<" F: "<<tLNode::FloodName(morenode->getFloodStatus())
		<<" ce-slope: "<<cee->getSlope()<<std::endl;
	  }
	}

      }
    std::cout<<"   "<<std::endl;
  }
}
