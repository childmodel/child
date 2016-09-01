/**************************************************************************/
/**
**  @file tVegetation.cpp
**  @brief Functions for tVegetation and tVegCover classes
**
**  Classes tVegetation and tVegCover represents the vegetation cover
**  across a terrain. Class tVegCover represents the properties of
**  vegetation at a point on the landscape (e.g., the %cover), while
**  class tVegetation represents the landscape-wide ("global")
**  properties (e.g., regrowth coefficient). Class tVegCover is designed
**  to be embedded in a node (in other words, each node "has a"
**  vegetation cover).
**
**  Created January, 2000, GT
**  State read back from file if needed - November 2003, AD
**
**  $Id: tVegetation.cpp,v 1.16 2004-05-10 10:52:52 childcvs Exp $
*/
/**************************************************************************/

#include <assert.h>
#include <algorithm>  // replaces max and min macros
  using std::min;  // ditto

#include "tVegetation.h"
#include "../tMesh/tMesh.h"
#include "../globalFns.h"
#include "../tRunTimer/tRunTimer.h"
#include "../tStorm/tStorm.h"

/*
**  Functions for tFire objects.
**  A tFire object generates random fires assuming an exponential
**  distribution of time to the next fire. 
**  Inspired by a conversation with Mike Wimberly, 4/2000, it assumes
**  that, when a fire occurs, everything is "burned" by simply killing
**  the vegetation.
**  Its services include
**  reading the necessary parameters from a tInputFile, generating a new      
**  fire, and reporting its various values.
**    The GammaDev() function is provided for future reference; it is not
**  actually used in version 1.0.
**    tFire objects could be easily modified (or inherited from) to use
**  different distributions. 
**
**  Stephen Lancaster, April, 2000. Added to CHILD (CU) 8/2010
*/

/*
**  tFire::tFire:  Constructor for fires. The default constructor
**                   assigns a value of unity to interfire duration. (Note:
**                   this constructor does not allow option for sinusoidal
**                   variation in means).
*/
tFire::tFire()
  : rand(0)
{
   ifrdurMean = 1.0;
   ifrdur = 1.0;
}

tFire::tFire( const tFire& orig )
  : optRandom(orig.optRandom), ifrdurMean(orig.ifrdurMean), 
    ifrdur(orig.ifrdur), time_to_burn(orig.time_to_burn) 
{
  rand = new tRand(*orig.rand);
}
/*
**  tFire::tFire:  The constructor that's really used assigns values for
**                   mean interfire duration,
**                   initializes current duration to the mean values,
**                   and initializes the random number generator. (Note:
**                   this constructor does not allow option for sinusoidal
**                   variation in means).
*/
tFire::tFire( double mif, unsigned sd, tRunTimer* ptr )
  : rand(0)
{
   ifrdurMean = mif;
   ifrdur = ifrdurMean;
   long seed = sd;
   rand = new tRand( seed );
   timePtr = ptr;
   ifrdur = ifrdurMean * rand->ExpDev();
   time_to_burn = timePtr->getCurrentTime() + ifrdur;
}


/*
**  tFire::tFire
**
**  Alternative constructor that reads parameters directly from a tInputFile
**  object. Reads mean value for interfire period, and a random seed to 
**  initialize the random number generator.
**  
*/
tFire::tFire( const tInputFile &infile, tRunTimer* tPtr
	      , bool no_write_mode /*=false*/ )
  : rand(0)
{
   ifrdurMean = infile.ReadItem( ifrdurMean, "IFRDUR" );
   ifrdur = ifrdurMean;

   //srand( seed );
   // I/O stuff:
   timePtr = tPtr;
   optRandom = infile.ReadBool( "OPTRANDOMFIRES", false );
   if( optRandom )
     {
       // Read and initialize seed for random number generation
       long seed = infile.ReadItem( seed, "FSEED" );
       rand = new tRand( seed );
       ifrdur = ifrdurMean * rand->ExpDev();
       // If random fires used, create a file for writing them
       if( !no_write_mode )
	 {
	   char fname[87];
#define THEEXT ".fir"
	   infile.ReadItem( fname, sizeof(fname)-sizeof(THEEXT), "OUTFILENAME" );
	   strcat( fname, THEEXT );
#undef THEEXT
	   firefs.open( fname );
	   if( !firefs.good() )
	     std::cerr << "Warning: unable to create fire data file '"
		       << fname << "'\n";
	 }
     }
   time_to_burn = timePtr->getCurrentTime() + ifrdur;
}

/**************************************************************************\
**
**  tFire::TurnOnOutput, TurnOffOutput
**
**  Open output file so output will be directed to it, 
**  or close output file so there won't be output.
**
**  SL, 10/2010
**
\**************************************************************************/
void tFire::TurnOnOutput( const tInputFile& infile )
{
     // If variable storms used, create a file for writing them
  if( !firefs.good() && optRandom )
   {
     char fname[87];
#define THEEXT ".fir"
     infile.ReadItem( fname, sizeof(fname)-sizeof(THEEXT), "OUTFILENAME" );
     strcat( fname, THEEXT );
#undef THEEXT
     firefs.open( fname );
     if( !firefs.good() )
       std::cerr << "Warning: unable to create fire data file '"
		 << fname << "'\n";
   }
}

void tFire::TurnOffOutput()
{
  if( firefs.good() )
    firefs.close();
}


/*
**  tFire::Burn_WhenItsTime()
**
**  Generates a new fire by drawing new value of ifrdur from
**  an exponential distribution.
**  Keeps track of time till next burn, and burns when it reaches zero.
**  Upon burning, finds new ifrdur and resets time_to_burn.
**
**  Members updated:  ifrdur, time_to_burn
**                    
*/
bool tFire::Burn_WhenItsTime()
{
   if( timePtr->getCurrentTime() >= time_to_burn )
   {
      // find new time for fire:
     if( optRandom )
       {
	 ifrdur = ifrdurMean * rand->ExpDev();
	 if( firefs.good() )
	   firefs << time_to_burn << std::endl;
       }
      time_to_burn += ifrdur;
      // burn now:
      return true;
   }
   return false;
}





/************************************************************************\
**
**  tForest: data and functions to enable basic forest evolution
**
**  Created: 11/98 SL; Added to CHILD (CU) 8/10
**
\************************************************************************/
tForest::tForest( const tForest& orig )
  : mesh(0), storm(0), rand(0),
    rootdecayK(orig.rootdecayK), rootdecayN(orig.rootdecayN), 
    rootgrowthA(orig.rootgrowthA), rootgrowthB(orig.rootgrowthB), 
    rootgrowthC(orig.rootgrowthC), rootgrowthF(orig.rootgrowthF), 
    rootstrengthJ(orig.rootstrengthJ), MVRC(orig.MVRC), MLRC(orig.MLRC), 
    RSPartition(orig.RSPartition), heightindex(orig.heightindex), 
    weightMax(orig.weightMax), weightA(orig.weightA), weightB(orig.weightB), 
    weightC(orig.weightC), weightK(orig.weightK), wooddensity(orig.wooddensity), 
    blowparam(orig.blowparam), diamB0(orig.diamB0), diamB1(orig.diamB1), 
    diamB2(orig.diamB2), wooddecayK(orig.wooddecayK)
{
  if( orig.rand )
    rand = new tRand( *orig.rand );  
}

void tForest::setMeshPtr( tMesh<tLNode>* ptr ) 
{
  mesh = ptr;
  tMesh<tLNode>::nodeListIter_t nIter( mesh->getNodeList() );
  tLNode *cn=0;
  for( cn = nIter.FirstP(); nIter.IsActive(); cn = nIter.NextP() )
    {
      cn->getVegCover().getTrees()->setForestPtr( this );
      cn->getVegCover().getTrees()->setNodePtr( cn );
    }
}

tForest::tForest( const tInputFile& infile, tMesh<tLNode>* mPtr, tStorm *sPtr )
{
  ForestInitialize( infile, mPtr, sPtr );
}

/***************************************************************************\
**
**  ForestInitialize(  ): function to find and set rootstrength
**    according to Benda & Dunne 97a, Sidle 92, Sidle 91
**
\***************************************************************************/
void tForest::ForestInitialize( const tInputFile& infile, tMesh<tLNode>* mPtr, 
				tStorm *sPtr )
{
  if( &infile != 0 )
    {
      rootdecayK = infile.ReadItem( rootdecayK, "ROOTDECAY_K" );
      rootdecayN = infile.ReadItem( rootdecayN, "ROOTDECAY_N" );
      rootgrowthA = infile.ReadItem( rootgrowthA, "ROOTGROWTH_A" );
      rootgrowthB = infile.ReadItem( rootgrowthB, "ROOTGROWTH_B" );
      rootgrowthC = infile.ReadItem( rootgrowthC, "ROOTGROWTH_C" );
      rootgrowthF = infile.ReadItem( rootgrowthF, "ROOTGROWTH_F" );
      rootstrengthJ = infile.ReadItem( rootstrengthJ, "ROOTSTRENGTH_J" );
      MVRC = infile.ReadItem( MVRC, "MAXVERTROOTCOHESION" );
      MLRC = infile.ReadItem( MLRC, "MAXLATROOTCOHESION" );
      RSPartition = MLRC / MVRC;
      heightindex = infile.ReadItem( heightindex, "TREEHEIGHTINDEX" );
      weightMax = infile.ReadItem( weightMax, "VEGWEIGHT_MAX" );
      weightA = infile.ReadItem( weightA, "VEGWEIGHT_A" );
      weightB = infile.ReadItem( weightB, "VEGWEIGHT_B" );
      weightC = infile.ReadItem( weightC, "VEGWEIGHT_C" );
      weightK = infile.ReadItem( weightK, "VEGWEIGHT_K" );
      wooddensity = infile.ReadItem( wooddensity, "WOODDENSITY" );
      blowparam = infile.ReadItem( blowparam, "BLOWDOWNPARAM" );
      rand = 0;
      if( blowparam > 0.0 )
	{
	  long seed = infile.ReadItem( seed, "BLOW_SEED" );
	  rand = new tRand( seed );
	}
      diamB0 = infile.ReadItem( diamB0, "TREEDIAM_B0" );
      diamB1 = infile.ReadItem( diamB1, "TREEDIAM_B1" );
      diamB2 = infile.ReadItem( diamB2, "TREEDIAM_B2" );
      wooddecayK = infile.ReadItem( wooddecayK, "WOODDECAY_K" );
    }
  mesh = mPtr;
  storm = sPtr;
  // initialize nodes
  int opt;
  if ( (opt = infile.ReadItem( opt, "OPTREADINPUT" ))
       == OPTREADINPUT_PREVIOUS) 
    {
      // Start from a previous computation
      tListInputDataForest inputForestData( infile );
      const int nnodes = mesh->getNodeList()->getSize();
      if( inputForestData.rootstrength.getSize() != static_cast<size_t>(nnodes) )
	ReportFatalError( "tForest(): invalid number of records"
			  " in input file." );
      // Rely on the fact that the IDs have not been re-numbered.
      // for fast lookup per ID
      const tMesh< tLNode >::tIdArrayNode_t NodeTable(*(mesh->getNodeList()));
      tTrees tmp;
      tmp.setForestPtr( this );
      for( int id=0; id < nnodes; ++id )
	{
	  tmp.setNodePtr( NodeTable[id] );
	  tmp.setStandDeathTime( inputForestData.standdeathtime[id] );
	  tmp.setMaxRootStrength( inputForestData.maxrootstrength[id] );
	  // root strength is complicated; use RootStrengthInit after
	  // standdeathtime and maxrootstrength are set:
	  tmp.setRootStrength( RootStrengthInit( &tmp ) );	  
// 	  tmp.setRootStrength( inputForestData.rootstrength[id] );
	  tmp.setMaxHeightStand( inputForestData.maxheightstand[id] );
	  tmp.setBioMassStand( inputForestData.biomassstand[id] );
	  tmp.setBioMassDown( inputForestData.biomassdown[id] );
	  tmp.setMaxDiamStand( Tree_diam_from_height( &tmp ) );
	  tTrees *tPtr = new tTrees( tmp );
	  NodeTable[id]->getVegCover().setTrees( tPtr );
	}
    }
  else
    {
      // initialize from stand age input:
      double standage = infile.ReadItem( standage, "INITSTANDAGE" );
      tMesh<tLNode>::nodeListIter_t nIter( mesh->getNodeList() );
      tLNode *cn=0;
      for( cn = nIter.FirstP(); nIter.IsActive(); cn = nIter.NextP() )
	{
	  // this constructor initializes data based on stand age 
	  tTrees *ptr = new tTrees( cn, this, standage );
	  cn->getVegCover().setTrees( ptr );
	}
      // need to initialize non-active nodes for obscure reasons:
      for( ; !nIter.AtEnd(); cn = nIter.NextP() )
	{
	  tTrees *ptr = new tTrees();
	  cn->getVegCover().setTrees( ptr );
	}
    }
}

// invert function for determination of height from diameter at breast height:
double tForest::Tree_diam_from_height( tTrees *tree )
{
   if( tree->getMaxHeightStand() <= 1.37 ) return 0.0;
   const double logquant = 
     1.0 - pow( ( tree->getMaxHeightStand() - 1.37 ) / diamB0, 1.0 / diamB2 );
   assert( logquant > 0.0 );
   // divide by 100 to convert cm to m:
   return 
     ( 1.0 + 1.37 / ( tree->getMaxHeightStand() - 1.37 ) ) * log( logquant ) 
     / diamB1 / 100.0;
}

// use to initialize tTrees::maxheightstand
double tForest::TreeHeightInit( tTrees *tree )
{
  if( tree->getStandDeathTime() < 0.0 )
    return heightindex;
  if( tree->getStandDeathTime() == 0.0 ) return 0.0;
  return 
    Richards_Chapman_equ( tree->getStandDeathTime(), 
			  37.57 + 0.71698 * heightindex,
			  0.00055019 * heightindex,
			  0.95516 + 0.0072776 * heightindex );
}

// use to initialize tTrees::biomassstand
double tForest::WoodVolumeInit( tTrees *tree )
{
   if( tree->getStandDeathTime() < 0.0 )
      return weightMax / wooddensity / GRAV;
   if( tree->getStandDeathTime() == 0.0 ) return 0.0;
   return
     weightMax / wooddensity / GRAV
     * min( 1.0,
	    1.0 / 
	    ( weightA + 
	      weightB * exp( -weightK * tree->getStandDeathTime() ) ) 
	    + weightC );
}

// use to initialize tTrees::rootstrength; sets rootgrowth and rootdecay
double tForest::RootStrengthInit( tTrees *tree )
{
  if( tree->getStandDeathTime() < 0.0 )
    {
      tree->setRootGrowth( 1.0 );
      tree->setRootDecay( 0.0 );
    }
  else
    {
      tree->setRootGrowth( 1.0 /
			   ( rootgrowthA + 
			     rootgrowthB * 
			     exp( -rootgrowthF * tree->getStandDeathTime() ) )
			   + rootgrowthC );
      tree->setRootDecay( exp( -rootdecayK * 
			       pow( tree->getStandDeathTime(), rootdecayN ) ) );
    }
  // if haven't already set maxrootstrength, set it to maximum; 
  // negative if not set:
  if( tree->getMaxRootStrength() < 0.0 )
    tree->setMaxRootStrength( MVRC + MLRC );
  // find new root strength with new soil depth:
//   const double h = tree->getNodePtr()->getRegolithDepth();
//   const double costheta = cos( atan( tree->getNodePtr()->calcSlope() ) );
//   assert( h >= 0.0 );
//   const double expnegjZ = exp( -rootstrengthJ * h * costheta );
  const double Cv0 = tree->getMaxRootStrength() / ( 1.0 + RSPartition );
  const double Cl0 = RSPartition * Cv0;
  tree->setRootStrengthLat( tree->getRootGrowth() * MLRC + 
			    tree->getRootDecay() * Cl0 );
  tree->setRootStrengthVert( tree->getRootGrowth() * MVRC +
			     tree->getRootDecay() * Cv0 );
//   tree->setRootStrengthVert( ( tree->getRootGrowth() * MVRC +
// 			       tree->getRootDecay() * Cv0 ) * expnegjZ );
  return tree->getRootStrengthLat() + tree->getRootStrengthVert();
}


// returns the height to add to maxheightstand
double tForest::TreeHeightEvol( tTrees *tree )
{
  double height = 
    Richards_Chapman_equ( tree->getStandDeathTime(), 
			  37.57 + 0.71698 * heightindex,
			  0.00055019 * heightindex,
			  0.95516 + 0.0072776 * heightindex );
  if( height > diamB0 ) height = diamB0;
  return height;
}

// returns the biomass to add to biomassstand
double tForest::WoodVolumeEvol( tTrees *tree )
{
   return weightMax / wooddensity / GRAV
     * min( 1.0,
	    1.0 / 
	    ( weightA + 
	      weightB * exp( -weightK * tree->getStandDeathTime() ) ) 
	    + weightC );
}

// returns the new root strength
double tForest::RootStrengthEvol( tTrees *tree )
{
  tree->setRootGrowth( 1.0 /
		       ( rootgrowthA + 
			 rootgrowthB * 
			 exp( -rootgrowthF * tree->getStandDeathTime() ) )
		       + rootgrowthC );
  tree->setRootDecay( exp( -rootdecayK * 
			   pow( tree->getStandDeathTime(), rootdecayN ) ) );
  // find new root strength with new soil depth:
//   const double h = tree->getNodePtr()->getRegolithDepth();
//   const double costheta = cos( atan( tree->getNodePtr()->calcSlope() ) );
//   assert( h >= 0.0 );
//   const double expnegjZ = exp( -rootstrengthJ * h * costheta );
  const double Cv0 = tree->getMaxRootStrength() / ( 1.0 + RSPartition );
  const double Cl0 = RSPartition * Cv0;
  tree->setRootStrengthLat( tree->getRootGrowth() * MLRC + 
			    tree->getRootDecay() * Cl0 );
  tree->setRootStrengthVert( tree->getRootGrowth() * MVRC +
			     tree->getRootDecay() * Cv0 );
//   tree->setRootStrengthVert( ( tree->getRootGrowth() * MVRC +
// 			       tree->getRootDecay() * Cv0 ) * expnegjZ );
  return tree->getRootStrengthLat() + tree->getRootStrengthVert();
}

// set root decay constants such that root strength decays from present value
double tForest::RootDeath( tTrees *tree ) 
{
  // finds current root strength for zero soil depth to initialize decay
  return 
    ( MVRC + MLRC ) * tree->getRootGrowth() 
    + tree->getMaxRootStrength() * tree->getRootDecay();
}

// Next 3 functions are to do stochastic tree mortality via "blowdown"
// related to rainfall intensity. These functions are the reason tForest
// needs a pointer to a tStorm object. That said, blowdown has not been
// fully implemented in the CU version of CHILD, partly because the simulations
// envisioned at this point will use a discretization of ~1 m, which is quite
// a bit less than the expected spacing between the trees originally 
// simulated, i.e., Douglas-fir.
// -SL, 8/10

// exponentially distributed number of trees blown down at a node,
// dependent on ratio of wind drag force to roots' resisting force;
// wind speed assumed linearly proportional to rain rate.
#define PiValue 3.14159265358979323846
int tForest::NumBlowDown( tTrees *tree )
{
  if( rand == 0 ) return 0;
  const double diameter = tree->getMaxDiamStand();
  if( diameter == 0.0 ) return 0;
  const double height = tree->getMaxHeightStand();
   // volume of tree is 1/3 volume of cylinder with maxdiamstand, i.e., a cone:
   const double falldepth = 
     height * PiValue * diameter * diameter / 12.0 
     / tree->getNodePtr()->getVArea();
   const double wooddepth = tree->getBioMassStand();
   if( falldepth > wooddepth ) 
       return 0;
   const double rain = storm->getRainrate();
   const double rootstrength = tree->getRootStrength();
   const int n = 
     ROUND( rain * rain / rootstrength * blowparam * rand->ExpDev() );
   int i=0;
   double totfalldepth = 0.0;
   while( i<n && totfalldepth < wooddepth )
   {
      ++i;
      totfalldepth += falldepth;
   }
   if( totfalldepth > wooddepth ) 
       --i;
   return i;
}

int tForest::MaxNumBlowDown( tTrees *tree )
{
  const double diameter = tree->getMaxDiamStand();
  if( diameter == 0.0 ) return 0;
  const double height = tree->getMaxHeightStand();
   if( height == 0.0 ) return 0;
   // volume of tree is 1/3 volume of cylinder with maxdiamstand, i.e., a cone:
   const double falldepth = 
     height * PiValue * diameter * diameter / 12.0
     / tree->getNodePtr()->getVArea();
   const double wooddepth = tree->getBioMassStand();
   if( falldepth > wooddepth ) 
       return 0;
   int i=0;
   double totfalldepth = 0.0;
   while( totfalldepth < wooddepth )
   {
      ++i;
      totfalldepth += falldepth;
   }
   --i;
   return i;
}

// carries out the mechanics of felling the tree already determined to fall
void tForest::TreeFall( tTrees *tree )
{
  if( rand == 0 ) return;
  // choose sine and cosine of fall angle separately by choosing
  // two uniformly distributed random numbers between -1 and 1:
  const double sinangle = ( rand->ran3() - 0.5 ) * 2.0;
  const double cosangle = ( rand->ran3() - 0.5 ) * 2.0;
  int inode=0;
  tArray< double > uv0(2);
  tArray< double > uv1(2);
  tArray< double > uvtree(2);
  uvtree[0] = cosangle;
  uvtree[1] = sinangle;
  // find triangle falling into:
  tSpkIter sI( tree->getNodePtr() );
  tEdge* ce;
  for( ce = sI.FirstP(); !sI.AtEnd(); ce = sI.NextP() )
    {
      uv0 = UnitVector( ce );
      uv1 = UnitVector( sI.ReportNextP() );
      if( PointsCCW( uv0, uvtree, uv1 ) )
	break;
    }
  ce = sI.NextP();
  tTriangle* ct = ce->TriWithEdgePtr();
  tTriangle* nt = 0;
  // find coords of end of log by multiplying the tree height by the
  // sine and cosine and adding that to the node coords
  const double totlen = tree->getMaxHeightStand();
  const double diameter = tree->getMaxDiamStand();
  double dist = 0.0;
  while( dist < totlen )
    {
      ++dist; // do log one meter at a time
      const double segx = tree->getNodePtr()->getX() + dist * cosangle;
      const double segy = tree->getNodePtr()->getY() + dist * sinangle;
      nt = ct->NbrToward( segx, segy );
      if( nt > 0 )
	ct = nt;
      // find which node to give wood:
      int i;
      double mindist2 = 10000.;
      for( i=0; i<3; ++i )
	{
	  const double xsep = ct->pPtr(i)->getX() - segx;
	  const double ysep = ct->pPtr(i)->getY() - segy;
	  const double dist2 = xsep * xsep + ysep * ysep;
	  if( dist2 < mindist2 )
	    {
	      mindist2 = dist2;
	      inode = i;
	    }
	}
      const double segdiam = diameter * ( 1.0 - ( dist - 0.5 ) / totlen );
      const double addbiomass = PiValue * segdiam * segdiam / 4.0;
      // take wood away:
      tree->addBioMassStand( -addbiomass / tree->getNodePtr()->getVArea() );
      if( ct->pPtr( inode )->isNonBoundary() )
	{
	  tLNode *rLNode = static_cast<tLNode*>( ct->pPtr( inode ) );
	  // give wood away:
	  rLNode->getVegCover().getTrees()->
	    addBioMassStand( addbiomass / rLNode->getVArea() );
	}
    }
}
#undef PiValue

/**************************************************************************\
**                 FUNCTIONS FOR CLASS tTrees
\**************************************************************************/
// negative value of maxrootstrength is flag to set it later:
tTrees::tTrees( tLNode* nPtr, tForest* fPtr, double val )
  : maxrootstrength(-1.0), maxheightdown(0.0), biomassdown(0.0)
{
  // Initialize stand:
  TreesInitialize( nPtr, fPtr, val );
}



void tTrees::TreesInitialize( tLNode* nPtr, tForest* fPtr, double val )
{
  node = nPtr;
  forest = fPtr;
  standdeathtime = val;
  maxheightstand = forest->TreeHeightInit( this );
  maxdiamstand = forest->Tree_diam_from_height( this );
  biomassstand = forest->WoodVolumeInit( this );
  // also sets rootdecay, rootgrowth, maxrootstrength:
  rootstrength = forest->RootStrengthInit( this ); 
}


/***************************************************************************\
**
**  TreesEvolve( duration ): function to find and set rootstrength
**    according to Benda & Dunne 97a, Sidle 92, Sidle 91
**
**  11/98 SL
**  Added to CHILD (CU), 8/10, SL
\***************************************************************************/
void tTrees::TreesEvolve( double duration )
{
  if( rand == 0 ) 
    {
      TreesEvolveSimple( duration );
      return;
    }
  standdeathtime += duration;
  maxheightstand = 
    forest->TreeHeightEvol( this );
  biomassstand = 
    forest->WoodVolumeEvol( this );
  biomassdown *= forest->WoodDecayFactor( duration );
  rootstrength = 
    forest->RootStrengthEvol( this );
  maxdiamstand = forest->Tree_diam_from_height( this );
  // for tree fall (tree height must be greater than breast height, 1.37m,
  // for diameter formula, o.w., NaN):
  if( biomassstand > 0.0 && maxdiamstand > 0.0 )
    {
      int n = forest->NumBlowDown( this );
      int i;
      for( i=0; ( i<n && biomassstand > 0.0 ); ++i )
	forest->TreeFall( this );
      if( biomassstand <= 0.0 )
	// all trees fell down, roots removed, root decay constant 
	// and standdeathtime zeroed:
	biomassstand = maxheightstand = rootstrength = maxrootstrength 
	  = standdeathtime = 0.0;
    }
}

/***************************************************************************\
**
**  TreesEvolveSimple( duration ): function to find and set rootstrength
**    according to Benda & Dunne 97a, Sidle 92, Sidle 91. 
**    "Simple" in that it doesn't do treefall.
**
**  11/98 SL
**  Added to CHILD (CU), 8/10, SL
\***************************************************************************/
void tTrees::TreesEvolveSimple( double duration )
{
  standdeathtime += duration;
  maxheightstand = 
    forest->TreeHeightEvol( this );
  biomassstand = 
    forest->WoodVolumeEvol( this );
  biomassdown *= forest->WoodDecayFactor( duration );
  rootstrength = 
    forest->RootStrengthEvol( this );
  maxdiamstand = forest->Tree_diam_from_height( this );
}

void tTrees::TreesKill()
{
  maxrootstrength = 
    forest->RootDeath( this ); // kill roots
   standdeathtime = 0.0;// set standdeathtime to zero
   // knock down trees but keep root strength:
   int n = forest->MaxNumBlowDown( this );
   double i=0;
   while( i<n )
   {
      ++i;
      forest->TreeFall( this );
   }
   // if any left over, put it at "this" node:
   biomassdown += biomassstand;
   biomassstand = maxheightstand = 0.0;
}

void tTrees::TreesKillSimple()
{
   maxrootstrength = 
     forest->RootDeath( this ); // kill roots
   standdeathtime = 0.0;// set standdeathtime to zero
   // don't use TreeFall; keep standing biomass at "this" node:
   biomassdown += biomassstand;
   biomassstand = maxheightstand = 0.0;
}

void tTrees::TreesCut()
{
   maxrootstrength = 
     forest->RootDeath( this ); // kill roots
   // zero out stand and leave logs lying
   biomassstand = maxheightstand = standdeathtime = 0.0;
}

void tTrees::TreesAllFallDown()
{
   // zero out roots and stand:
   rootstrength = maxrootstrength = standdeathtime = 0.0;
   // knock over all the trees:
   int n = forest->MaxNumBlowDown( this );
   double i=0;
   while( i<n )
   {
      ++i;
      forest->TreeFall( this );
   }
   // if any left over, put it at "this" node:
   biomassdown += biomassstand;
   biomassstand = maxheightstand = 0.0;
}

/**************************************************************************\
**                 FUNCTIONS FOR CLASS tVegetation
\**************************************************************************/

/**************************************************************************\
**
**  tVegetation constructors:
**    1. Default
**    2. Input file
**
\**************************************************************************/

tVegetation::tVegetation()
  :
  mdKvd(0),
  mdTVeg(1),
  mdTauCritBare(0), mdTauCritVeg(0),
  fire(0), forest(0)
{}

tVegetation::tVegetation( tMesh<class tLNode> * meshPtr, const tInputFile &infile,
			  bool no_write_mode /*=false*/, 
			  tRunTimer *tPtr /*=0*/, tStorm *stormPtr /*=0*/ )
  :
  mdKvd(0),
  mdTVeg(1),
  mdTauCritBare(0), mdTauCritVeg(0),
  fire(0), forest(0)
{
  optGrassSimple = infile.ReadBool( "OPTGRASS_SIMPLE", true );
  if( optGrassSimple )
    {
      mdKvd = infile.ReadItem( mdKvd, "VEG_KVD" );
      mdTVeg = infile.ReadItem( mdTVeg, "VEG_TV" );
      mdTauCritBare = infile.ReadItem( mdTauCritBare, "TAUC" );
      mdTauCritVeg = infile.ReadItem( mdTVeg, "VEG_TAUCVEG" );
    }
  bool optForest = infile.ReadBool( "OPTFOREST", false );
  if( optForest )
    forest = new tForest( infile, meshPtr, stormPtr );
  bool optFire = infile.ReadBool( "OPTFIRE", false );
  if( optFire )
    fire = new tFire( infile, tPtr, no_write_mode );

  // Loop through nodes and set initial vegetation cover & threshold
  // (for now, assume constant; later need to add restart capability)
  // Note: assumes initially 100% cover.
  tMesh<tLNode>::nodeListIter_t niter( meshPtr->getNodeList() );

  // unused
  // intlVegCover = infile.ReadItem( intlVegCover, "INTLVEGCOV" );

  // Initialise nodes
  int opt;
  if ( (opt = infile.ReadItem( opt, "OPTREADINPUT" ))
       == OPTREADINPUT_PREVIOUS) {
    // Start from a previous computation
    tListInputDataVegetation inputVegData( infile );
    const int nnodes = meshPtr->getNodeList()->getSize();
    if (inputVegData.vegCov.getSize() != static_cast<size_t>(nnodes))
      ReportFatalError( "tVegetation(): invalid number of records"
			" in input file." );
    // Rely on the fact that the IDs have not been re-numbered.
    // for fast lookup per ID
    const tMesh< tLNode >::tIdArrayNode_t NodeTable(*(meshPtr->getNodeList()));
    for( int id=0; id < nnodes; ++id )
      NodeTable[id]->getVegCover().mdVeg = inputVegData.vegCov[id];
    for( tLNode *cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
      cn->setTauCrit( mdTauCritBare + mdTauCritVeg );
  } else {
    // Start from scratch
    for( tLNode *cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
      {
	cn->getVegCover().mdVeg = 1.0;
	cn->setTauCrit( mdTauCritBare + mdTauCritVeg );
      }
  }
}

tVegetation::tVegetation( const tVegetation& orig )
  : optGrassSimple(orig.optGrassSimple), mdKvd(orig.mdKvd), mdTVeg(orig.mdTVeg), 
    mdTauCritBare(orig.mdTauCritBare), mdTauCritVeg(orig.mdTauCritVeg),
    fire(0), forest(0)
{
  if( orig.fire )
    fire = new tFire( *orig.fire );
  if( orig.forest )
    forest = new tForest( *orig.forest );
}

tVegetation::~tVegetation()
{ 
  if( fire )
    {
      delete fire;
      fire = NULL;
    }
  if( forest )
    {
      delete forest;
      forest = NULL;
    }
}
/**************************************************************************\
**
**  tVegetation::UpdateVegetation
**
**  This routine updates the % vegetation cover at each node in response
**  to a storm event and its ensuing interstorm period. It also updates
**  the critical shear stress at each node according to its present
**  vegetation cover.
**
**  Erosion of vegetation during a storm is computed as:
**    dV/dt = -Kvd V ( tau - tauc )
**  Where V represents the proportional cover (0 to 1), tau is shear
**  stress, tauc is critical shear stress, and Kvd is a vegetation
**  erodibility coefficient. The equation is solved using an analytical
**  solution, using the expedient approximation that tauc is constant
**  during a given storm.
**
**  Regrowth of vegetation following a storm is computed as:
**    dV/dt = 1/Tv ( 1 - V )
**  where Tv is the timescale of vegetation regrowth. An analytical
**  solution is used to update the vegetation cover.
**
**  Finally, the critical shear stress is updated at each node using:
**    Tc = Tcb + V Tcv
**  where Tcb is critical shear stress in the absence of cover (bare)
**  and Tcv is critical shear stress under 100% cover.
**
**  Created January 2000, GT
**
\**************************************************************************/

void tVegetation::UpdateVegetation( tMesh<class tLNode> *meshPtr,
				    double dt,
                                    double interstormdur ) const
{
  ErodeVegetation( meshPtr, dt );
  GrowVegetation( meshPtr, interstormdur );
}


void tVegetation::ErodeVegetation(  tMesh<class tLNode> *meshPtr,
				    double dt ) const
{
  if( !optGrassSimple ) return;
  
   tMesh<tLNode>::nodeListIter_t niter( meshPtr->getNodeList() ); // Node iterator
   tLNode * cn;   // Ptr to current node
   double tauex,  // Excess shear stress
       veg;       // Fractional vegetation cover
   
   // Loop on active nodes, computing erosion during the time step of
   // duration dt.
   //   For both erosion and regrowth, we use an analytical solution for
   // veg cover given its initial value and duration of erosion/regrowth.
   for( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
   {
      // Erosion of vegetation during storm (if any)
      tauex = cn->getTau() - cn->getTauCrit();
      veg = cn->getVegCover().getVeg();
      if( tauex>0.0 )
      {
          veg = veg * exp( -mdKvd * tauex * dt );
          cn->getVegCover().mdVeg = veg;
	  assert( veg >= 0.0 );
	  assert( veg <= 1.0 );
          //cout << "veg after erosion: " << veg << endl;
      }

      // Update critical shear stress
      cn->setTauCrit( mdTauCritBare + veg*mdTauCritVeg );
      //cout << "tau crit: " << mdTauCritBare + veg*mdTauCritVeg << endl;
   }   
}


void tVegetation::GrowVegetation(  tMesh<class tLNode> *meshPtr,
				    double duration ) const
{
  tMesh<tLNode>::nodeListIter_t niter( meshPtr->getNodeList() ); // Node iterator
  tLNode * cn;   // Ptr to current node
  if( optGrassSimple )
    {
      double veg;       // Fractional vegetation cover
      // Loop on active nodes, computing regrowth during the interstorm period.
      //   For both erosion and regrowth, we use an analytical solution for
      // veg cover given its initial value and duration of erosion/regrowth.
      for( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
	{
	  // Regrowth following storm
	  //cout << "veg before regrowth: " << veg << endl;
	  veg = cn->getVegCover().getVeg();
	  veg = 1.0 - (1.0 - veg) * exp( -duration / mdTVeg );
	  cn->getVegCover().mdVeg = veg;
	  //cout << "veg after regrowth: " << veg << endl;
	  assert( veg >= 0.0 );
	  assert( veg <= 1.0 );

	  // Update critical shear stress
	  cn->setTauCrit( mdTauCritBare + veg*mdTauCritVeg );
	  //cout << "tau crit: " << mdTauCritBare + veg*mdTauCritVeg << endl;
	}
    }
  // check whether the tForest pointer is non-zero; if so, evolve the forest:
  if( forest )
    for( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
      cn->getVegCover().getTrees()->TreesEvolveSimple( duration );
  // check whether the tFire pointer is non-zero:
  if( fire )
    // find out whether there's a fire this time around:
    if( fire->Burn_WhenItsTime() )
      {
	if( optGrassSimple )
	  // adding an effect of fire on the existing simple grass-like
	  // veg cover for now; any and all should feel free to do 
	  // something else (SL, 8/10).
	  // remove vegetation cover:
	  for( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
	    {
	      cn->getVegCover().mdVeg = 0.0;
	      cn->setTauCrit( mdTauCritBare );
	    }
	// forest fires kill trees and leave behind the wood and roots,
	// which then decay:
	if( forest )
	  for( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
	    cn->getVegCover().getTrees()->TreesKillSimple();
      }
}

