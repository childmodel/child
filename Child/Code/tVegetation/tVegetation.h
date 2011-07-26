//-*-c++-*- 

/**************************************************************************/
/**
**  @file tVegetation.h
**  @brief Header file for the tVegetation and tVegCover classes
**
**  Classes tVegetation and tVegCover represents the vegetation cover
**  across a terrain. Class tVegCover represents the properties of
**  vegetation at a point on the landscape (e.g., the %cover), while
**  class tVegetation represents the landscape-wide ("global")
**  properties (e.g., regrowth coefficient). Class tVegCover is designed
**  to be embedded in a node (in other words, each node "has a"
**  vegetation cover).
**
**  Design notes:
**  Classes tVegetation and tVegCover work together, with one being
**  "local" (at-a-node) and the other "global" (common properties across
**  the terrain and functions for terrain-wide updating). Class tVegCover
**  holds tVegetation as a friend, so tVegetation methods can directly
**  access its members. However, tVegCover provides no public "set"
**  methods, so that only the routines in tVegetation can modify the
**  vegetation properties.
**
**  Created January, 2000, GT
**
**  STL, 8/2010: Added forest and fires adapted from the OSU version of
**  CHILD. New classes: tFire, tForest, tTrees.
**  The idea is that tForest and tTrees mirror the divisions between
**  tVegetation and tVegCover: tForest contains the parameters and rules,
**  and tTrees are contained by each node and therefore contain local data 
**  and a few functions (mainly for different sorts of killing of trees).
**  The forest functions, and the parameters that have been used, are 
**  based on the work of several authors, Duan, Sidle, Benda and Dunne,
**  primarily to model Douglas-fir forests of the Pacific Northwest west
**  of the Cascades crest. Similarly, tFire has been used to model the
**  severe, infrequent fires of that region.
**  
**  $Id: tVegetation.h,v 1.12 2003-10-22 13:04:31 childcvs Exp $
*/
/**************************************************************************/

#ifndef TVEGETATION_H
#define TVEGETATION_H

#include <math.h>
#include "../tInputFile/tInputFile.h"
#include "../Mathutil/mathutil.h"
#include <iosfwd>
//#include "../tMesh/tMesh.h"
class tLNode;
template<class tLNode> class tMesh;
class tVegetation;
class tTrees;
class tStorm;
class tRunTimer;

class tFire
{
public:
  tFire();
  tFire( const tFire& );
  tFire( double, unsigned, tRunTimer* );
  tFire( const tInputFile&, tRunTimer*, bool no_write_mode = false );
  bool Burn_WhenItsTime();
  double interfireDur();
  double getMeanInterfireDur() const;
  void setMeanInterfireDur( double );
  void setTimePtr( tRunTimer* ptr ) {timePtr = ptr;}
  void TurnOnOutput( const tInputFile& );
  void TurnOffOutput();
 
private:
  bool optRandom;   
  double ifrdurMean;  // Mean time between fires
  double ifrdur;      // Actual time between fires
  double time_to_burn; // time remaining until burn happens
  tRunTimer* timePtr;
  std::ofstream firefs;  // fstream for writing out fire sequence
  tRand* rand; // give fires their own random number generator
};

class tForest
{
  friend class tTrees;
  friend class tVegetation;
public:
  tForest( const tForest& );
  tForest( const tInputFile&, tMesh<tLNode>*, tStorm* );
  // "get" and "set" functions:
  double getRootDecayK() const {return rootdecayK;}
  void setRootDecayK( double val ) {rootdecayK = val;}
  double getRootDecayN() const {return rootdecayN;}
  void setRootDecayN( double val ) {rootdecayN = val;}
  double getRootGrowthA() const {return rootgrowthA;}
  void setRootGrowthA( double val ) {rootgrowthA = val;}
  double getRootGrowthB() const {return rootgrowthB;}
  void setRootGrowthB( double val ) {rootgrowthB = val;}
  double getRootGrowthC() const {return rootgrowthC;}
  void setRootGrowthC( double val ) {rootgrowthC = val;}
  double getRootGrowthF() const {return rootgrowthF;}
  void setRootGrowthF( double val ) {rootgrowthF = val;}
  double getRootStrengthJ() const {return rootstrengthJ;}
  void setRootStrengthJ( double val ) {rootstrengthJ = val;}
  double getMVRC() const {return MVRC;}
  void setMVRC( double val ) {MVRC = val;}
  double getMLRC() const {return MLRC;}
  void setMLRC( double val ) {MLRC = val;}   
  double getRSPartition() {return RSPartition;}
  double getHeightIndex() const {return heightindex;}
  void setHeightIndex( double val ) {heightindex = val;}
  double getWeightMax() const {return weightMax;}
  void setWeightMax( double val ) {weightMax = val;}
  double getWeightA() const {return weightA;}
  void setWeightA( double val ) {weightA = val;}
  double getWeightB() const {return weightB;}
  void setWeightB( double val ) {weightB = val;}
  double getWeightC() const {return weightC;}
  void setWeightC( double val ) {weightC = val;}
  double getWeightK() const {return weightK;}
  void setWeightK( double val ) {weightK = val;}
  double getWoodDensity() const {return wooddensity;}
  void setWoodDensity( double val ) {wooddensity = val;}
  double getWoodDecayK() const {return wooddecayK;}
  void setWoodDecayK( double val ) {wooddecayK = val;}

  void setMeshPtr( tMesh<tLNode>* ptr );
  void setStormPtr( tStorm* ptr ) {storm = ptr;}

  // functions that do things:
  double Tree_diam_from_height( tTrees* ); // simple calc
  // initialization:
  double TreeHeightInit( tTrees* );
  double WoodVolumeInit( tTrees* );
  double RootStrengthInit( tTrees* );
  void ForestInitialize( const tInputFile&, tMesh<tLNode>*, tStorm* );
  // evolution:
  double TreeHeightEvol( tTrees* );
  double WoodVolumeEvol( tTrees* );
  double RootStrengthEvol( tTrees* );
  double WoodDecayFactor( double );
  double RootDeath( tTrees* );

  // Treefall functions
  int NumBlowDown( tTrees* );
  int MaxNumBlowDown( tTrees* );
  void TreeFall( tTrees* );

private:
  tMesh<tLNode> *mesh;
  tStorm *storm;
  tRand *rand;
  double rootdecayK;
  double rootdecayN;
  double rootgrowthA;
  double rootgrowthB;
  double rootgrowthC;
  double rootgrowthF;
  double rootstrengthJ;
  double MVRC;
  double MLRC;
  double RSPartition;
  double heightindex;
  double weightMax;
  double weightA;
  double weightB;
  double weightC;
  double weightK;
  double wooddensity;
  double blowparam;
  double diamB0;
  double diamB1;
  double diamB2;
  double wooddecayK;
};


class tVegetation
{
  public:
   tVegetation();
  tVegetation( const tVegetation& );
  tVegetation( tMesh<class tLNode> *meshPtr, const tInputFile &infile, 
	       bool no_write_mode = false, tRunTimer *tPtr = 0, tStorm *stormPtr = 0 );
  ~tVegetation();
   void UpdateVegetation( tMesh<class tLNode> *, double, double ) const;
   void GrowVegetation( tMesh<class tLNode> *, double ) const;
   void ErodeVegetation( tMesh<class tLNode> *, double ) const;
  tFire* FirePtr() {return fire;}
  tForest* ForestPtr() {return forest;}

  private:
  bool optGrassSimple; // option for simple grass
   double mdKvd;   // Vegetation erosion coefficient (dim's LT/M)
   double mdTVeg;  // Vegetation regrowth time scale (years)
   double mdTauCritBare;  // Erosion threshold on bare soil
   double mdTauCritVeg;   // Erosion threshold under 100% cover
   // unused
   //double intlVegCover;   // Initial vegetation cover
  tFire *fire; // pointer to fire object
  tForest *forest; // pointer to forest object
};

class tTrees
{
public:
  tTrees() : forest(0), node(0), maxrootstrength(-1.0) {}
  tTrees( const tTrees& );
  tTrees( tLNode*, tForest*, double );
  ~tTrees() {}
  tTrees& operator=( const tTrees& );
  // "get" and "set":
  void setForestPtr( tForest* ptr ) {forest = ptr;}
  tLNode* getNodePtr() {return node;}
  void setNodePtr( tLNode* ptr ) {node = ptr;}
  double getRootStrength() const {return rootstrength;}
  void setRootStrength( double val ) {rootstrength = val;}
  void addRootStrength( double val ) {rootstrength += val;}
  double getRootStrengthLat() const {return rootstrengthLat;}
  void setRootStrengthLat( double val ) {rootstrengthLat = val;}
  void addRootStrengthLat( double val ) {rootstrengthLat += val;}
  double getRootStrengthVert() const {return rootstrengthVert;}
  void setRootStrengthVert( double val ) {rootstrengthVert = val;}
  void addRootStrengthVert( double val ) {rootstrengthVert += val;}
  double getRootGrowth() const {return rootgrowth;}
  void setRootGrowth( double val ) {rootgrowth = val;}
  void addRootGrowth( double val ) {rootgrowth += val;}
  double getRootDecay() const {return rootdecay;}
  void setRootDecay( double val ) {rootdecay = val;}
  void addRootDecay( double val ) {rootdecay += val;}
  double getMaxRootStrength() const {return maxrootstrength;}
  void setMaxRootStrength( double val ) {maxrootstrength = val;}
  double getMaxHeightStand() const {return maxheightstand;}
  void setMaxHeightStand( double val ) {maxheightstand = val;}
  void addMaxHeightStand( double val ) {maxheightstand += val;}
  double getMaxDiamStand() const {return maxdiamstand;}
  void setMaxDiamStand( double val ) {maxdiamstand = val;}
  void addMaxDiamStand( double val ) {maxdiamstand += val;}
  double getBioMassStand() const {return biomassstand;}
  void setBioMassStand( double val ) {biomassstand = val;}
  void addBioMassStand( double val ) {biomassstand += val;}
  double getMaxHeightDown() const {return maxheightdown;}
  void setMaxHeightDown( double val ) {maxheightdown = val;}
  void addMaxHeightDown( double val ) {maxheightdown += val;}
  double getBioMassDown() const {return biomassdown;}
  void setBioMassDown( double val ) {biomassdown = val;}
  void addBioMassDown( double val ) {biomassstand += val;}
  double getStandDeathTime() const {return standdeathtime;}
  void setStandDeathTime( double val ) {standdeathtime = val;}
  void addStandDeathTime( double val ) {standdeathtime += val;}
  // functions that do stuff:
  // initialization:
  void TreesInitialize( tLNode*, tForest*, double );
  // evolution:
  void TreesEvolve( double );
  void TreesEvolveSimple( double );
  // various ways of killing and/or removing
  void TreesKill();
  void TreesKillSimple();
  void TreesCut();
  void TreesRemove();
  void TreesRemove( double );
  void TreesAllFallDown();

private:
  tForest *forest; // pointer to tForest object
  tLNode *node; // pointer to container node
  double rootstrength;
  double rootstrengthLat; // lateral cohesion
  double rootstrengthVert; // vertical cohesion
  double rootgrowth;
  double rootdecay;
  double maxrootstrength; // root strength at time of death
  double maxheightstand;
  double maxdiamstand;
  double biomassstand; // mass or volume per unit area
  double maxheightdown;
  double biomassdown; // mass or volume per unit area
  double standdeathtime;
};

class tVegCover
{
   friend class tVegetation;

  public:
   tVegCover();
   tVegCover( const tVegCover & );
   double getVeg() const;
  tTrees* getTrees() const {return trees;}
  void setTrees( tTrees* ptr ) {trees = ptr;}
   
  private:
   double mdVeg;
  tTrees *trees;
};

//
// inline functions for tFire
//

/*
**  InterfireDur
**
**  Returns the interfire duration.
*/
inline double tFire::interfireDur()
{
   return ifrdur;
}

inline double tFire::getMeanInterfireDur() const {return ifrdurMean;}

inline void tFire::setMeanInterfireDur( double tm ) {ifrdurMean = (tm>0.0) ? tm : 0.0;}

//
// inline functions for tForest
//

// decay downed wood:
// formula from Harmon et al., Ecology of  CWD in Temperate Ecosystems, _Adv. 
// in Ecol. Res._, 15, pp. 133-302, 1986.
// find factor by which to multiply biomassdown:
inline double tForest::WoodDecayFactor( double evoltime )
{
   assert( evoltime >= 0.0 );
   return exp( -wooddecayK * evoltime );
}


//
// inline functions for tTrees
//
inline tTrees::tTrees( const tTrees& orig )
  : forest(0),
    node(0),
    rootstrength(orig.rootstrength),
    rootstrengthLat(orig.rootstrengthLat),
    rootstrengthVert(orig.rootstrengthVert),
    rootgrowth(orig.rootgrowth),
    rootdecay(orig.rootdecay),
    maxrootstrength(orig.maxrootstrength),
    maxheightstand(orig.maxheightstand),
    maxdiamstand(orig.maxdiamstand),
    biomassstand(orig.biomassstand),
    maxheightdown(orig.maxheightdown),
    biomassdown(orig.biomassdown),
    standdeathtime(orig.standdeathtime)
{}


inline tTrees& tTrees::operator=( const tTrees& right )
{
   if( &right != this )
   {
      node = 0;
      forest = 0;
      rootstrength = right.rootstrength;
      rootstrengthLat = right.rootstrengthLat;
      rootstrengthVert = right.rootstrengthVert;
      rootgrowth = right.rootgrowth;
      rootdecay = right.rootdecay;
      maxrootstrength = right.maxrootstrength;
      maxheightstand = right.maxheightstand;
      maxdiamstand = right.maxdiamstand;
      biomassstand = right.biomassstand;
      maxheightdown = right.maxheightdown;
      biomassdown = right.biomassdown;
      standdeathtime = right.standdeathtime;
   }
   return *this;
}  
 
inline void tTrees::TreesRemove()
{
   // zero out roots and stand and logs:
   rootstrength = maxrootstrength = biomassstand = maxheightstand = 
     biomassdown = standdeathtime = 0.0;
}

// modified 5/22/02 SL: if vegetation removal results in live biomass
// that is less than 90% of maximum, do equivalent of TreesAllFallDown, 
// but put any leftover biomass at "this" node
#define TOL 1e-8
inline void tTrees::TreesRemove( double frac )
{
   if( frac < TOL ) return;
   if( frac > 1.0-TOL )
   {
      TreesRemove();
      return;
   }
   double remain = 1.0 - frac;
   // remove roots:
   rootstrength *= remain;
   // remove standing trees:
   biomassstand *= remain;
   // remove down logs:
   biomassdown *= remain;
   // if vegetation removal results in live biomass that is less than
   // 90% of maximum, do equivalent of TreesAllFallDown, but put 
   // any leftover biomass at "this" node:
   if( biomassstand < 0.9 * forest->weightMax / forest->wooddensity / GRAV )
   {
      // zero out roots and stand:
      rootstrength = maxrootstrength = standdeathtime = 0.0;
      // put any left over biomass at "this" node:
      biomassdown += biomassstand;
      biomassstand = maxheightstand = 0.0;
   }
}
#undef TOL

// inline functions for tVegCover

inline tVegCover::tVegCover() :
  mdVeg(1.), trees(0)
{}

inline tVegCover::tVegCover( const tVegCover &orig ) :
  mdVeg(orig.mdVeg), trees(0)
{
  if( orig.trees )
    trees = new tTrees( *orig.trees );
}

inline double tVegCover::getVeg() const {return mdVeg;}


#endif
