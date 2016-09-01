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
**  $Id: tVegetation.h,v 1.12 2003/10/22 13:04:31 childcvs Exp $
*/
/**************************************************************************/

#ifndef TVEGETATION_H
#define TVEGETATION_H

#include <math.h>
#include "../tInputFile/tInputFile.h"
//#include "../tMesh/tMesh.h"
template<class tLNode> class tMesh;


class tVegetation
{
  public:
   tVegetation();
   tVegetation( tMesh<class tLNode> *, const tInputFile & );
   void UpdateVegetation( tMesh<class tLNode> *, double, double ) const;
   void GrowVegetation( tMesh<class tLNode> *, double ) const;
   void ErodeVegetation( tMesh<class tLNode> *, double ) const;

  private:
   double mdKvd;   // Vegetation erosion coefficient (dim's LT/M)
   double mdTVeg;  // Vegetation regrowth time scale (years)
   double mdTauCritBare;  // Erosion threshold on bare soil
   double mdTauCritVeg;   // Erosion threshold under 100% cover
   // unused
   //double intlVegCover;   // Initial vegetation cover
};

class tVegCover
{
   friend class tVegetation;

  public:
   tVegCover();
   tVegCover( const tVegCover & );
   double getVeg() const;
   
  private:
   double mdVeg;
};

inline tVegCover::tVegCover() :
  mdVeg(1.)
{}

inline tVegCover::tVegCover( const tVegCover &orig ) :
   mdVeg(orig.mdVeg)
{}

inline double tVegCover::getVeg() const
{
   return mdVeg;
}


#endif
