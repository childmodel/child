
/***************************************************************************/
/**
**  @file
**  @brief definition of ParamMMFS_t
*/
/***************************************************************************/

#include <math.h>
#include "../tListInputData/tListInputData.h"
#include "ParamMesh_t.h"

ParamMMFS_t::ParamMMFS_t(tInputFile &infile) :
  kSloped(0),
  xout(0), yout(0)
{
  //reads in size of mesh (meters)
  xGrid = infile.ReadItem( xGrid, "X_GRID_SIZE" );
  yGrid = infile.ReadItem( yGrid, "Y_GRID_SIZE" );
  //read type of open boundary:
  //  0 = corner outlet (lower left)
  //  1 = one open side (lower)
  //  2 = two opposite sides (upper and lower)
  //  3 = all sides
  //  4 = specify outlet coordinates
  boundType = infile.ReadItem( boundType, "TYP_BOUND" );
   //ng 12/99 added so that the initial surface could be sloped
   //with one side open bndry cndtn if kSloped = 1.
   if(boundType == kOpenSide){
      kSloped = infile.ReadItem( kSloped, "SLOPED_SURF" );
      if(kSloped)
          upperZ = infile.ReadItem( upperZ, "UPPER_BOUND_Z" );
   }
   //read mean elevation & amplitude for random variation in elevation
   mElev = infile.ReadItem( mElev, "MEAN_ELEV" );
   randElev = infile.ReadItem( randElev, "RAND_ELEV" );
   //reads method of point placement:
   //  0 = uniform grid;
   //  1 = perturbed grid;
   //  2 = random scatter;
   ptPlace = infile.ReadItem( ptPlace, "OPT_PT_PLACE" );
   //read point spacing or number of points (excluding four boundary points)
   if( ptPlace == kUniformMesh || ptPlace == kPerturbedMesh )
     {
       delGrid = infile.ReadItem( delGrid, "GRID_SPACING" );
       if( delGrid >= xGrid || delGrid >= yGrid )
	 ReportFatalError(
			  "Mesh point spacing must be smaller than total mesh width." );
     }
   else
     {
       numPts = infile.ReadItem( numPts, "NUM_PTS" );
       delGrid = sqrt( xGrid * yGrid / numPts );
     }
   if( boundType == kOppositeSidesOpen )
     {
       upperZ = infile.ReadItem( upperZ, "UPPER_BOUND_Z" );
     }
   else if( boundType == kSpecifyOutlet )
     {
       // Read outlet coordinates from input file
       xout = infile.ReadItem( xout, "OUTLET_X_COORD" );
       yout = infile.ReadItem( yout, "OUTLET_Y_COORD" );
     }
}

