//-*-c++-*-

/***************************************************************************/
/**
**  @file
**  @brief declaration of ParamMMFS_t
*/
/***************************************************************************/

#include "../Definitions.h"
class tInputFile;

/**
   @class ParamMMFS_t
   @brief Parameters read in Input File for tMesh::MakeMeshFromScratch
*/
class ParamMMFS_t {
public:
  double xGrid, yGrid;
  tOpenBoundary_t boundType;
  int kSloped;
  double upperZ;
  double mElev;
  double randElev;              // amplitude factor for random var in elev
  tMeshType_t ptPlace;
  double delGrid;               // average spacing between nodes
  double xout, yout;            // coordinates of user-specified outlet
  int numPts;                   // total no. of interior pts (if random)

  ParamMMFS_t(const tInputFile & );
private:
  ParamMMFS_t();
};

