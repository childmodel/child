//-*-c++-*- 

/***************************************************************************/
/**
**  @file
**  @brief declaration of ParamMMFS_t
*/
/***************************************************************************/

/**
   @class ParamMMFS_t
   @brief Parameters read in Input File for MakeMeshFromScratch
*/
class ParamMMFS_t {
public:
  double xGrid, yGrid;
  int boundType;
  int kSloped;
  double upperZ;
  double mElev;
  int ptPlace;
  double delGrid;               // average spacing between nodes
  double xout, yout;            // coordinates of user-specified outlet
  int numPts;                   // total no. of interior pts (if random)

  ParamMMFS_t(tInputFile & );
};

