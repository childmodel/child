//-*-c++-*-

/***************************************************************************/
/**
**  @file
**  @brief declaration of ParamMMFS_t
*/
/***************************************************************************/

class tInputFile;

/**
   @class ParamMMFS_t
   @brief Parameters read in Input File for tMesh::MakeMeshFromScratch
*/
class ParamMMFS_t {
public:
  typedef enum { // type of open boundary (used in tMesh::MakeMeshFromScratch)
    kCornerOutlet = 0,       // corner outlet (lower left)
    kOpenSide = 1,           // one open side (lower)
    kOppositeSidesOpen = 2,  // two opposite sides (upper and lower)
    kAllSidesOpen = 3,       // all sides
    kSpecifyOutlet = 4 ,     // specify outlet coordinates
    kAllSideClosed = -1      // all sides closed
  } tOpenBoundary_t;

  typedef enum {      // method of grid construction
    kUniformMesh = 0,        // uniform grid;
    kPerturbedMesh = 1,      // perturbed grid
    kRandomMesh = 2          // random scatter
  } tMeshType_t;

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

