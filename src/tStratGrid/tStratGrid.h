//-*-c++-*-

/**************************************************************************/
/**
 **  @file  tStratGrid.h
 **  @brief header file for class tStratGrid and its helper class
 **         tStratNode.
 **
 **  Class tStratNode takes care of the functions belonging to the rectangular grid
 **  used to attach the stratigraphy, Class tStratNode represents the functions
 **  and data members associated with an individual rect. grid node.
 **
 **  Created 5/2003 (QC)
 **
 **  $Id: tStratGrid.h,v 1.7 2005-03-15 17:17:30 childcvs Exp $
 */
/**************************************************************************/

#ifndef  STRATGRID_H
#define  STRATGRID_H

#include "../tInputFile/tInputFile.h"
#include "../tMatrix/tMatrix.h"
#include "../tList/tList.h"

template< class T > class tMesh;
class tTriangle;

class tLayer;
class tLNode;
// necessary to get tLayer::tSed_t
#include "../tLNode/tLNode.h"

//---------------------STRATNODE---------------------------------------------


/**************************************************************************/
/**							         STRATNODE
 **  @class tStratNode
 **
 **  @brief Data record used in constructing a list of flood nodes.
 **  Contains a selection of data members and functions copied from tLNode
 **  ,but is only used for stratigraphy operations on a rectangular grid
 **
 **  Created 5/2003 (QC & GT)
 */
/**************************************************************************/

class tStratNode
{
public:
  tStratNode();
  tStratNode( int ); // for tMatrix
  tStratNode( tInputFile const &infile );
  tStratNode( double , double );
  tStratNode( double , double , const tStratNode &);
  tStratNode( const tStratNode & );
  tStratNode &operator=( const tStratNode & );
  ~tStratNode();
  // Simple Set and Get functions

  void setX( double );            // sets x coord
  void setY( double );            // sets y coord
  void setZ( double );            // sets z value
  void setI( int );
  void setJ( int );
  void setSectionBase(double);	  // sets the base of the stratigraphic section
  void setNewZ( double );	  // sets new z value, used for interpolation

  double getX() const {           // returns x coord
    return x;
  }
  double getY() const {            // returns y coord
    return y;
  }
  double getZ() const {            // returns z value
    return z;
  }
  double getNewZ() const {	  // returns a new z value
    return newz;
  }
  double getSectionBase() const { // returns the original base of the section
    return sectionZ;
  }
  int getI() const {
    return i;
  }
  int getJ() const {
    return j;
  }



  double getTotalLayerDepth() const;
  void setGrade( int, double ) const;
  double getGrade( int ) const;
  const tArray< double >& getGrade( ) const;
  int getNumg() const;
  void setNumg( int ) const;
  double getMaxregdep() const;
  double getLayerCtime(int) const;		// Creation time
  double getLayerRtime(int) const;		// Recent activity time
  double getLayerEtime(int) const;		// exposure time
  double getLayerDepth(int) const;		// thickness
  double getLayerErody(int) const;		// eodability
  double getLayerPaleoCurrent(int) const;	// paleocurrent direction
  double getAgeAtDepth(double) const;
  tLayer::tSed_t getLayerSed(int) const;
  double getLayerDgrade(int, int) const;  // first int is layer index
  int getNumLayer() const;
  void setLayerCtime(int, double);
  void setLayerRtime(int, double);
  void setLayerEtime(int, double);
  void addLayerEtime(int, double);
  void setLayerDepth(int, double);
  void setLayerErody(int, double);
  void setLayerPaleoCurrent(int, double);
  void setLayerSed(int, tLayer::tSed_t);
  void setLayerDgrade(int, int, double);
  tArray<double> EroDep(int, tArray<double>, double);
  void EroDepSimple(int, tArray<double>, double, double);         // layer, dh[i], time, paleocurrent

  tArray<double> addtoLayer(int, double,double);  // layer, g, fill[g], tt
  void addtoLayer(int, int, double, double, double);

  void makeNewLayerBelow(int, tLayer::tSed_t, double, tArray<double>const&,
			 double, double);
  double getPaleoCurrent(int) const;
  void removeLayer(int);
  void InsertLayerBack( tLayer const & );
  void CopyLayerList( tStratNode const * );
  double FindCurrentSectionBase();
  double AlluvialColumnThickness();
  void setClosestNode( tLNode *n_ ){
    ClosestNode = n_;
  }
  tLNode* getClosestNode();

#ifndef NDEBUG
  void TellAll();
#endif

protected:
  tList< tLayer > layerlist;    // Stratigraphic layerlist at this tStratGrid node location
  tLNode *ClosestNode;		// ptr to the closest tLNode (unused?)
  double x;         		  // x coordinate
  double y;         		  // y coordinate
  double z;         		  // z value (representing height)
  double sectionZ;		  // original base of the stratigraphic column
  double newz;
  int i;
  int j;

  static int numg;
  static tArray< double > grade;
  static double maxregdep;      // Active layer, and stratigraphic layer thickness.
  static double KRnew;          // erodability of a new layer
};

//--------------------------END OF STRATNODE--------------------------------


//--------------------------STRATGRID----------------------------------------

/**************************************************************************/
/**
 **  @class tStratGrid
 **
 **  Class tStratGrid contains the values and functions for a equidistant
 **  rectangular grid, which carries the stratigraphy.
 **
 */
/**************************************************************************/
class tStratGrid
{
  tStratGrid();

public:
  tStratGrid( tInputFile const &infile, tMesh<tLNode> *mp );
  tStratGrid( const tStratGrid& );
  tStratGrid& operator=(const tStratGrid&);
  ~tStratGrid();

  int getImax() const
  {
    return imax;
  }

  int getJmax() const
  {
    return jmax;
  }

  int getSectionLocation(int) const;
  tMatrix<tStratNode> const *getStratNodeMatrix() const {
    return StratNodeMatrix;
  }
  tMatrix<tTriangle*> const *getStratConnect() const{
    return StratConnect;
  }
  typedef enum {
    k0,
    k1,
    k2,
    k3,
    k4
  } tUpdate_t;

  void UpdateStratGrid(tUpdate_t, double);
  void ResetAccummulatedDh();
  void InterpolateElevations();
  void setSectionBase();
  void InterpolateErodep(double);
  void SweepChannelThroughRectGrid(double);      //
  void InterpolateErodepFromElevations(double);  // does the same as sweep, but is called in the beginning of a timestep/main loop
  void CheckSectionBase(int);
  void setSurface(int, double);
  void setSubsurface(int, double );
  void setSubsurface_mbelt(int, double);
  void setOutputTime(int, double );
  double getSurface(int) const;
  double getSubsurface(int) const;
  double getSubsurface_mbelt(int) const;
  double getOutputTime(int) const;
  int getnWrite() const;
  void setMesh( tMesh<tLNode>* ptr ){ mp = ptr;}
  void updateConnect();
  double CalculateMeanderCurrent(tTriangle *, double, double) const;
  double CompassAngle(tLNode *,tLNode *) const;

protected:			          // can be accessed by friend classes



private:

  int xcorner, ycorner;	                 // starting point for spanning StratGrid
  double griddx;		 	 // inter-node distance of tStratGrid
  double grwidth;              	         // x-width of tStratGrid
  double grlength;             	         // y-width of tStratGrid
  tMesh<tLNode> *mp;              	 // ptr to triangular mesh
  tMatrix<tStratNode> *StratNodeMatrix;  // Matrix of StratNodes
  tMatrix<tTriangle*> *StratConnect;  // Connectivity to tMesh
  int imax, jmax;
  int optSurferFiles;
  int nWrite;
  tArray<int>section;		          // array with 10 section locations, 5x, and 5y
  tArray<double>surface;		 // timeslice specific surface area
  tArray<double>subsurface;              // timeslice specific subsurface cummulative height in the entire floodplain
  tArray<double>subsurface_mbelt;        // timeslice specific subsurface cummulative height in meander belt
  tArray<double>outputTime;              // array with all the exact timings of the output /writesteps
};

//----------------------END OF STRATGRID-------------------------------------


#endif
