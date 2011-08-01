//-*-c++-*-

/**************************************************************************/
/**
**  @file tLithologyManager.h
**
**  @brief Header file for the tLithologyManager class.
**
**  A tLithologyManager is a class that reads, sets up, and/or modifies 
**  CHILD's lithology database.
**
**  The class was first written in April 2011 by GT for a project that
**  explores dynamic changes in rock strength as a function of deformation,
**  based on an idea proposed by Peter Koons. It is also designed to make
**  the process of configuring and modifying CHILD's representation of
**  lithology easier than it has been in the past (which has required
**  directly modifying the .lay output files, then re-starting).
**
**  In CHILD, lithology data is attached directly to the mesh, in the form
**  of sequences of (one or more) "layers" below each node. The idea
**  behind tLithologyManager is to provide methods (functions) that 
**  modify this layer information to suit the needs of different 
**  applications. For example, to explore the strain softening and rock
**  strength hypothesis, we want to have a way to signal to CHILD that
**  it should alter the strength of rock layers at particular (x,y) 
**  (locations) (or perhaps (x,y,z) locations). Equally, for applications
**  that explore how variations in rock strength in a "layercake" 
**  sequence influence patterns of landscape development, we would like
**  to have a simple way to set up variations with depth.
**
**  For information regarding this program, please contact Greg Tucker at:
**
**     Cooperative Institute for Research in Environmental Sciences (CIRES)
**     and Department of Geological Sciences
**     University of Colorado
**     2200 Colorado Avenue, Campus Box 399
**     Boulder, CO 80309-0399
*/
/**************************************************************************/

#ifndef TLITHOLOGYMANAGER_H
#define TLITHOLOGYMANAGER_H

#include <vector>
#include <string>
#include "../tInputFile/tInputFile.h"
#include "../tMesh/tMesh.h"
#include "../tLNode/tLNode.h"

#define kDefaultBulkDensity 2000.0

class tLithologyManager
{
public:

  // Constructor for basic initialization
  tLithologyManager();
  
  void InitializeFromInputFile( tInputFile &infile, tMesh<tLNode> *mesh);
  
  void SetLithologyFromChildOutputFile( string output_filename );
	
  void SetPropertiesAtNode( int nodeID, double erodibility, int sediment_flag, 
						    vector<double> grain_sizes = vector<double>(), 
						    double bulk_density = kDefaultBulkDensity,
						    double creation_time = 0.0, double rtime = 0.0 
						    double etime = 0.0, double paleocurrent = 0.0 );
  void SetErodibilitiesAtNodes( vector<double> erodibilities );
  void SetErodibilityWithinPolygon();
	
  // Destructor
  ~tLithologyManager();
  
private:
  tMesh<tLNode> *meshPtr_;
	
};


#endif
