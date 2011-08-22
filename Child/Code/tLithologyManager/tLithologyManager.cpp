//-*-c++-*-

/**************************************************************************/
/**
 **  @file tLithologyManager.cpp
 **
 **  @brief Implementation of the tLithologyManager class.
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

#include <vector>
#include <sstream>
#include "tLithologyManager.h"

using namespace std;

/**************************************************************************/
/**
 **  Basic constructor
 **
 **  (does nothing at the moment)
 */
/**************************************************************************/
tLithologyManager::
tLithologyManager() :
  meshPtr_(0)
{
  if(1) cout << "tLithologyManager default constructor" << endl;
}


/**************************************************************************/
/**
 **  Destructor
 **
 **  Closes and deletes any remaining ofstream objects.
 */
/**************************************************************************/
tLithologyManager::
~tLithologyManager()
{
  if(1) cout << "tLithologyManager destructor" << endl;
}


/**************************************************************************/
/**
 **  InitializeFromInputFile
 **
 */
/**************************************************************************/
void tLithologyManager::
InitializeFromInputFile( tInputFile &inputFile, tMesh<tLNode> *meshPtr )
{
  if(1) cout << "tLithologyManager::InitializeFromInputFile" << endl;
  
  meshPtr_ = meshPtr;
  
  // Read from the input file and do what the user requests
  
  /* Notes:
  
  Things we ought to handle include:
  - read from a restart run
  - read from another CHILD-format lith file with given name
  - set lithology according to a strat sequence of uniform thicknesses
    - with or without a uniform or varying regolith layer on top
  - set lithology according to a map pattern (uniform in depth)
    - with or without a uniform or varying regolith layer on top
  - set regolith thickness according to a map pattern
  - set to simple "bedrock and regolith" (current default)
  - ?? set according to a set of "contact surfaces" that represent the lower
    boundary of a unit
  
  */
	
	// See if the user wants to start with an existing .lay file. If so, read it
	// and set layers accordingly.
	bool user_wants_to_read_from_layfile = inputFile.ReadBool( "OPT_READ_LAYFILE", false );
	if( user_wants_to_read_from_layfile )
	{
		SetLithologyFromChildLayFile( inputFile );
	}
  
	// See if the user wants to modify layers according to an Etch File.
  // An Etch File specifies one or more layers, with given properties, to be
  // "etched in" to the current topography and lithology.
	bool user_wants_to_read_from_etchfile = inputFile.ReadBool( "OPT_READ_ETCHFILE", false );
	if( user_wants_to_read_from_etchfile )
	{
		SetLithologyFromEtchFile( inputFile );
	}
}


/**************************************************************************/
/**
 **  SetLithologyFromChildLayFile
 **
 **  This function sets the layers and properties based on data in a
 **  specified CHILD-format .lay file. The previous layers are removed.
 */
/**************************************************************************/
void tLithologyManager::
SetLithologyFromChildLayFile( const tInputFile &infile )
{
  int i, item, numl;
  int numlayernodes;
  double time;
  double ditem;
  char thestring[80], inname[80];
  std::ifstream layerinfile;
  infile.ReadItem( thestring, sizeof(thestring), "INPUT_LAY_FILE" );
  
  if (1) //DEBUG
    std::cout<<"in SetLithologyFromChildLayFile..."<<std::endl;
  
  strcpy( inname, thestring );
  layerinfile.open(inname); /* Layer input file pointer */
  assert( layerinfile.good() );
  
  // Read first line, which should contain the time corresponding to the run
  layerinfile >> time;

  if(1) //DEBUG
    std::cout<<"Time="<<time<<std::endl;
  
  // Read second line, which specified number of interior nodes in file.
  layerinfile >> numlayernodes;
  if(1) //DEBUG
    std::cout<<"nnodes in file="<<numlayernodes<<std::endl;
  if( numlayernodes != meshPtr_->getNodeList()->getActiveSize() )
    ReportFatalError( "Number of nodes in layer file doesn't match number in mesh" );
  
  // Create and initialize a tLayer object that will serve as the basic template
  // for all layers. Set its # of grain size classes.
  tLayer layer_template;
  int numg;
  numg = infile.ReadInt( "NUMGRNSIZE" );
  layer_template.setDgradesize(numg);
  
  // Get an iterator for the node list
  tLNode * cn;
  tMesh<tLNode>::nodeListIter_t ni( meshPtr_->getNodeList() );
    // Hack: make layers input backwards compatible for simulations
    // without bulk density; :
  bool optNewLayersInput = infile.ReadBool( "OPT_NEW_LAYERSINPUT", false );
	double rock_density, sed_density;
	if( !optNewLayersInput )
	{
		rock_density = infile.ReadDouble( "ROCKDENSITYINIT", false );
		if( rock_density==0.0 ) rock_density = kDefaultRockBulkDensity;
		sed_density = infile.ReadDouble( "SOILBULKDENSITY", false );
		if( sed_density==0.0 ) sed_density = kDefaultSoilBulkDensity;
	}
  
  for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
  {
    // Remove pre-existing layers
    for(i=cn->getNumLayer()-1; i>=0; i-- )
      cn->removeLayer(i);
    
    // Read the number of layers at the current node
    layerinfile >> numl;
    if( numl<1 )
    {
      std::cout << "In SetLithologyFromChildLayFile: node " 
      << cn->getPermID() << " has " << numl << " layers." 
      << std::endl;
      ReportFatalError( "Each interior point must have at least one layer." );
    }
    
    // For each layer, read its properties, assign them to our template layer,
    // and insert a copy of that layer to the bottom of the layer list for this
    // node.
    for(i = 1; i<=numl; i++)
    {
      // Read layer properties
      layerinfile >> ditem;     // layer creation time
      layer_template.setCtime(ditem);
      layerinfile >> ditem;     // layer "recent activity" time
      layer_template.setRtime(ditem);
      layerinfile >> ditem;     // layer surface exposure duration
      layer_template.setEtime(ditem);
      layerinfile >> ditem;     // layer thickness (m)
      if( ditem <= 0. )
      {
        std::cout << "MakeLayersFromInputData: Layer " << i 
        << " at node " << cn->getID()
        << " has thickness " << ditem << std::endl;
        ReportFatalError("Layers must have positive thickness.");
      }
      layer_template.setDepth(ditem);
      layerinfile >> ditem;     // layer erodibility
      layer_template.setErody(ditem);

			// layer alluvial-bedrock flag
			tLayer::tSed_t bedrock_alluvial_flag;
			int item;
      layerinfile >> item;
			if( item==0 )
				bedrock_alluvial_flag = tLayer::kBedRock;
			else 
			{
				if( item==1 )
					bedrock_alluvial_flag = tLayer::kSed;      
				else
				{
					std::cout << "MakeLayersFromInputData: Layer " << i 
					<< " at node " << cn->getID()
					<< " has flag value " << bedrock_alluvial_flag << std::endl;
					ReportFatalError("Flag value must be either zero or one.");
				}
			}
		  layer_template.setSed(bedrock_alluvial_flag);
				
			// Hack: to make this backwards compatible,
			// read bulk density if optNewLayersInput == true; otherwise, use values
			// specified in input file or defaults, as above.
      if( optNewLayersInput )
      {
        layerinfile >> ditem;   // layer bulk density
        layer_template.setBulkDensity(ditem);
      }
			else 
			{
				if( bedrock_alluvial_flag==tLayer::kSed ) layer_template.setBulkDensity( sed_density );
				else layer_template.setBulkDensity( rock_density );
			}
			
			for(int g=0; g<numg; g++){
        layerinfile >> ditem;   // equivalent thickness of grain size g
        layer_template.setDgrade(g, ditem);
      }
      
      // Insert a copy of the layer on the bottom of the layer stack
      cn->InsertLayerBack( layer_template );
    }
    
  }
    
}


/**************************************************************************/
/**
 **  SetLithologyFromEtchFile
 **
 **  This function sets the layers and properties based on data in a
 **  specified Etchfile.
 */
/**************************************************************************/
void tLithologyManager::
SetLithologyFromEtchFile( const tInputFile &infile )
{

}
    

/**************************************************************************/
/**
 **  SetRockErodibilityValuesAtAllDepths
 **
 **  This function sets the value of bedrock erodibility at each node
 **  according to the values in the parameter "erody". This parameter is
 **  assumed to be a vector of size equal to the number of nodes, and
 **  ordered according to node (permanent) ID. All rock layers at a node
 **  are given the same erodibility values. Sediment layers are not
 **  affected.
 */
/**************************************************************************/
void tLithologyManager::
SetRockErodibilityValuesAtAllDepths( vector<double> erody )
{
  tLNode * cn;
  tMesh<tLNode>::nodeListIter_t ni( meshPtr_->getNodeList() );
  
  for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
  {
    int node_id = cn->getPermID();
    double local_erody = erody[node_id];
    std::cout << "There are " << cn->getNumLayer() << " layers at node " << node_id << std::endl;
    for( int i=0; i<cn->getNumLayer(); i++ )
    {
      std::cout << "  Working on layer " << i << std:: endl;
      if( cn->getLayerSed(i)==tLayer::kBedRock )
        cn->setLayerErody(i, local_erody);
    }
  }
    
}

/**************************************************************************/
/**
 */
/**************************************************************************/
void tLithologyManager::
EtchLayerAboveHeightAtNode( double layer_base_height, tLNode * node,
                           tLayer &layer_properties, bool keep_regolith )
{
}

/**************************************************************************/
/**
 **  PointInPolygon
 **
 **  Determines whether the point (point_x, point_y) lies inside the
 **  polygon described by vertx, verty. Ray casting algorithm,
 **  re-written in C++ from web-published C code by W. Randolph Franklin.
 */
/**************************************************************************/
bool tLithologyManager::
PointInPolygon( std::vector<double> vertx, 
               std::vector<double> verty,
               double point_x, double point_y )
{
  unsigned i, j;
  unsigned n_vertices = vertx.size();
  assert( n_vertices==verty.size() );
  bool point_is_in = false;
  for( i=0, j = n_vertices-1; i<n_vertices; j = i++ )
  {
    if ( ((verty[i]>point_y) != (verty[j]>point_y)) &&
        (point_x < (vertx[j]-vertx[i]) * (point_y-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
      point_is_in = !point_is_in;
  }
  return point_is_in;
}

/**************************************************************************/
/**
 **  EtchLayerAbove2DSurface
 **
 **  This function inserts a new layer with a given set of properties
 **  above a 2D surface described by a polynomial, and within the map
 **  area described by a second (optional) polynomial.
 */
/**************************************************************************/
void tLithologyManager::
EtchLayerAbove2DSurface(  vector<double> &poly_coefs_x,
                          vector<double> &poly_coefs_y,
                          tLayer &layer_properties,
                          bool keep_regolith,
                          bool use_bounding_polygon,
                          vector<double> bounding_poly_x,
                          vector<double> bounding_poly_y )
{
  tLNode * cn;
  tMesh<tLNode>::nodeListIter_t ni( meshPtr_->getNodeList() );
  
  for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
  {
    bool layer_is_present_here = true;
    if( use_bounding_polygon )
      if( !PointInPolygon( bounding_poly_x, bounding_poly_y, cn->getX(), cn->getY() ) )
         layer_is_present_here = false;
    if( layer_is_present_here )
    {
      double x = cn->getX();
      double y = cn->getY();
      double new_layer_base_height = poly_coefs_x[0]*x*x*x +
                              poly_coefs_x[1]*x*x +
                              poly_coefs_x[2]*x +
                              poly_coefs_x[3] +
                              poly_coefs_y[0]*y*y*y +
                              poly_coefs_y[1]*y*y +
                              poly_coefs_y[2]*y +
                              poly_coefs_y[3];
      if( new_layer_base_height < cn->getZ() )
        EtchLayerAboveHeightAtNode( new_layer_base_height, cn, layer_properties, keep_regolith );
    }
  }
    
}

/**************************************************************************/
/**
 */
/**************************************************************************/

