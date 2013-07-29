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
 **  Etchlayer's TellData() method
 **
 **  (used for debugging)
 */
/**************************************************************************/
void Etchlayer::
TellData()
{
  cout << "Data in Etchlayer\n";
  cout << "  Erody = " << layer_properties_.getErody() << endl;
  cout << "  Sed type = " << layer_properties_.getSed() << endl;
  cout << "  ax = " << ax << endl;
  cout << "  bx = " << bx << endl;
  cout << "  cx = " << cx << endl;
  cout << "  ay = " << ay << endl;
  cout << "  by = " << by << endl;
  cout << "  cy = " << cy << endl;
  cout << "  d = " << d << endl;
  cout << "  keep_regolith = " << keep_regolith_ << endl;
  cout << "  use poly = " << use_bounding_polygon_ << endl;
  cout << "  in poly = " << layer_is_inside_poly_ << endl;
  if( layer_is_inside_poly_ )
  {
    cout << "  poly coords:\n";
    for( unsigned i=0; i<px.size(); i++ )
      cout << "    x: " << px[i] << " y: " << py[i] << endl;
  }
  
}


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

	// See if the user wants to set initial rock erodibility values at all depths
    // based on values in a file.
	bool user_wants_erody_from_file = inputFile.ReadBool( "OPT_SET_ERODY_FROM_FILE", false );
	if( user_wants_erody_from_file )
	{
		SetRockErodyFromFile( inputFile );
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
  int i, numl;
  int numlayernodes;
  double time;
  double ditem;
  char thestring[80], inname[80];
  std::ifstream layerinfile;
  infile.ReadItem( thestring, sizeof(thestring), "INPUT_LAY_FILE" );
  
  if (0) //DEBUG
    std::cout<<"in SetLithologyFromChildLayFile..."<<std::endl;
  
  meshPtr_->RenumberIDCanonically();
  
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
  
  // Loop over node IDs from 0 to Na-1, where Na is the # of active nodes.
  // We assume the layer file is ordered by ID!
  for( int k=0; k<meshPtr_->getNodeList()->getSize(); k++ )
  {
    // Find the node with the current ID #
    cn = ni.GetP( k );
    
    if (0) std::cout << "read lays at node " << cn->getID() << " pid=" << cn->getPermID() << "(x,y)=" << cn->getX() << "," << cn->getY() << std::endl;
    
    // Remove pre-existing layers
    for(i=cn->getNumLayer()-1; i>=0; i-- )
      cn->removeLayer(i);
    
    // Read the number of layers at the current node
    layerinfile >> numl;
    if (0) std::cout << " " << numl << "lays\n";
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
      if (0) std::cout << "  lay" << i << " thickness=" << ditem << std::endl;
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
        if (0) std::cout << "   thick of " << g << "=" << ditem << std::endl;
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
  if(1) std::cout << "tLithologyManager::SetLithologyFromEtchFile here\n";
  
  // Open the etch file
  ifstream etchfile;
  string etchfile_name = infile.ReadString( "ETCHFILE_NAME" );
  etchfile.open( etchfile_name.c_str() );
  if( !etchfile.good() )
    ReportFatalError( "Unable to find or open etchfile" );
  
  // Read through the etch file, assigning new layers accordingly
  // First, find out the number of layers to etch.
  unsigned num_grain_sizes = infile.ReadInt( "NUMG" );
  unsigned num_layers_to_etch;
  etchfile >> num_layers_to_etch;
  if( num_layers_to_etch>1000 )
  {
    cout << "WARNING: The etch file is asking for " << num_layers_to_etch << " layers.\n";
    cout << "Are you sure there isn't an error in the etch file?\n";
  }
  
  // Next, read through layer by layer.
  Etchlayer new_etch_layer;
  new_etch_layer.layer_properties_.setDgradesize( num_grain_sizes );
  std::cout << " new etch seems to allow for " << new_etch_layer.layer_properties_.getDgradesize() << std::endl;
  for( unsigned i=0; i<num_layers_to_etch; i++ )
  {
    double erody, bulk_density, grain_size_proportion, cum_fraction;
    unsigned sedrockflag;
    //tLayer layer_properties;
    
    // Read properties for the current layer: general properties
    if(1) cout << "Reading properties for layer " << i << endl;
    etchfile >> erody;
    new_etch_layer.layer_properties_.setErody( erody );
    etchfile >> bulk_density;
    new_etch_layer.layer_properties_.setBulkDensity( bulk_density );
    etchfile >> sedrockflag;
    tLayer::tSed_t flag = ( sedrockflag>0 ) ? tLayer::kSed : tLayer::kBedRock;
    if(1) cout << "  Sed/rock flag = " << sedrockflag << endl;
    new_etch_layer.layer_properties_.setSed( flag );
    if(1) cout << "  Confirming: sed/rock flag = " << new_etch_layer.layer_properties_.getSed() << endl;
    cum_fraction = 0.0;
    
    // Read properties: grain-size related
    for( unsigned j=1; j<num_grain_sizes; j++ )
    {
      if(1) std::cout << " Reading dgrade info for size " << j << " of " << num_grain_sizes << std::endl;
      etchfile >> grain_size_proportion;   // Note: dgrades are meant to be proportion, not
      if( grain_size_proportion < 0.0 || grain_size_proportion > 1.0 )
        ReportFatalError( "Grain size proportions must be 0 to 1." );
      cum_fraction += grain_size_proportion;
      new_etch_layer.layer_properties_.setDgrade( j, grain_size_proportion );
    }
    new_etch_layer.layer_properties_.setDgrade( 0, 1.0 - cum_fraction );
    
    // Read coefficients for 2D polynomial lower surface of layer, and also
    // bounding polygon if applicable
    etchfile >> new_etch_layer.ax >> new_etch_layer.bx >> new_etch_layer.cx; 
    etchfile >> new_etch_layer.ay >> new_etch_layer.by >> new_etch_layer.cy;
    etchfile >> new_etch_layer.d;
    etchfile >> new_etch_layer.keep_regolith_ >> new_etch_layer.use_bounding_polygon_;
    
    if(1) new_etch_layer.TellData();

    if( new_etch_layer.use_bounding_polygon_ )
    {
      etchfile >> new_etch_layer.layer_is_inside_poly_;
      if(1) std::cout << "layer in poly = " << new_etch_layer.layer_is_inside_poly_ << std::endl;
      int npoints;
      etchfile >> npoints;
      if( npoints < 3 || npoints>1000 )
      {
        std::cout << "Your etchfile asks for a polygon with " << npoints 
        << " points, which seems unlikely.\n";
        ReportFatalError( "Error in etchfile format" );
      }
      if(1) std::cout << "Using bounding polygon with " << npoints << " points.\n";
      new_etch_layer.px.resize( npoints );
      new_etch_layer.py.resize( npoints );
      for( unsigned j=0; j<npoints; j++ )
      {
        etchfile >> new_etch_layer.px[j] >> new_etch_layer.py[j];
      }
    }
    if(1) new_etch_layer.TellData();
    
    // Do some error checking.
    // Some potential signs of a problem with the etch file format:
    //  - if the depth parameter d is below the center of the earth (6400 km)
    if( new_etch_layer.d < -6400000.0 )
    {
      std::cout << "I'm reading the depth to the base of the etch layer as " << new_etch_layer.d << ", which seems improbable!\n";
      ReportFatalError( "Error in etch file format." );
    }
    
    // "etch" this layer into the existing stratigraphy
    EtchLayerAbove2DSurface( new_etch_layer );
  }
  
}

    
/**************************************************************************/
/**
 **  SetRockErodyFromFile
 **
  */
/**************************************************************************/
void tLithologyManager::
SetRockErodyFromFile( const tInputFile &infile )
{
  // Open the etch file
  ifstream erodyfile;
  string erodyfile_name = infile.ReadString( "ERODYFILE_NAME" );
  erodyfile.open( erodyfile_name.c_str() );
  if( !erodyfile.good() )
    ReportFatalError( "Unable to find or open erodyfile" );
  
  // Read the contents
  unsigned nnodes = meshPtr_->getNodeList()->getSize();
  vector<double> erody( nnodes );
  for( unsigned i=0; i<nnodes; i++ )
  {
	int node_id;
	double node_erody;
    if( erodyfile.eof() )
      ReportFatalError( "Reached EOF while reading erody file" );
	erodyfile >> node_id >> node_erody;
	erody[node_id] = node_erody;
	  //erodyfile >> erody[i];
  }
  
  // Set the erodibility field
  SetRockErodibilityValuesAtAllDepths( erody );
  
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
EtchLayerAboveHeightAtNode( double new_layer_base_height, tLNode * node,
                           tLayer &layer_properties, bool keep_regolith )
{
  bool layer_found = false;
  //int numlay = node->getNumLayer();
  tListIter<tLayer> li ( node->getLayersRefNC() );
  tLayer * curlay;
  double current_layer_base_height = node->getZ(); // Elevation of layer's base
  int layer_number = 0;
  curlay = li.FirstP();
  
  // Find which layer the new layer cuts through, if any.
  // To do this, we'll start at the top layer and work downward until we
  // either find an existing layer whose base is lower than that of the new
  // layer to be added, or we hit the bottom of the layer stack.
  while( !layer_found && !li.AtEnd() )
  {
    if(0) std::cout << "Checking layer " << layer_number << std::endl;
    
    // To get the height of the current layer's base, we take its current
    // value (which starts at the land surface) and subtract the thickness
    // ("depth") of the current layer.
    current_layer_base_height -= curlay->getDepth();
    
    // If the the base of the new layer to be added is at a higher
    // elevation than the base of the current layer, then the new layer
    // will cut through the current layer, so we flag that we've found the
    // layer we were looking for. Otherwise, we move down to the next
    // layer in the stack and continue the loop.
    if( new_layer_base_height > current_layer_base_height )
    {
      layer_found = true;
    }
    else
    {
      curlay = li.NextP();
      layer_number++;
    }
    
  }
  
  // TODO: IMPLEMENT PRESERVATION OF REGOLITH
  
  // If the new layer doesn't cut through an existing one, then it must cut 
  // *below* (or right at the base of) the layers, so we'll need to insert a
  // a layer at the bottom of the stack
  if( !layer_found )
  {
    if(0) std::cout << "Did not find a layer at depth "
      << new_layer_base_height << " at node " << node->getPermID() << std::endl;
    
    // Create a new layer to be copied and added at the bottom of the current
    // node's layer list. We'll initialize it with the properties of the 
    // current bottom layer, then modify these below.
    tLayer new_bottom_layer( *curlay );

    // Set the layer's thickness such that its base lies 
    // at new_layer_base_height.
    new_bottom_layer.setDepth( current_layer_base_height - new_layer_base_height );
    
    // Add a copy of the new layer to the bottom of the layer list
    node->InsertLayerBack( new_bottom_layer );
    
    // Set the new bottom layer as the "current" layer, 
    // and set its base height.
    curlay = li.LastP();  // our new layer is the last on the list
    current_layer_base_height = new_layer_base_height;
  }

  // At this point, the base of the current layer could be either equal to or 
  // deeper (smaller) than the base of the new etch layer. If it is *deeper*, 
  // then we need to split it in two, adding a new layer above. (We only call
  // the base "deeper" if it is deeper by more than a small tolerance value 
  // equal to one tenth of the standard layer thickness "maxregdepth").
  if( (new_layer_base_height-current_layer_base_height) > 0.1*node->getMaxregdep() )
  {
    if(0) std::cout << "Now splitting layer ...\n";
    // Our layer splitting function will use the following algorithm:
    //  - Create a new layer, copying from the current layer
    //  - Assign the bottom (new) layer the top layer's original thickness 
    //    minus its new thickness
    //  - Insert the new layer in the list after the current layer
    //  - Assign the top layer its new thickness
    
    // Create a new layer, with default properties identical to those of the 
    // current layer. This layer represents the "bottom" portion of the layer
    // we are slicing in two.
    tLayer new_layer( *curlay );
    
    // Set its thickness: its top is at new_layer_base_height, and its bottom
    // is at current_layer_base_height, so its thickness is the difference 
    // between these two.
    double layer_thickness = new_layer_base_height - current_layer_base_height;
    new_layer.setDepth( layer_thickness );
    if(0) std::cout << "Setting thickness of new layer to " << layer_thickness << std::endl;
    
    // Insert this layer below the current one.
    node->getLayersRefNC().insertAtNext( new_layer, li.NodePtr() );
    
    // Reduce the thickness of the "top" portion (which is pointed to by 
    // (curlay).
    layer_thickness = curlay->getDepth() - layer_thickness;
    assert( layer_thickness > 0.0 );
    curlay->setDepth( layer_thickness );
    if(0) {
      std::cout << "The remaining upper layer has been reduced to " << layer_thickness << std::endl;
      std::cout << "To confirm: " << curlay->getDepth() << std::endl;
      std::cout << "Top layer thickness: " << node->getLayerDepth(0) << std::endl;
    }
  }
  
  // Set the properties of the current layer and all layers above it to those
  // specified in layer_properties. To do this, we'll walk down from the top 
  // layer. We know that the lowest layer that will have our new "etch" 
  // properties is number "layer_number", so that's our count.
  curlay = li.FirstP();
  for( int i=0; i<=layer_number; i++, curlay=li.NextP() )
  {
    if(0) {
      std::cout << "Layer " << i << " has thickness " << curlay->getDepth() << std::endl;
      std::cout << "Via the node: " << node->getLayerDepth(i) << std::endl;
    }
    
    // Here we apply the properties of this "etch layer" to the current
    // layer. We don't use the overloaded assignment operator, because that
    // would also apply the thickness.
    curlay->setErody( layer_properties.getErody() );
    curlay->setCtime( 0.0 );
    curlay->setRtime( 0.0 );
    curlay->setEtime( 0.0 );
    curlay->setBulkDensity( layer_properties.getBulkDensity() );
    curlay->setSed( layer_properties.getSed() );
    
    // Now we're applying the grain size information. We assume that
    // layer_properties stores the PROPORTION not the THICKNESS of each
    // size class, so we have to translate.
    double thickness = curlay->getDepth();
    curlay->setDgradesize( layer_properties.getDgradesize() );
    for( size_t j=0; j<layer_properties.getDgradesize(); j++ )
    {
      double proportion_of_this_size = layer_properties.getDgrade(j);
      if(0) std::cout << "Setting dgrade " << j << " to " << proportion_of_this_size << " times " << thickness << " = " << proportion_of_this_size*thickness << std::endl;
      curlay->setDgrade( j, proportion_of_this_size*thickness );
    }
    if(0) std::cout << "Now at end of loop we have layer thickness " << curlay->getDepth() << " and " << node->getLayerDepth(i) << std::endl;
  }

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
EtchLayerAbove2DSurface( Etchlayer &etchlay )
{
  if(0) std::cout << "tLithologyManager::EtchLayerAbove2DSurface here\n";
  
  tLNode * cn;
  tMesh<tLNode>::nodeListIter_t ni( meshPtr_->getNodeList() );
  
  // Loop over all active (interior) nodes
  for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
  {
    // Find out whether the node falls within the layer's bounding polygon,
    // if there is one.
    bool layer_is_present_here = true;
    if( etchlay.use_bounding_polygon_ )
      if( !PointInPolygon( etchlay.px, etchlay.py, cn->getX(), cn->getY() ) )
         layer_is_present_here = false;
    
    // If indeed the layer should be added here, go ahead and add it
    if( layer_is_present_here )
    {
      // Find out the height of the base of the layer at this node by
      // applying the polynomial coefficients to this particular node
      // location.
      double x = cn->getX();
      double y = cn->getY();
      double new_layer_base_height = etchlay.ax*x*x*x +
                              etchlay.bx*x*x +
                              etchlay.cx*x +
                              etchlay.ay*y*y*y +
                              etchlay.by*y*y +
                              etchlay.cy*y +
                              etchlay.d;

      // If the base of the layer falls below the ground surface, then
      // go ahead and "etch" it in.
      if(0) std::cout << "About to etch at " << new_layer_base_height << std::endl;
      if( new_layer_base_height < cn->getZ() )
        EtchLayerAboveHeightAtNode( new_layer_base_height, cn, 
                                    etchlay.layer_properties_, 
                                    etchlay.keep_regolith_ );
      if(0) std::cout << "Our top layer is now " << cn->getLayerDepth(0) << std::endl;
      
    }
  }
    
}

/**************************************************************************/
/**
 */
/**************************************************************************/

