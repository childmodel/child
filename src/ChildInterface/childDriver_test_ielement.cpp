/**************************************************************************/
/**
**  childDriver.cpp: This provides a test and example of the CHILD
**  interface.
**
**  July 2008
**
**  For information regarding this program, please contact Greg Tucker at:
**
**     Cooperative Institute for Research in Environmental Sciences (CIRES)
**     and Department of Geological Sciences
**     University of Colorado
**     2200 Colorado Avenue, Campus Box 399
**     Boulder, CO 80309-0399
**
*/
/**************************************************************************/

#include "childInterface.h"

int main( int argc, char **argv )
{
	childInterface myChildInterface;
	
	myChildInterface.Initialize( argc, argv );

	if(1) // make this zero to use "example 2" below
	{
        std::cout << "\nInitialization done. The time remaining for the run is: ";
        std::cout << myChildInterface.GetRemainingRunTime() << std::endl;
        std::cout << "This number should be the same as RUNTIME in the input file.\n";

		// Example 1: using "Run" method, and setting run duration to zero so model reads duration from input file
		myChildInterface.Run( 0 );
		
        std::cout << "Run done. The time remaining for the run is: ";
        std::cout << myChildInterface.GetRemainingRunTime() << std::endl;
        std::cout << "This number should be zero.\n";

		int num_nodes;
		int i, j;
		string result_string;
		vector<int> face_vertex_indices;
		
		std::cout << "TEST OF CHILD'S IElement INTERFACE:\n\n";
		result_string = myChildInterface.getID();
		std::cout << "The model ID is " << result_string << std::endl;
		result_string = myChildInterface.getDescription();
		std::cout << "The model description is " << result_string << std::endl;
		std::cout << "The ElementType is " << myChildInterface.getElementType() << std::endl;
		std::cout << "The version is " << myChildInterface.getVersion() << std::endl;
		std::cout << "The 'element ID' (string) for element 0 is " << myChildInterface.GetElementID( 0 ) << std::endl;
		num_nodes = myChildInterface.getElementCount();
		std::cout << "There are " << num_nodes << " nodes in the mesh.\n";
		for( i=0; i<num_nodes; i++ )
		{
		    int num_faces = myChildInterface.GetFaceCount( i );
		    std::cout << " Node " << i << " has " << myChildInterface.GetVertexCount( i ) << " vertices and ";
			std::cout << num_faces << " faces.\n";
			for( j=0; j<num_faces; j++ )
			{
			   face_vertex_indices = myChildInterface.GetFaceVertexIndices( i, j );
			   std::cout << "  Face " << j << " vertex indices are: right=" << face_vertex_indices[0] << ", left=" << face_vertex_indices[1] << std::endl;
			   std::cout << "   The coords of vertex " << j << " are: x=" << myChildInterface.GetXCoordinate( i, j );
			   std::cout << " y=" << myChildInterface.GetYCoordinate( i, j );
			   std::cout << " z=" << myChildInterface.GetZCoordinate( i, j ) << std::endl;
		    }
		}
		
		// Next we test the custom interface function "GetNodeCoords"
		std::vector<double> node_coords;
		std::vector<double> elevations;
		node_coords = myChildInterface.GetNodeCoords();
		elevations = myChildInterface.GetValueSet( "elevations" );
		std::cout << "\nThe coordinates of the nodes are as follows:\n";
		std::cout << "ID\tX\tY\tZ\tZ from GetValueSet\n";
		for( i=0; i<num_nodes; i++ )
		{
		    std::cout << i << "\t" << node_coords[3*i] << "\t"
		              << node_coords[3*i+1] << "\t"
		              << node_coords[3*i+2] << "\t"
		              << elevations[i] << "\n";
		}
		
		// Now we'll use SetValueSet to re-set the elevations according 
		// to their x,y position
		for( i=0; i<num_nodes; i++ )
		{
		    elevations[i] = 0.01*node_coords[i+num_nodes] + 0.02*node_coords[i+2*num_nodes];
		}
		myChildInterface.SetValueSet( "elevations", elevations );
		
		// Ask CHILD what its new coords are and display
		node_coords = myChildInterface.GetNodeCoords();
		elevations = myChildInterface.GetValueSet( "elevations" );
		std::cout << "\nThe NEW coordinates of the nodes are as follows:\n";
		std::cout << "ID\tX\tY\tZ\tZ from GetValueSet\n";
		for( i=0; i<num_nodes; i++ )
		{
      std::cout << i << "\t" << node_coords[3*i] << "\t"
      << node_coords[3*i+1] << "\t"
      << node_coords[3*i+2] << "\t"
      << elevations[i] << "\n";
		}
    
    std::cout << "\nTest of AdjustElevations() and AdjustInteriorElevations:\n\n";
    std::vector<double> dz( nn, 1.0 );
    myChildInterface.AdjustInteriorElevations( dz );
    std::vector<double> new_elevs = myChildInterface.GetValueSet( "elevation" );
    std::cout << "The interior nodes should now be 1m higher, but not the boundary nodes:\n";
    std::cout << "ID\tOld z\tNew z:\n";
    for( long i=0; i<nt; i++ )
    {
      std::cout << i << "\t";
      std::cout << elevs[i] << "\t" << new_elevs[i] << std::endl;
    }
    myChildInterface.AdjustElevations( dz );
    new_elevs = myChildInterface.GetValueSet( "elevation" );
    std::cout << "Now all nodes should be another 1m higher:\n";
    std::cout << "ID\tOriginal z\tNew z:\n";
    for( long i=0; i<nt; i++ )
    {
      std::cout << i << "\t";
      std::cout << elevs[i] << "\t" << new_elevs[i] << std::endl;
    }
    
		
    
	}
	else
	{
		// Example 2: using "RunOneStorm" 
		double mytime = 0;
		double myrunduration = 100000;
		
		while( mytime<myrunduration )
		{
			mytime = myChildInterface.RunOneStorm();
		}
	}
	
	// Note that calling CleanUp() isn't strictly necessary, as the destructor will automatically clean it
	// up when myChildInterface is deleted ... but it's nice to be able to do this at will (and free up
	// memory)
	myChildInterface.CleanUp();
	
	return 0;
}
