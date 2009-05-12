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
		// Example 1: using "Run" method, and setting run duration to zero so model reads duration from input file
		myChildInterface.Run( 0 );
		
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
