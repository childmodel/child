//-*-c++-*-

/**************************************************************************/
/**
**  @file TravelDis.cpp
**
**  @brief Implementation of the TravelDis class
**
**  TravelDis calculates the probability that the sediment depositing at
**  target cells at the end of a storm was sourced from each of the cells
**  upstream and the minimum travel time since the material was sourced
**
**  MP, Jun 2013
**
*/
/**************************************************************************/

#include <vector>
#include <sstream>
#include <ctime>
#include <list>

#include "TravelDis.h"

using namespace std;

/**************************************************************************/
/**
**  Basic constructor - runs at the start of the storm to find the tracking nodes
**
*/
/**************************************************************************/
TravelDis::TravelDis()
: active(false)
{

}

TravelDis::TravelDis( tInputFile &inputFile, tMesh<tLNode>* nPtr)
{
    
    stormcounter=0;
    outputcounterTD=1;
    
    GridSpacing = inputFile.ReadInt( "GRID_SPACING" );
    
    // Remember the base name for the sed/water flux output files
    output_file_base_name_ = inputFile.ReadString( "TRAVELDIS_OUTFILENAME" );
    
    int number_of_nodes_to_track = inputFile.ReadInt( "NUMBER_OF_TARGET_NODES" );
    tracking_node_list_.resize( number_of_nodes_to_track, NULL );
    input_string = inputFile.ReadString( "COORDS_OF_TARGET_NODES" );
    
    CreateAndOpenWaterAndSedimentOutputFiles();
//    
//    kb = inputFile.ReadItem( kb, "KB" );
//    mb = inputFile.ReadItem( mb, "MB" );
//    nb = inputFile.ReadItem( nb, "NB" );
    
    outputperiod=inputFile.ReadItem( outputperiod, "OPINTRVL");
    
}



///**************************************************************************/
///**
//**  Destructor
//**
//**  Closes and deletes any remaining ofstream objects.
//*/
///**************************************************************************/
TravelDis::
~TravelDis()
{
  if(0) cout << "TravelDis destructor" << endl;
  
  for( unsigned i=0; i<output_file_list_.size(); i++ )
  {
    output_file_list_[i]->close();
    delete output_file_list_[i];
  }
}


/**************************************************************************/
/**
 **  TravelDis_InsideStorm
 **
 **  This method splits the work that has to be done at every subtimestep in a storm
 **
 **  Called by erosion.cpp
 */
/**************************************************************************/
void TravelDis::TravelDis_InsideStorm( tMesh<tLNode> *mesh, double tot_time, double storm_time, double t_left, double flux_duration)
{
    
    double time_sofar = storm_time-t_left;
    
    if (tot_time<outputperiod) outputcounterTD=1;
    
    
    if (tot_time>=outputcounterTD*outputperiod && tot_time>0.0){
        
        
    vector<tLNode *> List_of_ptrs_to_All_nodes; // list of pointers to elements of tlNode - all target nodes in sequence
    List_of_ptrs_to_All_nodes.reserve(500);
    
    vector<int> List_numUpstreamNodes; // list of number of stream nodes, one per target node
    List_numUpstreamNodes.reserve(tracking_node_list_.size());
    
    DefineWshed(mesh, List_of_ptrs_to_All_nodes,List_numUpstreamNodes);
    
    if(t_left==storm_time){
        
        WriteToNodeFile(tot_time,List_of_ptrs_to_All_nodes,List_numUpstreamNodes);
        
    }
    
    
    WriteToProbDataFile(tot_time, time_sofar, flux_duration, List_of_ptrs_to_All_nodes,List_numUpstreamNodes);
    
        if(t_left<=flux_duration){
            outputcounterTD++;
        }
        
        
    }
}


/**************************************************************************/
/**
 **  DefineWshed
 **
 **  This method finds the nodes upstream of each of the target nodes
 ** List_of_ptrs_to_All_nodes  includes the points upstream of all target nodes
 ** Use List_numUpstreamNodes to separate
 */
/**************************************************************************/
#define kMaxSpokes 100
#define kBigNumber 1e10
 void TravelDis::DefineWshed( tMesh<tLNode> *mesh, vector<tLNode *> &List_of_ptrs_to_All_nodes, vector<int> &List_numUpstreamNodes)
    {
        
    stringstream input_stringstream;
        
        double tempx, tempy;
        int number_of_nodes_to_track=tracking_node_list_.size();
        
        input_stringstream.str( input_string );   // use the stringstream to get #s from the input line

        vector<double> x; x.reserve(number_of_nodes_to_track);
        vector<double> y; y.reserve(number_of_nodes_to_track);
        
        
        for( int i=0; i<number_of_nodes_to_track; i++ )
          {
           if( input_stringstream.eof() )
             ReportFatalError( "Reached end of line in reading COORDS_OF_TARGET_NODES:\n make sure there is at least one x and y value for each target node." );
              input_stringstream >> tempx;
              x.push_back(tempx);
           if( input_stringstream.eof() )
              ReportFatalError( "Reached end of line in reading COORDS_OF_TARGET_NODES:\n make sure there is at least one x and y value for each target node." );
              
            input_stringstream >> tempy;
              y.push_back(tempy);
              
          }


        // Find the closest high DA node to each (x,y) pair and put it on the tracking list.
        // Isolate all nodes within 5 GRID_SPACING of the node of interest. Pick the one with the highest DA
        vector<double> min_distance_squared( number_of_nodes_to_track, kBigNumber );   // Min dist^2 found so far
        tLNode *current_node;
        tMesh< tLNode >::nodeListIter_t ni( mesh->getNodeList() );
        
        vector<tLNode *> closepoints;
        double DA;
        double maxDA=0;
        
        double searchdistance;
        
        searchdistance=15*GridSpacing;
        
        //if(searchdistance>200) searchdistance=200;
        
        
        for( int i=0; i<number_of_nodes_to_track; i++ ) {
            
            
            closepoints.clear();
            maxDA=0;
            
            for( current_node=ni.FirstP(); ni.IsActive(); current_node=ni.NextP() ) {
                double dist_squared = pow(x[i]-current_node->getX(),2)+pow(y[i]-current_node->getY(),2);
                
            
                
                if( dist_squared <= searchdistance )
                {
                    closepoints.push_back(current_node);   // ... and set this as the closest node.
                }
            }
            
            
            for(int j=0; j<closepoints.size(); j++){
                DA=closepoints[j]->getDrArea();
                if(DA>=maxDA){
                    maxDA=DA;
                    tracking_node_list_[i]=closepoints[j];
                }
                
            }
            
            
        }
        
    
    
    double slp=0;               // slope of edge
    tLNode *currentnode;
    tLNode *newnode;           // ptr to new downstream node
    
    tEdge *firstedg;   // ptr to first edg
    tEdge *currentedg;     // pointer to current edge
        
    int ctr;
        
        double maxslope;
    
    vector<tLNode *> List_of_ptrs_to_nodes; // list of pointers to elements of tlNode - one target node
    vector<tLNode *> List_of_ptrs_to_nodes_done; // list of pointers to elements of tlNode that have already been searched
    vector<tLNode>::iterator it;
    vector<int> List_of_IDs; // List of node IDs for checking
    
    vector<int>::iterator CheckIter;
        
        
        for( int j=0; j<number_of_nodes_to_track; j++){ // for each target node
            currentnode = tracking_node_list_[j];  //  pointer to target node
            
            maxslope=0.0;
            
            // add ptr of target node to the list
            List_of_ptrs_to_nodes_done.clear();
            List_of_ptrs_to_nodes.clear();
            List_of_IDs.clear();
            
            List_of_ptrs_to_nodes.push_back(currentnode); // List of nodes to search
            List_of_IDs.push_back(currentnode->getPermID());
            
            unsigned i=0;
                while(i<List_of_ptrs_to_nodes.size())
                {
                    currentnode=List_of_ptrs_to_nodes[i];
                    
                    if(static_cast<tNode *>(currentnode)->getBoundaryFlag() == 0 )
                    {
                        
                    
                            // get the first edge
                        firstedg =  currentnode->getEdg();
            
                        // if it's a boundary, move to the next edge
                        if( unlikely(firstedg == 0) ) {
                            currentedg = firstedg->getCCWEdg();
                        }
                        
                        else { // check the slope
                            
                            slp = firstedg->getSlope();
                        
                            if(slp<maxslope){
                                
                                newnode = static_cast<tLNode *>(firstedg->getDestinationPtrNC());
                                    
                                CheckIter=find(List_of_IDs.begin(),List_of_IDs.end(),newnode->getPermID());
                        
                                if(CheckIter!=List_of_IDs.end()){
                                }
                                else {
                                    List_of_ptrs_to_nodes.push_back(newnode);
                                    List_of_IDs.push_back(newnode->getPermID());
                                }
                            
                                }
                                
                            
                            currentedg = firstedg->getCCWEdg();	// Go to the next counter clockwise edge
            
                        }
                        
                        ctr = 0;
                        
            
                        while( currentedg!=firstedg ) // while we haven't returned to the current edge
                        {
                           
                           // if it's a boundary, move to the next edge
                           if( unlikely(currentedg == 0) ) {
                               currentedg = currentedg->getCCWEdg();
                           }
                           else { // check the slope
                               
                               slp = currentedg->getSlope();
                               
                               if(slp<maxslope){
                                   newnode = static_cast<tLNode *>(currentedg->getDestinationPtrNC());
                                   CheckIter=find(List_of_IDs.begin(),List_of_IDs.end(),newnode->getPermID());
                                   
                                   if(CheckIter!=List_of_IDs.end()){
                                   }
                                   else {
                                       List_of_ptrs_to_nodes.push_back(newnode);
                                       List_of_IDs.push_back(newnode->getPermID());
                                   }
                                   
                                   
                               }
                               
                               currentedg = currentedg->getCCWEdg();	// Go to the next counter clockwise edge
                               
                           }
                            
                            ctr++;
            
                            if( unlikely(ctr>kMaxSpokes) ) // Make sure to prevent endless loops
                                
                            {
                                std::cerr << "Mesh error: node " << currentnode->getID()
                                << " going round and round"
                                << std::endl;
                                ReportFatalError( "Bailing out of DefineWshed()" );
                            }
                           
                        } // end of loop through spokes
                        
                        
                        // add the point just searched to the list of points done
                        List_of_ptrs_to_nodes_done.push_back(currentnode);
                        
                    }
                        
                        i++;
                    } // end of loop through list of nodes
                    // save the list of all pointers
            
            
                    List_of_ptrs_to_All_nodes.insert(List_of_ptrs_to_All_nodes.end(), List_of_ptrs_to_nodes_done.begin(),List_of_ptrs_to_nodes_done.end());
                    List_numUpstreamNodes.push_back(List_of_ptrs_to_nodes_done.size());
            
            
            
        } // end of loop through target nodes
        
}
#undef kMaxSpokes
#undef kBigNumber

/**************************************************************************/
/**
 ** CreateAndOpenWaterAndSedimentOutputFiles - creates and opens the text files
 ** one for node information and one for data used for probability
*/
/**************************************************************************/
void TravelDis::
CreateAndOpenWaterAndSedimentOutputFiles()
{
  
    stringstream ss1;
    ss1 << output_file_base_name_;
    ss1 << "_probabilities.probnodes";
    output_file_list_.push_back( new ofstream( ss1.str().c_str() ) );
    
    stringstream ss2;
    ss2 << output_file_base_name_;
    ss2 << "_probabilities.probdata";
    output_file_list_.push_back( new ofstream( ss2.str().c_str() ) );
      
    if( !output_file_list_[0]->good() )
    {
      cout << "When trying to create output file '" << ss1.str() << endl;
      ReportFatalError( "Unable to create file." );
    }
    // X Y Z B ID IDdownstream Lspoke Qflowedge Sspoke
    //*output_file_list_[0] << "X Y Z B ID IDdownstreamnode Lspoke Qflowedge Slope" << endl;
    
    if( !output_file_list_[1]->good() )
    {
        cout << "When trying to create output file '" << ss2.str() << endl;
        ReportFatalError( "Unable to create file." );
    }
    //*output_file_list_[1] << "Header of prob data file" << endl;
  
}

/**************************************************************************/
/**
 ** WriteToNodeFile - writes the information about the nodes into the .nodes file
 */
/**************************************************************************/
void TravelDis::WriteToNodeFile(double tot_time, vector<tLNode *> &List_of_ptrs, vector<int> &List_numpts)
{
    
    if( !output_file_list_[0]->good() )
    {
        cout << "Node File disappeared!" << endl;
        ReportFatalError( "Unable to write to file." );
    }
    
    
    int numTargetNodes = List_numpts.size();
    int cumsumNodes=0;
    
    *output_file_list_[0] << "STORMINFO" << endl;
    *output_file_list_[0] <<  tot_time << " " << stormcounter << " " << numTargetNodes << endl;
    
    tEdge *firstedg, *currentedg, *DSspoke;
    double edgslope;
    double minslope;
    
    
    tNode *DSnode;
    
    
    for(int i=0; i<numTargetNodes; i++){ // For each target node
        
        
        // Subset the ptrs to just the wshed of this target node
        int first=cumsumNodes;
        cumsumNodes+=List_numpts[i];
        int last = cumsumNodes-1;
        
        // Target node header
        *output_file_list_[0] << "TARGETNODE" << endl;
        *output_file_list_[0] << i << " " << List_numpts[i] << " " << tracking_node_list_[i]->getPermID() << endl;
        
        //vector<tLNode *>::iterator it = Target_ptrs.begin();
        
        for(int j=first; j<last+1; j++){
            
            minslope=0.0;
            
            firstedg=List_of_ptrs[j]->getEdg();
            
//            if (unlikely(firstedg==0)){
//                currentedg=firstedg->getCCWEdg();
//            }
//            else {
                edgslope=firstedg->getSlope();
                
                if(edgslope>minslope){
                    minslope=edgslope;
                    DSnode=firstedg->getDestinationPtrNC();
                    DSspoke=firstedg;
                }
                
                
                currentedg=firstedg->getCCWEdg();
                
                
//            }
            
            
            while(currentedg!=firstedg){
                
//                if( unlikely(currentedg == 0) ) {
//                    currentedg = currentedg->getCCWEdg();
//                }
//                else {
                
                    edgslope = currentedg->getSlope();
                    
                    if(edgslope>minslope){
                        minslope=edgslope;
                        DSnode=currentedg->getDestinationPtrNC();
                        DSspoke=currentedg;
                    }
                    
                    currentedg = currentedg->getCCWEdg();	// Go to the next counter clockwise edge
                    
//                }
                
    
            }

            
            *output_file_list_[0] << List_of_ptrs[j]->getPermID() << " " << List_of_ptrs[j]->getFlowEdg()->getDestinationPtrNC()/**DSnode**/->getPermID() << " " << static_cast<tNode *>(List_of_ptrs[j])->getX() << " " << static_cast<tNode *>(List_of_ptrs[j])->getY() << " " << static_cast<tNode *>(List_of_ptrs[j])->getZ() << " " << List_of_ptrs[j]->getFlowEdg()/**DSspoke**/ ->getVEdgLen() << endl;

            
        }
        
    }
    
}


/**************************************************************************/
/**
 ** WriteToProbDataFile - writes the information about the nodes into the .probdata file
 */
/**************************************************************************/
void TravelDis::WriteToProbDataFile(double tot_time, double time_sofar,  double flux_duration, vector<tLNode *> &List_of_ptrs, vector<int> &List_numpts)
{
    
    if( !output_file_list_[1]->good() )
    {
        cout << "Probability Data file disappeared!" << endl;
        ReportFatalError( "Unable to write to file." );
    }
    
    int numTargetNodes = List_numpts.size();
    int cumsumNodes=0;
    double Qsin, Qs, depth, DrDt, dt, Area;
    
    double slp, MaterialAvailable, MaterialDeposited, Q;
    
    
    
    
    
    
    *output_file_list_[1] << "STORMINFO" << endl;
    *output_file_list_[1] <<  tot_time << " " << stormcounter << " " << time_sofar << " " << numTargetNodes << endl;
    
    for(int i=0; i<numTargetNodes; i++){ // For each target node
        
        // Subset the ptrs to just the wshed of this target node
        int first=cumsumNodes;
        cumsumNodes+=List_numpts[i];
        int last = cumsumNodes-1;
        
        // Target node header
        *output_file_list_[1] << "TARGETNODE" << endl;
        *output_file_list_[1] << i << " " << List_numpts[i] << " " << tracking_node_list_[i]->getPermID() << endl;
        
        //vector<tLNode *>::iterator it = Target_ptrs.begin();
        
        for(int j=first; j<last+1; j++){
            
            
            DrDt=List_of_ptrs[j]->getDrDt();
            Qsin=List_of_ptrs[j]->getQsin();
            Qs=List_of_ptrs[j]->getQs();
            
            depth= List_of_ptrs[j]->getChanDepth();
            
            dt=flux_duration;
            Area=List_of_ptrs[j]->getVArea();
            
            slp = List_of_ptrs[j]->calcSlope();
            Q=List_of_ptrs[j]->getQ();
            
            
            MaterialAvailable=(Qsin+(-DrDt)) * dt;
            
            //MaterialDeposited=(MaterialAvailable-Qs) * dt / Area;
            MaterialDeposited = List_of_ptrs[j]->getDzDt() * dt;
            
            
            *output_file_list_[1] << List_of_ptrs[j]->getPermID() << " " << MaterialAvailable << " " << MaterialDeposited << " " << Q << " " << slp <<  endl;
            
        }
    }
    
    
}


