//-*-c++-*-

/**************************************************************************/
/**
**  @file TravelDis.h
**
**  @brief Header file for the TravelDis class.
**
*/
/**************************************************************************/

#ifndef TravelDis_H
#define TravelDis_H

#include <vector>
//#include <array>
#include <string>
#include <algorithm>
#include "../tInputFile/tInputFile.h"
#include "../tMesh/tMesh.h"
#include "../tLNode/tLNode.h"
#include "../tStreamNet/tStreamNet.h"



class TravelDis
{
public:
    
    TravelDis();

  // Constructor for basic initialization
  TravelDis( tInputFile& , tMesh<tLNode>* );

  // Destructor
  ~TravelDis();
  

//    
    void TravelDis_InsideStorm(tMesh<tLNode> *mesh, double tot_time, double storm_time, double t_left, double flux_duration);
    
    void CreateAndOpenWaterAndSedimentOutputFiles();
    
    void DefineWshed(tMesh<tLNode> *mesh, vector<tLNode *> &List_of_ptrs_to_All_nodes, vector<int> &List_numUpstreamNodes);
    
    void WriteToProbDataFile(double tot_time, double time_sofar, double flux_duration, vector<tLNode *> &List_of_ptrs, vector<int> &List_numpts);
    
    void WriteToNodeFile(double tot_time, vector<tLNode *> &List_of_ptrs_to_All_nodes, vector<int> &List_numUpstreamNodes);


private:

    
  std::vector<tLNode *> tracking_node_list_;  // List of ptrs to nodes to track
  std::vector<ofstream *> output_file_list_;  // List of ptrs to file I/O
  std::string output_file_base_name_;         // Base name for output files
    int stormcounter;                           // storm counter
    int GridSpacing;
    string input_string;
    double /**kb, mb, nb,**/ outputperiod;
    int outputcounterTD;
  bool active;

 
};


#endif
