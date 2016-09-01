#ifndef CHILD_CHILDINTERFACE_CHILD_H_
#define CHILD_CHILDINTERFACE_CHILD_H_

#include <string>
#include <sstream>
#include <vector>
#include "../trapfpe.h"
#include "../Inclusions.h"
#include "../tFloodplain/tFloodplain.h"
#include "../tStratGrid/tStratGrid.h"
#include "../tEolian/tEolian.h"
#include "../tOption/tOption.h"
#include "../tWaterSedTracker/tWaterSedTracker.h"
#include "../tMeshList/tMeshList.h"
#include "../tLithologyManager/tLithologyManager.h"

class Child {
  public:
    void Initialize (std::string);
    double RunOneStorm();
    void Run(double run_duration);
    void CleanUp();
    void MaskNodesBelowElevation(double elev);  // Mask out nodes below elev (eg, sea level)
    void CopyNodeElevations (double * const dest);
    void CopyNodeErosion (double * const dest);
    void CopyNodeDischarge (double * const dest);
    void CopyNodeSedimentFlux (double * const dest);
    void SetNodeElevations (const double * src);
    void SetNodeUplift (const double * src);
  
    bool initialized;      // Flag indicated whether model has been initialized
    bool optNoDiffusion,   // Option to turn off diffusive processes (default to false)
         optNoFluvial,        // Option to turn off fluvial processes (default to false)
         optNoUplift,       // Option to turn off uplift (default to false) 
         optDetachLim,      // Option for detachment-limited erosion only
         optFloodplainDep,  // Option for floodplain (overbank) deposition
         optLoessDep,       // Option for eolian deposition
         optVegetation,     // Option for dynamic vegetation cover
         optFire,           // Option for fire object within tVegetation
         optForest,         // Option for forest object within tVegetation
         optMeander,        // Option for stream meandering
         optDiffuseDepo,    // Option for deposition / no deposition by diff'n
         optStratGrid,      // Option to enable stratigraphy grid
         optTrackWaterSedTimeSeries,  // Option to record timeseries Q and Qs
         optNonlinearDiffusion, // Option for nonlinear creep transport
         optDepthDependentDiffusion, // Option for depth dependent creep transport
         optLandslides, // Option for landsliding
         opt3DLandslides, // Option for determining which landslide function to use
         optChemicalWeathering, // Option for chemical weathering
         optPhysicalWeathering; // Option for physical weathering

    bool optStreamLineBoundary; // Option for converting streamlines to open boundaries
         tRand *rand;             // -> random number generator
         tMesh<tLNode> *mesh;        // -> mesh object
         tLOutput<tLNode> *output;   // -> output handler
         tStorm *storm;              // -> storm generator
         tStreamNet *strmNet;	    // -> stream network module
         tErosion *erosion;          // -> erosion module
         tUplift *uplift;            // -> uplift/baselevel module
         tRunTimer *time;             // -> run timer
         tWaterSedTracker water_sed_tracker_;   // Water and sediment tracker
	       tLithologyManager lithology_manager_;  // Lithology manager
         tVegetation *vegetation;  // -> vegetation object
         tFloodplain *floodplain;  // -> floodplain object
         tStratGrid *stratGrid;     // -> Stratigraphy Grid object
         tEolian *loess;           // -> eolian deposition object
         tStreamMeander *strmMeander; // -> stream meander object
};

#endif

