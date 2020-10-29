#ifndef CHILD_CHILDINTERFACE_CHILD_H_
#define CHILD_CHILDINTERFACE_CHILD_H_

#include <string>

#include "../Erosion/erosion.h"
#include "../Mathutil/mathutil.h"
#include "../tEolian/tEolian.h"
#include "../tFloodplain/tFloodplain.h"
#include "../tLithologyManager/tLithologyManager.h"
#include "../tMesh/tMesh.h"
#include "../tOutput/tOutput.h"
#include "../tRunTimer/tRunTimer.h"
#include "../tStorm/tStorm.h"
#include "../tStratGrid/tStratGrid.h"
#include "../tStreamMeander/tStreamMeander.h"
#include "../tStreamNet/tStreamNet.h"
#include "../tUplift/tUplift.h"
#include "../tVegetation/tVegetation.h"
#include "../tWaterSedTracker/tWaterSedTracker.h"


class Child {
  public:
    void Initialize (std::string);
    double RunOneStorm();
    void Run(double run_duration);
    void CleanUp();
    void MaskNodesBelowElevation(double elev);  // Mask out nodes below elev (eg. sea level)
    void CopyNodeElevations (double * const dest);
    void CopyNodeErosion (double * const dest);
    void CopyNodeDischarge (double * const dest);
    void CopyNodeSedimentFlux (double * const dest);
    void SetNodeElevations (const double * src);
    void SetNodeUplift (const double * src);
  
    bool initialized;  // Flag indicated whether model has been initialized
    bool optNoDiffusion;  // Option to turn off diffusive processes (default to false)
    bool optNoFluvial;  // Option to turn off fluvial processes (default to false)
    bool optNoUplift;  // Option to turn off uplift (default to false) 
    bool optDetachLim;  // Option for detachment-limited erosion only
    bool optFloodplainDep;  // Option for floodplain (overbank) deposition
    bool optLoessDep;  // Option for eolian deposition
    bool optVegetation;  // Option for dynamic vegetation cover
    bool optFire;  // Option for fire object within tVegetation
    bool optForest;  // Option for forest object within tVegetation
    bool optMeander;  // Option for stream meandering
    bool optDiffuseDepo;  // Option for deposition / no deposition by diff'n
    bool optStratGrid;  // Option to enable stratigraphy grid
    bool optTrackWaterSedTimeSeries;  // Option to record timeseries Q and Qs
    bool optNonlinearDiffusion;  // Option for nonlinear creep transport
    bool optDepthDependentDiffusion;  // Option for depth dependent creep transport
    bool optLandslides;  // Option for landsliding
    bool opt3DLandslides;  // Option for determining which landslide function to use
    bool optChemicalWeathering;  // Option for chemical weathering
    bool optPhysicalWeathering;  // Option for physical weathering
    bool optStreamLineBoundary;  // Option for converting streamlines to open boundaries

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

