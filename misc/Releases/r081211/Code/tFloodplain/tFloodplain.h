//-*-c++-*-

/**************************************************************************/
/**
**  @file tFloodplain.h
**  @brief header file for class tFloodplain and its helper class
**         tFloodNode.
**
**  Class tFloodplain simulates overbank deposition in the CHILD model,
**  using a modified version of Howard's (1996) diffusive overbank
**  deposition model. Implementation of the model is described in
**  tFloodplain.cpp. The class consists of parameters for the model,
**  a pointer to the mesh, a constructor to read in parameters, and
**  a DepositOverbank function to implement deposition during a given
**  storm event. Helper class tFloodNode stores info for a "flood node"
**  (one containing a channel big enough to generate a sedimentologically
**  significant flood), and is used to construct a list of flood nodes.
**
**  (Created 1/99 by GT)
**
**  $Id: tFloodplain.h,v 1.21 2004/06/16 13:37:30 childcvs Exp $
*/
/**************************************************************************/

#ifndef TFLOODPLAIN_H
#define TFLOODPLAIN_H

#include <assert.h>
#include "../tMesh/tMesh.h"
#include "../tLNode/tLNode.h"
#include "../tInputFile/tInputFile.h"
#include "../tTimeSeries/tTimeSeries.h"

#define kVeryFar 1.0e12


/**************************************************************************/
/**
**  @class tMainChannelDriver
**
**  Class tMainChannelDriver controls the altitude behavior of the main
**  channel within the floodplain. It is used, at the user's option, to
**  control the gradient and relative altitude of the main channel as an
**  imposed boundary condition. It is only used when the user switches
**  on the option to impose the main channel height/slope as a boundary
**  condition, as opposed to allowing it to evolve as a function of
**  erosion, deposition, and baselevel change.
**    A tMainChannelDriver object sits within a tFloodplain object.
**
**    Created May 2003 by GT based on earlier implementation within a
**    a special version of the main file (gamain1.cpp) designed by
**    GT and NMG.
*/
/**************************************************************************/
class tMainChannelDriver
{
  tMainChannelDriver(const tMainChannelDriver&);
  tMainChannelDriver& operator=(const tMainChannelDriver&);
  tMainChannelDriver();

public:
  tMainChannelDriver( const tInputFile &infile );
  void UpdateMainChannelElevation( double tm, tLNode * inletNode );
  void RaiseBanks( double, tLNode *, tLNode*,double);

private:
  tTimeSeries InletElevationVariation; // Elev oscillation [L]
  double drop;      // Channel elevation drop from top to bottom of valley [L]

  int num_grnsize_fractions;  // No. of grain size fractions
  double kdb;             // depth-disch coeff (lumped; see tFloodplain.cpp)
  double mqs;             // depth-disch at-a-station exponent
  double mqbmqs;          // bankfull minus at-a-station exponents
  std::ofstream meanderfile;   // file in which every timestep some information of the channel is written
};


/**************************************************************************/
/**
**  @class tFloodplain
**
**  Class tFloodplain contains the parameters for the modified Howard
**  floodplain deposition model, as well as a constructor to read in
**  parameters from an input file and a DepositOverbank function to
**  implement the model.
**
**  Modifications:
**    - 3/13/99 GT added deparr variable which is used by DepositOverbank
**              (could be a local var but this saves having to reinit
**              it with every call to tFloodplain)
**    - May 2003 GT added code to implement control of main channel
**              elevations as a boundary condition.
**
*/
/**************************************************************************/
class tFloodplain
{
  tFloodplain(const tFloodplain&);
  tFloodplain& operator=(const tFloodplain&);
  tFloodplain();

public:
  tFloodplain( const tInputFile &infile, tMesh<tLNode> *mp );
  ~tFloodplain();
  void DepositOverbank( double precip, double delt, double ctime );
  bool OptControlMainChan() const;
  void UpdateMainChannelHeight( double tm, tLNode * inletNode );
  double FloodplainDh(double, double, tLNode * ) const;		// Howard, 1992
  double FloodplainDh2(double , double, double , double ) const;// Gross and Small, 1998
  double getSuspendedConcentration(double) const;
  double ConcentrationToHeight(double, tLNode *, double ) const;

private:
  tTimeSeries fpmuVariation;   // "mu" parameter of Howard model, value dependent of time
  int fpmode;                  // 1) Howard   2) modified form based on suspension in column and where C~kQcha.
  double fplamda;         // "lamda" (distance-decay) parameter
  double kdb;             // depth-disch coeff (lumped; see tFloodplain.cpp)
  double event_min;       // bankfull event precip rate
  double drarea_min;      // min drainage area for a "flood node"
  double mqs;             // depth-disch at-a-station exponent
  double mqbmqs;          // bankfull minus at-a-station exponents
  tMainChannelDriver *chanDriver;   // if user wants to control chan elev
  tMesh<tLNode> *meshPtr; // ptr to mesh
  tArray<double> deparr;  // depth deposited (# grn size; all but 1st=0)
  tArray<double> deparrRect; // same but for rectangular grid
  bool optControlMainChan;  // option to treat chan elev as boundary cond
};


/**************************************************************************/
/**
**  @class tFloodNode
**
**  @brief Data record used in constructing a list of flood nodes.
**  Contains a pointer to the node and the water surface height.
**
*/
/**************************************************************************/
class tFloodNode
{
    tFloodNode& operator=(const tFloodNode&);

    friend class tFloodplain;
public:
    tFloodNode() : nodePtr(0), wsh(-1.) {}
    tFloodNode(tLNode *p_, double w_) : nodePtr(p_), wsh(w_) {}
    tFloodNode(const tFloodNode& c) : nodePtr(c.nodePtr), wsh(c.wsh) {}
private:
    tLNode *nodePtr; // ptr to flood node
    double wsh;      // water surface height at flood node
};

#endif
