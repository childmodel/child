/**************************************************************************\
**
**  tStreamMeander.h: Header file for class tStreamMeander.
**
**  tStreamMeander objects contain data and functions related to lateral
**  erosion and channel migration across the landscape surface.
**
**    * Contains functions for constructing reaches and calculating lateral
**      channel migration for nodes in those reaches.
**
**    * Contains pointers to a tGrid, a tStreamNet, and a tInputFile.
**
**    * Kicks ass and takes names.
**
**  $Id: tStreamMeander.h,v 1.10 1998-02-12 01:46:28 stlancas Exp $
\**************************************************************************/
#ifndef TSTREAMMEANDER_H
#define TSTREAMMEANDER_H

#include <math.h>
#include <assert.h>
#include <string.h>
#include "../Classes.h"
#include "../Definitions.h"
#include "../Mathutil/mathutil.h"
#include "../tArray/tArray.h"
#include "../tPtrList/tPtrList.h"
#include "../tGridList/tGridList.h"
#include "../tList/tList.h"
#include "../tStorm/tStorm.h"
#include "../tStreamNet/tStreamNet.h"
#include "../tLNode/tLNode.h"
#include "../GridElements/gridElements.h"
#include "../tGrid/tGrid.h"
#include "../tInputFile/tInputFile.h"
#include "../GlobalFns.h"
/*#include "../Inclusions.h"*/

/**Class tStreamMeander****/
class tStreamMeander
{
public:
   //sets everything to zero:
   tStreamMeander();
   //what should be called, after tStreamNet obj. is done
   //(will crash if ptr to tStreamNet obj. is zero);
   //does not call MakeReaches():
   tStreamMeander( tStreamNet &, tGrid< tLNode > &, tInputFile & );
   ~tStreamMeander();
   //"get" ready, "set"...go!
   const tList< tPtrList< tLNode > > &getReachList();
   tList< tPtrList< tLNode > > &getReachListNC();
   void setReachList( const tList< tPtrList< tLNode > > & );
   //finds "meanderability" and calls tLNode::setMeanderStatus
   //for each active node:
   void FindMeander();
   //sets up reach lists; does not add points to interpolate channel:
   void FindReaches();
   //calls FindMeander, FindReaches, and InterpChannel;
   //loops while points added: 
   void MakeReaches();
   //find hydraulic and channel geometries, respectively;
   //FindHydrGeom is contingent upon current storm conditions
   //and storm variability;
   //FindChanGeom is based on the 1.5-yr storm event,
   //or the mean rainrate if no variability:
   void FindHydrGeom();
   void FindChanGeom();
   //interpolates/adds points along channel; returns 1 if points added, else 0
   int InterpChannel();
   //routines that add nodes on the floodplain;
   //Make... called before nodes are actually moved; Add... called after:
    void MakeChanBorder( tList<double> &, tList<double> &, tList<double> & );
      //tList< tArray< double > > & );
    void AddChanBorder( tList<double> &, tList<double> &, tList<double> & );
      //tList< tArray< double > > & );
   //finds the erodibility of each bank, returns an array [right, left]:
   tArray< double > FindBankErody( tLNode * );
   //CheckBanksTooClose, CheckFlowedgCross, and CheckBrokenFlowedg are to check
   //for 'violations' peculiar to moving channels.
   void CheckBanksTooClose();
   void CheckFlowedgCross();
   void CheckBrokenFlowedg();
   //calls meander_; pass it the running time, storm duration, and cum. mvmt.:
   void CalcMigration( double &, double &, double & ); 
   //calls CalcMigration, routines to keep the mesh and net happy,
   //routines to add nodes;
   //only routine you need to make meandering happen for you:
    void Migrate();
    
protected:
      //ptrs and list stuff:
    tGrid< tLNode > *gridPtr;//ptr to tGrid obj. containing meandering reaches
    tStreamNet *netPtr;      //ptr to tStreamNet obj., just to make sure it exists
    tInputFile *infilePtr;   //ptr to tInputFile obj. containing parameters
    tList< tPtrList< tLNode > > reachList; //list of tPtrLists of reach node ptrs
    tListIter< tPtrList< tLNode > > rlIter;//iterator for reachList
      //tList< tArray< double > > bList;  //list of channel border coords
      //tListIter< tArray< double > > blIter; //iter. for bList
      //data items/parameters
    double critflow;  //minimum discharge for meandering
    int optdiamvar;  //flag w/ 1=>multiple grain sizes
    int optrainvar;  //flag w/ 1=>varying storms=>hydraulic geom != chan. geom.
    double meddiam;   //median grain diameter, if optdiamvar = 0
    double kwds, ewds, ewstn;//coeff's & exp's for dwnstrm & at-a-stn hydr. width
    double knds, ends, enstn;//coeff's & exp's for dwnstrm & at-a-stn hydr. roughness
    double dscrtwids; //nominal channel discretization in channel widths
    double allowfrac; //maximum channel point mvmt. in channel widths
    double leavefrac; //distance in widths to add new node
    double vegerod;  //erodibility of vegetated surface or bank
    double rockerod; //erodibility of bedrock
    tArray< int > nrnodes; //array w/ #elements=#reaches, # "active" reach nodes
    tArray< double > reachlen; // " , length of active reach in reach list
    tArray< double > taillen;  // " , length of inactive "tail" in reach list
    long seed;
};

#endif
