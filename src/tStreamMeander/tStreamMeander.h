//-*-c++-*-

/**************************************************************************/
/**
**  @file tStreamMeander.h
**  @brief Header file for class tStreamMeander.
**
**  tStreamMeander objects contain data and functions related to lateral
**  erosion and channel migration across the landscape surface.
**
**    * Contains functions for constructing reaches and calculating lateral
**      channel migration for nodes in those reaches.
**
**    * Contains pointers to a tMesh, a tStreamNet, and a tInputFile.
**
**    * Kicks ass and takes names.
**
**  Modifications:
**   - 6/99: GT commented out latadjust and vegerod in favor of using
**           a single parameter, rockerod, to describe the rate of bank
**           erosion per unit bank shear stress.
**
**  $Id: tStreamMeander.h,v 1.39 2005-03-15 17:17:30 childcvs Exp $
*/
/**************************************************************************/

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
#include "../tList/tList.h"
#include "../tStorm/tStorm.h"
#include "../tStreamNet/tStreamNet.h"
#include "../tLNode/tLNode.h"
#include "../MeshElements/meshElements.h"
#include "../tMesh/tMesh.h"
#include "../tInputFile/tInputFile.h"
#include "../globalFns.h"
/*#include "../Inclusions.h"*/

/** @class tStreamMeander
*/
class tStreamMeander
{
  tStreamMeander& operator=(const tStreamMeander&);
public:
   //sets everything to zero:
   tStreamMeander();
   //what should be called, after tStreamNet obj. is done
   //(will crash if ptr to tStreamNet obj. is zero);
   //does not call MakeReaches():
   tStreamMeander( tStreamNet &, tMesh< tLNode > &, tInputFile &, tRand *rand );
  tStreamMeander(const tStreamMeander&);
   ~tStreamMeander();
   //finds "meanderability" and calls tLNode::setMeanderStatus
   //for each active node:
   void FindMeander();
   //sets up reach lists; does not add points to interpolate channel:
   void FindReaches();
   //calls FindMeander, FindReaches, and InterpChannel;
   //loops while points added:
   //MakeReaches is sent the current time
   void MakeReaches( double );
   //interpolates/adds points along channel; returns 1 if points added, else 0
   //time is passed to it
   const tArray< double > FindInterpCoords( tLNode*, tLNode* );
   tLNode* FindShortcutNode( tLNode*, tLNode* );
   tLNode* BlockShortcut( tLNode*, tLNode*, tLNode&, const tArray<double>&, double );
   int InterpChannel( double );
   //routines that add nodes on the floodplain;
   //Make... called before nodes are actually moved; Add... called after:
   void MakeChanBorder();
   void AddChanBorder( double );
   void FindBankCoords( tLNode*, tArray< double >& ); // carved from Make...
   void ResetEffNbrCoords( tLNode * );
   //finds the erodibility of each bank, returns an array [right, left]:
   tArray< double > FindBankErody( tLNode * ) const;
   //CheckBanksTooClose, CheckFlowedgCross, and CheckBrokenFlowedg are to check
   //for 'violations' peculiar to moving channels.
   int InChannel( tLNode *, tLNode const * ); //called by CheckBanksTooClose()
   tLNode* getUpstreamMeander(tLNode*);
   void CheckBndyTooClose();
   void CheckBanksTooClose();
   void CheckFlowedgCross();
   tLNode* FixBrokenFlowedg( tLNode*, tLNode*, double );
   void CheckBrokenFlowedg( double );
   //calls meander_; pass it the running time, storm duration, and cum. mvmt.:
   void CalcMigration( double &, double const &, double & );
   //calls CalcMigration, routines to keep the mesh and net happy,
   //routines to add nodes;
   //only routine you need to make meandering happen for you:
   //Migrate is sent the current time
   void Migrate( double );
  inline void setMeshPtr( tMesh<tLNode>* Ptr ) {meshPtr = Ptr;}
  inline void setNetPtr( tStreamNet* Ptr ) {netPtr = Ptr;}
  inline void setInfilePtr( tInputFile* Ptr ) {infilePtr = Ptr;}
  inline void setRandPtr( tRand* Ptr ) {rand = Ptr;}

protected:
      //ptrs and list stuff:
    tMesh< tLNode > *meshPtr;//ptr to tMesh obj. containing meandering reaches
    tStreamNet *netPtr;      //ptr to tStreamNet obj., just to make sure it exists
    tInputFile *infilePtr;   //ptr to tInputFile obj. containing parameters

    typedef tListNodeBasic< tPtrList< tLNode > >  rlListNode_t;
    tList< tPtrList< tLNode > > reachList; //list of tPtrLists of reach node ptrs
    tListIter< tPtrList< tLNode > > rlIter;//iterator for reachList
      //data items/parameters
#define FIXCRITFLOWBUG 1
#if FIXCRITFLOWBUG
    double critarea;  //minimum drainage area for meandering
#else
    double critflow;  //minimum discharge for meandering
#endif
#undef FIXCRITFLOWBUG
    bool optdiamvar;  //flag w/ 1=>multiple grain sizes
    bool optrainvar;  //flag w/ 1=>varying storms=>hydraulic geom != chan. geom.
    double meddiam;   //median grain diameter, if optdiamvar = 0
    double kwds, ewds, ewstn;//coeff's & exp's for dwnstrm & at-a-stn hydr. width
    double knds, ends, enstn;//coeff's & exp's for dwnstrm & at-a-stn hydr. roughness
    double kdds, edds, edstn;//coeff's & exp's for dwnstrm & at-a-stn hydr. width
   double klambda, elambda; //coeff & exp for downstrm bank roughness length
    double dscrtwids; //nominal channel discretization in channel widths
    double allowfrac; //maximum channel point mvmt. in channel widths
    double leavefrac; //distance in widths to add new node
    //double vegerod;  //erodibility of vegetated surface or bank
    double rockerod; //erodibility of bedrock
    //double latadjust; //ratio bank erody:bed erody; lets us speed up lat. erosion
   double Pdz; //dependence of bank erody on bank height, 0=none, 1=all
    tArray< int > nrnodes; //array w/ #elements=#reaches, # "active" reach nodes
    tArray< double > reachlen; // " , length of active reach in reach list
    tArray< double > taillen;  // " , length of inactive "tail" in reach list
    tRand *rand;
};

#endif






