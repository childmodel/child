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
**  $Id: tStreamMeander.h,v 1.5 1998-01-21 23:24:53 stlancas Exp $
\**************************************************************************/
#ifndef TSTREAMMEANDER_H
#define TSTREAMMEANDER_H

/**Class tStreamMeander****/
class tStreamMeander
{
public:
   //sets everything to zero:
   tStreamMeander();
   //what should be called, after tStreamNet obj. is done
   //(will crash if ptr to tStreamNet obj. is zero):
   tStreamMeander( tStreamNet &, tGrid< tLNode > &, tInputFile & );
   //tStreamMeander( tStreamNet &, tInputFile & );
   ~tStreamMeander();
   //"get" ready, "set"...go!
   const tList< tReach > &getReachList();
   tList< tReach > &getReachListNC();
   void setReachList( const tList< tReach > & );
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
   void MakeChanBorder( tList< tArray< float > > & );
   void AddChanBorder( tList< tArray< float > > & );
   //finds the erodibility of each bank, returns an array [right, left]:
   tArray< float > FindBankErody( tLNode * );
     //CheckBanksTooClose, CheckFlowedgCross, and CheckBrokenFlowedg are to check
     //for 'violations' peculiar to moving channels.
   int CheckBanksTooClose();
   int CheckFlowedgCross();
   int CheckBrokenFlowedg();
   //calls meander_; pass it the running time, storm duration, and cum. mvmt.:
   void CalcMigration( float &, float &, float & ); 
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
   //data items/parameters
   float critflow;  //minimum discharge for meandering
   int optdiamvar;  //flag w/ 1=>multiple grain sizes
   int optrainvar;  //flag w/ 1=>varying storms=>hydraulic geom != chan. geom.
   float meddiam;   //median grain diameter, if optdiamvar = 0
   float kwds, ewds, ewstn;//coeff's & exp's for dwnstrm & at-a-stn hydr. width
   float knds, ends, enstn;//coeff's & exp's for dwnstrm & at-a-stn hydr. roughness
   float dscrtwids; //nominal channel discretization in channel widths
   float allowfrac; //maximum channel point mvmt. in channel widths
   float leavefrac; //distance in widths to add new node
   tArray< int > nrnodes; //array w/ #elements=#reaches, # "active" reach nodes
   tArray< float > reachlen; // " , length of active reach in reach list
   tArray< float > taillen;  // " , length of inactive "tail" in reach list
   long seed;
};

#endif
