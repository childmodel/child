/**************************************************************************\
**
**  tStreamMeander.h: Header file for class tStreamMeander.
**
**  tStreamMeander objects contain objects and functions related to lateral
**  erosion and channel migration across the landscape surface.
**  Contains functions for constructing reaches and calculating lateral
**  channel migration for nodes in those reaches.
**
**
**  $Id: tStreamMeander.h,v 1.4 1998-01-20 20:33:04 stlancas Exp $
\**************************************************************************/
#ifndef TSTREAMMEANDER_H
#define TSTREAMMEANDER_H

/**Class tStreamMeander****/
class tStreamMeander
{
public:
   tStreamMeander();
   tStreamMeander( tStreamNet &, tGrid< tLNode > &, tInputFile & );
   //tStreamMeander( tStreamNet &, tInputFile & );
   ~tStreamMeander();
   const tList< tReach > &getReachList();
   tList< tReach > &getReachListNC();
   void setReachList( const tList< tReach > & );
   void FindMeander();
   void MakeReaches();
   void FindHydrGeom();
   void FindChanGeom();
   void InterpChannel();
   void MakeChanBorder();
   tArray< float > FindBankErody( tLNode * );
     //CheckBanksTooClose, CheckFlowedgCross, and CheckBrokenFlowedg are to check
     //for 'violations' peculiar to moving channels.
   int CheckBanksTooClose();
   int CheckFlowedgCross();
   int CheckBrokenFlowedg();
   void CalcMigration( tStorm & ); //calls meander_; pass it a ref to a storm
   void Migrate();
   void AddChanBorder();

protected:
   //ptrs and ptr lists
   tGrid< tLNode > *gridPtr;
   tStreamNet *netPtr;
   //tList< tReach > reachList;
   tList< tPtrList< tLNode > > reachList;
   //tListIter< tReach > rlIter;
   tListIter< tPtrList< tLNode > > rlIter;
   //data items/parameters
   float critflow;
   int optdiamvar;
   float meddiam;
   float kwds, ewds, ewstn;//coeff's & exp's for dwnstrm & at-a-stn hydr. width
   float knds, ends, enstn;//coeff's & exp's for dwnstrm & at-a-stn hydr. roughness
   tArray< int > nrnodes;
   tArray< float > reachlen;
   tArray< float > taillen;
};

#endif
