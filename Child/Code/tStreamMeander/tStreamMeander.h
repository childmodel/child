/**************************************************************************\
**
**  tStreamMeander.h: Header file for class tStreamMeander.
**
**  tStreamMeander objects contain objects and functions related to lateral
**  erosion and channel migration across the landscape surface.
**  Contains functions for constructing reaches and calculating lateral
**  channel migration for nodes in those reaches.
**
**  $Id: tStreamMeander.h,v 1.1 1998-01-16 19:16:23 stlancas Exp $
\**************************************************************************/
#ifndef TSTREAMMEANDER_H
#define TSTREAMMEANDER_H

/**Class tStreamMeander****/
class tStreamMeander
{
public:
   tStreamMeander();
   tStreamMeander( tGrid< tLNode > &, tInputFile & );
   tStreamMeander( tStreamNet &, tInputFile & );
   ~tStreamMeander();
   const tList< tReach > &getReachList();
   tList< tReach > &getReachListNC();
   void setReachList( const tList< tReach > & );
   void FindMeander();
   void MakeReaches();
   void InterpChannel();
   void FindHydrGeom();
   void MakeChanBorder();
   void FindBankErody();
     //CheckBanksTooClose, CheckFlowedgCross, and CheckBrokenFlowedg are to check
     //for 'violations' peculiar to moving channels.
   int CheckBanksTooClose();
   int CheckFlowedgCross();
   int CheckBrokenFlowedg();
   void CalcMigration();
   void Migrate();
   void LeavePoints();

protected:
   tGrid< tLNode > *gridPtr;
   tStreamNet *netPtr;
   tList< tReach > reachList;
   tListIter< tReach > rlIter;
};

#endif
