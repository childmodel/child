/**************************************************************************\
**
**  tStreamMeander.h: Header file for class tStreamMeander.
**
**  tStreamMeander objects contain objects and functions related to lateral
**  erosion and channel migration across the landscape surface.
**  Contains functions for constructing reaches and calculating lateral
**  channel migration for nodes in those reaches.
**
how about this: tStreamMeander contains a pointer to a tStreamNet object;
it has a constructor that takes a ref. to a tGrid object as an argument;
this constructor automatically creates a tStreamNet object if the pointer
is NULL; grid data will then be accessed through the grid ptr in tStreamNet;
calling this constructor guarantees that the necessary data items calc'ed
by tStreamNet are indeed done. tStreamNet can still stand alone, but you can't
meander without a tStreamNet. Could do something similar for erosion?
Since the tStreamNet doesn't contain its own data to speak of, it wouldn't
matter too much if more than one net were constructed. On the other hand, also 
have a meander constructor which takes a reference to a net object;
could also have a gridPtr member of meander, but that would need to assume
that the net stuff had been done...
**
**  $Id: tStreamMeander.h,v 1.3 1998-01-17 23:29:23 stlancas Exp $
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
   void AddChanBorder();

protected:
   tGrid< tLNode > *gridPtr;
   tStreamNet *netPtr;
   //tList< tReach > reachList;
   tList< tPtrList< tLNode > > reachList;
   //tListIter< tReach > rlIter;
   tListIter< tPtrList< tLNode > > rlIter;
   float critflow;
   tArray< int > nrnodes;
   tArray< float > reachlen;
   tArray< float > taillen;
};

#endif
