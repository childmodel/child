/**************************************************************************\
**
**  tLNode.cpp
**
**  Functions for derived class tLNode and its member classes
**
**  $Id: tLNode.cpp,v 1.4 1998-01-20 20:31:12 stlancas Exp $
\**************************************************************************/

#include <iostream.h>
#include <fstream.h>
#include <assert.h>
#include <math.h>
#include "../Definitions.h"
#include "../Classes.h"
#include "../tListNode/tListNode.h"
#include "../tList/tList.h"
#include "../tPtrListNode/tPtrListNode.h"
#include "../tPtrList/tPtrList.h"
#include "../tArray/tArray.h"
#include "../tNode/tNode.h"
#include "../tEdge/tEdge.h"

#include "tLNode.h"

/*************************************************************************
**  tDeposit::tDeposit : Constructor function for tDeposit
**                       should be fed pointer to layer below it.
**
**  created 7/10/97  NG
*************************************************************************/
tDeposit::tDeposit ()
        : dgrade( NUMG )
{
   dpth=0;
     //cout << "tDeposit()" << endl;
}

tDeposit::tDeposit ( int num )
        : dgrade( num )
{
   dpth=0;
     //cout << "tDeposit( num )" << endl;
}

tDeposit::tDeposit( const tDeposit &orig )                                //tDeposit
        :dgrade( orig.dgrade )
{
   if( &orig != 0 )
   {
      dpth = orig.dpth;
   }
}

tDeposit::~tDeposit()
{
     //cout << "    ~tDeposit()" << endl;
}

const tDeposit &tDeposit::operator=( const tDeposit &right )     //tDeposit
{
   if( &right != this )
   {
      dgrade = right.dgrade;
      dpth = right.dpth;
   }
   return *this;
}

tErode::tErode()                                                    //tErode
{
     //erodtype = 0;
   sedinput = totdz = zp = qs = qsp = qsin = qsinp = tau = dz = 0.0;
   nsmpts = 0;
     //cout << "  tErode()" << endl;
}

tErode::tErode( const tErode &orig )                                //tErode
{
   if( &orig != 0 )
   {
        //erodtype = orig.erodtype;
      sedinput = orig.sedinput;
      totdz = orig.totdz;
      zp = orig.zp;
      qs = orig.qs;
      qsp = orig.qsp;
      qsin = orig.qsin;
      qsinp = orig.qsinp;
      tau = orig.tau;
      nsmpts = orig.nsmpts;
      dz = orig.dz;
      newdz = orig.newdz;
      smooth = orig.smooth;
   }
     //cout << "  tErode( orig )" << endl;
}

tErode::tErode( int numg, int nums )                                //tErode
        : newdz( numg ), smooth( nums )
{
   nsmpts = nums;
   sedinput = totdz = zp = qs = qsp = qsin = qsinp = tau = dz = 0.0;
     //cout << "  tErode( numg, nums )" << endl;
}

tErode::~tErode()                                                   //tErode
{
     //cout << "    ~tErode()" << endl;
}

//assignment
const tErode &tErode::operator=( const tErode &right )     //tErode
{
   if( &right != this )
   {
      sedinput = right.sedinput;
      totdz = right.totdz;
      zp = right.zp;
      qs = right.qs;
      qsp = right.qsp;
      qsin = right.qsin;
      qsinp = right.qsinp;
      tau = right.tau;
      nsmpts = right.nsmpts;
      dz = right.dz;
      newdz = right.newdz;
      smooth = right.smooth;
   }
   return *this;
}

tMeander::tMeander()                                              //tMeander
{
   meander = head = reachmember = 0;
   newx = newy = deltax = deltay = zoldright = zoldleft = 0.0;
     //cout << "  tMeander()" << endl;
}

tMeander::tMeander( const tMeander &orig )                        //tMeander
{
   if( &orig != 0 )
   {
      meander = orig.meander;
      head = orig.head;
      reachmember = orig.reachmember;
      newx = orig.newx;
      newy = orig.newy;
      deltax = orig.deltax;
      deltay = orig.deltay;
      zoldright = orig.zoldright;
      zoldleft = orig.zoldleft;
   }
     //cout << "  tMeander( orig )" << endl;
}

tMeander::tMeander( int state, float x, float y )                //tMeander
{
     //chanptr = 0;
   meander = state;
   newx = x;
   newy = y;
   head = reachmember = 0;
   deltax = deltay = zoldright = zoldleft = 0.0;
     //cout << "  tMeander( state, x, y )" << endl;
}

tMeander::~tMeander()                                             //tMeander
{
     //cout << "    ~tMeander()" << endl;
}

//assignment
const tMeander &tMeander::operator=( const tMeander &right )     //tMeander
{
   if( &right != this )
   {
      meander = right.meander;
      head = right.head;
      reachmember = right.reachmember;
      newx = right.newx;
      newy = right.newy;
      deltax = right.deltax;
      deltay = right.deltay;
      zoldright = right.zoldright;
      zoldleft = right.zoldleft;
   }
   return *this;
}

tBedrock::tBedrock()                                             //tBedrock
{
   erodibility = 0;
     //cout << "  tBedrock()" << endl;
}

tBedrock::tBedrock( const tBedrock &orig )                       //tBedrock
{
   if( &orig != 0 )
   {
      erodibility = orig.erodibility;
   }
     //cout << "  tBedrock( orig )" << endl;
}

tBedrock::~tBedrock()                                            //tBedrock
{
     //cout << "    ~tBedrock()" << endl;
}

//assignment
const tBedrock &tBedrock::operator=( const tBedrock &right )     //tBedrock
{
   if( &right != this )
   {
      erodibility = right.erodibility;
   }
   return *this;
}

tSurface::tSurface()                                             //tSurface
{
     //vegresistfactor = 0;
   veg = tauc = 0.0;
     //cout << "  tSurface()" << endl;
}

tSurface::tSurface( const tSurface &orig )                       //tSurface
{
   if( &orig != 0 )
   {
        //vegresistfactor = orig.vegresistfactor;
      veg = orig.veg;
      tauc = orig.tauc;
   }
     //cout << "  tSurface( orig )" << endl;
}

tSurface::~tSurface()                                            //tSurface
{
     //cout << "    ~tSurface()" << endl;
}

//assignment
const tSurface &tSurface::operator=( const tSurface &right )     //tSurface
{
   if( &right != this )
   {
      veg = right.veg;
      tauc = right.tauc;
   }
   return *this;
}

tRegolith::tRegolith()                                          //tRegolith
        : dgrade( NUMG ), dpth( ACTDEPTH )
{
   thickness = 0;
   numal = 0;
     //dpth = ACTDEPTH;
     //cout << "  tRegolith()" << endl;
}

tRegolith::tRegolith( const tRegolith &orig )                   //tRegolith
        : dgrade( orig.dgrade ), depositList( orig.depositList )
{
   if( &orig != 0 )
   {
      thickness = orig.thickness;
      numal = orig.numal;
      dpth = orig.dpth;
   }
     //cout << "  tRegolith( orig )" << endl;
}

tRegolith::tRegolith( int numg, float vald )                    //tRegolith
        : dgrade( numg ), dpth( vald )
{
   assert( vald > 0.0 );
   thickness = 0.0;
   numal = 0;
     //cout << "  tRegolith( numg, vald )" << endl;
}

tRegolith::~tRegolith()                                         //tRegolith
{
     //cout << "    ~tRegolith()" << endl;
}

//assignment
const tRegolith &tRegolith::operator=( const tRegolith &right )     //tRegolith
{
   if( &right != this )
   {
      dgrade = right.dgrade;
      depositList = right.depositList;
      dpth = right.dpth;
      thickness = right.thickness;
      numal = right.numal;
   }
   return *this;
}

tChannel::tChannel()                                             //tChannel
{
   drarea = q = chanwidth = hydrwidth = channrough = hydrnrough =
       chandepth = hydrdepth = diam = 0;
   diam = kVeryHigh;
     //cout << "  tChannel()" << endl;
}

tChannel::tChannel( const tChannel &orig )                       //tChannel
        : erosion( orig.erosion ), migration( orig.migration )
{
   if( &orig != 0 )
   {
      drarea = orig.drarea;
      q = orig.q;
      chanwidth = orig.chanwidth;
      chandepth = orig.chandepth;
      channrough = orig.channrough;
      hydrwidth = orig.hydrwidth;
      hydrdepth = orig.hydrdepth;
      hydrnrough = orig.hydrnrough;
      diam = orig.diam;
   }
     //cout << "  tChannel( orig )" << endl;
}

tChannel::~tChannel()                                            //tChannel
{
     //cout << "    ~tChannel()" << endl;
}

//assignment
const tChannel &tChannel::operator=( const tChannel &right )     //tChannel
{
   if( &right != this )
   {
      erosion = right.erosion;
      migration = right.migration;
      drarea = right.drarea;
      q = right.q;
      chanwidth = right.chanwidth;
      chandepth = right.chandepth;
      channrough = right.channrough;
      hydrwidth = right.hydrwidth;
      hydrdepth = right.hydrdepth;
      hydrnrough = right.hydrnrough;
      diam = right.diam;
   }
   return *this;
}

tLNode::tLNode()                                                   //tLNode
        : tNode(), rock(), surf(), reg(), chan()
{
     //cout << "=>tLNode()" << endl;
   flood = 0;
   flowedge = 0;
}

tLNode::tLNode( const tLNode &orig )                               //tLNode
        : tNode( orig ),
          rock( orig.rock ), surf( orig.surf ),
          reg( orig.reg ), chan( orig.chan )
{
   flowedge = orig.flowedge;
   flood = orig.flood;
     //cout << "=>tLNode( orig )" << endl;
}

tLNode::~tLNode()                                                  //tLNode
{
     //cout << "    ~tLNode()" << endl;
}

//assignment
const tLNode &tLNode::operator=( const tLNode &right )                  //tNode
{
   if( &right != this )
   {
      tNode::operator=( right );
      rock = right.rock;
      surf = right.surf;
      reg = right.reg;
      chan = right.chan;
      flowedge = right.flowedge;
      flood = right.flood;
   }
   return *this;
}


int tLNode::GetFloodStatus() {   return flood; }

void tLNode::SetFloodStatus( int status )
{
   flood = status;
}

tEdge * tLNode::GetFlowEdg() 
{
   return flowedge;
}

void tLNode::SetFlowEdg( tEdge * edg )
{
   assert( edg > 0 );  // Fails when passed an invalid edge
   //cout << "Setting flow edge to edge # " << edg->getID() << endl;
   
   flowedge = edg;
}

void tLNode::SetDrArea( float val ) {chan.drarea = ( val >= 0 ) ? val : 0;}

void tLNode::AddDrArea( float val ) {chan.drarea += ( val >= 0 ) ? val : 0;}

tLNode * tLNode::GetDownstrmNbr()
{
   return (tLNode *)flowedge->getDestinationPtrNC();     
}

int tLNode::Meanders() const {return chan.migration.meander;}
void tLNode::SetMeanderStatus( int val )
{chan.migration.meander = ( val == 0 || val == 1 ) ? val : 0;}


float tLNode::getHydrWidth() const {return chan.hydrwidth;}
float tLNode::getChanWidth() const {return chan.chanwidth;}
float tLNode::getHydrDepth() const {return chan.hydrdepth;}
float tLNode::getChanDepth() const {return chan.chandepth;}
float tLNode::getHydrRough() const {return chan.hydrnrough;}
float tLNode::getChanRough() const {return chan.channrough;}

void tLNode::setHydrWidth( float val )  {chan.hydrwidth = ( val > 0 ) ? val : 0;}
void tLNode::setChanWidth( float val )  {chan.chanwidth = ( val > 0 ) ? val : 0;}
void tLNode::setHydrDepth( float val )  {chan.hydrdepth = ( val > 0 ) ? val : 0;}
void tLNode::setChanDepth( float val )  {chan.chandepth = ( val > 0 ) ? val : 0;}
void tLNode::setHydrRough( float val )  {chan.hydrnrough = ( val > 0 ) ? val : 0;}
void tLNode::setChanRough( float val )  {chan.channrough = ( val > 0 ) ? val : 0;}
float tLNode::getDrArea() const {return chan.drarea;}

tArray< float >
tLNode::getZOld() const
{
   tArray< float > rl(2);
   rl[0] = chan.migration.zoldright;
   rl[1] = chan.migration.zoldleft;
   return rl;
}

void tLNode::setZOld( float right, float left )
{
   chan.migration.zoldright = right;
   chan.migration.zoldleft = left;
}

tArray< float >                                                   //tNode
tLNode::getNew2DCoords() const
{
   tArray< float > xy(2);
   xy[0] = chan.migration.newx;
   xy[1] = chan.migration.newy;
   return xy;
}

tArray< float >                                                   //tNode
tLNode::getNew3DCoords() const
{
   tArray< float > xyz(3);
   xyz[0] = chan.migration.newx;
   xyz[1] = chan.migration.newy;
   xyz[2] = z;
   return xyz;
}

void tLNode::setNew2DCoords( float val1, float val2 )                    //tNode
{
   chan.migration.newx = val1;
   chan.migration.newy = val2;
}

tArray< float > tLNode::
getLatDisplace() const
{
   tArray< float > xy(2);
   xy[0] = chan.migration.deltax;
   xy[1] = chan.migration.deltay;
   return xy;
}

void tLNode::setLatDisplace( float dx, float dy )
{
   chan.migration.deltax = dx;
   chan.migration.deltay = dy;
}

void tLNode::addLatDisplace( float dx, float dy )
{
   chan.migration.deltax += dx;
   chan.migration.deltay += dy;
}



void tLNode::setDischarge( float val ) {chan.q = ( val > 0 ) ? val : 0;}

void tLNode::RevertToOldCoords()
{
   chan.migration.newx = x;
   chan.migration.newy = y;
}

void tLNode::UpdateCoords()
{
   x = chan.migration.newx;
   y = chan.migration.newy;
}

// nb: if channel is integrated into node, change this
float tLNode::GetQ()
{
   return chan.q;
}

/************************************************************************\
**  GetSlope: Computes and returns the slope of the node's flowedg, or
**  zero if the slope is less than zero.
**
**  Assumptions: edge lengths up to date and nonzero, flowedg's up to
**    date.
\************************************************************************/
float tLNode::GetSlope()
{
   assert( GetFlowEdg()->getLength()>0 ); // failure means lengths not init'd
   float slp = (z - GetDownstrmNbr()->getZ() ) / GetFlowEdg()->getLength();
   if( slp>=0.0 ) return slp;
   else return 0.0;
}

   


//void tLNode::insertFrontSpokeList( tEdge *eptr )             //tNode
//{spokeList.insertAtFront( eptr );}

//void tLNode::insertBackSpokeList( tEdge *eptr )              //tNode
//{spokeList.insertAtBack( eptr );}

/*****************************************************************************\
**
**
**      DistNew: replaces disnew; the distance of the node's (newx, newy)
**              from the line 
**              formed by points  p0->(newx, newy) and p1->(newx, newy)
**      Data members updated: 
**      Called by: 
**      Calls:  
**        
**
\*****************************************************************************/
float tLNode::DistNew(tLNode * p0,tLNode * p1 )
{
  float a,b,c,res;

  a = (p1->chan.migration.newy) - (p0->chan.migration.newy);  
  b = -((p1->chan.migration.newx) - (p0->chan.migration.newx));
  c = -((a*(p0->chan.migration.newx)) + (b*(p0->chan.migration.newy)));
  res = (a*chan.migration.newx + b*chan.migration.newy + c) / sqrt(a*a+b*b);
  if (res<0) res=-res;
  return(res);
}



/**************************************************************************\
**
**  TellAll
**
**  ...is a debugging routine that prints out a bunch of info about a node.
**
\**************************************************************************/
#ifndef NDEBUG
void tLNode::TellAll()
{
   tLNode * nbr;
   
   cout << " NODE " << id << ":\n";
   if( edg ) {
      cout << "  x=" << x << " y=" << y << " z=" << z;
      cout << "  points to edg #" << edg->getID() << endl;
      cout << "  dr area: " << getDrArea() << "  discharge: " << GetQ()
           << "  boundary status: " << boundary << "  flood status: "
           << flood << endl;
      
      if( flowedge ) {
         nbr = (tLNode *)flowedge->getDestinationPtrNC();
         cout << "  Flows along edg " << flowedge->getID() << " to node "
              << nbr->getID() << " at (" << nbr->getX() << ","
              << nbr->getY() << "," << nbr->getZ() << ")\n";
      }
      else cout << "  Flowedg is undefined\n";
   }
   else cout << "  edg is undefined!\n";
   
}
#endif

/**************************************************************************\
**
**  Tracer-sorting routines:
**
**  These routines are utilities that are used in sorting the nodes
**  according to their position within the drainage network. The main
**  sorting algorithm is implemented in tStreamNet::SortNodesByNetOrder().
**  The sorting method works by introducing a "tracer" at each point,
**  then allowing the tracers to iteratively cascade downstream. At each
**  step any nodes not containing tracers are moved to the back of the
**  list. The result is a list sorted in upstream-to-downstream order.
**
**  These utilities do the following:
**    ActivateSortTracer -- injects a single tracer at a node
**    AddTracer -- adds a tracer to a node (ignored if node is a bdy)
**    MoveSortTracerDownstream -- removes a tracer and sends it to the
**                                downstream neighbor (unless the node is
**                                a sink; then the tracer just vanishes)
**    NoMoreTracers -- reports whether there are any tracers left here
**
**  Created by GT 12/97.
**
\**************************************************************************/

void tLNode::ActivateSortTracer()
{ tracer = 1; }

void tLNode::MoveSortTracerDownstream()
{
   tracer--;
   if( !flood ) GetDownstrmNbr()->AddTracer();
}

void tLNode::AddTracer()
{
   if( !boundary ) tracer++;
}

int tLNode::NoMoreTracers()
{
   assert( tracer>=0 );
   return( tracer==0 );
}


/**************************************************************************\
**
**  EroDep
**
**  Erodes or deposits material of depth dz.
**
**  Note: as of now, just changes elev. Should also change alluvial
**  thickness and set to zero if negative (bedrock erosion).
**
**  1/15/98 gt
\**************************************************************************/
void tLNode::EroDep( float dz )
{
   z += dz;
   cout << "  eroding " << id << " by " << dz << endl;
   
   //sed += dz;
   //if( sed<0 ) sed=0;
}
