/**************************************************************************\
**
**  tStreamMeander.cpp
**
**  Functions for class tStreamMeander.
**
**  $Id: tStreamMeander.cpp,v 1.4 1998-01-20 20:30:25 stlancas Exp $
\**************************************************************************/

#include "Inclusions.h"

extern "C" 
{
   void meander_( int *, int *, float *, float *, float *, float *, 
                  float *, float *, float *, float *, float *, 
                  float *, float *, float *, float *, float *, float * ); 
}


/**************************************************************************\
**
**  Constructors
**
**  1) takes grid and inputfile references and constructs a tStreamNet
**     object if ptr to tStreamNet is zero
**
**  2) takes grid ptr and inputfile reference; assumes tStreamNet object
**     has already done its thing
**
**  2) takes streamnet and inputfile references
**
\**************************************************************************/
tStreamMeander::tStreamMeander()
        : reachList(), rlIter( reachList )
{
   gridPtr = 0;
   netPtr = 0;
   critflow = meddiam = = dwds = ewds = ewstn = dnds = ends = enstn = 0;
}

tStreamMeander::tStreamMeander( tStreamNet &netRef, tGrid< tLNode > &gRef,
                                tInputFile &infile )
        : reachList(), rlIter( reachList )
{
   //if( netPtr != 0 ) netPtr = new tStreamNet( gRef );
   netPtr = &netRef;
   assert( netPtr != 0 );
   gridPtr = &gRef;
   assert( gridPtr != 0 );
   critflow = infile.ReadItem( critflow, "CRITICAL_FLOW" );
   optdiamvar = infile.ReadItem( optdiamvar, "OPT_VAR_SIZE" );
   if( !optdiamvar )
       meddiam = infile.ReadItem( meddiam, "MEDIAN_DIAMETER" );
   kwds = infile.ReadItem( kwds, "HYDR_WID_COEFF_DS" );
   ewds = infile.ReadItem( ewds, "HYDR_WID_EXP_DS" );
   ewstn = infile.ReadItem( ewstn, "HYDR_WID_EXP_STN" );
   knds = infile.ReadItem( knds, "HYDR_ROUGH_COEFF_DS" );
   ends = infile.ReadItem( ends, "HYDR_ROUGH_EXP_DS" );
   enstn = infile.ReadItem( enstn, "HYDR_ROUGH_EXP_STN" );
   FindMeander();
   MakeReaches();
   assert( &reachList != 0 );
}

tStreamMeander::~tStreamMeander()
{
   //if( netPtr != 0 ) delete netPtr;
   gridPtr = 0;
   netPtr = 0;
}


/*****************************************************************************\
**
**       
**                
**
**      Data members updated: 
**      Called by: 
**      Calls: 
**
**      Created:  YC
**      Added:   YC
**      Modified:          
**
**
\*****************************************************************************/

void tStreamMeander::FindMeander()
{
   tLNode * cn;
   tGridListIter< tLNode > nodIter( gridPtr->GetNodeList() );
   for( cn = nodIter.FirstP(); cn->isActive(); cn = nodIter.NextP() )
   {
      if( cn->GetQ() >= critflow )
          cn->SetMeanderStatus( kMeanderNode );
      else
          cn->SetMeanderStatus( kNonMeanderNode );
   }
}

void tStreamMeander::FindHydrGeom()
{
   int i, j, num;
   float hradius, kwdspow, kndspow, widpow, npow, radfactor, qpsec;
   float width, depth, rough, slope;
   tLNode *cn;

   kwdspow = pow(kwds, ewstn / ewds);
   kndspow = pow(knds, enstn / ends);
   widpow = 1.0 - ewstn / ewds;
   npow = 1.0 - enstn / ends;
   //timeadjust = 86400 * days;  /* 86400 = seconds in a day */
   tPtrListIter< tLNode > nIter;
   tPtrList< tLNode > *plPtr;
   for( plPtr = rlIter.FirstP(), i=0; !(rlIter.AtEnd()); plPtr = rlIter.NextP(), i++ )
   {
      nIter.Reset( *plPtr );
      num = nrnodes[i];
      for( cn = nIter.FirstP(), j=0; j<num; cn = nIter.NextP(), j++ )
      {
         qpsec = cn->GetQ();
         width = pow(cn->getChanWidth(), widpow) * kwdspow * pow(qpsec, ewstn);
         cn->setHydrWidth( width );
         rough = pow(cn->getChanRough(), npow) * kndspow * pow(qpsec, enstn);
         cn->setHydrRough( rough );
         slope = cn->GetSlope();
         assert( slope > 0 );
         radfactor = qpsec * rough / width / sqrt(slope);
         hradius = pow(radfactor, 0.6);
         depth = width / ( width / hradius - 2.0 );
         cn->setHydrDepth( depth );
      }
   }
}

   

void tStreamMeander::FindChanGeom()
{
   int i, j, num;
   float qbf, hradius, qbffactor=0, radfactor, width, depth, rough, slope;
   tLNode *cn;
   tPtrListIter< tLNode > nIter;
   tPtrList< tLNode > *plPtr;
   //timeadjust = 86400 * days;  /* 86400 = seconds in a day */
   if (isdmn)  qbffactor = pmn * log(1.5 / isdmn);
// qbffactor is now in m^3/s
   for( plPtr = rlIter.FirstP(), i=0; !(rlIter.AtEnd()); plPtr = rlIter.NextP(), i++ )
   {
      nIter.Reset( *plPtr );
      num = nrnodes[i];
      for( cn = nIter.FirstP(), j=0; j<num; cn = nIter.NextP(), j++ )
      {
         qbf = cn->getDrArea() * qbffactor;
         if (! qbf) qbf=curnode->GetQ();  // q is now in m^3/s
         width = kwds * pow(qbf, ewds);
         rough = knds * pow(qbf, ends);
         slope = cn->GetSlope();
         assert( slope > 0 );
         cn->setChanWidth( width );
         cn->setChanRough( rough );
         radfactor = qbf * rough / width / sqrt(slope);
         hradius = pow(radfactor, 0.6); 
         depth = width / (width / hradius - 2.0);
         cn->setChanDepth( depth );
      }
   }
}



/*****************************************************************************\
**
**      MakeReaches : constructs the reach objects:
**         1) finds channel head nodes, i.e., furthest upstream extent of
**            channel that meanders, and puts each head at the beginning of
**            a ptr list of meandering nodes; these lists are the reaches;
**            adds each one-member list to a list of reaches;
**         2) goes downstream from heads and adds a node to a reach list if it
**            is not already a reach member, i.e., has not been added to
**            another reach list; when a node is added to the list, its
**            'reachmember' flag is set to keep it from being added to another
**            reach;
**         3) because the meandering "forces" are distal, i.e., they act on
**            points downstream, we tack a "tail" onto each reach; tails have
**            arbitrary length, some number * the channel width at the last
**            reach node;
**         4) during the above, keep track of the number of nodes in each
**            reach, the reach length, and the tail length; these are kept
**            in arrays which are dimensioned after we know how many reaches
**            there are.
**                
**
**      Data members updated: 
**      Called by: 
**      Calls: 
**
**      Created:  1/98 SL
**      Modified:          
**
**
\*****************************************************************************/

void tStreamMeander::MakeReaches()
{
   float curwidth, ctaillen;
   int i, j, nmndrnbrs;
   tLNode *cn;
   tEdge *ce;
   tGridListIter< tLNode > nodIter( gridPtr->GetNodeList() );
   tPtrListIter< tEdge > spokIter;
   tPtrListIter< tLNode > rnIter;
   tPtrList< tLNode > rnodList, *plPtr;
   tArray< int > *iArrPtr;
   tArray< float > *fArrPtr;
   if( !(reachList.isEmpty()) ) reachList.Flush();
   for( cn = nodIter.FirstP(); !(nodIter.isActive()); cn = nodIter.NextP() )
   {
      if( cn->Meanders() )
      {
         spokIter.Reset( cn->getSpokeListNC() );
         nmndrnbrs = 0;
         for( ce = spokIter.FirstP(); !(spokIter.AtEnd()); ce = spokIter.NextP() )
         {
            if( ce->getDestinationPtrNC()->GetDownstrmNbr() == cn 
                && ce->getDestinationPtrNC()->Meanders() )
            {
               nmndrnbrs++;
               break;
            }
         }
         if (nmndrnbrs==0)
         {
            rnodList.Flush();
            rnodList.insertAtFront( cn );
            reachList.insertAtBack( rnodList );
         }
      }
   }
   iArrPtr = new tArray< int >( reachList.getSize() );
   nrnodes = *iArrPtr;
   delete iArrPtr;
   fArrPtr = new tArray< float >( reachList.getSize() );
   reachlen = *fArrPtr;
   taillen = *fArrPtr;
   delete fArrPtr;
   for( plPtr = rlIter.FirstP(), i=0; !(rlIter.AtEnd()); plPtr = rlIter.NextP(), i++ )
   {
      rnIter.Reset( *plPtr );
      cn = rnIter.FirstP();
      do
      {
         nrnodes[i]++;
         reachlen[i] += cn->GetFlowEdg()->getLength();
         cn->setReachMember( 1 );
         plPtr->insertAtBack( cn );
         cn = cn->GetDownstrmNbr();
      }
      while (!cn->getReachMember() && cn->getBoundFlag() == kNonBoundary);
      curwidth = rnIter.LastP()->getChanWidth();
      taillen[i] = 10.0*curwidth;
      ctaillen = 0.0;
      while (ctaillen <= taillen[i] && cn->getBoundFlag() == kNonBoundary)
      {
         ctaillen+= cn->flowedg->len;
         plPtr->insertAtBack( cn );
         cn = cn->GetDownstrmNbr();
      }
      taillen[i] = ctaillen;
   }
}


/****************************************************************\
**	CalcMigrate: 	This is the routine that makes the 
**				arrays to pass to the fortran
**				routine 'meander_', and finally 
**				takes care of some of the 
** 				necessary post-processing;
**        makes a bunch of tArrays, and passes the array to
**        fortran by way of the tArray member ptr, gotten
**        with getArrayPtr().
**
**		Parameters:	allowfrac -- fraction of chanwidth 
**				a node is allowed to move in a 
**				given meander iteration
**				duration -- storm duration, copy 
**				of 'storm.stdur'
**				dtm -- meandering time step, 
**				subdivision of 'duration'
**		Called by:	Main
**		Created: 5/1/97  SL
**
\***************************************************************/
void tStreamMeander::CalcMigrate( tStorm &storm )
{
   int i, j, *stations, *stnserod, nttlnodes;
   tArray< float > xa, ya, xsa, qa, rerodya, lerodya, delsa,
       slopea, widtha, deptha, diama, deltaxa, deltaya,
       rdeptha, ldeptha;
   tArray< float > *dumArrPtr, delta(2), newxy(2);
   float rerody, lerody, allowfrac, rz, lz, width;
   float maxfrac, displcmt, a, b, time, dtm, tmptim, frac, xs;
   long seed;
   float num;
   float duration = storm.GetStormDuration();
   //tReach * creach;
   tPtrList< tLNode > *creach;
   tPtrListIter< tLNode > rnIter;
   tLNode *curnode;
   tEdge *fedg;
   tArray< float > bankerody;
   time = 0.0;
   allowfrac = 0.1;
   //timeadjust = 86400. * pr->days;
   while( time < duration)
   {
      //loop through reaches:
      for( creach = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
           creach = rlIter.NextP(), i++ )
      {
         rnIter.Reset( *creach );
         //initialize deltax, deltay, newx, newy:
         for( curnode = rnIter.FirstP(); !(rnIter.AtEnd()); curnode = rnIter.NextP() )
         {
            curnode->setLatDisplace( 0.0, 0.0 );
            curnode->setNew2DCoords( curnode->getX(), curnode->getY() );
         }
         stations = &(nrnodes[i]);
         nttlnodes = creach->getSize();
         stnserod = &nttlnodes;
         dumArrPtr = new tArray< float >( nttlnodes );
         xa = ya = xsa = qa = rerodya = lerodya = delsa =
             slopea = widtha = deptha = diama = deltaxa =
             deltaya = rdeptha = ldeptha = *dumArrPtr;
         xs = 0.0;
         for( curnode = rnIter.FirstP(), j=0; !(rnIter.AtEnd());
              curnode = rnIter.NextP(), j++ )
         {
            fedg = curnode->GetFlowEdg();
            xa[j] = curnode->getX();
            ya[j] = curnode->getY();
            xsa[j] = xs;
            xs += fedg->getLength();
            qa[j] = curnode->GetQ();
            delsa[j] = fedg->getLength();
            xs += delsa[j];
            bankerody = FindBankErody( curnode );
            //curnode->GetErodibility( rerody, lerody, pr->kf);
            rerodya[j] = bankerody[0];
            lerodya[j] = bankerody[1];
            slopea[j] = fedg->getSlope();
            widtha[j] = curnode->getHydrWidth();
            deptha[j] = curnode->getHydrDepth();
            /*diama[j] = curnode->diam;*/
            diama[j] = ( optdiamvar ) ? curnode->getDiam() : meddiam;
         }
         cout << "stations, stnserod: " << *stations <<" "<< *stnserod
              << endl << flush;
         //this looks horrible, but we need to pass the pointer to the array
         //itself, not the tArray object, which contains the array pointer.
         meander_( stations,
                   stnserod,
                   xa->getArrayPtr(),
                   ya->getArrayPtr(),
                   xsa->getArrayPtr(), 
                   delsa->getArrayPtr(),
                   qa->getArrayPtr(),
                   rerodya->getArrayPtr(),
                   lerodya->getArrayPtr(), 
                   slopea->getArrayPtr(),
                   widtha->getArrayPtr(),
                   deptha->getArrayPtr(),
                   diama->getArrayPtr(),
                   deltaxa->getArrayPtr(),
                   deltaya->getArrayPtr(),
                   rdeptha->getArrayPtr(),
                   ldeptha->getArrayPtr() );
         //Now reset the node values according to the arrays:
         for( curnode = rnIter.FirstP(), j=0; !(rnIter.AtEnd());
              curnode = rnIter.NextP(), j++ )
             curnode->addLatDisplace( deltaxa[j], deltaya[j] );
//( 0.1, 0.01*(ran3(&seed)-.5)
         num = nrnodes[i];
         for( curnode = rnIter.FirstP(), j=0; j<num; curnode = rnIter.NextP(), j++ )
         {
            rz = curnode->getZ() + deptha[j] - rdeptha[j];
            lz = curnode->getZ() + deptha[j] - ldeptha[j];
            curnode->setZOld( rz, lz );
         }
         delete dumArrPtr;
      }
      //if( bugtime > 89 ) DumpEdges();
      //calculate ratio of total displacement to length of flow edge,
      //find maximum of ratios and make sure it is less than or equal,
      //to the allowed fraction (e.g., 1/10) and scale displacements
      //as necessary.
      maxfrac = 0.0;
      for( creach = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
           creach = rlIter.NextP(), i++ )
      {
         rnIter.Reset( *creach );
         num = nrnodes[i];
         for( curnode = rnIter.FirstP(), j=0; j<num; curnode = rnIter.NextP(), j++ )
         {
            delta = curnode->getLatDisplace();
            displcmt = sqrt( delta[0] * delta[0] + delta[1] * delta[1] );
            width = curnode->getChanWidth();
            if( width > 0.0 ) frac = displcmt / width;
            else frac = 0.0;
            if( frac > maxfrac ) maxfrac = frac;
         }
      }
      if( maxfrac < allowfrac ) 
      {
         dtm = 1.0;
         //cumvmt += maxfrac;
      }
      else
      {
         dtm = allowfrac / maxfrac;
         //cumvmt += allowfrac;
      }
      for( creach = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
           creach = rlIter.NextP(), i++ )
      {
         rnIter.Reset( *creach );
         num = nrnodes[i];
         for( curnode = rnIter.FirstP(), j=0; j<num; curnode = rnIter.NextP(), j++ )
         {
            //dtm = 1.0;
            tmptim = time + dtm;
            if( tmptim > duration ) dtm = duration - time;
            delta = curnode->getLatDisplace();
            delta[0] *= dtm;
            delta[1] *= dtm;
            curnode->setLatDisplace( delta[0], delta[1] );
            newxy = curnode->getNew2DCoords();
            newxy[0] += delta[0];
            newxy[1] += delta[1];
            curnode->setNew2DCoords( newxy[0], newxy[1] );
         }
      }
      time += dtm;
      //if( bugtime > 89 ) DumpEdges();
   }
}

/******************************************************************************\
**
**      MakeChanBorder(): Make stack of point coords for left and right banks
**                        of meandering channels; called after _meander but
**                        before points are actually moved on the grid.
**
**              Parameters:     
**              Called by:      StreamMeander
**              Created:        8/18/97 SL
**
\*****************************************************************************/
void tGrid::MakeChanBorder()
{
   int j;
   float x, y, z, delx, dely, phi;
   tReach * cr;
   tNode * cn;
   for( cr = firstreach; cr; cr = cr->next )
   {
      for( j = 0; j < cr->nrnodes; j++ )
      {
         cn = cr->reachnode[j];
         delx = cn->flowedg->dest->x - cn->x;
         dely = cn->flowedg->dest->y - cn->y;
         phi = atan2( dely, delx );
         x = cn->x + 0.5 * cn->hydrwidth * sin(phi);
         y = cn->y - 0.5 * cn->hydrwidth * cos(phi);
         z = cn->zoldright;
         AddPointStack( x, y, z );
         x = cn->x - 0.5 * cn->hydrwidth * sin(phi);
         y = cn->y + 0.5 * cn->hydrwidth * cos(phi);
         z = cn->zoldleft;
         AddPointStack( x, y, z );
      }
   }
}

/******************************************************************************\
**
**	GetErodibility:	This is the routine that finds the effective erodibility 
**                              of each bank
**				based on a reach node's neighbor's erodibility 
**                    and relative height 
**				above the channel.
**
**		Parameters:	erody -- copy of pr->kf, the fluvial erodibility
**		Called by:	StreamMeander
**		Created: 5/1/97 SL
**
\*****************************************************************************/
   tArray< float > FindBankErody( tLNode * );
void tNode::GetErodibility( float rerody, float lerody, float erody )
{
   float x1, y1, x2, y2, dx, dy, dx1, dy1, a, b, c, d, dmin, dlast,
       dzleft, dzright;
   tEdge * curedg, * ledg, * redg;
//simplest thing: find point on left which has smallest distance 
//to the perpendicular; once that distance starts getting larger 
//again, find other (right bank) point:
   x1=x;
   y1=y;
   dx = flowedg->dest->x - x;
   dy = flowedg->dest->y - y;
   dx1 = -dy;
   dy1 = dx;
   a = dy1;
   b = -dx1;
   c = dx1 * y1 - dy1 * x1;
   dmin = 1000.0;
   dlast = dmin;
   d = dmin - 1.;
   curedg = flowedg;
   do 
   {
      curedg = curedg->nextedg;
      dlast = d;
      x2 = curedg->dest->x;
      y2 = curedg->dest->y;
      d = fabs( a * x2 + b * y2 + c );
      if( d < dmin )
      {
         dmin = d;
         ledg = curedg;
      }
   } while( d < dlast );
   dmin = 1000.0;
   do
   {
      x2 = curedg->dest->x;
      y2 = curedg->dest->y;
      d = fabs( a * x2 + b * y2 + c );
      if( d < dmin )
      {
         dmin = d;
         redg = curedg;
      }
      curedg = curedg->nextedg;
   } while( curedg != flowedg );
   dzleft = ledg->dest->z - z;
   dzright = redg->dest->z - z;
   if( dzleft > hydrdepth ) lerody = erody * hydrdepth / dzleft;
   else lerody = erody;
   if( dzright > hydrdepth ) rerody = erody * hydrdepth / dzright;
   else rerody = erody;
}

/******************************************************************************\
**
**	CheckPointStack: After meandering points have been moved and the
**                       triangulation adjusted, check points in stack to 
**                       find whether they should be discarded or added.
**
**              Parameters:     uses node->hydrwidth and, indirectly, 
**                              parameters which determine hydraulic geometry
**              Called by:      StreamMeander
**              Created:        8/18/97 SL
**
\*****************************************************************************/
void tGrid::CheckPointStack()
{
   int i;
   float halfwid, dist, mindist = 10000000.;
   tPointStack * cps, * cpshold;
   tTriangle * ct;
   tNode * cn, channode;
   for( cps = firststackpoint; cps; cps = cps->next )
   {
      ct = LocateTriangle( cps->x, cps->y );
      if( ct )
      {
         halfwid = NULL;
         for( i=0; i<3; i++ )
         {
            cn = ct->p[i];
            dist = VectorLength( cps->x, cps->y, cn->x, cn->y );
            if( dist < mindist ) mindist = dist;
            if( cn->meander ) halfwid = cn->hydrwidth * 0.5;
         }
         if( halfwid ) if( mindist > halfwid )
             CreateAdd( cps->x, cps->y, cps->z, 0 );
      }
   }
   for( cps = firststackpoint; cps; )
   {
      cpshold = cps->next;
      delete cps;
      cps = cpshold;
   }
}


