/**************************************************************************\
**
**  tStreamMeander.cpp
**
**  Functions for class tStreamMeander.
**
**  $Id: tStreamMeander.cpp,v 1.2 1998-01-16 22:05:59 stlancas Exp $
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
   critflow = 0;
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
/*****************************************************************************\
**
**      GetOnlyReaches : constructs the reach objects       
**                
**
**      Data members updated: 
**      Called by: 
**      Calls: 
**
**      Created:  
**      Added:   YC
**      Modified:          
**
**
\*****************************************************************************/

void tGrid::GetOnlyReaches()
{
  float curwidth, ctaillen;
  int i, j, nmndrnbrs;
  tNode *curn = firstnode;
  tNode *crn, *frnode;
  tEdge *ce;
  nheads = 0;
  if( firstreach )
  {
     delete [] firstreach;
     firstreach = lastreach = NULL;
  }

  for (curn = firstnode , i=0 ; i<nActiveNodes ; i++)
  {
    if (curn->meander)
    {
      nmndrnbrs = 0;
      curn->CountNbrs();
      for (j=0, ce=curn->edg; j<curn->nnbrs; j++, ce=ce->nextedg)
      {
        if (ce->dest->flowedg->dest == curn 
	      && ce->dest->meander) {nmndrnbrs++; break;}
      }
      if (nmndrnbrs==0)
      {
        nheads++;
        curn->head = TRUE;
      }
    }
    curn = curn->next;
  }
  firstreach = new tReach[nheads];

  for (curn=firstnode , j=0 , i=0; i<nActiveNodes; i++)
  {
    if (curn->head)
    {
      firstreach[j].reachlen=0.0;
      firstreach[j].nrnodes = 0;
      crn = curn;
      frnode = curn;
      do
      {
        firstreach[j].nrnodes++;
        curwidth = crn->chanwidth;
        firstreach[j].reachlen+= crn->flowedg->len;
        crn->reachmember = TRUE;
        crn = crn->flowedg->dest;
      }
      while (!crn->reachmember && crn->boundary != kClosedBoundary);
      firstreach[j].taillen = 10.0*curwidth;
      ctaillen = 0.0;
      firstreach[j].ntnodes=0;
      do
      {
        ctaillen+= crn->flowedg->len;
        crn = crn->flowedg->dest;
        firstreach[j].ntnodes++;
      }
      while (ctaillen <= firstreach[j].taillen &&
             crn->boundary != kClosedBoundary);
        //firstreach[j].ClearNodes();
      firstreach[j].GetReachNodes(frnode);
      if (j<nheads-1) firstreach[j].next = &firstreach[j+1];
      else
      {
         firstreach[j].next = NULL;
         lastreach = &firstreach[j];
      }
      j++;
    }
    curn=curn->next;
  }
}


/****************************************************************\
**	StreamMeander: 	This is the routine that converts the 
**				reaches to arrays, passes the 
**				array pointers to the fortran
**				routine 'meander_', and finally 
**				takes care of some of the 
** 				necessary post-processing.
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
void tGrid::StreamMeander( tParameters * pr, float duration )
{
   int i, j, *stations, *stnserod, nttlnodes;
   tfloatArray * xa, * ya, * xsa, * qa, * rerodya, * lerodya, * delsa,
       * slopea, * widtha, * deptha, * diama, * deltaxa, * deltaya,
       * rdeptha, * ldeptha;
     /*float * xa, * ya, * xsa, * qa, * rerodya, * lerodya, * delsa;
   float * slopea, * widtha, * deptha, * diama, * deltaxa, * deltaya;
   float * rdeptha, * ldeptha;*/
   float rerody, lerody, allowfrac, timeadjust;
   float maxfrac, displcmt, a, b, time, dtm, frac, xs;
   long seed;
   float num;
   tReach * creach;
   tNode * curnode;
   time = 0.0;
   allowfrac = 0.1;
   timeadjust = 86400. * pr->days;
   do
   {
      creach = firstreach;
      curnode = firstnode;
      for (i=0; i<nnodes; i++)
      {
         curnode->deltax = 0.0;
         curnode->deltay = 0.0;
         curnode->newx = curnode->x;
         curnode->newy = curnode->y;
      }
      for (i=0; i<nheads; i++)
      {
         stations = &creach->nrnodes;
         nttlnodes = creach->nrnodes + creach->ntnodes;
         stnserod = &nttlnodes;
         xa = new tfloatArray( nttlnodes );
         ya = new tfloatArray( nttlnodes );
         xsa = new tfloatArray( nttlnodes );
         qa = new tfloatArray( nttlnodes );
         rerodya = new tfloatArray( nttlnodes );
         lerodya = new tfloatArray( nttlnodes );
         delsa = new tfloatArray( nttlnodes );
         slopea = new tfloatArray( nttlnodes );
         if( bugtime > 89 ) DumpEdges();
         widtha = new tfloatArray( nttlnodes );
         deptha = new tfloatArray( nttlnodes );
         diama = new tfloatArray( nttlnodes );
         deltaxa = new tfloatArray( nttlnodes );
         deltaya = new tfloatArray( nttlnodes );
         rdeptha = new tfloatArray( nttlnodes );
         ldeptha = new tfloatArray( nttlnodes );
         xs = 0.0;
         for (j=0; j<nttlnodes; j++)
         {
            xa->fvalue[j] = creach->reachnode[j]->x;
            ya->fvalue[j] = creach->reachnode[j]->y;
            xs += creach->reachnode[j]->flowedg->len;
            xsa->fvalue[j] = xs;
            qa->fvalue[j] = creach->reachnode[j]->q / timeadjust;
            delsa->fvalue[j] = creach->reachnode[j]->flowedg->len;
            creach->reachnode[j]->GetErodibility( rerody, lerody, pr->kf);
            rerodya->fvalue[j] = rerody;
            lerodya->fvalue[j] = lerody;
            slopea->fvalue[j] = creach->reachnode[j]->flowedg->slope;
            widtha->fvalue[j] = creach->reachnode[j]->hydrwidth;
            deptha->fvalue[j] = creach->reachnode[j]->hydrdepth;
              /*diama[j] = creach->reachnode[j]->diam;*/
            diama->fvalue[j] = pr->meddiam;
         }
         cout << "stations, stnserod: " << *stations <<" "<< *stnserod
              << endl << flush;
           /*meander_(stations, stnserod, xa, ya, xsa, 
             delsa, qa, rerodya, lerodya, 
             slopea, widtha, deptha, diama,
             deltaxa, deltaya, rdeptha, ldeptha);*/
           //Now reset the node values according to the arrays:
         for (j=0; j<nttlnodes; j++)
         {
            creach->reachnode[j]->deltax += 0.1;              //deltaxa[j];
            num = ran3(&seed) - 0.5;
            creach->reachnode[j]->deltay += 0.01 * num;       //deltaya[j];
         }
         for (j=0; j<creach->nrnodes; j++)
         {
            creach->reachnode[j]->zoldright = creach->reachnode[j]->z 
                + deptha->fvalue[j] - rdeptha->fvalue[j];
            creach->reachnode[j]->zoldleft = creach->reachnode[j]->z
                + deptha->fvalue[j] - ldeptha->fvalue[j];
         }
         delete xa;
         delete ya;
         delete xsa;
         delete delsa;
         delete qa;
         delete rerodya;
         delete lerodya;
         delete slopea;
         delete widtha;
         delete deptha;
         delete diama;
         delete deltaxa;
         delete deltaya;
         delete rdeptha;
         delete ldeptha;
         creach = creach->next;
      }
      if( bugtime > 89 ) DumpEdges();
     //calculate ratio of total displacement to length of flow edge,
        //find maximum of ratios and make sure it is less than or equal,
        //to the allowed fraction (e.g., 1/10) and scale displacements
        //as necessary.
      maxfrac = 0.0;
      creach = firstreach;
      for (i=0; i<nheads; i++)
      {
         for (j=0; j<creach->nrnodes; j++)
         {
            a = creach->reachnode[j]->deltax;
            b = creach->reachnode[j]->deltay;
            displcmt = sqrt(a * a + b * b);
            if( creach->reachnode[j]->chanwidth != 0.0)
                frac = displcmt / creach->reachnode[j]->chanwidth;
            else frac = 0.0;
            if( frac > maxfrac ) maxfrac = frac;
         }
         creach = creach->next;
      }
      if( maxfrac < allowfrac ) 
      {
         dtm = 1.0;
         cumvmt += maxfrac;
      }
      else
      {
         dtm = allowfrac / maxfrac;
         cumvmt += allowfrac;
      }
      creach = firstreach;
      for (i=0; i<nheads; i++)
      {
         for (j=0; j<creach->nrnodes; j++)
         {
            dtm = 1.0;
            if( dtm > duration ) dtm = duration;

            creach->reachnode[j]->deltax *= dtm;
            creach->reachnode[j]->deltay *= dtm;
            creach->reachnode[j]->newx += creach->reachnode[j]->deltax;
            creach->reachnode[j]->newy += creach->reachnode[j]->deltay;
         }
         creach = creach->next;
      }
      time += dtm;
        //MakeChanBorder(); //puts point coords for right and left banks in stack
      //ChanCloseDel();
      if( bugtime > 89 ) DumpEdges();
      MoveNodes();	// (frmr PreApply) changes 'x' to 'newx', etc. and calls
        //ApplyChanges which adjusts grid triangulation
        //CheckPointStack(); // adds points which should be added, deletes stack
      UpdateMesh();
   }
   while( time < duration);
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


