/**************************************************************************\
**
**  tStreamMeander.cpp
**
**  Functions for class tStreamMeander.
**
**  $Id: tStreamMeander.cpp,v 1.5 1998-01-21 23:23:39 stlancas Exp $
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
   infilePtr = 0;
   optdiamvar = optrainvar = 0;
   critflow = meddiam = kwds = ewds = ewstn = knds = ends = enstn =
       dscrtwids = leavefrac = 0;
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
   infilePtr = &infile;
   assert( infilePtr != 0 );
   critflow = infilePtr->ReadItem( critflow, "CRITICAL_FLOW" );
   optdiamvar = infilePtr->ReadItem( optdiamvar, "OPT_VAR_SIZE" );
   if( !optdiamvar )
   {
      meddiam = infilePtr->ReadItem( meddiam, "MEDIAN_DIAMETER" );
      assert( meddiam > 0 );
   }
   optrainvar = infilePtr->ReadItem( optrainvar, "OPTVAR" );
   kwds = infilePtr->ReadItem( kwds, "HYDR_WID_COEFF_DS" );
   assert( kwds > 0 );
   ewds = infilePtr->ReadItem( ewds, "HYDR_WID_EXP_DS" );
   ewstn = infilePtr->ReadItem( ewstn, "HYDR_WID_EXP_STN" );
   knds = infilePtr->ReadItem( knds, "HYDR_ROUGH_COEFF_DS" );
   assert( knds > 0 );
   ends = infilePtr->ReadItem( ends, "HYDR_ROUGH_EXP_DS" );
   enstn = infilePtr->ReadItem( enstn, "HYDR_ROUGH_EXP_STN" );
   dscrtwids = infilePtr->ReadItem( dscrtwids, "DEF_CHAN_DISCR" );
   assert( dscrtwids > 0 );
   allowfrac = infilePtr->ReadItem( allowfrac, "FRAC_WID_MOVE" );
   assert( allowfrac > 0 );
   leavefrac = infilePtr->ReadItem( leavefrac, "FRAC_WID_ADD" );
   assert( leavefrac > 0 );
   MakeReaches();
   assert( &reachList != 0 );
}

tStreamMeander::~tStreamMeander()
{
   //if( netPtr != 0 ) delete netPtr;
   gridPtr = 0;
   netPtr = 0;
   infilePtr = 0;
}


/*****************************************************************************\
**
**       FindMeander: simply goes through active grid nodes and assigns
**                    meander status flag based on a threshold discharge
**                
**
**		  Parameters: critflow
**      Data members updated: tLNode->chan.migration.meander
**      Called by: MakeReaches
**      Calls: no "major" functions
**
**      Created: SL 1/98         
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
**       FindHydrGeom: goes through reach nodes and calculates/assigns
**                     values for hydraulic geometry.
**                
**
**		  Parameters: kwds, ewds, ewstn, knds, ends, enstn, optrainvar
**      Data members updated: tLNode->chan.hydrwidth, hydrnrough, hydrdepth
**      Called by: FindReaches (needs to know how long to make reach "tails"
**      Calls: no "major" functions
**
**      Created: SL 1/98         
**
**
\*****************************************************************************/
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
   for( plPtr = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
        plPtr = rlIter.NextP(), i++ )
   {
      nIter.Reset( *plPtr );
      num = nrnodes[i];
      for( cn = nIter.FirstP(), j=0; j<num; cn = nIter.NextP(), j++ )
      {
         //if rainfall varies, find hydraulic width "at-a-station"
         //based on the channel width "downstream":
         if( optrainvar)
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
         //if rainfall does not vary, set hydraulic geom. = channel geom.
         else
         {
            width = cn->getChanWidth();
            rough = cn->getChanRough();
            depth = cn->getChanDepth();
            cn->setHydrWidth( width );
            cn->setHydrRough( rough );
            cn->setHydrDepth( depth );
         }
      }
   }
}

   

/*****************************************************************************\
**
**       FindChanGeom: goes through reach nodes and calculates/assigns
**                     values for channel geometry.
**                
**
**		  Parameters: kwds, ewds, ewstn, knds, ends, enstn
**      Data members updated: tLNode->chan.hydrwidth, hydrnrough, hydrdepth
**      Called by: FindReaches (needs to know how long to make reach "tails"
**      Calls: no "major" functions
**
**      Created: SL 1/98         
**
**
\*****************************************************************************/
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
   for( plPtr = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
        plPtr = rlIter.NextP(), i++ )
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


/****************************************************************\
**
**	InterpChannel: Fill in points between widely spaced
**			channel nodes with noisy interpolation
**			
**	
**		Parameters: dscrtwids -- default spacing in # hydr widths
**		Called by:  MakeReaches
**    Calls: tGrid::AddNode, tStreamNet::UpdateNet
**		Created:  5/16/97 SL
**    Updated: 1/98 SL
**
\****************************************************************/

int tStreamMeander::InterpChannel()
{
   int npts, num;
   float curwidth;
   float curseglen, defseglen, maxseglen, bigseglen;
   float x, y, z, val, phi, x0, y0, z0, x1, y1, slope;
   tPtrListIter< tLNode > rnIter;
   tArray< float > xp, yp, zp, *arrPtr;
   tLNode *crn, nn, *nPtr;
   int change = 0; //haven't added any nodes yet
   //loop through reaches:
   for( creach = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
        creach = rlIter.NextP(), i++ )
   {
      rnIter.Reset( *creach );
      num = nrnodes[i];
      //loop through reach nodes:
      for( crn = rnIter.FirstP(), j=0; j<num;
           crn = rnIter.NextP(), j++ )
      {
         curwidth = crn->getHydrWidth();
         nPtr = crn->GetDownstreamNbr();
         slope = crn->GetSlope();
         curseglen = fedg->getLength();
         defseglen = dscrtwids * curwidth;
         maxseglen = 2.0 * defseglen;
         //if current flowedg length is greater than
         //a given proportion of hydraulic width
         if( curseglen > maxseglen )
         {
            change = 1;//flag so we know that we need to update the network
            nn = *crn; //added nodes are copies of the upstream node except xyz.
            x0 = crn->getX();
            y0 = crn->getY();
            z0 = crn->getZ();
            x1 = nPtr->getX();
            y1 = nPtr->getY();
            y = y1 - y0;
            x = x1 - x0;
            phi = atan2( y, x );
            //cout << "pref'd seg length is " << prefseglen << endl << flush;
            bigseglen = 3.0 * defseglen;
            //if we can get the desired spacing by adding a single node
            //use linear interpolation with noise added
            //(for Delaunay mesh, avoid precisely co-linear points):
            if( curseglen <= bigseglen )
            {
               val = ran3(seed) - 0.5;
               x = x0 + curseglen * cos(phi) /
                   2.0 + 0.01 * val * pow(defseglen, 2.0);
               val = ran3(seed) - 0.5;
               y = y0 + curseglen * sin(phi) /
                   2.0 + 0.01 * val * pow(defseglen, 3.0); 
               z = z0 - curseglen / 2.0 * slope;
               nn.set3DCoords( x, y, z );
               gridPtr->AddNode( nn );
            }
            //otherwise, if we need to add more than one point,
            //generate a random walk with uniform spacing in x
            //and exponentially weighted steps in y
            //(weighting to keep last step, back to the orig. downstream node,
            //from being too much of a doozy):
            else
            {
               npts = ROUND( curseglen / prefseglen );
               arrPtr = new tArray< float >( npts );
               xp = yp = zp = *arrPtr;
               for(int i=1; i<npts; i++ )
               {
                  xp[i] = (float)i * defseglen;
                  val = ran3(seed) - 0.5;
                  yp[i] = yp[i-1] + 0.01 * val * exp(-yp[i-1] * val) 
                      * pow(defseglen, 3.0);
                  zp[i] = xp[i] * slope; 
                  x = sqrt(xp[i] * xp[i] + yp[i] * yp[i]) *
                      cos(atan2(yp[i], xp[i]) + phi);
                  y = sqrt(xp[i] * xp[i] + yp[i] * yp[i]) *
                      sin(atan2(yp[i], xp[i]) + phi);
                  x += x0;
                  y += y0;
                  z = z0 - zp[i];
                  cout<<"InterpChannel: call AddNode"<<endl<<flush;
                  nn.set3DCoords( x, y, z );
                  gridPtr->AddNode( nn );
               }
               delete arrPtr;
            }
         }
      }
   }
   if( change )
   {
      netPtr->UpdateNet();
      return 1;
   }
   else return 0;
}

/*****************************************************************************\
**
**      MakeReaches : takes care of the loop of
**                    1) find meandering nodes
**                    2) make reaches out of meandering nodes
**                    3) go through reaches and add points
**                       between nodes which are too far apart
**                    4) Update the areas, stream net, etc.
**                    5) go back to #1 if any points were added
**
**      Data members updated: reachList, nrnodes, reachlen, taillen
**      Called by: tStreamMeander(...) constructor
**      Calls: FindMeander, FindReaches,
**             InterpChannel=>calls netPtr->UpdateNet() if points are added
**
**      Created:  1/98 SL
**      Modified:          
**
**
\*****************************************************************************/

void tStreamMeander::MakeReaches()
{
   do
   {
      FindMeander();
      FindReaches();
   }
   while( InterpChannel() );
}


/*****************************************************************************\
**
**      FindReaches : constructs the reach objects; does
**                    not interpolate channel:
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
**       **Might be good to put in something that will excise     **
**       **very sharp bends, i.e., bank undulations, from reaches.**
**                
**
**      Data members updated: reachList, nrnodes, reachlen, taillen
**      Called by: MakeReaches
**      Calls: FindChanGeom, FindHydrGeom
**
**      Created:  1/98 SL
**      Modified:          
**
**
\*****************************************************************************/
void tStreamMeander::FindReaches()
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
   //loop through active nodes
   for( cn = nodIter.FirstP(); !(nodIter.isActive()); cn = nodIter.NextP() )
   {
      //if node meanders
      if( cn->Meanders() )
      {
         spokIter.Reset( cn->getSpokeListNC() );
         nmndrnbrs = 0;
         //loop through spokes to find upstream meandering nbr
         for( ce = spokIter.FirstP(); !(spokIter.AtEnd()); ce = spokIter.NextP() )
         {
            if( ce->getDestinationPtrNC()->GetDownstrmNbr() == cn 
                && ce->getDestinationPtrNC()->Meanders() )
            {
               nmndrnbrs++;
               break;
            }
         }
         //if no upstream meandering nbr, then it's a reach "head";
         //put it at the beginning of an otherwise empty reach
         //and add the reach to the list of reaches:
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
   //loop through reaches
   for( plPtr = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
        plPtr = rlIter.NextP(), i++ )
   {
      //go downstream from reach head
      //and add nodes to the reach
      //if they're not already reach members:
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
      //now we'll need to know the
      //channel and hydraulic geometry:
      FindChanGeom();
      FindHydrGeom();
      //construct a tail for the reach
      //which is some number of hydr. widths long:
      curwidth = rnIter.LastP()->getHydrWidth();
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
**	CalcMigration: 	This is the routine that makes the 
**				arrays to pass to the fortran
**				routine 'meander_' and sets "new" x and y values;
**        makes a bunch of tArrays, and passes the arrays to
**        fortran by way of the tArray member ptr, gotten
**        with getArrayPtr().
**
**		Parameters:	allowfrac -- fraction of chanwidth 
**				          a node is allowed to move in a 
**				          given meander iteration
**				        duration -- storm duration,
**                  copy from netPtr->stormPtr->stdur
**                time -- running tab on how long we've meandered
**
**		Called by:	Migrate
**    Calls: FindBankErody, external fortran routine _meander
**		Created: 5/1/97  SL
**
\***************************************************************/
void tStreamMeander::CalcMigration( float &time, float &duration )
{
   int i, j, *stations, *stnserod, nttlnodes;
   tArray< float > xa, ya, xsa, qa, rerodya, lerodya, delsa,
       slopea, widtha, deptha, diama, deltaxa, deltaya,
       rdeptha, ldeptha;
   tArray< float > *dumArrPtr, delta(2), newxy(2);
   float rerody, lerody, rz, lz, width;
   float maxfrac, displcmt, a, b, dtm, tmptim, frac, xs;
   long seed;
   float num;
   tPtrList< tLNode > *creach;
   tPtrListIter< tLNode > rnIter;
   tLNode *curnode;
   tEdge *fedg;
   tArray< float > bankerody;
   //loop through reaches:
   for( creach = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
        creach = rlIter.NextP(), i++ )
   {
      rnIter.Reset( *creach );
      num = nrnodes[i];
      for( curnode = rnIter.FirstP(), j=0; j<num;
           curnode = rnIter.NextP(), j++ )
          //initialize deltax, deltay, newx, newy:
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
          //arbitrary change for simple debugging:
          //( 0.1, 0.01*(ran3(&seed)-.5)
      //go through and set the elevations at the right and left banks:
      num = nrnodes[i];
      for( curnode = rnIter.FirstP(), j=0; j<num;
           curnode = rnIter.NextP(), j++ )
      {
         rz = curnode->getZ() + deptha[j] - rdeptha[j];
         lz = curnode->getZ() + deptha[j] - ldeptha[j];
         curnode->setZOld( rz, lz );
      }
      delete dumArrPtr;
   }
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
      for( curnode = rnIter.FirstP(), j=0; j<num;
           curnode = rnIter.NextP(), j++ )
      {
         delta = curnode->getLatDisplace();
         displcmt = sqrt( delta[0] * delta[0] + delta[1] * delta[1] );
         width = curnode->getChanWidth();
         if( width > 0.0 ) frac = displcmt / width;
         else frac = 0.0;
         if( frac > maxfrac ) maxfrac = frac;
      }
   }
   //if( maxfrac < allowfrac )
   //    dtm = ( maxfrac > 0 ) ? allowfrac / maxfrac : 1.0;
   //else dtm = allowfrac / maxfrac;
   
   //maximize the time step:
   //the smaller of the time to move the allowed distance...
   dtm = ( maxfrac > 0 ) ? allowfrac / maxfrac : 1.0;
   for( creach = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
        creach = rlIter.NextP(), i++ )
   {
      rnIter.Reset( *creach );
      num = nrnodes[i];
      for( curnode = rnIter.FirstP(), j=0; j<num;
           curnode = rnIter.NextP(), j++ )
      {
         tmptim = time + dtm;
         //...and the time remaining in the storm
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
   cummvmt += maxfrac * dtm;
}

/****************************************************************\
**	
**  Migrate: 
**		Parameters:
**		Called by: Main
**    Calls: CalcMigration, MakeChanBorder, CheckBanksTooClose,
**           CheckFlowedgCross, CheckBrokenFlowedg,
**           tGrid::MoveNodes, AddChanBorder
**		Created: 1/98 SL
**
\***************************************************************/
void tStreamMeander::Migrate()
{
   tList< tArray< float > > bList;
   float duration = netPtr->getStormPtrNC()->GetStormDuration();
   float time = 0.0;
   float cummvmt = 0.0;
   //timeadjust = 86400. * pr->days;
   while( time < duration)
   {
      CalcMigration( time, duration, cummvt ); //incremented time
      MakeChanBorder( bList ); //bList of coordinate arrays made
      CheckBanksTooClose();
      CheckFlowedgCross();
      CheckBrokenFlowedg();
      gridPtr->MoveNodes();
      AddChanBorder( bList );
      //after the channel migrates a certain amount
      //(here, maximum migration distances, in units of hydr. width,
      //at each iteration are summed and compared to 1.0)
      //update the stream net, redo the reaches,
      //interpolate if necessary, etc.
      if( cummvmt > 1.0 )
      {
         netPtr->UpdateNet();
         MakeReaches();
         cummvmt = 0.0;
      }
   }
}


/******************************************************************************\
**
**      MakeChanBorder(): Make stack of point coords for left and right banks
**                        of meandering channels; called after _meander but
**                        before points are actually moved on the grid.
**
**              Parameters:     
**              Called by:      Migrate
**              Created:        8/18/97 SL
**              Updated: 1/98 SL
**
\*****************************************************************************/

void tStreamMeander::MakeChanBorder( tList< tArray< float > > &bList )
{
   int i, j, num;
   float x0, y0, x1, y1, x, y, z, delx, dely, phi, width, xdisp, ydisp;
   tPtrList< tLNode > *cr;
   tPtrListIter< tLNode > rnIter;
   tLNode *cn, *cnbr;
   tArray< float > xyz(3), rl;
   for( cr = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
        cr = rlIter.NextP(), i++ )
   {
      rnIter.Reset( *cr );
      num = nrnodes[i];
      for( cn = rnIter.FirstP(), j=0; j<num;
           cn = rnIter.NextP(), j++ )
      {
         width = cn->getHydrWidth();
         rl = cn->getZOld();
         cnbr = cn->GetDownstreamNbr();
         x0 = cn->getX();
         y0 = cn->getY();
         x1 = cnbr->getX();
         y1 = cnbr->getY()
         delx = x1 - x0;
         dely = y1 - y0;
         phi = atan2( dely, delx );
         xdisp = 0.5 * width * sin(phi);
         ydisp = 0.5 * width * cos(phi);
         xyz[0] = x0 + xdisp;
         xyz[1] = y0 - ydisp;
         xyz[2] = rl[0];
         bList.insertAtBack( xyz );
         xyz[0] = x0 - xdisp;
         xyz[1] = y0 + ydisp;
         xyz[2] = rl[1];
         bList.insertAtBack( xyz );
      }
   }
}

/******************************************************************************\
**
**	AddChanBorder: After meandering points have been moved and the
**                       triangulation adjusted, check points in list to 
**                       find whether they should be discarded or added.
**
**              Parameters:     uses node->hydrwidth and, indirectly, 
**                              parameters which determine hydraulic geometry
**              Called by:      Migrate
**              Created:        8/18/97 SL
**              Updated:        1/98 SL
**
\*****************************************************************************/
void tStreamMeander::AddChanBorder( tList< tArray< float > > &bList )
{
   int i;
   float x, y, halfwid, dist, mindist = 10000000.;
   tArray< float > cp, xy;
   tTriangle *ct;
   tLNode *cn, *channodePtr, channode;
   //go through list of coordinates made by MakeChanBorder:
   while( bList.removeFromFront( cp ) ) //copies the tArray in the list to cp
   {
      //find triangle in which coords lie:
      ct = gridPtr->LocateTriangle( cp[0], cp[1] );
      if( ct != 0 )
      {
         halfwid = 0;
         //find shortest distance between coords and tri vertices:
         for( i=0; i<3; i++ )
         {
            cn = (tLNode *) ct->pPtr(i);
            xy = cn->get2DCoords();
            dist = sqrt( (cp[0]-xy[0])*(cp[0]-xy[0]) +
                         (cp[1]-xy[1])*(cp[1]-xy[1]) );
            if( dist < mindist ) mindist = dist;
            if( cn->Meanders() )
            {
               halfwid = cn->getHydrWidth() * 0.5;
               channodePtr = cn;
            }
         }
         //if one of the vertices is a meandering node,
         //and the smallest distance to a vertex is less than
         //half a channel width (any way to do a different distance?),
         //copy the meand. node, reset
         //its coords., and add it to the grid:
         if( halfwid > 0 && mindist > halfwid )
         {
            channode = *channodePtr;
            channode.set3DCoords( cp[0], cp[1], cp[2] );
            gridPtr->AddNode( channode );
         }
      }
   }
   assert( bList.isEmpty() ); //prob. don't need this check, but...
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


/****************************************************************\
**	
**
**		Parameters:	
**		Called by:	
**		Created: 1/98 SL
**
\***************************************************************/
