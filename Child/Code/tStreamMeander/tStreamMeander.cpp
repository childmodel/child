/**************************************************************************\
**
**  tStreamMeander.cpp
**
**  Functions for class tStreamMeander.
**
**  $Id: tStreamMeander.cpp,v 1.8 1998-01-29 19:49:29 stlancas Exp $
\**************************************************************************/

#include <math.h>
#include <assert.h>
#include "tStreamMeander.h"

extern "C" 
{
   void meander_( int *, int *, double *, double *, double *, double *, 
                  double *, double *, double *, double *, double *, 
                  double *, double *, double *, double *, double *, double * ); 
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
       dscrtwids = leavefrac = vegerod = rockerod = 0;
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
   vegerod = infile.ReadItem( vegerod, "VEG_ERODY" );
   rockerod = infile.ReadItem( rockerod, "KB" );
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
   for( cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP() )
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
   double hradius, kwdspow, kndspow, widpow, npow, radfactor, qpsec;
   double width, depth, rough, slope;
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
   double qbf, hradius, qbffactor=0, radfactor, width, depth, rough, slope;
   tLNode *cn;
   tPtrListIter< tLNode > nIter;
   tPtrList< tLNode > *plPtr;
   //timeadjust = 86400 * days;  /* 86400 = seconds in a day */
   tStorm *sPtr = netPtr->getStormPtrNC();
   double isdmn = sPtr->getMeanInterstormDur();
   double pmn = sPtr->getMeanPrecip();
   if (isdmn > 0 )  qbffactor = pmn * log(1.5 / isdmn);
// qbffactor is now in m^3/s
   for( plPtr = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
        plPtr = rlIter.NextP(), i++ )
   {
      nIter.Reset( *plPtr );
      num = nrnodes[i];
      for( cn = nIter.FirstP(), j=0; j<num; cn = nIter.NextP(), j++ )
      {
         qbf = cn->getDrArea() * qbffactor;
         if( !qbf ) qbf = cn->GetQ();  // q is now in m^3/s
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
   int i, j, npts, num;
   double curwidth;
   double curseglen, defseglen, maxseglen, bigseglen;
   double x, y, z, val, phi, x0, y0, z0, x1, y1, slope;
   tPtrListIter< tLNode > rnIter;
   tPtrList< tLNode > *creach;
   tArray< double > xp, yp, zp, *arrPtr;
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
         nPtr = crn->GetDownstrmNbr();
         slope = crn->GetSlope();
         curseglen = crn->GetFlowEdg()->getLength();
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
               val = ran3(&seed) - 0.5;
               x = x0 + curseglen * cos(phi) /
                   2.0 + 0.01 * val * pow(defseglen, 2.0);
               val = ran3(&seed) - 0.5;
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
               npts = ROUND( curseglen / defseglen );
               arrPtr = new tArray< double >( npts );
               xp = yp = zp = *arrPtr;
               for(int i=1; i<npts; i++ )
               {
                  xp[i] = (double)i * defseglen;
                  val = ran3(&seed) - 0.5;
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
   double curwidth, ctaillen;
   int i, j, nmndrnbrs;
   tLNode *cn, *lnPtr;
   tEdge *ce;
   tGridListIter< tLNode > nodIter( gridPtr->GetNodeList() );
   tPtrListIter< tEdge > spokIter;
   tPtrListIter< tLNode > rnIter;
   tPtrList< tLNode > rnodList, *plPtr;
   tArray< int > *iArrPtr;
   tArray< double > *fArrPtr;
   if( !(reachList.isEmpty()) ) reachList.Flush();
   //loop through active nodes
   for( cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP() )
   {
      //if node meanders
      if( cn->Meanders() )
      {
         spokIter.Reset( cn->getSpokeListNC() );
         nmndrnbrs = 0;
         //loop through spokes to find upstream meandering nbr
         for( ce = spokIter.FirstP(); !(spokIter.AtEnd()); ce = spokIter.NextP() )
         {
            lnPtr = (tLNode *) ce->getDestinationPtrNC();
            if( lnPtr->GetDownstrmNbr() == cn && lnPtr->Meanders() )
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
   fArrPtr = new tArray< double >( reachList.getSize() );
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
      while (!cn->getReachMember() && cn->getBoundaryFlag() == kNonBoundary);
      //now we'll need to know the
      //channel and hydraulic geometry:
      FindChanGeom();
      FindHydrGeom();
      //construct a tail for the reach
      //which is some number of hydr. widths long:
      curwidth = rnIter.LastP()->getHydrWidth();
      taillen[i] = 10.0*curwidth;
      ctaillen = 0.0;
      while (ctaillen <= taillen[i] && cn->getBoundaryFlag() == kNonBoundary)
      {
         ctaillen+= cn->GetFlowEdg()->getLength();
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
void tStreamMeander::CalcMigration( double &time, double &duration,
                                    double &cummvmt )
{
   int i, j, *stations, *stnserod, nttlnodes;
   tArray< double > xa, ya, xsa, qa, rerodya, lerodya, delsa,
       slopea, widtha, deptha, diama, deltaxa, deltaya,
       rdeptha, ldeptha;
   tArray< double > *dumArrPtr, delta(2), newxy(2);
   double rerody, lerody, rz, lz, width;
   double maxfrac, displcmt, a, b, dtm, tmptim, frac, xs;
   double num;
   tPtrList< tLNode > *creach;
   tPtrListIter< tLNode > rnIter;
   tLNode *curnode;
   tEdge *fedg;
   tArray< double > bankerody;
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
      dumArrPtr = new tArray< double >( nttlnodes );
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
                xa.getArrayPtr(),
                ya.getArrayPtr(),
                xsa.getArrayPtr(), 
                delsa.getArrayPtr(),
                qa.getArrayPtr(),
                rerodya.getArrayPtr(),
                lerodya.getArrayPtr(), 
                slopea.getArrayPtr(),
                widtha.getArrayPtr(),
                deptha.getArrayPtr(),
                diama.getArrayPtr(),
                deltaxa.getArrayPtr(),
                deltaya.getArrayPtr(),
                rdeptha.getArrayPtr(),
                ldeptha.getArrayPtr() );
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
   tList< tArray< double > > bList;
   double duration = netPtr->getStormPtrNC()->GetStormDuration();
   double time = 0.0;
   double cummvmt = 0.0;
   //timeadjust = 86400. * pr->days;
   while( time < duration)
   {
      CalcMigration( time, duration, cummvmt ); //incremented time
      MakeChanBorder( bList ); //bList of coordinate arrays made
      CheckBanksTooClose();
      CheckFlowedgCross();
      //CheckBrokenFlowedg();
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

void tStreamMeander::MakeChanBorder( tList< tArray< double > > &bList )
{
   int i, j, num;
   double x0, y0, x1, y1, x, y, z, delx, dely, phi, width, xdisp, ydisp;
   tPtrList< tLNode > *cr;
   tPtrListIter< tLNode > rnIter;
   tLNode *cn, *cnbr;
   tArray< double > xyz(3), rl;
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
         cnbr = cn->GetDownstrmNbr();
         x0 = cn->getX();
         y0 = cn->getY();
         x1 = cnbr->getX();
         y1 = cnbr->getY();
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
void tStreamMeander::AddChanBorder( tList< tArray< double > > &bList )
{
   int i;
   double x, y, halfwid, dist, mindist = 10000000.;
   tArray< double > cp, xy;
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
**	FindBankErody :	This is the routine that finds the effective erodibility 
**                              of each bank
**				based on a reach node's neighbor's erodibility 
**                    and relative height 
**				above the channel.
**
**		Parameters:	tSurface::vegerody; tBedrock::erodibility
**		Called by: CalcMigration
**		Created: 5/1/97 SL
**
\*****************************************************************************/
tArray< double > tStreamMeander::FindBankErody( tLNode *nPtr )
{
   double x1, y1, x2, y2, dx, dy, dx1, dy1, a, b, c, d, dmin, dlast,
       dzleft, dzright, depth, dfactor, dtotal, erody, sed;
   tArray< double > xy, xyz1, xy2, dxy, rlerody(2);
   tEdge * curedg, * ledg, * redg, *ce;
   tLNode *cn, *dn, *rn, *ln;
   tPtrListIter< tEdge > spokIter( nPtr->getSpokeListNC() );
   tPtrList< tLNode > rList, lList;
   tPtrListIter< tLNode > rIter( rList ), lIter( lList );
//simplest thing: find points on right and left which have smallest distance 
//to the perpendicular; those determine right and left erodibility, resp.
//(another way: weighted average according to distance from perpendicular.)
   xyz1 = nPtr->get3DCoords();
   dn = nPtr->GetDownstrmNbr();
   xy2 = dn->get2DCoords();
   for( ce = spokIter.FirstP(); !(spokIter.AtEnd()); ce = spokIter.NextP() )
   {
      cn = (tLNode *) ce->getDestinationPtrNC();
      if( cn != dn && cn->GetDownstrmNbr() != nPtr )
      {
         xy = cn->get2DCoords();
         a = (xyz1[1] - xy[1]) * (xy2[0] - xy[0]);
         b = (xyz1[0] - xy[0]) * (xy2[1] - xy[1]);
         c = a - b;
         if( c > 0.0 ) rList.insertAtBack( cn );
         else lList.insertAtBack( cn );
      }
   }
   dxy[0] = xy2[0] - xyz1[0];
   dxy[1] = xy2[1] - xyz1[1];
   a = dxy[0];
   b = dxy[1];
   c = -dxy[1] * xyz1[1] - dxy[0] * xyz1[0];
   dmin = 100000.0;
   for( cn = rIter.FirstP(); !(rIter.AtEnd()); cn = rIter.NextP() )
   {
      xy = cn->get2DCoords();
      d = fabs( a * xy[0] + b * xy[1] + c );
      if( d < dmin )
      {
         dmin = d;
         rn = cn;
      }
   }
   dmin = 100000.0;
   for( cn = lIter.FirstP(); !(lIter.AtEnd()); cn = lIter.NextP() )
   {
      xy = cn->get2DCoords();
      d = fabs( a * xy[0] + b * xy[1] + c );
      if( d < dmin )
      {
         dmin = d;
         ln = cn;
      }
   }
   dzright = rn->getZ() - xyz1[2];
   dzleft = ln->getZ() - xyz1[2];
     //Xbank erodibility is the nominal erodibility X depth / dz
   //make erodibility depend only on whether alluvium or bedrock:
   depth = nPtr->getHydrDepth();
   if( dzright > 0.0 )
   {
      sed = rn->getAlluvThickness();
      if( dzright > sed ) rlerody[0] = rockerod;
      else rlerody[0] = vegerod;
   }
   if( dzleft > 0.0 )
   {
      sed = ln->getAlluvThickness();
      if( dzleft > sed ) rlerody[1] = rockerod;
      else rlerody[1] = vegerod;
   }
   
   /*
   if( dzright <= 0.0 ) dfactor = 2.0;
   else dfactor = depth / dzright;
   rlerody[0] = dfactor * erody;
   if( dzleft <= 0.0 ) dfactor = 2.0;
   else dfactor = depth / dzleft;
   rlerody[1] = dfactor * erody;
   */
   return rlerody;
}

/*****************************************************************************\
**
**        CheckBanksTooClose()
**
**		Parameters:	
**		Called by:	
**		Created: 1/98 SL
**
**
\*****************************************************************************/
void tStreamMeander::CheckBanksTooClose()
{
   cout << "CheckBanksTooClose()..." << flush << endl;
   int tooclose;
   tGridListIter< tLNode > nodIter( gridPtr->GetNodeList() );
   tPtrListIter< tEdge > spokIter;
   tLNode * cn, *pointtodelete;
   tEdge *ce;
     //check for proximity to channel:
   tooclose = TRUE;
   do
   {
      tooclose = FALSE;
        //nodIter.First();
      for( cn = nodIter.FirstP(); !(nodIter.AtEnd()); cn = nodIter.NextP() )
      {
         pointtodelete = 0;
         if( cn->Meanders() )
         {
            spokIter.Reset( cn->getSpokeListNC() );
            for( ce = spokIter.FirstP(); !( spokIter.AtEnd() );
                 ce = spokIter.NextP() )
            {
               if( ce->getLength() < cn->getHydrWidth()/2.0 )
               {
                  cout<<"too close: cn, cn->hydrwidth "<<cn->getID()<<" "
                      <<cn->getHydrWidth()<<endl<<flush;
                  tooclose = TRUE;
                  pointtodelete = (tLNode *) ce->getDestinationPtrNC();
                  if( pointtodelete->getDrArea() > cn->getDrArea() )
                      pointtodelete = cn;
                  break;
               }
            }
         }
         if( pointtodelete != 0 ) gridPtr->DeleteNode( pointtodelete );
      }
   } while( tooclose );
   cout << "finished" << endl;
}

/*****************************************************************************\
**
**       CheckFlowedgCross(): checks to see whether any nodes have been
**          crossed/swept over by flowedges; i.e., whether the channel has
**          eroded the node away.
**          -finds tri's with !NewTriCCW;
**          -are two vtcs. connected by a flowedg?
**          -if so, is >1 nbr tri !NewTriCCW, or does a spoke of the 3rd vtx.
**           cross the flowedg?
**          -if so, mark the 3rd node for deletion;
**          -if a point is to be deleted, is it a meandering node with greater
**           drainage area than the upstream node of the found flowedg?
**          -if so, mark the upstream node of the flowedg for deletion instead.
**          -if a point is to be deleted, is it a boundary node?
**          -if so, if either/both of the other nodes have moved outside the
**           boundaries, revert it/them to the old coords., and assert that
**           one or the other is indeed outside the boundary, o.w., exit w/err.
**          -if the point to be deleted is not a boundary node, delete it. 
**
**		Parameters:	
**		Called by:	Migrate
**		Created: 1/98 SL; moved/modified from tGrid routine
**
\*****************************************************************************/
void tStreamMeander::CheckFlowedgCross()
{
   cout << "CheckFlowedgCross()..." << flush << endl;
   int i, j, nv, nvopp, id0, id1;
   int ft;
   int flipped;
   int crossed;
   tArray< double > p0, p1, p2, xy0, xy1;
   tLNode *pointtodelete, *lnodePtr, *nod, *cn, *dscn;  
   tEdge * fedg, * cedg, * ccedg, *ce;
   tTriangle * ct, *nt;
   tListIter< tTriangle > triIter( gridPtr->GetTriList() );
   tPtrListIter< tEdge > spokIter;
     //delete node crossed by flowedg:
     //if a new triangle is !CCW and two vtcs. are connected by a flowedg AND
     //  a spoke of the third vtx. intersects the flowedg OR
     //  more than one nbr. tri. is also !CCW
     //then delete third vtx. because it has been "crossed" by the flowedg
   crossed = 1;
   do
   {
      crossed = 0;
      for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
      {
         if( !NewTriCCW( ct ) )
         {
            ft = 0;
            pointtodelete = 0;
            for( i=0; i<3; i++ )
            {
               cn = (tLNode *) ct->pPtr(i);
               dscn = cn->GetDownstrmNbr();
               for( j=1; j<3; j++ )
               {
                  nod = (tLNode *) ct->pPtr( (i+j)%3 );
                    //if another vtx is cn's downstream nbr
                  if( nod == dscn )
                  {
                     ft = 1;
                     if( j == 1 ) nod = (tLNode *) ct->pPtr( (i+2)%3 );
                     else nod = (tLNode *) ct->pPtr( (i+1)%3 );
                  }
                  if( ft ) break;
               }
               if( ft ) break;
            }
            for( i=0, j=0; i<3; i++ )
            {
               nt = ct->tPtr(i);
               if( nt != 0 )
                   if( !NewTriCCW( nt ) )
                       j++;
            }
              //if >1 nbr tri !CCW
            if( j>1 )
            {
               pointtodelete = nod;
            }
            else
            {
               spokIter.Reset( nod->getSpokeListNC() );
               fedg = cn->GetFlowEdg();
               for( ce = spokIter.FirstP(); !( spokIter.AtEnd() );
                    ce = spokIter.NextP() )
               {
                    //if flowedg crossed by spoke(s) of the third vtx.
                  if( Intersect( ce, fedg ) )
                  {
                     pointtodelete = nod;
                     break;
                  }
               }
            }
            if( pointtodelete != 0 )
            {
               crossed = 1;
               if( pointtodelete->Meanders() &&
                   cn->getDrArea() < pointtodelete->getDrArea() )
               {
                  pointtodelete = cn;
               }
               else if( pointtodelete->getBoundaryFlag() )
               {
                  pointtodelete = 0;
               }
               if( pointtodelete != 0 )
               {
                  gridPtr->DeleteNode( pointtodelete );
               }
               else
               {
                  xy0 = cn->getNew2DCoords();
                  if( gridPtr->LocateTriangle( xy0[0], xy0[1] ) == 0 )
                      cn->RevertToOldCoords();
                  xy1 = dscn->getNew2DCoords();
                  if( gridPtr->LocateTriangle( xy1[0], xy1[1] ) == 0 )
                      dscn->RevertToOldCoords();
                  assert( gridPtr->LocateTriangle( xy0[0], xy0[1] ) == 0 ||
                          gridPtr->LocateTriangle( xy1[0], xy1[1] ) == 0 );
               }
            }
         }
      }
   } while( flipped );
   cout << "finished" << endl;
}
      
/*****************************************************************************\
**
**         CheckBrokenFlowedg()
**
\*****************************************************************************/
void tStreamMeander::CheckBrokenFlowedg()
{
   cout << "CheckBrokenFlowedg()..." << flush << endl;
   int i;
   int breakedge;
   int change = 0;
   int flipped;
   int flip = FALSE;
   double dis0, dis1;
   tArray< int > npop(3);
   tLNode *pointtodelete, *lnodePtr, * pedg[2];  
   tEdge * fedg, *tedg[2];
   tTriangle * ct, * trop[3];
   tListIter< tTriangle > triIter( gridPtr->GetTriList() );
     //check for broken flow edges that we don't want broken
   breakedge = TRUE;
   tPtrList< tTriangle > triPtrList;
   do
   {
      breakedge = FALSE;
        //cout << "breaking edge loop: form tri stack" << endl << flush;
      for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
      {
         change = FALSE;
         for( i = 0; i < 3; i++ )
         {
            lnodePtr = (tLNode *) ct->pPtr(i);
            if( lnodePtr->Meanders() )
            {
               change = TRUE;
                 //cout << " add triangle " << ct->getID() << " to ptrlist" << endl;
               break;
            }
         }
         if( change ) triPtrList.insertAtBack( ct );
      }
        //for each triangle in the stack
      tPtrListIter< tTriangle > triPtrIter( triPtrList );
      for( ct = triPtrIter.FirstP(); !( triPtrIter.AtEnd() );
           ct = triPtrIter.NextP() )
      {
         assert( ct != 0 );
           //cout << " check triangle " << ct->getID() << endl;
         for( i=0; i<3; i++ )
         {
            trop[i] = ct->tPtr(i);
            if( trop[i] ) npop[i] = trop[i]->nVOp( ct );
            else npop[i] = NULL;
         }
         for( i=0; i<3; i++ )
         {
            if( gridPtr->CheckForFlip( trop[i], npop[i], flip ) )
            {
               pedg[0] = (tLNode *) ct->pPtr( (i+1)%3 );
               pedg[1] = (tLNode *) ct->pPtr( (i+2)%3 );
               tedg[0] = ct->ePtr( (i+1)%3 );
               tedg[1] = gridPtr->getEdgeCompliment( tedg[0] );
               fedg = 0;
               if( pedg[0]->GetFlowEdg() )
                   if( tedg[0]->getID() == pedg[0]->GetFlowEdg()->getID() )
                       fedg = tedg[0];
               if( pedg[1]->GetFlowEdg() )
                   if( tedg[1]->getID() == pedg[1]->GetFlowEdg()->getID() )
                       fedg = tedg[1];
               if( fedg )
               {
                  lnodePtr = (tLNode *) ct->pPtr(i);
                  dis0 = lnodePtr->DistNew( pedg[0], pedg[1] );
                  lnodePtr = (tLNode *) trop[i]->pPtr( npop[i] );
                  dis1 = lnodePtr->DistNew( pedg[0], pedg[1] );
                  if( dis0 < dis1 ) pointtodelete = (tLNode *) ct->pPtr(i);
                  else pointtodelete = (tLNode *) trop[i]->pPtr( npop[i] );
                  lnodePtr = (tLNode *) fedg->getOriginPtr();
                  if( pointtodelete->getDrArea() > lnodePtr->getDrArea() )
                      pointtodelete = NULL;
                  else if( pointtodelete->getBoundaryFlag() )
                  {
                     for( i=0; i<2; i++ )
                     {
                        pedg[i]->RevertToOldCoords();
                     }
                  }
                  else
                  {
                       //DumpTriangles();
                       //DumpNodes();
                       //DumpEdges();
                     cout<<"delete flow edge breaker"<<endl<<flush;
                     gridPtr->DeleteNode( pointtodelete );
                     breakedge = TRUE;
                     break;
                  }
               }  
            }
         }
         if( breakedge ) break;
      }
      triPtrList.Flush();
   } while( breakedge );
   cout << "finished" << endl;
}



/****************************************************************\
**	
**
**		Parameters:	
**		Called by:	
**		Created: 1/98 SL
**
\***************************************************************/
