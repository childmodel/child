/***************************************************************************\
**
**  erosion.cpp
**
**  Functions for equilibrium checking, sediment transport and bed erosion
**    (detachment) objects.
**  Equilibrium check objects:
**    tEquilibCheck
**  Transport objects:
**    tSedTransPwrLaw
**  Detachment objects:
**    tBedErodePwrLaw
**
**    Created 1/98 gt; add tEqChk 5/98 sl
**
**  $Id: erosion.cpp,v 1.28 1998-05-14 14:21:31 gtucker Exp $
\***************************************************************************/

#include <math.h>
#include <assert.h>
#include "erosion.h"

/***************************************************************************\
**  tEquilibCheck
**
**  Functions to find and report system mass change rate.
**
**  5/98 SL
\***************************************************************************/
/***************************************************************************\
**  Constructors: default (no args) and given grid and timer as args
\***************************************************************************/
tEquilibCheck::tEquilibCheck()
        : massList()
{
   gridPtr = 0;
   timePtr = 0;
   longTime = longRate = shortRate = 0.0;
}

tEquilibCheck::tEquilibCheck( tGrid< tLNode > &gridRef, tRunTimer &timeRef )
        : massList()
{
   gridPtr = &gridRef;
   timePtr = &timeRef;
   longTime = longRate = shortRate = 0.0;
   FindIterChngRate();
}

tEquilibCheck::tEquilibCheck( tGrid< tLNode > &gridRef, tRunTimer &timeRef,
                              tInputFile &fileRef )
        : massList()
{
   gridPtr = &gridRef;
   timePtr = &timeRef;
   longRate = shortRate = 0.0;
   longTime = fileRef.ReadItem( longTime, "EQUITIME" );
   FindIterChngRate();
}

tEquilibCheck::~tEquilibCheck()
{
   gridPtr = 0;
   timePtr = 0;
}

/***************************************************************************\
**  'get' and 'set' functions for tEquilibCheck:
\***************************************************************************/
double tEquilibCheck::getLongTime() const {return longTime;}

void tEquilibCheck::setLongTime( double val )
{longTime = ( val > 0 ) ? val : 0.0;}

const tGrid< tLNode > *tEquilibCheck::getGridPtr() const {return gridPtr;}

tGrid< tLNode > *tEquilibCheck::getGridPtrNC() {return gridPtr;}

void tEquilibCheck::setGridPtr( tGrid< tLNode > &Ref )
{gridPtr = ( &Ref > 0 ) ? &Ref : 0;}

void tEquilibCheck::setGridPtr( tGrid< tLNode > *Ptr )
{gridPtr = ( Ptr > 0 ) ? Ptr : 0;}

const tRunTimer *tEquilibCheck::getTimePtr() const {return timePtr;}

tRunTimer *tEquilibCheck::getTimePtrNC() {return timePtr;}

void tEquilibCheck::setTimePtr( tRunTimer &Ref )
{timePtr = ( &Ref > 0 ) ? &Ref : 0;}

void tEquilibCheck::setTimePtr( tRunTimer *Ptr )
{timePtr = ( Ptr > 0 ) ? Ptr : 0;}

double tEquilibCheck::getLongRate() const {return longRate;}

double tEquilibCheck::getShortRate() const {return shortRate;}

/***************************************************************************\
**  tEquilibCheck::FindIterChngRate()
**
**  Find average rate of elevation change since the last time the object was
**  called, as short a time as one model iteration.
**  Finds new average elevation, new time, and calculates rate based on last
**  avg. elev. and time.
\***************************************************************************/
double tEquilibCheck::FindIterChngRate()
{
   assert( timePtr > 0 && gridPtr > 0 );
   tArray< double > tmp(2), last;
   tmp[0] = timePtr->GetCurrentTime();
   tGridListIter< tLNode > nI( gridPtr->GetNodeList() );
   tListIter< tArray< double > > mI( massList );
   tLNode *cn;
   double mass = 0.0;
   double area = 0.0;
   for( cn = nI.FirstP(); nI.IsActive(); cn = nI.NextP() )
   {
      mass += cn->getZ() * cn->getVArea();
      area += cn->getVArea();
   }
   tmp[1] = mass / area;
   if( !(massList.isEmpty()) )
   {
      last = *(mI.LastP());
      double dt = (tmp[0] - last[0]);
      assert( dt > 0.0 );
      shortRate = (tmp[1] - last[1]) / dt;
   }
   else
   {
        //cout << "tEquilibCheck::FindIterChngRate(), Warning: empty massList\n";
      assert( tmp[0] > 0 );
      shortRate = tmp[1] / tmp[0];
   }
   massList.insertAtBack( tmp );
   return shortRate;
}

/***************************************************************************\
**  tEquilibCheck::FindLongTermChngRate()
**
**  Find average rate of elevation change over a set time, as short as one
**  model iteration.
**  Calls FindIterChngRate() to update massList, then searches massList for
**  time to provide rate over a time >= longTime.
\***************************************************************************/
double tEquilibCheck::FindLongTermChngRate()
{
   FindIterChngRate();
   tListIter< tArray< double > > mI( massList );
   tArray< double > last = *(mI.LastP());
   tArray< double > ca, na;
   double dt, targetTime = last[0] - longTime;
   if( longTime == 0.0 || mI.FirstP() == mI.LastP() ) longRate = shortRate;
   else
   {
      ca = *(mI.FirstP());
      na = *(mI.NextP());
      while( na[0] < targetTime && !(mI.AtEnd()) )
      {
         ca = na;
         na = *(mI.NextP());
      }
      dt = last[0] - ca[0];
      assert( dt > 0 );
      longRate = (last[1] - ca[1]) / dt;
   }
   return longRate;
}

/***************************************************************************\
**  tEquilibCheck::FindLongTermChngRate( double newtime )
**
**  Set longTime = newtime and call FindLongTermChngRate()
\***************************************************************************/
double tEquilibCheck::FindLongTermChngRate( double newtime )
{
   setLongTime( newtime );
   return FindLongTermChngRate();
}

   
/***************************************************************************\
**  Constructors:
**    Given a tInputFile as an argument, will read relevant parameters from
**  the input file.
\***************************************************************************/
tSedTransPwrLaw::tSedTransPwrLaw( tInputFile &infile )
{
   kf = infile.ReadItem( kf, "KF" );
   mf = infile.ReadItem( mf, "MF" );
   nf = infile.ReadItem( nf, "NF" );
}

//tBedErodePwrLaw::tBedErodePwrLaw()
//{ kb = 0; mb=0; nb=0; }

tBedErodePwrLaw::tBedErodePwrLaw( tInputFile &infile )
{
   kb = infile.ReadItem( kb, "KB" );
   mb = infile.ReadItem( mb, "MB" );
   nb = infile.ReadItem( nb, "NB" );
}



/***************************************************************************\
**  tBedErodePwrLaw functions
\***************************************************************************/

/***************************************************************************\
**  tBedErode::DetachCapacity
**
**  Computes the depth of erosion over a time interval dt assuming the
**  erosion rate = kb Q^mb S^nb
**
**  Input: n -- node at which to compute detachment capacity
**         dt -- time interval
**  Returns: the detachment depth
**  Assumptions: n->GetSlope() does not return a negative value (returns neg.
**               only if infinite loop in GetSlope()); kb, mb,
**               and nb all >=0.
\***************************************************************************/
double tBedErodePwrLaw::DetachCapacity( tLNode * n, double dt )
{
   double slp = n->GetSlope();
   if( slp < 0.0 )
       ReportFatalError("neg. slope in tBedErodePwrLaw::DetachCapacity(tLNode*,double)");
   return( kb*pow( n->GetQ(), mb )*pow( slp, nb )*dt );
}

/***************************************************************************\
**  tBedErode::DetachCapacity
**
**  Computes the rate of erosion  = kb Q^mb S^nb
**
**  Input: n -- node at which to compute detachment capacity
** 
**  Returns: the detachment rate
**  Assumptions: n->GetSlope() does not return a negative value (returns neg.
**               only if infinite loop in GetSlope()); kb, mb,
**               and nb all >=0.
\***************************************************************************/
double tBedErodePwrLaw::DetachCapacity( tLNode * n )
{
   double slp = n->GetSlope();
   if( slp < 0.0 )
       ReportFatalError("neg. slope in tBedErodePwrLaw::DetachCapacity(tLNode*)");
   double erorate =  kb*pow( n->GetQ(), mb )*pow( slp, nb );
   n->setDrDt( -erorate );
   return erorate;
}

/***************************************************************************\
**  tBedErode::SetTimeStep
**
**  Estimates a maximum time step as a function of the Courant stability
**  criterion: dt <= dx/v, where dx is the node spacing and v is the
**  wave velocity. If nb=1, v = kb Q^mb. If the equation is nonlinear
**  with nb!=1, the wave velocity is estimated as if the equation were
**  linear with a gradient term in the coefficient, ie
**      v S = [kb Q^mb S^(nb-1)] S  (recall S = -dz/dx)
**  The estimated time step limitation is therefore given by:
**      dt = 0.2 * ( dx / (kb Q^mb S^(nb-1) ) )
**  (the 0.2 is a factor to keep us comfortably below the Courant #).
**  If the denominator is zero, an arbirarily large number is returned.
**
**  Input: n -- node for which to estimate time step
**  Returns: the estimated maximum time step size
**  Assumptions: GetSlope() returns a value >=0, edge length>0.
**
\***************************************************************************/
double tBedErodePwrLaw::SetTimeStep( tLNode * n )
{
   double slp = n->GetSlope();
   if( slp < 0.0 )
       ReportFatalError("neg. slope in tBedErodePwrLaw::SetTimeStep(tLNode*)");
   assert( n->GetQ()>=0 );
   double eroterm = kb * pow( n->GetQ(), mb ) * pow( slp, nb-1.0 );
   if( eroterm==0 ) return 100000;
   return( 0.2*n->GetFlowEdg()->getLength() / eroterm );

}

/***************************************************************************\
**  tSedTransPwrLaw functions
\***************************************************************************/

/***************************************************************************\
**  tSedTransPwrLaw::TransCapacity
**
**  Computes sediment transport capacity using the simple power law
**  Qs = kf Q^mf S^nf
\***************************************************************************/
double tSedTransPwrLaw::TransCapacity( tLNode *node )
{
   double slp = node->GetSlope();
   if( slp < 0.0 )
       ReportFatalError("neg. slope in tBedErodePwrLaw::TransCapacity(tLNode*)");
   double cap = 0;
   if( !node->GetFloodStatus() )
       cap = kf * pow( node->GetQ(), mf ) * pow( slp, nf );
   node->setQs( cap );
   return cap;
}
   
/*************************************************************************\
** tSedTransWilcock functions
\**************************************************************************/

tSedTransWilcock::tSedTransWilcock( tInputFile &infile )
        : grade()
{
   int i;
   char add, name[20];
   double help, sum;

   cout << "tSedTransWilcock(infile)\n";
   add='1';
   for(i=1; i<=2; i++){
      strcpy( name, "GRAINDIAM");
      strcat( name, &add );
      help = infile.ReadItem( help, name);
      grade[i] = help;
   }

   taudim= RHO*GRAV;
   refs = (RHOSED-RHO)*9.81*grade[0];
   refg = (RHOSED-RHO)*9.81*grade[1];
   lowtaucs = 0.8*(grade[1]/grade[0])*0.040*refs*0.8531;
   lowtaucg = 0.04*refg*0.8531;;
   sandb = (lowtaucs-(0.04*10./40.))/(1-(10./40.));
   hightaucs = 0.04*refs*0.8531;
   hightaucg = 0.01*refg*0.8531; 
   sands = (lowtaucs-0.04)/(-30.);
}

double tSedTransWilcock::TransCapacity( tLNode *nd )
{
   double tau;
   double taucrit;
   double persand=grade[0]/(grade[0]+grade[1]);
   double timeadjust=31536000.00; /* number of seconds in a year */

   if( nd->GetSlope() < 0 ){
      nd->setQs(0, 0);
      nd->setQs(1, 0);
      nd->setQs(0);
      return 0.0;
   }

   taudim *= pow(nd->getHydrWidth(), 0.6);
   // nic, make sure about the units here ie units on Q
   tau = taudim*pow(nd->GetQ()*timeadjust,0.3)*pow( nd->GetSlope(), 0.7);

   
   if(persand<.10)
       taucrit=lowtaucs;
   else if(persand<.40)
       taucrit=((sands*persand*100)+sandb)*refs*0.8531;
   else
       taucrit=hightaucs;
   
   if(tau>taucrit)
       nd->setQs(0, ( 1829088*persand*pow(tau,1.5)*pow((1-sqrt(taucrit/tau)),4.5) ));
   else 
       nd->setQs( 0, 0.0 ) ;
   
       
}


/***************************************************************************\
**  tErosion functions
\***************************************************************************/
tErosion::tErosion( tGrid<tLNode> *gptr, tInputFile &infile )
        : bedErode( infile ), sedTrans( infile )
{
   assert( gptr!=0 );
   gridPtr = gptr;

   // Read parameters needed from input file
   kd = infile.ReadItem( kd, "KD" );  // Hillslope diffusivity coefficient
   
}


/*****************************************************************************\
**
**  ErodeDetachLim
**
**  Solves for erosion and deposition during a time interval dtg, for the
**  case in which any sediment detached from the landscape is assumed to
**  be carried away. This might be appropriate, for example, in a bedrock
**  channel system in which deposition is negligible. This case is handled
**  separately from the more general case in which transport capacity and
**  deposition are taken into account, because the numerical solutions
**  to detachment-limited erosion equations tend to be considerably more
**  stable than the general case.
**
**  The function will solve the erosion equation(s) over one or more time
**  intervals within the total "global" time period dtg. First, the
**  maximum time step size dt is computed. Then erosion is computed
**  iteratively in increments of dt until the total time dtg has been used.
**
**  Modified, 2/20/98, SL & GT: maximize time step such that slope in down-
**   stream direction does not reverse; don't use modified Courant condition.
**   This has the advantage that it does not grind to a halt when n<1; but, it
**   is not completely satisfactory; it seems to be good at finding out how
**   small the time step needs to be but not how large it can be; i.e., if
**   dtg is quite large, one can still run into problems. E.g., for m=0.3,
**   n=0.7, kb=3.16e-02 m^3/N, U=1e-3 m/yr, Q=100m^3/s, kwds=3, knds=0.03
**   => dtg=10 yrs seems to work well; dtg=100 yrs does not work well.
\*****************************************************************************/
void tErosion::ErodeDetachLim( double dtg )
{
   double dt,
       dtmax = 1000000.0; // time increment: initialize to arbitrary large val
   double frac = 0.9; //fraction of time to zero slope
   int i;
   tLNode * cn, *dn;
   int nActNodes = gridPtr->GetNodeList()->getActiveSize();
   tGridListIter<tLNode> ni( gridPtr->GetNodeList() );
   tArray<double> //dz( nActNodes ), // Erosion depth @ each node
       dzdt( nActNodes ); //Erosion rate @ ea. node
   double ratediff;

   // Iterate until total time dtg has been consumed
   do
   {
      //first find erosion rate:
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
          cn->setQs( -bedErode.DetachCapacity( cn ) );
      dtmax = dtg;
      //find max. time step s.t. slope does not reverse:
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      {
         dn = cn->GetDownstrmNbr();
         ratediff = dn->getQs() - cn->getQs();
         if( ratediff > 0 )
         {
            dt = ( cn->getZ() - dn->getZ() ) / ratediff * frac;
            if( dt > 0 && dt < dtmax ) dtmax = dt;
         }
      }
      //assert( dtmax > 0 );
      //apply erosion:
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
          cn->EroDep( cn->getQs() * dtmax );
      //update time:
      dtg -= dtmax;
   } while( dtg>0 );
   
}

   
void tErosion::ErodeDetachLim( double dtg, tUplift *UPtr )
{
   double dt,
       dtmax = 1000000.0; // time increment: initialize to arbitrary large val
   double frac = 0.1; //fraction of time to zero slope
   int i;
   tLNode * cn, *dn;
   int nActNodes = gridPtr->GetNodeList()->getActiveSize();
   tGridListIter<tLNode> ni( gridPtr->GetNodeList() );
   tArray<double> //dz( nActNodes ), // Erosion depth @ each node
       dzdt( nActNodes ); //Erosion rate @ ea. node
   double ratediff;
   double slp, dslpdt;
   double dtmin = dtg * 0.0001;

   // Iterate until total time dtg has been consumed
   do
   {
      //first find erosion rate:
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
          cn->setDzDt( -bedErode.DetachCapacity( cn ) );
      dtmax = dtg;
      //find max. time step s.t. slope does not reverse:
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      {
         /*slp = cn->GetSlope();
         dslpdt = cn->GetDSlopeDt();
         if( slp > 0.0 )
         {
            if( dslpdt < 0.0 )
            {
               dt = slp / (-dslpdt - UPtr->GetRate() ) * frac;
               if( dt > 0 && dt < dtmax ) dtmax = dt;
            }
         }*/
         
         dn = cn->GetDownstrmNbr();
         if( dn->getBoundaryFlag() == kNonBoundary )
             ratediff = dn->getDzDt() - cn->getDzDt();
         else
             ratediff = dn->getDzDt() - cn->getDzDt() - UPtr->GetRate();
         if( ratediff > 0 && cn->getZ() > dn->getZ() )
         {
            dt = ( cn->getZ() - dn->getZ() ) / ratediff * frac;
            if( dt > dtmin && dt < dtmax )
            {
               dtmax = dt;
            }
            else
            {
               dtmax = dtmin;
               //cout << "time step too small because of node at x,y,z "
               //     << cn->getX() << " " << cn->getY() << " " << cn->getZ()
               //     << endl;
            }
         }
      }
      //assert( dtmax > 0 );
      //apply erosion:
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
          cn->EroDep( cn->getDzDt() * dtmax );
      //update time:
      dtg -= dtmax;
   } while( dtg>0 );
   
}


#define kSmallTimeStep 1e-6
void tErosion::StreamErode( double dtg, tStreamNet *strmNet )
{
   double dt,
       dtmax;         // time increment: initialize to arbitrary large val
   double frac = 0.3; // fraction of time to zero slope
   int i;
   tLNode * cn, *dn;
   int nActNodes = gridPtr->GetNodeList()->getActiveSize();
   tGridListIter<tLNode> ni( gridPtr->GetNodeList() );
   double ratediff,  // Difference in ero/dep rate btwn node & its downstrm nbr
       cap,          // Transport capacity
       pedr,         // Potential erosion/deposition rate
       dcap,         // Bedrock detachment capacity
       dz,           // Depth of deposition/erosion (erosion = negative)
       dzr;          // Potential depth of bedrock erosion

   //cout << "tErosion::StreamErode\n";

   // Sort so that we always work in upstream to downstream order
   strmNet->SortNodesByNetOrder();

   // Compute erosion and/or deposition until all of the elapsed time (dtg)
   // is used up
   do
   {
      // Zero out sed influx
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
          cn->setQsin( 0.0 );

      // Compute erosion rates: when this block is done, the transport rate
      // (qs), influx (qsin), and deposition/erosion rate (dzdt) values are
      // set for each active node.
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      {
         // Transport capacity and potential erosion/deposition rate
         // (this also sets the node's Qs value)
         cap = sedTrans.TransCapacity( cn );
         pedr = (cn->getQsin() - cap ) / cn->getVArea();
        //sediment input:
         if( cn == strmNet->getInletNodePtr() )
             pedr += strmNet->getInSedLoad() / cn->getVArea();
         
         // If we're on bedrock, adjust accordingly
         if( cn->OnBedrock() && pedr<0 )
         {
            // Get detachment capacity (this also sets node's drdt)
            dcap = -bedErode.DetachCapacity( cn );
            if( dcap > pedr )
                pedr = dcap;
         }
         // Set the erosion (deposition) rate and send the corresponding
         // sediment influx downstream
         cn->setDzDt( pedr );
         cn->GetDownstrmNbr()->AddQsin( cn->getQsin() - pedr*cn->getVArea() );
         //cout << "RATE STEP:\n";
         //cn->TellAll();
      }

      // Given these rates, figure out how big a time-step we can get away with
      //  (Note: the division and subsequent multiplication by frac are done
      //   for performance reasons: we avoid having to multiply every dt value
      //   by frac)
      dtmax = dtg/frac;
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      {
         dn = cn->GetDownstrmNbr();
         ratediff = dn->getDzDt() - cn->getDzDt(); // Are the pts converging?
         if( ratediff > 0 && cn->getZ() > dn->getZ() )  // if yes, get time
         {                                              //  to zero slope
            dt = ( cn->getZ() - dn->getZ() ) / ratediff;
            if( dt < dtmax ) dtmax = dt;
            /*if( dt < kSmallTimeStep )
            {
               cout << "Very small dt " << dt << " at:\n";
                cn->TellAll();
                dn->TellAll();
            }*/
            
         }
      }
      dtmax *= frac;  // Take a fraction of time-to-flattening
      if( dtmax < kSmallTimeStep ) dtmax = kSmallTimeStep;
      //cout << "  dt " << dtmax << endl;

      // Zero out sed influx again, because depending on bedrock-alluvial
      // interaction it may be modified; if inlet, give it strmNet->inlet.inSedLoad
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
          cn->setQsin( 0.0 );
        //sediment input:
      if(strmNet->getInletNodePtrNC() != NULL)
          strmNet->getInletNodePtrNC()->setQsin( strmNet->getInSedLoad() );

      // Notes for multi-size adaptation:
      // qs, qsin, dz, etc could be arrays with dimensions (1..NUMG+1) with
      // the extra entry storing the total. "on bedrock" might be defined as
      // active layer depth less than its "normal" depth (even zero).
      // For bedrock, crit shear might become maximum because of protrusion
      // over the bed --- how to handle?
      // Bedrock scour could automatically generate a given distribution of
      // sizes (specified as a parameter)
      // The basic rule here, as in golem, is don't erode more bedrock than
      // you have capacity to carry. but what happens when you have plenty of
      // surplus capacity in one size and none (or deposition) in another?
      // Could use total capacity: limit TOTAL br erosion to that allowed by
      // TOTAL excess capacity. If some material is generated that can't be
      // carried, just put it in the active layer.
      // One way to proceed would be to _always_ have scour when the bed is
      // exposed, and count the scoured material as part of the sediment
      // influx. If the influx is too large for the avail capacity, just
      // leave it in the active layer. Probably ok as long as time steps
      // remain small and br erosion rate isn't too large relative to capacity.
      // In this case, algo might be:
      // - if on bedrock, erode br and "inject" resulting sed into Qsin,
      //   remembering depth of br erosion
      // - compute surplus/deficit capacity in each size frac. If erosion
      //   amount in a size exceed AL thickness, limit it.
      // - do mass balance on each size to get dz for each size and total
      // - apply br erosion depth
      // this would mean you could have simultaneous deposition _and_ br
      // erosion, but that's not totally unreasonable and would be unlikely
      // to last very long.
      
      // Now do erosion/deposition by integrating rates over dtmax
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      {
         // Depth of potential erosion due to excess transport capacity
         // Note: for multiple sizes, dz could be an array (1..NUMG),
         // maybe w/ an extra field for the total dz.
         dz = ( (cn->getQsin() - cn->getQs() ) / cn->getVArea() ) * dtmax;
         
         // If we're on bedrock, scour the bedrock
         if( cn->OnBedrock() && dz<0.0 ) {
            dzr = cn->getDrDt()*dtmax;  // potential bedrock erosion depth
            // If excess capacity depth-equivalent is greater than the depth
            // of sediment available on bed plus eroded from bedrock, limit
            // depth of erosion to alluvium plus bedrock erosion depth
            if( -dz > -dzr+cn->getAlluvThickness() ) {
               dz = dzr-cn->getAlluvThickness(); 
            }
         }
         
         //cout << "** THIS node has dz " << dz << endl << flush;
         //cn->TellAll();

         // Update alluvium thickness and node elevation
         /*if( cn->getID()==3214 )
         {
             cn->TellAll();
             cout << "Applying dz=" << dz << endl;
         }
         if( cn->getZ() < -1 ) {
            cout << "The following node is going below baselevel:\n";
            cn->TellAll();
         }*/
         cn->EroDep( dz );
         dn = cn->GetDownstrmNbr();
         /*if( dn->getID()==3214 )
         {
            cout << "SENDing to 3214: " << cn->getQsin() << " - " << dz*cn->getVArea()/dtmax << endl;
            cn->TellAll();
         }*/

         // Send sediment downstream: sediment flux is equal to the flux in
         // plus/minus rate of erosion/deposition times node area
         assert( dtmax>0 );
         dn->AddQsin( cn->getQsin() - dz*cn->getVArea()/dtmax );
         if( (cn->getQsin() - dz*cn->getVArea()/dtmax) < -0.1 )
         {
            cout << "NEG OUTFLUX! (dz=" << dz << ")\n";
            cout << "outflux: " << cn->getQsin() - dz*cn->getVArea()/dtmax
                 << endl;
            cn->TellAll();
            dn->TellAll();
         }
      }

      // Update time remaining
      dtg -= dtmax;
      
   } while( dtg>1e-6 ); // Keep going until we've used up the whole time intrvl
   
}


/*****************************************************************************\
**
**  tErosion::Diffuse
**
**  Computes slope-dependent mass transfer by hillslope creep-related
**  processes (hillslope diffusion). Volumetric flux of sediment occurs
**  across each Voronoi cell face; the volume rate of sediment transfer
**  between two nodes that share a Voronoi face of length Lv is
**    Fv = Kd * S * Lv,
**  with transfer being in the downhill (positive S) direction. The total
**  transfer during a time step dt is Fv*dt.
**    Because transfer occurs along edges (across Voronoi faces), we can
**  compute the solution most efficiently by calculating the exchange along
**  each edge pair, accumulating the net influx/outflux at each node, and
**  then summing the influx/outflux for each node and dividing by the node's
**  Voronoi area to get the change in elevation.
**    To ensure numerical stability, the routine begins by estimating a
**  maximum stable time-step size for each pair of nodes (ie, each edge)
**  using the Courant condition Dt <= (Le^2)/(2KdLv). The minimum of these
**  is used as the maximum step size. If the estimated maximum is smaller
**  than the "global" time size size (storm+interstorm duration) rt, the
**  solution is computed in a series of sub-steps until all of rt has been
**  used up.
**    The routine includes an option to "turn off" deposition in areas of
**  concave topography (= net deposition), on the assumption that stream
**  erosion would quickly remove such material.
**
**  Parameters:  rt -- time duration over which to compute diffusion
**               noDepoFlag -- if true, material is only eroded, never
**                             deposited
**  Modifies:  node elevations (z); node Qsin (sed influx) is reset and
**             used in the computations
**  Notes:  as of 3/98, does not differentiate between rock and sediment
**
\*****************************************************************************/
#define kVerySmall 1e-6
#define kEpsOver2 0.1
void tErosion::Diffuse( double rt, int noDepoFlag )
{
   tLNode * cn;
   tEdge * ce;
   double volout,  // Sediment volume output from a node (neg=input)
       denom,      // Denominator in Courant number (=Kd*Lve)
       delt,       // Max local step size
       dtmax;      // Max global step size (initially equal to total time rt)
   tGridListIter<tLNode> nodIter( gridPtr->GetNodeList() );
   tGridListIter<tEdge> edgIter( gridPtr->GetEdgeList() );

#if TRACKFNS
   cout << "tErosion::Diffuse()" << endl << flush;
#endif
   
   // Compute maximum stable time-step size based on Courant condition
   // (Note: for a fixed mesh, this calculation only needs to be done once;
   // performance could be improved by having this block only called if
   // mesh has changed since last time through)
   dtmax = rt;  // Initialize dtmax to total time rt
   for( ce=edgIter.FirstP(); edgIter.IsActive(); ce=edgIter.NextP() )
       if( (denom=kd*ce->getVEdgLen() ) > kVerySmall )
       {
          delt = kEpsOver2*(ce->getLength()/denom);
          if( delt < dtmax )
          {
             dtmax = delt;
             /*cout << "TIME STEP CONSTRAINED TO " << dtmax << " AT EDGE:\n";
             ce->TellCoords();*/
          }
       }

   // Loop until we've used up the entire time interval rt
   do
   {
      // Reset sed input for each node for the new iteration
      for( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
          cn->setQsin( 0 );

      // Compute sediment volume transfer along each edge
      for( ce=edgIter.FirstP(); edgIter.IsActive(); ce=edgIter.NextP() )
      {
         volout = kd*ce->CalcSlope()*ce->getVEdgLen()*dtmax;
         // Record outgoing flux from origin
         cn = (tLNode *)ce->getOriginPtrNC();
         cn->AddQsin( -volout );
         // Record incoming flux to dest'n
         cn = (tLNode *)ce->getDestinationPtrNC();
         cn->AddQsin( volout );
         /*cout << volout << " mass exch. from " << ce->getOriginPtr()->getID()
              << " to "
              << ce->getDestinationPtr()->getID()
              << " on slp " << ce->getSlope() << " ve " << ce->getVEdgLen()
              << "\nvp " << ce->getRVtx()[0] << " " << ce->getRVtx()[1] << endl;
         ((tLNode*)ce->getOriginPtr())->TellAll();
         ((tLNode*)ce->getDestinationPtr())->TellAll();
         cout << endl;*/
         
         ce = edgIter.NextP();  // Skip complementary edge
      }
   
      // Compute erosion/deposition for each node
      for( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
      {
             /*cout << "Node " << cn->getID() << " Qsin: " << cn->getQsin() << " dz: " << cn->getQsin() / cn->getVArea() << endl;*/
         if( noDepoFlag )
             if( cn->getQsin() > 0.0 )
                 cn->setQsin( 0.0 );
         cn->EroDep( cn->getQsin() / cn->getVArea() );  // add or subtract net flux/area
         //cout << cn->z << "Q: " << cn->q << "VA " << cn->varea << endl;
         /*if( cn->id==700 ) {
           cn->TellAll();
           }*/
      }

      rt -= dtmax;
      if( dtmax>rt ) dtmax=rt;
      
   } while( rt>0.0 );
   
   
}




   
         
               
         
