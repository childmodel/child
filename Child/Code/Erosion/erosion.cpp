/***************************************************************************\
**
**  erosion.cpp
**
**  Functions for sediment transport and bed erosion (detachment) objects.
**  Transport objects:
**    tSedTransPwrLaw
**  Detachment objects:
**    tBedErodePwrLaw
**
**    Created 1/98 gt
**
**  $Id: erosion.cpp,v 1.8 1998-03-09 22:48:30 gtucker Exp $
\***************************************************************************/

#include <math.h>
#include <assert.h>
#include "erosion.h"

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
**  Assumptions: n->GetSlope() does not return a negative value; kb, mb,
**               and nb all >=0.
\***************************************************************************/
double tBedErodePwrLaw::DetachCapacity( tLNode * n, double dt )
{
   double slp = n->GetSlope();
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
**  Assumptions: n->GetSlope() does not return a negative value; kb, mb,
**               and nb all >=0.
\***************************************************************************/
double tBedErodePwrLaw::DetachCapacity( tLNode * n )
{
   double slp = n->GetSlope();
   double erorate =  kb*pow( n->GetQ(), mb )*pow( slp, nb );
   n->SetDrDt( -erorate );
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
   assert( n->GetQ()>=0 );
   double eroterm = kb * pow( n->GetQ(), mb ) * pow( n->GetSlope(), nb-1.0 );
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
   double cap = 0;
   if( !node->GetFloodStatus() )
       cap = kf * pow( node->GetQ(), mf ) * pow( node->GetSlope(), nf );
   node->SetQs( cap );
   return cap;
}
   

/***************************************************************************\
**  tErosion functions
\***************************************************************************/
tErosion::tErosion( tGrid<tLNode> *gptr, tInputFile &infile )
        : bedErode( infile ), sedTrans( infile )
{
   assert( gptr!=0 );
   gridPtr = gptr;
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
          cn->SetQs( -bedErode.DetachCapacity( cn ) );
      dtmax = dtg;
      //find max. time step s.t. slope does not reverse:
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      {
         dn = cn->GetDownstrmNbr();
         ratediff = dn->GetQs() - cn->GetQs();
         if( ratediff > 0 )
         {
            dt = ( cn->getZ() - dn->getZ() ) / ratediff * frac;
            if( dt > 0 && dt < dtmax ) dtmax = dt;
         }
      }
      //assert( dtmax > 0 );
      //apply erosion:
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
          cn->EroDep( cn->GetQs() * dtmax );
      //update time:
      dtg -= dtmax;
   } while( dtg>0 );
   
}

   
void tErosion::ErodeDetachLim( double dtg, tUplift *UPtr )
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
          cn->SetQs( -bedErode.DetachCapacity( cn ) + UPtr->GetRate() );
      dtmax = dtg;
      //find max. time step s.t. slope does not reverse:
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      {
         dn = cn->GetDownstrmNbr();
         ratediff = dn->GetQs() - cn->GetQs();
         if( ratediff > 0 && cn->getZ() > dn->getZ() )
         {
            dt = ( cn->getZ() - dn->getZ() ) / ratediff * frac;
            if( dt > 0 && dt < dtmax ) dtmax = dt;
         }
      }
      //assert( dtmax > 0 );
      //apply erosion:
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
          cn->EroDep( cn->GetQs() * dtmax );
      //update time:
      dtg -= dtmax;
   } while( dtg>0 );
   
}


void tErosion::StreamErode( double dtg, tStreamNet *strmNet )
{
   double dt,
       dtmax;         // time increment: initialize to arbitrary large val
   double frac = 0.3; //fraction of time to zero slope
   int i;
   tLNode * cn, *dn;
   int nActNodes = gridPtr->GetNodeList()->getActiveSize();
   tGridListIter<tLNode> ni( gridPtr->GetNodeList() );
   double ratediff,  // Difference in ero/dep rate btwn node & its downstrm nbr
       cap,          // Transport capacity
       pedr,         // Potential erosion/deposition rate
       dcap,         // Bedrock detachment capacity
       dzs,          // Rate of alluvium ero/dep
       dzr;          // Rate of rock ero/dep

   // Sort so that we always work in upstream to downstream order
   strmNet->SortNodesByNetOrder();

   // Compute erosion and/or deposition until all of the elapsed time (dtg)
   // is used up
   do
   {
      cout << "Remaing time: " << dtg << endl;
      
      // Zero out sed influx
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
          cn->SetQsin( 0.0 );

      // Compute erosion rates: when this block is done, the transport rate
      // (qs), influx (qsin), and deposition/erosion rate (dzdt) values are
      // set for each active node.
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      {
         // Transport capacity and potential erosion/deposition rate
         // (this also sets the node's Qs value)
         cap = sedTrans.TransCapacity( cn );
         pedr = (cn->GetQsin() - cap ) / cn->getVArea();
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
         cn->SetDzDt( pedr );
         cn->GetDownstrmNbr()->AddQsin( cn->GetQsin() - pedr*cn->getVArea() );
      }

      // Given these rates, figure out how big a time-step we can get away with
      //  (Note: the division and subsequent multiplication by frac are done
      //   for performance reasons: we avoid having to multiply every dt value
      //   by frac)
      dtmax = dtg/frac;
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      {
         dn = cn->GetDownstrmNbr();
         ratediff = dn->GetDzDt() - cn->GetDzDt(); // Are the pts converging?
         if( ratediff > 0 && cn->getZ() > dn->getZ() )  // if yes, get time
         {                                              //  to zero slope
            dt = ( cn->getZ() - dn->getZ() ) / ratediff;
            if( dt < dtmax ) dtmax = dt;
            if( dt < 1e-6 )
            {
                cout << "Very small dt " << dt << " at:\n";
                cn->TellAll();
                dn->TellAll();
            }
            
         }
      }
      dtmax *= frac;  // Take a fraction of time-to-flattening

      // Zero out sed influx again, because depending on bedrock-alluvial
      // interaction it may be modified
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
          cn->SetQsin( 0.0 );

      // Now do erosion/deposition by integrating rates over dtmax
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      {
         // If we're on bedrock, scour the bedrock and "inject" the resulting
         // sediment flux into the node's sed influx field (qsin)
         if( cn->OnBedrock() )
             dzr = cn->GetDrDt()*dtmax;
         dzs = ( (cn->GetQsin() - cn->GetQs() ) / cn->getVArea() ) * dtmax;
         
         // If potential erosion depth is greater than available depth of
         // alluvium, limit erosion to alluvium plus bedrock scour, and
         // adjust sediment outflux accordingly.
         if( -dzs > cn->getAlluvThickness() )
         {
            cout << "(Following node is br-limited)\n";
             dzs = -cn->getAlluvThickness();
             cn->SetQs( cn->GetQsin() - ((dzr+dzs)*cn->getVArea())/dtmax );
         }

         cout << "** THIS node has dzs " << dzs << " & dzr " << dzr << endl;
         cn->TellAll();

         // Update alluvium thickness and node elevation
         cn->setAlluvThickness( cn->getAlluvThickness() + dzs );
         cn->setZ( cn->getZ() + dzs + dzr );
         dn = cn->GetDownstrmNbr();
         dn->AddQsin( cn->GetQs() );
      }

      // Update time remaining
      dtg -= dtmax;
      
   } while( dtg>1e-6 ); // Keep going until we've used up the whole time intrvl
   
}

   
         
               
         
