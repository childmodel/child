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
**    tSedTransWilcock
**  Detachment objects:
**    tBedErodePwrLaw
**
**    Created 1/98 gt; add tEqChk 5/98 sl
**
**  $Id: erosion.cpp,v 1.32 1998-07-10 21:34:32 gtucker Exp $
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
   tmp[0] = timePtr->getCurrentTime();
   tGridListIter< tLNode > nI( gridPtr->getNodeList() );
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
**  set longTime = newtime and call FindLongTermChngRate()
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
**  Assumptions: n->getSlope() does not return a negative value (returns neg.
**               only if infinite loop in getSlope()); kb, mb,
**               and nb all >=0.
\***************************************************************************/
double tBedErodePwrLaw::DetachCapacity( tLNode * n, double dt )
{
   double slp = n->getSlope();
   if( slp < 0.0 )
       ReportFatalError("neg. slope in tBedErodePwrLaw::DetachCapacity(tLNode*,double)");
   return( kb*pow( n->getQ(), mb )*pow( slp, nb )*dt );
}

/***************************************************************************\
**  tBedErode::DetachCapacity
**
**  Computes the rate of erosion  = kb Q^mb S^nb
**
**  Input: n -- node at which to compute detachment capacity
** 
**  Returns: the detachment rate
**  Assumptions: n->getSlope() does not return a negative value (returns neg.
**               only if infinite loop in getSlope()); kb, mb,
**               and nb all >=0.
\***************************************************************************/
double tBedErodePwrLaw::DetachCapacity( tLNode * n )
{
   double slp = n->getSlope();
   if( slp < 0.0 )
       ReportFatalError("neg. slope in tBedErodePwrLaw::DetachCapacity(tLNode*)");
   double erorate =  kb*pow( n->getQ(), mb )*pow( slp, nb );
   n->setDrDt( -erorate );
   return erorate;
}

/***************************************************************************\
**  tBedErode::DetachCapacity
**
**  Computes the rate of erosion  = e*kb Q^mb S^nb
**
**  Input: n -- node at which to compute detachment capacity
**         i -- layer which you are computing detachment of
** 
**  Returns: the detachment rate
**  Assumptions: n->getSlope() does not return a negative value (returns neg.
**               only if infinite loop in getSlope()); kb, mb,
**               and nb all >=0.
\***************************************************************************/
double tBedErodePwrLaw::DetachCapacity( tLNode * n, int i )
{
   int g;
   double slp = n->getSlope();
   if( slp < 0.0 )
       ReportFatalError("neg. slope in tBedErodePwrLaw::DetachCapacity(tLNode*)");
   double erorate =n->getLayerErody(i)*kb*pow( n->getQ(), mb )*pow( slp, nb );
   n->setDrDt( -erorate );
//    if(n->getID() == 11 ){
//       cout<<"node 11 is in detach capacity"<<endl;
//       cout<<"layer is "<<i<<endl;
//       cout<<"numg is "<<n->getNumg()<<endl;
//    }
   
   for(g=0; g<n->getNumg(); g++){
       n->setLayerDrDt(i,g,-erorate*n->getLayerDgrade(i,g)/n->getLayerDepth(i));
//        if(n->getID() == 11 ){
//           cout<<"g is "<<g<<endl;
//           cout<<"x is "<<n->getX()<<" y is "<<n->getY()<<endl;
//           cout<<"layer "<<i<<" drdt is set at "<<n->getLayerDrDt(i,g)<<endl;
//        }
   }
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
**  Assumptions: getSlope() returns a value >=0, edge length>0.
**
\***************************************************************************/
double tBedErodePwrLaw::SetTimeStep( tLNode * n )
{
   double slp = n->getSlope();
   if( slp < 0.0 )
       ReportFatalError("neg. slope in tBedErodePwrLaw::setTimeStep(tLNode*)");
   assert( n->getQ()>=0 );
   double eroterm = kb * pow( n->getQ(), mb ) * pow( slp, nb-1.0 );
   if( eroterm==0 ) return 100000;
   return( 0.2*n->getFlowEdg()->getLength() / eroterm );

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
   double slp = node->getSlope();
   if( slp < 0.0 )
       ReportFatalError("neg. slope in tBedErodePwrLaw::TransCapacity(tLNode*)");
   double cap = 0;
   if( !node->getFloodStatus() )
       cap = kf * pow( node->getQ(), mf ) * pow( slp, nf );
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

   cout << "tSedTransWilcock(infile)\n" << endl;
   add='1';
   grade.setSize(2);
   for(i=0; i<=1; i++){
      strcpy( name, "GRAINDIAM");
      strcat( name, &add );
      help = infile.ReadItem( help, name);
      add++;
      grade[i] = help;
   }

   taudim= RHO*GRAV;
   refs = (RHOSED-RHO)*9.81*grade[0];
   refg = (RHOSED-RHO)*9.81*grade[1];
   lowtaucs = 0.8*(grade[1]/grade[0])*0.040*refs*0.8531;
   lowtaucg = 0.04*refg*0.8531;
   hightaucs = 0.04*refs*0.8531;
   hightaucg = 0.01*refg*0.8531;
   sands = (lowtaucs-hightaucs)/(-0.3);
   // slope = m = del y / del x
   sandb = lowtaucs-(sands*0.1);
   // intercept = y - mx
   gravs = (lowtaucg-hightaucg)/(-0.3);
   gravb = lowtaucg-(gravs*0.1);
}

/*********************************************************************\
  tSedTransWilcock::TransCapacity

  This function uses the sediment transport model developed by
  P. Wilcock to calculate sediment transport rates of the sand and
  gravel fraction individually.  This function should only be used with
  two grain sizes and it is assumed that the grain size one is in the
  sand range and grain size 2 is in the gravel range.  The sediment
  transport rate of both grain sizes is calculated, and the sum of
  these two rates is returned. (rate here is in m^3/yr)
/***********************************************************************/

double tSedTransWilcock::TransCapacity( tLNode *nd )
{
   double tau;
   double taucrit;
   double persand=nd->getLayerDgrade(0,0)/(nd->getLayerDepth(0));
   double timeadjust=31536000.00; /* number of seconds in a year */
   double factor=nd->getLayerDepth(0)/nd->getMaxregdep();

   //if(nd->getID() == 93){
   // cout<<"WIL you're there with 93"<<endl;
   //}
   
   
   if( nd->getSlope() < 0 ){
      nd->setQs(0, 0);
      nd->setQs(1, 0);
      nd->setQs(0);
      return 0.0;
   }

   //if(nd->getID() == 93){
   // cout<<"WIL node 93 has "<<nd->getNumLayer()<<" layers"<< endl;
   // cout<<"WIL texture of surface is "<<nd->getLayerSed(0)<<endl;
   //}
   
   
   // units of Q are m^3/sec
   tau = taudim*pow(0.03, 0.6)*pow(nd->getQ(),0.3)*pow( nd->getSlope(), 0.7);
   //cout << "hydrrough is " << nd->getChanRough() << endl;
   //cout << "q is " << nd->getQ() << endl;
   //cout << "slope is " << nd->getSlope() << endl;
   //cout << "taudim is " << taudim << endl;

   //Calculate Sand transport rates first
   
   if(persand<.10)
       taucrit=lowtaucs;
   else if(persand<=.40)
       taucrit=((sands*persand)+sandb);
   else
       taucrit=hightaucs;

   //cout<<"nic value of tau is "<<tau<<" value of taucsand is "<<taucrit<<endl;
   
   if(tau>taucrit){
       nd->setQs(0, ((0.058/RHOSED)*nd->getLayerErody(0)*factor*pow(nd->getQ(),0.5)*timeadjust*persand*pow(tau,1.5)*pow((1-sqrt(taucrit/tau)),4.5) ));
       //cout << "nic sand transport rate is " << nd->getQs(0) << endl;
   }
   else 
       nd->setQs( 0, 0.0 ) ;

   //Now calculate Gravel transport rates

   if(persand<.10)
       taucrit=lowtaucg;
   else if(persand<=.40)
       taucrit=((gravs*persand)+gravb);
   else
       taucrit=hightaucg;

   //cout<<"nic value of tau is "<<tau<<" value of taucgrav is "<<taucrit<<endl;

   if(tau>taucrit){
       nd->setQs(1, (0.058*timeadjust*nd->getLayerErody(0)*factor*pow(nd->getQ(),0.5)/(RHOSED))*
                 (1-persand)*pow(tau,1.5)*pow((1-(taucrit/tau)),4.5));
       //  cout << "nic nic nic gravel transport is happening" << endl;
   }
   else
       nd->setQs(1,0.0);

   nd->setQs(nd->getQs(0)+nd->getQs(1));
   return nd->getQs();
       
}

/*********************************************************************\
  tSedTransWilcock::DetachCapacity

  This function uses the sediment transport model developed by
  P. Wilcock to calculate sediment detachment rates of the sand and
  gravel fraction individually.  This function should only be used with
  two grain sizes and it is assumed that grain size one is in the
  sand range and grain size 2 is in the gravel range.  The detachment
  rate of both grain sizes is calculated, and the sum of
  these two rates is returned. (rate here is in m/yr)

  NOTE : nic you are not very consistent here.  Factor, which is a
  part of transcapacity is NOT in this because it is put into
  the erosion algorithm.  Also the units on the return rates are
  different.  Maybe you should change this?

  takes : nd - node at which you are calculating rate
          i - layer which you are detaching from
/***********************************************************************/

double tSedTransWilcock::DetachCapacity( tLNode *nd, int i )
{
   double tau;
   double taucrit;
   double persand=nd->getLayerDgrade(i,0)/(nd->getLayerDepth(i));
   double timeadjust=31536000.00; /* number of seconds in a year */

   //if(nd->getID() == 93){
   //cout<<"You're in wilcock detach capacity"<<endl;
   //}
   
   
   if( nd->getSlope() < 0 ){
      nd->setDrDt(0.0);
      nd->setLayerDrDt(i,0,0.0);
      nd->setLayerDrDt(i,1,0.0);
      return 0.0;
   }

   //if(nd->getID() == 93){
   //cout<<"Wil Number of layers is "<<nd->getNumLayer()<<" layers"<< endl;
   //cout<<"WIL texture of layer is "<<nd->getLayerSed(i)<<endl;
   //}
   
   
   // units of Q are m^3/sec
   tau = taudim*pow(0.03, 0.6)*pow(nd->getQ(),0.3)*pow( nd->getSlope(), 0.7);
   //NIC you need to put channel roughness function and width function
   //into erosion or maybe tlnode
   //cout << "hydrrough is " << nd->getChanRough() << endl;
   //cout << "q is " << nd->getQ() << endl;
   //cout << "slope is " << nd->getSlope() << endl;
   //cout << "taudim is " << taudim << endl;
   //cout <<"area is "<<nd->getVArea()<<endl;

   //Calculate Sand transport rates first
   
   if(persand<.10)
       taucrit=lowtaucs;
   else if(persand<=.40)
       taucrit=((sands*persand)+sandb);
   else
       taucrit=hightaucs;

   //cout<<"nic value of tau is "<<tau<<" value of taucsand is "<<taucrit<<endl;
   
   if(tau>taucrit){
       nd->setLayerDrDt(i,0, ((0.058/(RHOSED*nd->getVArea()))*nd->getLayerErody(i)*pow(nd->getQ(),0.5)*timeadjust*persand*pow(tau,1.5)*pow((1-sqrt(taucrit/tau)),4.5) ));
       //cout << "nic sand detach rate is " << nd->getLayerDrDt(i,0) << endl;
   }
   else 
       nd->setLayerDrDt(i, 0, 0.0 ) ;

   //Now calculate Gravel transport rates

   if(persand<.10)
       taucrit=lowtaucg;
   else if(persand<=.40)
       taucrit=((gravs*persand)+gravb);
   else
       taucrit=hightaucg;

   //cout<<"nic value of tau is "<<tau<<" value of taucgrav is "<<taucrit<<endl;

   if(tau>taucrit){
       nd->setLayerDrDt(i,1, (0.058*timeadjust*nd->getLayerErody(i)*pow(nd->getQ(),0.5)/(RHOSED*nd->getVArea()))*(1-persand)*pow(tau,1.5)*pow((1-(taucrit/tau)),4.5));
       //cout << "nic gravel detach rate is " << nd->getLayerDrDt(i,1)<<endl;
   }
   else
       nd->setLayerDrDt(i,1,0.0);

   nd->setDrDt(nd->getLayerDrDt(i,0)+nd->getLayerDrDt(i,1));
   return nd->getDrDt();
       
}


/***************************************************************************\
**  tErosion functions
\***************************************************************************/
tErosion::tErosion( tGrid<tLNode> *gptr, tInputFile &infile )
        : bedErode( infile ), sedTrans( infile ), totalTrans( infile )
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
   int nActNodes = gridPtr->getNodeList()->getActiveSize();
   tGridListIter<tLNode> ni( gridPtr->getNodeList() );
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
         dn = cn->getDownstrmNbr();
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
   int nActNodes = gridPtr->getNodeList()->getActiveSize();
   tGridListIter<tLNode> ni( gridPtr->getNodeList() );
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
         /*slp = cn->getSlope();
         dslpdt = cn->getDSlopeDt();
         if( slp > 0.0 )
         {
            if( dslpdt < 0.0 )
            {
               dt = slp / (-dslpdt - UPtr->getRate() ) * frac;
               if( dt > 0 && dt < dtmax ) dtmax = dt;
            }
         }*/
         
         dn = cn->getDownstrmNbr();
         if( dn->getBoundaryFlag() == kNonBoundary )
             ratediff = dn->getDzDt() - cn->getDzDt();
         else
             ratediff = dn->getDzDt() - cn->getDzDt() - UPtr->getRate();
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


#define kSmallTimeStep 1e-8
void tErosion::StreamErode( double dtg, tStreamNet *strmNet )
{
   double dt,
       dtmax;         // time increment: initialize to arbitrary large val
   double frac = 0.3; // fraction of time to zero slope
   int i;
   tLNode * cn, *dn;
   int nActNodes = gridPtr->getNodeList()->getActiveSize();
   tGridListIter<tLNode> ni( gridPtr->getNodeList() );
   double ratediff,  // Difference in ero/dep rate btwn node & its downstrm nbr
       cap,          // Transport capacity
       pedr,         // Potential erosion/deposition rate
       dcap,         // Bedrock detachment capacity
       dz,           // Depth of deposition/erosion (erosion = negative)
       dzr;          // Potential depth of bedrock erosion
   int smallflag=0, smallcount=0;

   cout << "tErosion::StreamErode\n";

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
            // get detachment capacity (this also sets node's drdt)
            dcap = -bedErode.DetachCapacity( cn );
            if( dcap > pedr )
                pedr = dcap;
         }
         // set the erosion (deposition) rate and send the corresponding
         // sediment influx downstream
         cn->setDzDt( pedr );
         cn->getDownstrmNbr()->AddQsin( cn->getQsin() - pedr*cn->getVArea() );
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
         dn = cn->getDownstrmNbr();
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
      if( dtmax <= 0.01 && smallflag==0 )
      {
         smallflag=1;
         cout << "SMALL STEP: " << dtmax << endl;
      }
      if( smallflag==1 )
      {
         smallcount++;
         if( smallcount==100 )
         {
            cout << "TIME REMAINING: " << dtg << endl;
            smallcount=0;
         }
      }
      
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
         dn = cn->getDownstrmNbr();
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

   cout << "Leaving StreamErode()\n";
}

/***************************************************************************\
  tErosion::StreamErodeMulti

  Eventually there should only be one streamerode (not one for multi
  and one for single size).  This work is still in progress so both
  remain for the time being.  This function calculates the erosion and
  deposition rate of each grain size assuming the system is transport
  limited.  If there is less than a critical amount of sediment on the
  bed, the bedrock is eroded and added to the load.

/***************************************************************************/

void tErosion::StreamErodeMulti( double dtg, tStreamNet *strmNet, double time )
{
   double dt,
       dtmax;         // time increment: initialize to arbitrary large val
   double frac = 0.3; //fraction of time to zero slope
   double timegb; //time gone by - for layering time purposes
   int i,n;
   tLNode * cn, *dn;
   int nActNodes = gridPtr->getNodeList()->getActiveSize();
   tGridListIter<tLNode> ni( gridPtr->getNodeList() );
   double ratediff,  // Difference in ero/dep rate btwn node & its downstrm nbr
       cap,
       pedr,
       dcap,
       dzt, // Total amount of sediment erosion
       dzrt; // Total amount of bedrock erosion
   
   cn = ni.FirstP();
   
   tArray <double> dz( cn->getNumg() );
   // Depth of deposition/erosion  of each size (erosion = negative)
   tArray <double> dzr( cn->getNumg() );
   // Potential depth of bedrock erosion of each grain size
   tArray <double> retbr( cn->getNumg() ); //br amt actually ero'd/dep'd
   tArray <double> retsed( cn->getNumg()  ); //sed amt actually ero'd/dep'd

   // Sort so that we always work in upstream to downstream order
   strmNet->SortNodesByNetOrder();

   // Compute erosion and/or deposition until all of the elapsed time (dtg)
   // is used up
   timegb=0;
   do
   {
      //cout<<"AT BEGINNING!"<<endl;
      // Zero out sed influx of all sizes
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() ){
         //if(cn->getID()==93)
         //  cout<<"93 is active"<<endl;
         cn->setQsin(0.0); //totals are for ts calculation
         cn->setQs(0.0);
         for( i=0; i<cn->getNumg(); i++ ){
             cn->setQsin( i, 0.0 );
             cn->setQs(i,0.0);
         }
         //if(cn->getID()==93)
         //    cout<<"93 is active"<<endl;
      }

      // Compute erosion rates: when this block is done, the transport rates
      // (qsm), influx (qsinm), and deposition/erosion rate (dzdt) values are
      // set for each active node.
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      {
         // Transport capacity and potential erosion/deposition rate
         // (this also sets the node's Qs value)
         // remember that qs and qsin store totals of qsm and qsinm, resp.
         if(cn->getLayerSed(0)>0)
         {
            // Sediment in first layer so transport cap
            cap = sedTrans.TransCapacity( cn );
            pedr = (cn->getQsin() - cap ) / cn->getVArea();
         
            //sediment input:
            if( cn == strmNet->getInletNodePtr() ){
               pedr += strmNet->getInSedLoad() / cn->getVArea();
            }
            // If we're on bedrock, adjust accordingly
            if( cn->getLayerSed(1)==0 && pedr<0 &&  fabs(cn->getLayerDepth(0)-cn->getMaxregdep())>0.001 )
            {
               // Bedrock right below sediment so erode this
               // Get detachment capacity (this also sets node's drdt)
               // Limit detachment because there is also sediment present
               dcap = -bedErode.DetachCapacity(cn)*(1-(cn->getLayerDepth(0)/cn->getMaxregdep()));
               pedr += dcap;
               // if(fabs(dcap)>1){
//                   cout << "huge bedrock erosion of "<< dcap << " at node " << cn->getID()<<endl;
//                }
               
            }
            
         }
         else
         {
            // Top layer is bedrock, so just detach and go...
            //sediment input:
            // ignoring sediment input for now
            // nic - NEED TO CHANGE THIS
            //if( cn == strmNet->getInletNodePtr() ){
            //   pedr += strmNet->getInSedLoad() / cn->getVArea();
            //}
            //if( pedr<0 )
            //{
               // Get detachment capacity (this also sets node's drdt)
             dcap = -bedErode.DetachCapacity( cn );
             // if(fabs(dcap)>1){
//                 cout << "huge bedrock erosion of "<< dcap <<" at node " << cn->getID()<< endl;
//              }
                          
             //if( dcap > pedr )
                   pedr = dcap;
                   //}
         }

         // if(cn->getX()==2.75 && cn->getY()==1){
//             cout<<cn->getNumLayer()<<" layers" <<endl;
//             n=0;
//             while(n<cn->getNumLayer()){
//                cout << "layer " << n+1 << endl;
//                cout << "layer creation time is " << cn->getLayerCtime(n) << endl;
//                cout << "layer recent time is " << cn->getLayerRtime(n) << endl;
//                cout << "layer depth is " << cn->getLayerDepth(n) << endl;
//                cout << "layer erodibility is " << cn->getLayerErody(n) << endl;
//                cout << "is layer sediment? " << cn->getLayerSed(n) << endl;
//                cout << "dgrade 1 is " << cn->getLayerDgrade(n,0) << endl;
//                cout << "dgrade 2 is " << cn->getLayerDgrade(n,1) << endl;
//                n++;  
//             }
//             cout<<"qs0 is "<<cn->getQs(0)<<" qs1 is "<<cn->getQs(1)<<endl;
//             cout<<"texture of surface is "<<cn->getLayerSed(0)<<endl;
//          }
         
         // Set the erosion (deposition) rate and send the corresponding
         // sediment influx downstream
         cn->setDzDt( pedr );
         cn->getDownstrmNbr()->AddQsin( cn->getQsin() - pedr*cn->getVArea() );
         // only doing totals here
         // sediment transport rates for each grn size have been calculated
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
         dn = cn->getDownstrmNbr();
         ratediff = dn->getDzDt() - cn->getDzDt(); // Are the pts converging?
         if( ratediff > 0 && cn->getZ() > dn->getZ() )  // if yes, get time
         {                                              //  to zero slope
            dt = ( cn->getZ() - dn->getZ() ) / ratediff;
            if( dt < dtmax ) dtmax = dt;
            if( dt < 1e-6 )
            {
               
               cout << "Very small dt " << dt << " at:\n" << endl;
               //cout << "rate dif is " << ratediff << endl;
               //cout << "elev dif is " << cn->getZ() - dn->getZ() << endl;
               //cout << "dzdt upstream is " << cn->getDzDt() << endl;
               //cout << "dzdt downstream is " << dn->getDzDt() << endl;
               // cn->TellAll();
               //dn->TellAll();
               //cout << "arbitrarily set dt to 0.0015" << endl;
               dtmax=0.005;
            }
         }
      }
      dtmax *= frac;  // Take a fraction of time-to-flattening

      // Zero out sed influx again, because depending on bedrock-alluvial
      // interaction it may be modified; if inlet, give it strmNet->inlet.inSedLoad
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() ){
         cn->setQsin(0.0);
         for( i=0; i<cn->getNumg(); i++ )
             cn->setQsin( i, 0.0 );
      }
      

      //sediment input:
      // nic- ignoring sediment input for now CHANGE THIS LATER
      //if(strmNet->getInletNodePtrNC() != NULL){
      //   for( i=0; i<=cn->getNumg(); i++ )
      //      strmNet->getInletNodePtrNC()->setQsin(i,strmNet->getInSedLoad(i) );
      //}
      

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

      timegb += dtmax;
      // Now do erosion/deposition by integrating rates over dtmax
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      {
         // Depth of potential erosion due to excess transport capacity
         dzt=0;
         for( i=0; i<cn->getNumg(); i++ ){
            dz[i] = ( (cn->getQsin(i)-cn->getQs(i)) / cn->getVArea() ) * dtmax;
            dzt += dz[i];
            retbr[i]=0;
            retsed[i]=0;
         }
        //  if(cn->getX()==2.75 && cn->getY()==1){
//             cout<<"qsin0 qs0 "<<cn->getQsin(0)<<" "<<cn->getQs(0)<<endl;
//             cout<<"qsin1 qs1 "<<cn->getQsin(1)<<" "<<cn->getQs(1)<<endl;
//             cout<<"area "<<cn->getVArea()<<" time "<<dtmax<<endl;
//          }
         
         
         // If we're on bedrock, scour the bedrock
         dzrt = 0;
         if(cn->getLayerSed(0)<1){
            // Bedrock at surface
            for( i=0; i<cn->getNumg(); i++ ){
               dzr[i] = cn->getDrDt()*cn->getLayerDgrade(0,i)/cn->getLayerDepth(0)*dtmax;
               dzrt += dzr[i];
               // if(cn->getX()==2.75 && cn->getY()==1){
//                   cout<<"no sed drdt is "<<dzrt<<" of node "<<cn->getID()<<endl;
//                   cout<<"should erode some bedrock"<<endl;
//                }
            }
            
            if(dzrt<0)
                retbr=cn->EroDep(0,dzr,time+timegb);
            //if(cn->getID()==93){
            // cout<<"done with erosion nic"<< endl;
            //}
            
         }
         else if( fabs(cn->getLayerDepth(0)-cn->getMaxregdep())>0.001 && dzt<0.0 && cn->getLayerSed(1)<1)
         {
            // Bedrock not at surface, but not enough sediment either
            // nic this should work if only regolith and bedrock
            // and also IF layering is working correctly
            for( i=0; i<cn->getNumg(); i++ ){
               dzr[i] = cn->getDrDt()*cn->getLayerDgrade(1,i)/cn->getLayerDepth(1)*dtmax*((cn->getMaxregdep()-cn->getLayerDepth(0))/cn->getMaxregdep());
               dzrt += dzr[i];
               // potential bedrock erosion depth
               // If excess capacity depth-equivalent is greater than the depth
               // of sediment available on bed plus eroded from bedrock, limit
               // depth of erosion to alluvium plus bedrock erosion depth
               //if( -dz[i] > -dzr[i]+cn->getAlluvThickness( i ) )
               //  dz[i] = dzr[i]+cn->getAlluvThickness( i );
               // nic - don't think you need this
            }
            // if(cn->getX()==2.75 && cn->getY()==1){
//                cout<<"sed on top drdt is "<<dzrt<<" of node "<<cn->getID()<<endl;
//                cout<<"should erode some bedrock"<<endl;
//             }
            
            if(dzrt<0)
                retbr=cn->EroDep(1,dzr,time+timegb);
         }
                  
         //cout << "** THIS node has dzs " << dzs << " & dzr " << dzr
         //     << " & net change " << dzr+dzs << endl;
         //cn->TellAll();

         // Update alluvium thickness and node elevation
         if(fabs(dzt)>0){
            // if(cn->getX()==2.75 && cn->getY()==1){
//                 cout<<"eroding sediment at specified node, total of "<<dzt<<endl;
//                 cout<<"Total number of layers is "<<cn->getNumLayer()<<endl;
//                 cout<<"dz0 is "<<dz[0]<<" dz1 is "<<dz[1]<<endl;
                
//             }
            
            //cout << "dzt is " << dzt << endl;
            retsed=cn->EroDep(0, dz, time+timegb );
         }
         
         dn = cn->getDownstrmNbr();

         // Send sediment downstream: sediment flux is equal to the flux in
         // plus/minus rate of erosion/deposition times node area
         for(i=0; i<cn->getNumg(); i++){
            dn->AddQsin( i, cn->getQsin(i) - (retbr[i]+retsed[i])*cn->getVArea()/dtmax );
         }
      }
      

      // Update time remaining
      dtg -= dtmax;
   } while( dtg>1e-6 ); // Keep going until we've used up the whole time intrvl
   
}

/***********************************************************************\
  tErosion::DetachErode

  Algorithm for eroding sediment and bedrock where everything is
  considered to be detachment limited.  Material is only detached
  if the stream has the capacity to carry it.

/************************************************************************/

void tErosion::DetachErode(double dtg, tStreamNet *strmNet, double time )
{
   double dt,
       dtmax;         // time increment: initialize to arbitrary large val
   double frac = 0.3; //fraction of time to zero slope
   double timegb=time; //time gone by - for layering time purposes
   int i,j,nodenum;
   tLNode * cn, *dn;
   int nActNodes = gridPtr->getNodeList()->getActiveSize();
   tGridListIter<tLNode> ni( gridPtr->getNodeList() );
   double ratediff,  // Difference in ero/dep rate btwn node & its downstrm nbr
       cap,
       pedr,
       dcap,
       dzdt, 
       dzrt,
       drdt,
       dz,
       depck,
       qs,
       excap,
       factor;
   
   cn = ni.FirstP();
   
   tArray <double> ret( cn->getNumg() ); //amt actually ero'd/dep'd
   tArray <double> erolist( cn->getNumg() );
   tArray <int> erolay( 1 ); // number of layers to erode through at each node

   // Sort so that we always work in upstream to downstream order
   strmNet->SortNodesByNetOrder();
   
   // Compute erosion and/or deposition until all of the elapsed time (dtg)
   // is used up
   timegb=0;
   do
   {   
      // Zero out sed influx of all sizes
      nodenum=0;
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() ){
         //if(cn->getID()==93)
         //  cout<<"93 is active"<<endl;
         nodenum++; // nodenum is used to find number of active nodes
         cn->setQsin(0.0); //totals are for ts calculation
         cn->setQs(0.0);
         for( i=0; i<cn->getNumg(); i++ ){
            cn->setQsin( i, 0.0 );
            cn->setQs(i,0.0);
         }
      }
      erolay.setSize(nodenum);
      //cout << "size of erolay is " << erolay.getSize() << endl;
      //Greg - I realize this array size will be set each time
      //the time do loop is gone through.  Is it better to just
      //loop through all the nodes to set the array size before 
      //entering the do loop? --> NICOLE: why can't you just use nActNodes?
      nodenum=-1;
      // Estimate erosion rates and time-step size
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      {
         //GREG - I don't think it makes sense to use only one
         //detachcapacity function, even if it does consider
         //the eordibility of a layer.  I think this would be
         //going backwards from the idea of hiding and protrusion
         //or grain sizes.  Am I missing something here?  For now
         //the code below works as follows.  A total drdt (meaning
         //sum of all the sizes) is returned by the detachcapacity functions.
         //For now there are two detachcapacity functions, one if the layer
         //is bedrock (sed value in layer is <1) and another if the layer
         //is sediment (sed value in layer is >=1).  The erodibility
         //is still used as a coefficient in the rate calculation
         //so different types of sediment can be more or less erodible,
         //but all sediment will consider critical shear stress to be
         //a function of what is present on the bed.  What I decided to do
         //here to get around the problem of what if the surface layer
         //is very thin - do we still do detachment based on that layer alone?
         //I do a weighted average of the detachrate of all the surface layers
         //which make up some perscribed depth (maxregdep for now).
         //Note that the bottom layer is not necessarily weighted by
         //it's actual depth, only the depth which is needed to fufill
         //the perscribed depth (maxregdep)
         //The total drdt is just the sum of all of these.
         //There are alot of if statements below which I know is a coding
         //no no.  I don't know what else to do.  Any suggestions?
         //Maybe this whole method is kind of hokie, feel free to change it.
         //Also, with the multiple layer thing I couldn't figure out a way
         //to do the limiting of the time step if a whole layer is removed.
         //It would become an issue if a layer was already eroded from
         //and then erosion from a lower layer was limited, but the upper
         //layer had already been eroded from according to the original time
         //step.  (see the eroalgo.txt to see what I am talking about)
         //I just scrapped the reduction of the time step if a whole layer
         //was eroded.  Maybe this is a problem.  Again, feel free to
         //take issue with it.
         // NICOLE: weighted avg makes sense to me. I think it makes sense
         // to use channel depth as the averaging (and active layer) depth.
         // Not sure I understand what you mean by limiting time step if
         // whole layer is eroded. As far as separate detachcap fns go, I
         // guess the issue is how to define tau crit. Seems to me tau crit
         // for detachment is different from transport tau crit --- for
         // completely disaggregated particles, tauc "det" would be zero even
         // if entrainment/xport tauc was large! what do you think?
         nodenum++;
         depck=0;
         i=0;
         drdt=0;
         while((cn->getMaxregdep()-depck)>0.0001)
         {
            // Total detachment Capacity is corrected for the
            // amount that the layer is exposed - layer depths
            // are used as a surrogate for exposure
            // detachcapacity returns positive values
            //cout<<"depck at begining is "<<depck<<" index i is "<<i<<endl;
            //cout<<"layer depth is "<<cn->getLayerDepth(i)<<endl;
            //cout<<"is layer sed? "<<cn->getLayerSed(i)<<endl;
            if((depck+cn->getLayerDepth(i))<=cn->getMaxregdep()){
               if(cn->getLayerSed(i)>0) //sediment 
                   drdt-=sedTrans.DetachCapacity(cn,i)*cn->getVArea()*cn->getLayerDepth(i)/cn->getMaxregdep();
               else //bedrock
                   drdt-=bedErode.DetachCapacity(cn,i)*cn->getLayerDepth(i)/cn->getMaxregdep();
            }
            else{
               if(cn->getLayerSed(i)<1) //bedrock
               {
                  drdt-=bedErode.DetachCapacity(cn,i)*(1-(depck/cn->getMaxregdep()));
                  // if(cn->getID() == 11)
//                       cout<<"value from detach capacity is "<<bedErode.DetachCapacity(cn,i)<<endl;
               }
               
               else //sediment
                   drdt-=sedTrans.DetachCapacity(cn,i)*cn->getVArea()*(1-(depck/cn->getMaxregdep()));
            }
            depck+=cn->getLayerDepth(i); //need to keep this here for drdt calc
            i++;
            //cout<<"at end of loop "<<i-1<<" drdt is "<<drdt<<endl;
            //cout<<"depck at end is "<<depck<<" index i is "<<i<<endl;
         }
         erolay[nodenum]=i;
         //cout << "number of layers to erode at node "<<nodenum<<" is "<<erolay[nodenum]<<endl;
         //cout<<"final drdt is "<<drdt<<endl;
         
         //although drdt is set individually when each detachcap
         //function is called, it needs to be reset to the total
         //thought about changing to adding but didn't want to
         //muck up the other functions and thought things should be consistent
         cn->setDrDt(drdt);
         cn->setDzDt(drdt);
         //Greg, Am I using the right function here?
         qs=totalTrans.TransCapacity( cn );
         excap=(qs - cn->getQsin())/cn->getVArea();
         //excap negative = deposition; positive = erosion
         //Note that signs are opposite to what one
         //might expect.  This works out for Qsin addition.
         //Limit erosion to capacity of flow or deposition
         if( -drdt > excap ){
             cn->setDzDt(excap);
             //cout<<"at node "<<nodenum<<" erosion is limited to capacity of "<<excap<<endl;
         }
         cn->getDownstrmNbr()->AddQsin(cn->getDzDt() * cn->getVArea());
      }//ends for( cn = ni.FirstP...
      
      //Find local time-step based on dzdt
      //GREG Your outline didn't have details for
      //finding dtmax so I just used same method
      //from streamerode.  If you had something else in mind
      //pls feel free to change this.
      dtmax = dtg/frac;
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      {
         //Not for time step calculations, just utilizing loop
         cn->setQsin(0.0);
         // NIC, part below may be redundant, Check this
         for( i=0; i<cn->getNumg(); i++ )
             cn->setQsin( i, 0.0 );
         
         dn = cn->getDownstrmNbr();
         ratediff = dn->getDzDt() - cn->getDzDt(); //Are the pts converging?
         if( ratediff > 0 && cn->getZ() > dn->getZ() )  // if yes, get time
         {                                              //  to zero slope
            dt = ( cn->getZ() - dn->getZ() ) / ratediff;
            if( dt < dtmax ) dtmax = dt;
            if( dt < 1e-6 )
            {
               
               cout << "Very small dt " << dt <<  endl;
               //cout << "rate dif is " << ratediff << endl;
               //cout << "elev dif is " << cn->getZ() - dn->getZ() << endl;
               //cout << "dzdt upstream is " << cn->getDzDt() << endl;
               //cout << "dzdt downstream is " << dn->getDzDt() << endl;
               // cn->TellAll();
               //dn->TellAll();
               //cout << "arbitrarily set dt to 0.0015" << endl;
               dtmax=0.005; //GREG I added this just because things
               // were taking forever.  I kept it for now just for
               // testing stuff.  Maybe we should discuss this.
            }
         }
      }// End for( cn = ni.FirstP()..
      dtmax *= frac;  // Take a fraction of time-to-flattening
      timegb+=dtmax;
      //cout<<"dtmax is set at "<<dtmax<<endl;
      //At this point: we have drdt and qs for each node, plus dtmax
      
      // Do erosion/deposition
      nodenum=-1;
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      {
         //cout<<"entered erosion/depo loop"<<endl;
         
         nodenum++;
         excap=(cn->getQs() - cn->getQsin())/cn->getVArea();
         //cout<<"passed excap calculation, excap = "<<excap<<endl;
         //again, excap pos if eroding, neg if depositing
         //nic here is where drdt comes in again
         if( -cn->getDrDt() < excap ){
            dz = cn->getDrDt()*dtmax; // detach-lim
            //cout<<"detach-lim - dz is set at "<<dz<<endl;
         }
         else if(excap != 0){
            factor=-cn->getDrDt()/excap; //only use if eroding so negative is OK
            dz = -excap*dtmax; // trans-lim
            //cout<<"trans-lim - dz is set at "<<dz<<endl;
         }
         else
             dz=0;
         
         //cout<<"dz for real is "<<dz<<endl;
         if( dz<0 ) //erosion
         {
            //NIC&GREG Issue that you need to make sure that you erode
            //the correct amounts from each layer.  Thought about doing
            //a weighted average by erodibility and layer depth.
            //I have decided to store dzdt info (calculated
            //when detachcapacity is called) in each layer
            //because there should be a dzdt for each grain size too.
            //This part of the code strays alot from eroalgo.txt because
            //I couldn't figure out how to do the time thing with
            //the different layers.  You might be able to erode the
            //whole time-steps worth from some layers but not from
            //others.  So I basically ignored the whole dt thing,
            //hoping that the timestep would get split up enough
            //by the time step predictor.
            i=0;
            while(i<erolay[nodenum])
            {
               //cout<<"nodenum is "<<nodenum<<" number of layers is "<<erolay[nodenum]<<endl;
               //cout<<"nodeid is "<<cn->getID()<<endl;
               //cout<<"x is "<<cn->getX()<<" y is "<<cn->getY()<<endl;
               //cout<<"numg is "<<cn->getNumg()<<endl;
               for(j=0;j<cn->getNumg();j++){
                  //cout<<"j is "<<j<<endl;
                  //cout<<"drdt is "<<cn->getLayerDrDt(i,j)<<endl;
                  erolist[j]=cn->getLayerDrDt(i,j)*factor;
               }
               ret=cn->EroDep(i,erolist,timegb);
               for(j=0;j<cn->getNumg();j++)
                   cn->getDownstrmNbr()->AddQsin(j,-ret[j]*cn->getVArea()/dtmax);
               i++;
            }
         }//end if( dz<0 )
         else if(dz>0) //deposition -> need if cause erodep chokes with 0
         {
            //Get texture of stuff to be deposited
            //Don't know yet how Qs[i] will be set
            for(j=0;j<cn->getNumg();j++)
                erolist[j]=(cn->getQsin(j)-cn->getQs(j))*dtmax/cn->getVArea();
            ret=cn->EroDep(0,erolist,timegb);
         }
      } // Ends for( cn = ni.FirstP()...
      // Update time remaining
      dtg -= dtmax;
   } while( dtg>1e-6 );  //Keep going until we've used up the whole time intrvl

}// End erosion algorithm



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
   tGridListIter<tLNode> nodIter( gridPtr->getNodeList() );
   tGridListIter<tEdge> edgIter( gridPtr->getEdgeList() );

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




   
         
               
         
