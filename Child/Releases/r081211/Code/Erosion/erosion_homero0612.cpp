/***************************************************************************/
/**
 **  @file erosion.cpp
 **  @brief Functions for equilibrium checking, sediment transport and
 **          bed erosion(detachment) objects.
 **
 **    tEquilibCheck
 **  Transport objects:
 **    tSedTransPwrLaw
 **    tSedTransPwrLaw2
 **    tSedTransBridgeDom
 **    tSedTransWilcock
 **    tSedTransMineTailings //added 4/00 ng
 **  Detachment objects:
 **    tBedErodePwrLaw
 **    tBedErodePwrLaw2
 **
 **    Created 1/98 gt; add tEqChk 5/98 sl
 **
 **    Modifications:
 **     - detach capacity functions modified to use spatially variable
 **       (ie node-based, not "global") critical shear stress, in
 **       conjunction w/ veg module. GT, Jan 2000
 **     - added DensifyMesh function to add new nodes in areas of high
 **       erosion/deposition flux. (gt, 2/2000)
 **     - modified detachment and power-law transport formulas to use
 **       kt in truly SI units, with a conversion factor from Q in m3/yr
 **       to m3/s added in the constructor -- so kt is READ in SI units,
 **       and multiplied by the unit conversion factor before being
 **       used in the detachment and transport equations. (GT 6/01)
 **     - added functions for new class tSedTransPwrLawMulti (GT 2/02)
 **     - added functions for tBedErodePwrLaw2 (GT 4/02)
 **     - added functions for tSedTransPwrLaw2 (GT 4/02)
 **     - added functions for tSedTransBridgeDom (GT 5/02)
 **     - added a check in tErosion constructor that tests to make sure
 **       user really wants the erosion & transport options for which the
 **       executable is compiled (GT 7/02)
 **     - fixed obscure bug in EroDep in which erosion was being double-
 **       counted when both ero and dep of different sizes was happening
 **       simultaneously. Also added assertions. (GT 8/02)
 **
 **    Known bugs:
 **     - ErodeDetachLim assumes 1 grain size. If multiple grain sizes
 **       are specified in the input file and the detachment limited
 **       option is used, a crash will result when tLNode::EroDep
 **       attempts to access array indices above 1. TODO (GT 3/00)
 **
 **  $Id: erosion.cpp,v 1.139 2005/03/15 17:17:29 childcvs Exp $
 */
/***************************************************************************/

#include <math.h>
#include <assert.h>
# include <iomanip>
//#include <string>
#include "erosion.h"

// Here follows a table for transport and detachment laws, which are
// chosen at compile time using #define switches.
//
// ("X()" trick exposed in:
// The New C: X Macros, Randy Meyers, C/C++ Users Journal,
// 19(5), May 2001)

#define TRANSPORT_LAW_TABLE \
X(PowerLaw1,"Power-law transport formula"), \
X(PowerLaw2,"Power-law transport formula, form 2"), \
X(BridgeDominic,"Bridge-Dominic form of Bagnold bedload formula"), \
X(Wilcock,"Wilcock sand-gravel formula"), \
X(PowerLawMulti,"Multi-size power-law formula"), \
X(MineTailings,"Willgoose/Riley mine tailings formula")

#define TRANSPORT_LAW_TABLE2 \
X(PowerLaw1,tSedTransPwrLaw) \
X(PowerLaw2,tSedTransPwrLaw2) \
X(BridgeDominic,tSedTransBridgeDom) \
X(Wilcock,tSedTransWilcock) \
X(PowerLawMulti,tSedTransPwrLawMulti) \
X(MineTailings,tSedTransMineTailings)

#define X(a,b) a
enum {
  TRANSPORT_LAW_TABLE
};
#undef X

#define X(a,b) b
char const * const TransportLaw[] =
  {
    TRANSPORT_LAW_TABLE
  };
#undef X

const int NUMBER_OF_TRANSPORT_LAWS =
sizeof(TransportLaw)/sizeof(TransportLaw[0]);

#define DETACHMENT_LAW_TABLE \
X(DetachPwrLaw1,"Power law, form 1"), \
X(DetachPwrLaw2,"Power law, form 2")

#define DETACHMENT_LAW_TABLE2 \
X(DetachPwrLaw1,tBedErodePwrLaw) \
X(DetachPwrLaw2,tBedErodePwrLaw2)

#define X(a,b) a
enum {
  DETACHMENT_LAW_TABLE
};
#undef X

#define X(a,b) b
char const * const DetachmentLaw[] =
  {
    DETACHMENT_LAW_TABLE
  };
#undef X

const int NUMBER_OF_DETACHMENT_LAWS =
sizeof(DetachmentLaw)/sizeof(DetachmentLaw[0]);

//#define channel
/***************************************************************************\
 **  FUNCTIONS FOR CLASS tEquilibCheck
 **
 **  Functions to find and report system mass change rate.
 **
 **  5/98 SL
\***************************************************************************/

/***************************************************************************\
 **  Constructors: default (no args) and given mesh and timer as args
\***************************************************************************/
tEquilibCheck::tEquilibCheck()
  :
  meshPtr(0),
  timePtr(0),
  longTime(0.),
  massList(),
  longRate(0.), shortRate(0.)
{}

tEquilibCheck::tEquilibCheck( tMesh< tLNode > &meshRef, tRunTimer &timeRef )
  :
  meshPtr(&meshRef),
  timePtr(&timeRef),
  longTime(0.),
  massList(),
  longRate(0.), shortRate(0.)
{
  FindIterChngRate();
}

tEquilibCheck::tEquilibCheck( tMesh< tLNode > &meshRef, tRunTimer &timeRef,
                              const tInputFile &fileRef )
  :
  meshPtr(&meshRef),
  timePtr(&timeRef),
  longTime(0.),
  massList(),
  longRate(0.), shortRate(0.)
{
  longTime = fileRef.ReadItem( longTime, "EQUITIME" );
  FindIterChngRate();
}

tEquilibCheck::~tEquilibCheck()
{
  meshPtr = 0;
  timePtr = 0;
}


/***************************************************************************\
 **  'get' and 'set' functions for tEquilibCheck:
\***************************************************************************/
double tEquilibCheck::getLongTime() const {return longTime;}

void tEquilibCheck::setLongTime( double val )
{longTime = ( val > 0 ) ? val : 0.0;}

const tMesh< tLNode > *tEquilibCheck::getMeshPtr() const {return meshPtr;}

tMesh< tLNode > *tEquilibCheck::getMeshPtrNC() {return meshPtr;}

void tEquilibCheck::setMeshPtr( tMesh< tLNode > &Ref )
{meshPtr = &Ref;}

void tEquilibCheck::setMeshPtr( tMesh< tLNode > *Ptr )
{meshPtr = Ptr;}

const tRunTimer *tEquilibCheck::getTimePtr() const {return timePtr;}

tRunTimer *tEquilibCheck::getTimePtrNC() {return timePtr;}

void tEquilibCheck::setTimePtr( tRunTimer &Ref )
{timePtr = &Ref;}

void tEquilibCheck::setTimePtr( tRunTimer *Ptr )
{timePtr = Ptr;}

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
  assert( timePtr != 0 && meshPtr != 0 );

  double ma;
  {
    tMesh< tLNode >::nodeListIter_t nI( meshPtr->getNodeList() );
    tLNode *cn;
    double mass = 0.0;
    double area = 0.0;
    for( cn = nI.FirstP(); nI.IsActive(); cn = nI.NextP() )
      {
	mass += cn->getZ() * cn->getVArea();
	area += cn->getVArea();
      }
    ma = mass / area;
  }
  const tArray2< double > tmp( timePtr->getCurrentTime(), ma);

  tListIter< tArray2< double > > mI( massList );
  if( !(massList.isEmpty()) )
    {
      tArray2< double > *last = mI.LastP();
      double dt = (tmp.at(0) - (*last).at(0));
      assert( dt > 0.0 );
      shortRate = (tmp.at(1) - (*last).at(1)) / dt;
    }
  else
    {
      //cout << "tEquilibCheck::FindIterChngRate(), Warning: empty massList\n";
      assert( tmp.at(0) > 0 );
      shortRate = tmp.at(1) / tmp.at(0);
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
  tListIter< tArray2< double > > mI( massList );
  tArray2< double > *last = mI.LastP();
  tArray2< double > *ca, *na;
  double dt, targetTime = (*last).at(0) - longTime;
  if( longTime == 0.0 || mI.FirstP() == mI.LastP() ) longRate = shortRate;
  else
    {
      ca = mI.FirstP();
      na = mI.NextP();
      while( (*na).at(0) < targetTime && !(mI.AtEnd()) )
	{
	  ca = na;
	  na = mI.NextP();
	}
      dt = (*last).at(0) - (*ca).at(0);
      assert( dt > 0 );
      longRate = ((*last).at(1) - (*ca).at(1)) / dt;
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
 **  FUNCTIONS FOR CLASS tBedErodePwrLaw
\***************************************************************************/

/***************************************************************************\
 **  tBedErodePwrLaw Constructor
 **
 **  Reads in coefficients and exponents for the power law. "mb" is
 **  assumed to be the exponent on specific discharge (L2/T), i.e.,
 **  Dc ~ q^mb. It is then converted to a total discharge exponent,
 **  i.e., Dc ~ Q^mb, by multiplying by (1-ws), where ws is the
 **  at-a-station width-discharge relation. Before converting, we
 **  define the drainage area exponent, ma = mb*(ws-wb). This
 **  essentially describes the difference between at-a-station and
 **  downstream width variations with discharge.
 **
 **  Modifications:
 **   - Erosion rate is now calculated from specific discharge, which is
 **     explicitly computed from Q/W, W = channel width computed in calls
 **     to tStreamNet::FindChanGeom and tStreamNet::FindHydrGeom. As a
 **     result, ma is now obsolete and mb is no longer modified after being
 **     read in from file. mb is still the specific-discharge exponent,
 **     NOT the total-discharge exponent. (GT 2/01)
 **   - kt now read in SI units, with conversion factor between discharge
 **     in m3/yr and shear stress (sic) in SI units factored in here in
 **     this constructor. The derivation is:
 **
 **       Tau (SI) = kt(SI) * Uconv * (Q(m/yr)/W)^mb * S^nb
 **
 **     where Tau is shear stress or power or whatever, depending on
 **     dimensions, and Uconv = SPY^-mb where SPY = # secs in a year
 **     (GT 6/01)
 **
\***************************************************************************/
//constructor: reads and sets the parameters
tBedErodePwrLaw::tBedErodePwrLaw( const tInputFile &infile )
{
  const double secPerYear = SECPERYEAR;  // # secs in one year

  kb = infile.ReadItem( kb, "KB" );
  kt = infile.ReadItem( kt, "KT" );
  ns = infile.ReadItem( ns, "NS" );
  mb = infile.ReadItem( mb, "MB" ); // Specific q exponent
  //wb = infile.ReadItem( wb, "HYDR_WID_EXP_DS" );
  //ws = infile.ReadItem( ws, "HYDR_WID_EXP_STN" );
  //ma = mb*(ws-wb);  // Drainage area exponent
  //mb = mb*(1-ws);   // Convert mb to total-discharge exponent
  nb = infile.ReadItem( nb, "NB" );
  pb = infile.ReadItem( pb, "PB" );
  taucd = infile.ReadItem( taucd, "TAUCD" );

  // Add unit conversion factor for kt -- this is required to convert
  // the quantity (Q/W)^mb from units of years to units of seconds.
  kt *= pow( secPerYear, -mb );
}


/***************************************************************************\
 **  tBedErode::DetachCapacity (1 of 3)
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
 **  Modifications:
 **   - replaced uniform erodibility coefficient kb with erodibility of
 **     topmost layer at node n (GT 8/98)
 **   - added threshold term and third exponent pb (GT 4/99)
 **   - added shear coefficient kt, width term, and moved erodibility coeff
 **     (GT 5/99) (to be ver 2.0.2)
 **   - replaced spatially uniform taucd parameter with node variable
 **     tauc (GT 1/00)
 **   - placed channel width explicitly in the denominator rather than
 **     have it be buried in exponents and coefficients, so as to be
 **     consistent with the transport equations.
 **     (GT 2/01)
\***************************************************************************/
double tBedErodePwrLaw::DetachCapacity( tLNode * n, double dt )
{
  if( n->getFloodStatus() != tLNode::kNotFlooded) return 0.0;
  double slp = n->calcSlope();
  /*This subrutine calculates the backslope of the heacuts-----------*/
  if(slp > 0.84 && atan(slp) > 2* atan((n->getDownstrmNbr())->calcSlope()))
    {
        tPtrList<tLNode > drainingNodes;
        searchDrainingNodes(n, drainingNodes);
        tLNode *ndi;
        tPtrListIter<tLNode> dNIter(drainingNodes);
        //int sizeDN = drainingNodes.getSize();
            double maxarea = 0.0;
        for(ndi=dNIter.FirstP(); !(dNIter.AtEnd()); ndi=dNIter.NextP() )
	  {  
	  if(ndi->getDrArea() > maxarea)
	    {
	    maxarea = ndi->getDrArea();
	    slp = ndi->calcSlope();
	    }
	  }   
    }
  /*-----------------------------------------------------------------*/
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw::DetachCapacity(tLNode*,double)");

#ifdef channel //option to have shear stress assuming channel or overland_flow
  const double tau =  kt*pow( n->getQ() / n->getHydrWidth(), mb ) * pow( slp, nb ) * pow ( ns, 1.5 )/ pow( (ns + n->getVegRough()), 0.9 );
#else
  const double tau = kt*pow( n->getQ() / (n->getFlowEdg())->CalcVEdgLen(), mb ) * pow( slp, nb ) * pow ( ns, 1.5 )/ pow( (ns + n->getVegRough()), 0.9 );
#endif

  n->setTau( tau );
  double tauex = tau - n->getTauCrit();
  tauex = (tauex>0.0) ? tauex : 0.0;
  return( n->getLayerErody(0)*pow(tauex,pb)*dt );
}


/***************************************************************************\
 **  tBedErode::DetachCapacity (2 of 3)
 **
 **  Computes the rate of erosion  = kb Q^mb S^nb
 **
 **  Input: n -- node at which to compute detachment capacity
 **
 **  Returns: the detachment rate
 **  Assumptions: n->calcSlope() does not return a negative value (returns neg.
 **               only if infinite loop in calcSlope()); kb, mb,
 **               and nb all >=0.
 **  Modifications:
 **   - replaced uniform erodibility coefficient kb with erodibility of
 **     topmost layer at node n (GT 8/98)
 **   - added threshold term and third exponent pb (GT 4/99)
 **   - added shear coefficient kt, width term, and moved erodibility coeff
 **     (GT 5/99) (ver 2.0.2)
 **   - replaced spatially uniform taucd parameter with node variable
 **     tauc (GT 1/00)
 **   - placed channel width explicitly in the denominator rather than
 **     have it be buried in exponents and coefficients, so as to be
 **     consistent with the transport equations.
 **     (GT 2/01)
\***************************************************************************/
double tBedErodePwrLaw::DetachCapacity( tLNode * n )
{
  if(0) //DEBUG
    std::cout<<"in detach capacity "<<std::endl;
  assert( n->getQ()>=0.0 );

  if( n->getFloodStatus() != tLNode::kNotFlooded) return 0.0;
  double slp = n->calcSlope();
  /*This subrutine calculates the backslope of the heacuts-----------*/
  /*This subrutine calculates the backslope of the heacuts-----------*/
  if(slp > 0.84 && atan(slp) > 2* atan((n->getDownstrmNbr())->calcSlope()))
    {
        tPtrList<tLNode > drainingNodes;
        searchDrainingNodes(n, drainingNodes);
        tLNode *ndi;
        tPtrListIter<tLNode> dNIter(drainingNodes);
        //int sizeDN = drainingNodes.getSize();
            double maxarea = 0.0;
        for(ndi=dNIter.FirstP(); !(dNIter.AtEnd()); ndi=dNIter.NextP() )
	  {  
	  if(ndi->getDrArea() > maxarea)
	    {
	    maxarea = ndi->getDrArea();
	    slp = ndi->calcSlope();
	    }
	  }   
    }
  /*-----------------------------------------------------------------*/
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw::DetachCapacity(tLNode*)");

#ifdef channel //option to have shear stress assuming channel or overland_flow
  const double tau =  kt*pow( n->getQ() / n->getHydrWidth(), mb ) * pow( slp, nb ) * pow ( ns, 1.5 )/ pow( (ns + n->getVegRough()), 0.9 );
#else
  const double tau = kt*pow( n->getQ() / (n->getFlowEdg())->CalcVEdgLen(), mb ) * pow( slp, nb ) * pow ( ns, 1.5 )/ pow( (ns + n->getVegRough()), 0.9 );
#endif

  if( n->getQ()<0.0 || n->getDrArea()<0.0 ) n->TellAll();
  assert( n->getQ()>=0.0 );
  assert( n->getDrArea()>=0.0 );
  n->setTau( tau );
  double erorate = tau - n->getTauCrit();
  if(0) { //DEBUG
    std::cout << "tau " << tau;
    std::cout << " tauc " << n->getTauCrit() << std::endl;
  }
  erorate = (erorate>0.0) ? erorate : 0.0;
  erorate = n->getLayerErody(0)*pow( erorate, pb );
  n->setDrDt( -erorate );
  return erorate;
}


/***************************************************************************\
 **  tBedErode::DetachCapacity (3 of 3)
 **
 **  Computes the rate of erosion  = e* Q^mb S^nb
 **  Here erodibility of layer is used as the coefficient for detach capacity
 **
 **  TODO: have this just call the other DetachCapacity and multiply by dt!
 **  Also: consolidate w/ 2 of 3 by using default 0 parameter for layer #.
 **
 **  Input: n -- node at which to compute detachment capacity
 **         i -- layer which you are computing detachment of
 **
 **  Returns: the detachment rate
 **  Assumptions: n->calcSlope() does not return a negative value (returns neg.
 **               only if infinite loop in calcSlope()); kb, mb,
 **               and nb all >=0.
 **  Modifications:
 **   - added shear coefficient kt, width term, and moved erodibility coeff
 **     (GT 5/99) (ver 2.0.2)
 **   - replaced spatially uniform taucd parameter with node variable
 **     tauc (GT 1/00)
 **   - placed channel width explicitly in the denominator rather than
 **     have it be buried in exponents and coefficients, so as to be
 **     consistent with the transport equations.
 **     (GT 2/01)
\***************************************************************************/
double tBedErodePwrLaw::DetachCapacity( tLNode * n, int i )
{
  if( n->getFloodStatus() != tLNode::kNotFlooded) return 0.0;
  double slp = n->calcSlope();
  /*This subrutine calculates the backslope of the heacuts-----------*/
  if(slp > 0.84 && atan(slp) > 2* atan((n->getDownstrmNbr())->calcSlope()))
    {
        tPtrList<tLNode > drainingNodes;
        searchDrainingNodes(n, drainingNodes);
        tLNode *ndi;
        tPtrListIter<tLNode> dNIter(drainingNodes);
        //int sizeDN = drainingNodes.getSize();
            double maxarea = 0.0;
        for(ndi=dNIter.FirstP(); !(dNIter.AtEnd()); ndi=dNIter.NextP() )
	  {  
	  if(ndi->getDrArea() > maxarea)
	    {
	    maxarea = ndi->getDrArea();
	    slp = ndi->calcSlope();
	    }
	  }   
    }
  /*-----------------------------------------------------------------*/
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw::DetachCapacity(tLNode*)");

#ifdef channel //option to have shear stress assuming channel or overland_flow
  const double tau =  kt*pow( n->getQ() / n->getHydrWidth(), mb ) * pow( slp, nb ) * pow ( ns, 1.5 )/ pow( (ns + n->getVegRough()), 0.9 );
  n->setTau( tau );

  double erorate = tau - n->getTauCrit();
  if(0) //DEBUG
    std::cout << "erorate: " << erorate << std::endl;
  erorate = (erorate>0.0) ? erorate : 0.0;
  erorate = n->getLayerErody(i)*pow( erorate, pb );
  n->setDrDt( -erorate );

#else
  const double tau = kt*pow( n->getQ() / (n->getFlowEdg())->CalcVEdgLen(), mb ) * pow( slp, nb ) * pow ( ns, 1.5 )/ pow( (ns + n->getVegRough()), 0.9 );

  n->setTau( tau );
  double erorate = tau - n->getTauCrit();
  if(0) //DEBUG
    std::cout << "erorate: " << erorate << std::endl;
  erorate = (erorate>0.0) ? erorate : 0.0;
  erorate = n->getLayerErody(0)*pow( erorate, pb ); //in this case 0 instead of i, because we assume overland flow, not channelized flow (hf, jul06).
  n->setDrDt( -erorate );
  if(0) //DEBUG
  if (n->getID() ==500)
    std::cout << "erody " << i << ": " << n->getLayerErody(i) << ", erody 0 : " << n->getLayerErody(0) << std::endl;

#endif

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
 **  Assumptions: calcSlope() returns a value >=0, edge length>0.
 **
 **  TODO: update this to handle threshold term taucd and pb
\***************************************************************************/
double tBedErodePwrLaw::SetTimeStep( tLNode * n )
{
  const double slp = n->calcSlope();
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw::setTimeStep(tLNode*)");
  assert( n->getQ()>=0 );
  double eroterm = kb * pow( n->getQ(), mb ) * pow( slp, nb-1.0 );
  if( eroterm==0 ) return 100000.;
  return( 0.2*n->getFlowEdg()->getLength() / eroterm );

}

/***************************************************************************\
**  tBedErodePwrLaw::searchDrainingNodes
\***************************************************************************/

void tBedErodePwrLaw::searchDrainingNodes(tLNode *cn,  tPtrList<tLNode> &dNList)
{
  //tPtrList<tEdge> edgeList(cn->getSpokeList());//(or nodes drainingt to cn);
  //tPtrListIter<tEdge> spokeIter(edgeList);
  tSpkIter spokeIter(cn);
  tEdge * edg;
  /*int pointerfound = 0;*/
  for (edg =spokeIter.FirstP(); !(spokeIter.AtEnd()); edg=spokeIter.NextP() )
    {
      tLNode *ndi;
      ndi =(tLNode *)edg->getDestinationPtrNC();  
      if(ndi->getDownstrmNbr() != 0)
	{
	  if((ndi->getDownstrmNbr()->getID()) == cn->getID())
	    {
	      dNList.insertAtFront(ndi);  //the address of a pointer
	    }
	}
    }       
}


/***************************************************************************\
 **  FUNCTIONS FOR CLASS tBedErodePwrLaw2
\***************************************************************************/

/***************************************************************************\
 **  tBedErodePwrLaw2 Constructor
 **
\***************************************************************************/
//constructor: reads and sets the parameters
tBedErodePwrLaw2::tBedErodePwrLaw2( const tInputFile &infile )
{
  const double secPerYear = SECPERYEAR;  // # secs in one year

  kb = infile.ReadItem( kb, "KB" );
  kt = infile.ReadItem( kt, "KT" );
  mb = infile.ReadItem( mb, "MB" ); // Specific q exponent
  //wb = infile.ReadItem( wb, "HYDR_WID_EXP_DS" );
  //ws = infile.ReadItem( ws, "HYDR_WID_EXP_STN" );
  //ma = mb*(ws-wb);  // Drainage area exponent
  //mb = mb*(1-ws);   // Convert mb to total-discharge exponent
  nb = infile.ReadItem( nb, "NB" );
  pb = infile.ReadItem( pb, "PB" );
  taucd = infile.ReadItem( taucd, "TAUCD" );

  // Add unit conversion factor for kt -- this is required to convert
  // the quantity (Q/W)^mb from units of years to units of seconds.
  kt *= pow( secPerYear, -mb );
}


/***************************************************************************\
 **  tBedErodePwrlaw2::DetachCapacity (1 of 3)
 **
 **  Computes the depth of erosion over a time interval dt assuming the
 **  erosion rate = kb ( tau^pb - taucrit^pb )
 **
 **  Input: n -- node at which to compute detachment capacity
 **         dt -- time interval
 **  Returns: the detachment depth
 **  Assumptions: n->calcSlope() does not return a negative value (returns neg.
 **               only if infinite loop in calcSlope()); kb, mb,
 **               and nb all >=0.
 **  Modifications:
 **   - replaced uniform erodibility coefficient kb with erodibility of
 **     topmost layer at node n (GT 8/98)
 **   - added threshold term and third exponent pb (GT 4/99)
 **   - added shear coefficient kt, width term, and moved erodibility coeff
 **     (GT 5/99) (to be ver 2.0.2)
 **   - replaced spatially uniform taucd parameter with node variable
 **     tauc (GT 1/00)
 **   - placed channel width explicitly in the denominator rather than
 **     have it be buried in exponents and coefficients, so as to be
 **     consistent with the transport equations.
 **     (GT 2/01)
\***************************************************************************/
double tBedErodePwrLaw2::DetachCapacity( tLNode * n, double dt )
{
  if( n->getFloodStatus() != tLNode::kNotFlooded) return 0.0;
  const double slp = n->calcSlope();
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw2::DetachCapacity(tLNode*,double)");
  const double tau =  kt*pow( n->getQ() / n->getHydrWidth(), mb ) * pow( slp, nb );
  n->setTau( tau );
  double tauexpb = pow( tau, pb ) - pow( n->getTauCrit(), pb );
  if(0) //DEBUG
    std::cout << "tauexpb: " << tauexpb << std::endl;
  tauexpb = (tauexpb>0.0) ? tauexpb : 0.0;
  return( n->getLayerErody(0)*tauexpb*dt );
}


/***************************************************************************\
 **  tBedErodePwrLaw2::DetachCapacity (2 of 3)
 **
 **  Computes the rate of erosion  = kb ( tau^pb - taucrit^pb )
 **
 **  Input: n -- node at which to compute detachment capacity
 **
 **  Returns: the detachment rate
 **  Assumptions: n->calcSlope() does not return a negative value (returns neg.
 **               only if infinite loop in calcSlope()); kb, mb,
 **               and nb all >=0.
\***************************************************************************/
double tBedErodePwrLaw2::DetachCapacity( tLNode * n )
{
  if(0) //DEBUG
    std::cout<<"in detach capacity "<<std::endl;
  assert( n->getQ()>=0.0 );

  if( n->getFloodStatus() != tLNode::kNotFlooded) return 0.0;
  const double slp = n->calcSlope();
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw2::DetachCapacity(tLNode*)");
  const double tau = kt*pow( n->getQ() / n->getHydrWidth(), mb )*pow( slp, nb );
  if( n->getQ()<0.0 || n->getDrArea()<0.0 ) n->TellAll();
  assert( n->getQ()>=0.0 );
  assert( n->getDrArea()>=0.0 );
  n->setTau( tau );
  double erorate = pow( tau, pb ) - pow( n->getTauCrit(), pb );
  if(0) { //DEBUG
    std::cout << "tau " << tau;
    std::cout << " tauc " << n->getTauCrit() << std::endl;
  }
  erorate = (erorate>0.0) ? erorate : 0.0;
  erorate = n->getLayerErody(0)*erorate;
  n->setDrDt( -erorate );
  return erorate;
}


/***************************************************************************\
 **  tBedErodePwrLaw2::DetachCapacity (3 of 3)
 **
 **  Computes the rate of erosion  = e* Q^mb S^nb
 **  Here erodibility of layer is used as the coefficient for detach capacity
 **
 **  TODO: have this just call the other DetachCapacity and multiply by dt!
 **  Also: consolidate w/ 2 of 3 by using default 0 parameter for layer #.
 **
 **  Input: n -- node at which to compute detachment capacity
 **         i -- layer which you are computing detachment of
 **
 **  Returns: the detachment rate
 **  Assumptions: n->calcSlope() does not return a negative value (returns neg.
 **               only if infinite loop in calcSlope()); kb, mb,
 **               and nb all >=0.
\***************************************************************************/
double tBedErodePwrLaw2::DetachCapacity( tLNode * n, int i )
{
  if( n->getFloodStatus() != tLNode::kNotFlooded) return 0.0;
  const double slp = n->calcSlope();
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw2::DetachCapacity(tLNode*)");
  const double tau = kt*pow( n->getQ() / n->getHydrWidth(), mb )*pow( slp, nb );
  n->setTau( tau );
  double erorate = pow( tau, pb ) - pow( n->getTauCrit(), pb );
  if(0) //DEBUG
    std::cout << "erorate: " << erorate << std::endl;
  erorate = (erorate>0.0) ? erorate : 0.0;
  erorate = n->getLayerErody(i)*erorate;
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
 **  Assumptions: calcSlope() returns a value >=0, edge length>0.
 **
 **  TODO: update this to handle threshold term taucd and pb
\***************************************************************************/
double tBedErodePwrLaw2::SetTimeStep( tLNode * n )
{
  const double slp = n->calcSlope();
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw::setTimeStep(tLNode*)");
  assert( n->getQ()>=0 );
  double eroterm = kb * pow( n->getQ(), mb ) * pow( slp, nb-1.0 );
  if( eroterm==0 ) return 100000.;
  return( 0.2*n->getFlowEdg()->getLength() / eroterm );

}



/***************************************************************************\
 **  FUNCTIONS FOR CLASS tSedTransPwrLaw
\***************************************************************************/

/***************************************************************************\
 **  tSedTransPwrLaw Constructor:
 **
 **    Given a tInputFile as an argument, will read relevant parameters from
 **  the input file.
 **
 **  Modifications:
 **   - kt adjusted to include unit conversion from Q in m3/yr to shear
 **     stress (sic) in SI units. See comment above for tBedErodePwrLaw
 **     constructor. (GT 06/01)
 **
\***************************************************************************/
tSedTransPwrLaw::tSedTransPwrLaw( const tInputFile &infile )
{
  const double secPerYear = SECPERYEAR;  // # secs in one year

  kf = infile.ReadItem( kf, "KF" );
  kt = infile.ReadItem( kt, "KT" );
  ns = infile.ReadItem( ns, "NS" );
  mf = infile.ReadItem( mf, "MF" );
  nf = infile.ReadItem( nf, "NF" );
  pf = infile.ReadItem( pf, "PF" );
  tauc = infile.ReadItem( tauc, "TAUCD" );

  // Add unit conversion factor for kt -- this is required to convert
  // the quantity (Q/W)^mb from units of years to units of seconds.
  kt *= pow( secPerYear, -mf );
}


/***************************************************************************\
 **  tSedTransPwrLaw::TransCapacity
 **
 **  Computes sediment transport capacity using the simple power law
 **  Qs = kf W ( kt (Q/W)^mf S^nf - tauc )^pf
 **
 **  Modifications:
 **   - Now incorporates threshold term; previously was Qs = kf Q^mf S^nf
 **     (GT 5/99) (ver 2.0.2)
\***************************************************************************/
double tSedTransPwrLaw::TransCapacity( tLNode *node )
{
  double slp = node->calcSlope();
  /*This subrutine calculates the backslope of the heacuts-----------*/
  if(slp > 0.84 && atan(slp) > 2* atan((node->getDownstrmNbr())->calcSlope()))
    {
        tPtrList<tLNode > drainingNodes;
        searchDrainingNodes(node, drainingNodes);
        tLNode *ndi;
        tPtrListIter<tLNode> dNIter(drainingNodes);
        //int sizeDN = drainingNodes.getSize();
            double maxarea = 0.0;
        for(ndi=dNIter.FirstP(); !(dNIter.AtEnd()); ndi=dNIter.NextP() )
	  {  
	  if(ndi->getDrArea() > maxarea)
	    {
	    maxarea = ndi->getDrArea();
	    slp = ndi->calcSlope();
	    }
	  }   
    }
  /*-----------------------------------------------------------------*/
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw::TransCapacity(tLNode*)");
  double tau, tauex, cap = 0.;
  if( node->getFloodStatus() == tLNode::kNotFlooded )
    {

#ifdef channel
      tau = kt * pow( node->getQ()/node->getHydrWidth(), mf ) * pow( slp, nf ) * pow ( ns, 1.5 )/ pow( (ns + node->getVegRough()), 0.9 );
#else
      tau = kt*pow( node->getQ()/(node->getFlowEdg())->CalcVEdgLen(), mf ) * pow( slp, nf ) * pow ( ns, 1.5 )/ pow( (ns + node->getVegRough()), 0.9 );
#endif

      node->setTau( tau );
      if(0) //DEBUG
	std::cout << "kt=" << kt << " Q=" << node->getQ() << " W="
	     << node->getHydrWidth() << " S=" << node->calcSlope() << std::endl;
      tauex = tau - tauc;
      tauex = (tauex>0.0) ? tauex : 0.0;

#ifdef channel
      cap = kf * node->getHydrWidth() * pow( tauex, pf );
#else
      cap = kf *  (node->getFlowEdg())->CalcVEdgLen() * pow( tauex, pf );
#endif

    }
  node->setQs( cap );
  return cap;
}


/***************************************************************************\
 **  tSedTransPwrLaw::TransCapacity
 **
 **  Computes sediment transport capacity using the simple power law
 **  Qs = weight kf W ( kt (Q/W)^mf S^nf - tauc )^pf
 **  This is a weighted version which is called from DetachErode.
 **  Weight is a weighting by depth of layer.
 **  Here, qsi is set by proportion in layer and threshold is constant
 **  The value returned should be in units of m^3/yr
 **
 **  Modifications:
 **   - Now incorporates threshold term; previously was Qs = kf Q^mf S^nf
 **     (GT 5/99) (ver 2.0.2)
\***************************************************************************/
double tSedTransPwrLaw::TransCapacity( tLNode *node, int lyr, double weight )
{
  double slp = node->calcSlope();
  /*This subrutine calculates the backslope of the heacuts-----------*/
  if(slp > 0.84 && atan(slp) > 2* atan((node->getDownstrmNbr())->calcSlope()))
    {
        tPtrList<tLNode > drainingNodes;
        searchDrainingNodes(node, drainingNodes);
        tLNode *ndi;
        tPtrListIter<tLNode> dNIter(drainingNodes);
        //int sizeDN = drainingNodes.getSize();
            double maxarea = 0.0;
        for(ndi=dNIter.FirstP(); !(dNIter.AtEnd()); ndi=dNIter.NextP() )
	  {  
	  if(ndi->getDrArea() > maxarea)
	    {
	    maxarea = ndi->getDrArea();
	    slp = ndi->calcSlope();
	    }
	  }   
    }
  /*-----------------------------------------------------------------*/
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tSedTransPwrLaw::TransCapacity(tLNode*)");
  double tau, tauex, cap = 0.;
  if( node->getFloodStatus() == tLNode::kNotFlooded )
    {

#ifdef channel
      tau = kt * pow( node->getQ()/node->getHydrWidth(), mf ) * pow( slp, nf ) * pow ( ns, 1.5 )/ pow( (ns + node->getVegRough()), 0.9 );
#else
      tau = kt*pow( node->getQ()/(node->getFlowEdg())->CalcVEdgLen(), mf ) * pow( slp, nf ) * pow ( ns, 1.5 )/ pow( (ns + node->getVegRough()), 0.9 );
#endif

      node->setTau( tau );
      if(0) //DEBUG
	std::cout << "kt=" << kt << " Q=" << node->getQ() << " W="
	     << node->getHydrWidth() << " S=" << node->calcSlope()
	     << " tau=" << tau << std::endl;
      tauex = tau - tauc;
      tauex = (tauex>0.0) ? tauex : 0.0;

#ifdef channel
      cap = weight * kf * node->getHydrWidth() * pow( tauex, pf );
#else
      cap = weight * kf *  (node->getFlowEdg())->CalcVEdgLen() * pow( tauex, pf );
#endif

    }
  for(size_t i=0; i<node->getNumg(); i++)
    node->addQs(i, cap*node->getLayerDgrade(lyr,i) / node->getLayerDepth(lyr));

  node->setQs( cap );
  return cap;
}

/*************************************************************************\
**   tSedTransPwrLaw::searchDrainingNodes
\**************************************************************************/
void tSedTransPwrLaw::searchDrainingNodes(tLNode *cn,  tPtrList<tLNode> &dNList)
{          
  //tPtrList<tEdge> edgeList(cn->getSpokeList());//(or nodes drainingt to cn);
  //tPtrListIter<tEdge> spokeIter(edgeList);
  tSpkIter spokeIter(cn);
  tEdge * edg;
  /*int pointerfound = 0;*/
  for (edg =spokeIter.FirstP(); !(spokeIter.AtEnd()); edg=spokeIter.NextP() )
    {
      tLNode *ndi;
      ndi =(tLNode *)edg->getDestinationPtrNC();  
      if(ndi->getDownstrmNbr() != 0)
	{
	  if((ndi->getDownstrmNbr()->getID()) == cn->getID())
	    {
	      dNList.insertAtFront(ndi);  //the address of a pointer
	    }
	}
    }       
}  

/***************************************************************************\
 **  FUNCTIONS FOR CLASS tSedTransPwrLaw2
\***************************************************************************/

/***************************************************************************\
 **  tSedTransPwrLaw2 Constructor:
 **
 **    Given a tInputFile as an argument, will read relevant parameters from
 **  the input file.
 **
\***************************************************************************/
tSedTransPwrLaw2::tSedTransPwrLaw2( const tInputFile &infile )
{
  const double secPerYear = SECPERYEAR;  // # secs in one year

  kf = infile.ReadItem( kf, "KF" );
  kt = infile.ReadItem( kt, "KT" );
  mf = infile.ReadItem( mf, "MF" );
  nf = infile.ReadItem( nf, "NF" );
  pf = infile.ReadItem( pf, "PF" );
  tauc = infile.ReadItem( tauc, "TAUCD" );

  // Add unit conversion factor for kt -- this is required to convert
  // the quantity (Q/W)^mb from units of years to units of seconds.
  kt *= pow( secPerYear, -mf );
}


/***************************************************************************\
 **  tSedTransPwrLaw2::TransCapacity
 **
 **  Computes sediment transport capacity using the simple power law
 **  Qs = kf W ( [kt (Q/W)^mf S^nf]^pf - tauc^pf )
 **
\***************************************************************************/
double tSedTransPwrLaw2::TransCapacity( tLNode *node )
{
  const double slp = node->calcSlope();
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw::TransCapacity(tLNode*)");
  double tau, tauexpf, cap = 0.;
  if( node->getFloodStatus() == tLNode::kNotFlooded )
    {
      tau = kt * pow( node->getQ()/node->getHydrWidth(), mf ) * pow( slp, nf );
      node->setTau( tau );
      if(0) //DEBUG
	std::cout << "kt=" << kt << " Q=" << node->getQ() << " W="
	     << node->getHydrWidth() << " S=" << node->calcSlope() << std::endl;
      tauexpf = pow( tau, pf ) - pow( tauc, pf );
      tauexpf = (tauexpf>0.0) ? tauexpf : 0.0;
      cap = kf * node->getHydrWidth() * tauexpf;
    }
  node->setQs( cap );
  return cap;
}


/***************************************************************************\
 **  tSedTransPwrLaw2::TransCapacity
 **
 **  Computes sediment transport capacity using the simple power law
 **  Qs = weight kf W ( [ kt (Q/W)^mf S^nf ]^pf - tauc^pf )
 **  This is a weighted version which is called from DetachErode.
 **  Weight is a weighting by depth of layer.
 **  Here, qsi is set by proportion in layer and threshold is constant
 **  The value returned should be in units of m^3/yr
 **
\***************************************************************************/
double tSedTransPwrLaw2::TransCapacity( tLNode *node, int lyr, double weight )
{
  const double slp = node->calcSlope();
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tSedTransPwrLaw::TransCapacity(tLNode*)");
  double tau, tauexpf, cap = 0.;

  if( node->getFloodStatus() == tLNode::kNotFlooded )
    {
      tau = kt * pow( node->getQ()/node->getHydrWidth(), mf ) * pow( slp, nf );
      node->setTau( tau );
      if(0) //DEBUG
	std::cout << "kt=" << kt << " Q=" << node->getQ() << " W="
	     << node->getHydrWidth() << " S=" << node->calcSlope()
	     << " tau=" << tau << std::endl;
      tauexpf = pow( tau, pf ) - pow( tauc, pf );
      tauexpf = (tauexpf>0.0) ? tauexpf : 0.0;
      cap = weight * kf * node->getHydrWidth() * tauexpf;
    }
  for(size_t i=0; i<node->getNumg(); i++)
    node->addQs(i, cap*node->getLayerDgrade(lyr,i)/node->getLayerDepth(lyr));

  node->setQs( cap );
  return cap;
}


/***************************************************************************\
 **  FUNCTIONS FOR CLASS tSedTransBridgeDom
\***************************************************************************/

/***************************************************************************\
 **  tSedTransBridgeDom Constructor:
 **
 **    Given a tInputFile as an argument, will read relevant parameters from
 **  the input file.
 **
\***************************************************************************/
tSedTransBridgeDom::tSedTransBridgeDom( const tInputFile &infile )
{
  const double secPerYear = SECPERYEAR;  // # secs in one year

  kf = infile.ReadItem( kf, "KF" );
  kt = infile.ReadItem( kt, "KT" );
  mf = infile.ReadItem( mf, "MF" );
  nf = infile.ReadItem( nf, "NF" );
  tauc = infile.ReadItem( tauc, "TAUCD" );
  sqrtTauc = sqrt( tauc );

  // Add unit conversion factor for kt -- this is required to convert
  // the quantity (Q/W)^mb from units of years to units of seconds.
  kt *= pow( secPerYear, -mf );
}


/***************************************************************************\
 **  tSedTransBridgeDom::TransCapacity
 **
 **  Computes sediment transport capacity using the Bridge-Dominic
 **  form of the Bagnold bedload equation:
 **
 **  Qs = kf' W ( tau - taucrit ) * ( ustar - ustarcrit )
 **    ustar = ( tau / rho ) ^ 0.5
 **  Qs = kf W ( tau - taucrit ) * ( sqrt( tau ) - sqrt( taucrit ) )
 **    kf = kf' / rho^0.5
 **  tau = kt * (Q/W)^mf * S^nf
 **
 **  Note: assumes fixed, uniform tauc
 **
\***************************************************************************/
double tSedTransBridgeDom::TransCapacity( tLNode *node )
{
  const double slp = node->calcSlope();
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw::TransCapacity(tLNode*)");
  double tau, tauex, ustarex, cap = 0.;
  if( node->getFloodStatus() == tLNode::kNotFlooded )
    {
      tau = kt * pow( node->getQ()/node->getHydrWidth(), mf ) * pow( slp, nf );
      node->setTau( tau );
      if(0) //DEBUG
	std::cout << "kt=" << kt << " Q=" << node->getQ() << " W="
	     << node->getHydrWidth() << " S=" << node->calcSlope() << std::endl;
      tauex = ( tau > tauc ) ? (tau - tauc) : 0.0;
      ustarex = ( tau > tauc ) ? (sqrt(tau) - sqrtTauc) : 0.0;
      cap = kf * node->getHydrWidth() * tauex * ustarex;
    }
  node->setQs( cap );
  return cap;
}


/***************************************************************************\
 **  tSedTransBridgeDom::TransCapacity
 **
 **  Computes sediment transport capacity using the Bridge-Dominic
 **  form of the Bagnold bedload equation:
 **
 **  Qs = kf' W ( tau - taucrit ) * ( ustar - ustarcrit )
 **    ustar = ( tau / rho ) ^ 0.5
 **  Qs = kf W ( tau - taucrit ) * ( sqrt( tau ) - sqrt( taucrit ) )
 **    kf = kf' / rho^0.5
 **  tau = kt * (Q/W)^mf * S^nf
 **
 **  Note: assumes fixed, uniform tauc
 **
 **  This is a weighted version which is called from DetachErode.
 **  Weight is a weighting by depth of layer.
 **  Here, qsi is set by proportion in layer and threshold is constant
 **  The value returned should be in units of m^3/yr
 **
\***************************************************************************/
double tSedTransBridgeDom::TransCapacity( tLNode *node, int lyr, double weight )
{
  const double cap = weight * TransCapacity( node );

  for(size_t i=0; i<node->getNumg(); i++)
    node->addQs(i, cap*node->getLayerDgrade(lyr,i)/node->getLayerDepth(lyr));

  //node->setQs( cap ); // <- this sets it to the weighted version ... correct?
  return cap;
}


/*************************************************************************\
 **  FUNCTIONS FOR CLASS tSedTransPwrLawMulti
\**************************************************************************/

/*************************************************************************\
 **
 **  tSedTransPwrLawMulti constructor
 **
\**************************************************************************/
tSedTransPwrLawMulti::tSedTransPwrLawMulti( const tInputFile &infile )
{
  const double secPerYear = SECPERYEAR;  // # secs in one year

  kf = infile.ReadItem( kf, "KF" );
  kt = infile.ReadItem( kt, "KT" );
  mf = infile.ReadItem( mf, "MF" );
  nf = infile.ReadItem( nf, "NF" );
  pf = infile.ReadItem( pf, "PF" );
  miNumgrnsizes = infile.ReadItem( miNumgrnsizes, "NUMGRNSIZE" );
  if( miNumgrnsizes>9 )
    {
      std::cout << "WARNING: maximum of 9 grain size classes exceeded.\n";
      std::cout << "Resetting to 9 size-fractions.\n";
      std::cout << "(That was a non-fatal warning, my friend!)\n";
      miNumgrnsizes = 9;
    }

  // Record diameter of each size-fraction
  mdGrndiam.setSize( miNumgrnsizes );
  mdTauc.setSize( miNumgrnsizes );
  //std::string taglinebase = "GRAINDIAM";
  //std::string digits = "123456789";
  //std::string tagline;
  char tagline[11], digit = '0';
  strcpy( tagline, "GRAINDIAM0");
  int i;
  const double thetac = 0.045,
    sig = RHOSED,
    rho = RHO,
    g = GRAV;
  for( i=0; i<miNumgrnsizes; i++ )
    {
      /*tagline = taglinebase + digits.substr(i,i);
	mdGrndiam[i] = infile.ReadItem( mdGrndiam[i], tagline.c_str() );*/
      digit++;
      tagline[9] = digit;
      mdGrndiam[i] = infile.ReadItem( mdGrndiam[i], tagline );
      mdTauc[i] = thetac * (sig-rho) * g * mdGrndiam[i];
      if(0) //DEBUG
	std::cout << "Diam " << i << " = " << mdGrndiam[i] << " tauc = "
	     << mdTauc[i] << std::endl;
    }

  // Add unit conversion factor for kt -- this is required to convert
  // the quantity (Q/W)^mb from units of years to units of seconds.
  kt *= pow( secPerYear, -mf );

  // Read hiding/protrusion exponent (should be btwn 0 to 1)
  mdHidingexp = infile.ReadItem( mdHidingexp, "HIDINGEXP" );

}



/***************************************************************************\
 **  tSedTransPwrLawMulti::TransCapacity
 **
 **
 **  Created: 02/02 GT & RAJR
 **
\***************************************************************************/
double tSedTransPwrLawMulti::TransCapacity( tLNode *node, int lyr, double weight )
{
  const double slp = node->calcSlope();
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tSedTransPwrLaw::TransCapacity(tLNode*)");
  double d50 = 0.0,   // Mean grain size
    tau,              // Shear stress
    tauex;            // Excess shear stress
  tArray<double> frac( miNumgrnsizes );
  int i;

  // Compute D50 and fraction of each size
  assert( node->getLayerDepth(lyr) > 0.0 );
  for( i=0; i<miNumgrnsizes; i++ )
    {
      frac[i] = node->getLayerDgrade(lyr,i) / node->getLayerDepth(lyr);
      assert( frac[i]>=0.0 );
      d50 += frac[i] * mdGrndiam[i];
      if( 0 )
	{
	  std::cout << "uh oh2: " << node->getLayerDgrade(lyr,0) << " "
	       << node->getLayerDgrade(lyr,1) << std::endl;
	  std::cout << frac[0] << " " << frac[1] << std::endl;
	}
      if(0) //DEBUG
	std::cout << "frac " << i << " = " << frac[i] << std::endl;
    }
  assert( d50>=0.0 );
  assert( d50<1e10 );
  if(0) //DEBUG
    std::cout << "D50 = " << d50 << std::endl;

  // Compute shear stress
  if( node->getFloodStatus() ) tau = 0.0;
  else {
    assert( node->getHydrWidth() > 0.0 );
    assert( slp >= 0.0 );
    tau = kt * pow( node->getQ()/node->getHydrWidth(), mf ) * pow( slp, nf );
  }
  node->setTau( tau );
  if(0) //DEBUG
    std::cout << "kt=" << kt << " Q=" << node->getQ() << " W="
	 << node->getHydrWidth() << " S=" << node->calcSlope()
	 << " tau=" << tau << std::endl;

  // Compute crit shear stress and xport capacity for each size fraction
  double totalcap = 0.0,
    tauc;
  for( i=0; i<miNumgrnsizes; i++ )
    {
      tauc = mdTauc[i] * pow( mdGrndiam[i] / d50, -mdHidingexp );
      if(0) //DEBUG
	std::cout << "tauc " << i << " = " << tauc << std::endl;
      tauex = tau - tauc;
      tauex = (tauex>0.0) ? tauex : 0.0;
      // Xport capacity for a given size
      const double cap = frac[i] * weight * kf * node->getHydrWidth() * pow( tauex, pf );
      totalcap += cap;
      node->addQs( i, cap );
    }

  node->setQs( totalcap );
  return totalcap;
}


double tSedTransPwrLawMulti::TransCapacity( tLNode * /* node */ )
{
  return 0.0;
}


/*************************************************************************\
 **  FUNCTIONS FOR CLASS tSedTransWilcock
\**************************************************************************/

/*************************************************************************\
 **
 **  tSedTransWilcock constructor
 **
\**************************************************************************/
tSedTransWilcock::tSedTransWilcock( const tInputFile &infile )
  : grade(2)
{
  if(0) //DEBUG
    std::cout << "tSedTransWilcock(infile)\n" << std::endl;
  //strcpy( add, "1" );  // GT changed from add = '1' to prevent glitch
  /*for(i=0; i<=1; i++){
    strcpy( name, "GRAINDIAM");
    strcat( name, add );
    help = infile.ReadItem( help, name);
    add[0]++;
    grade[i] = help;
    }*/
  grade[0] = infile.ReadItem( grade[0], "GRAINDIAM1" );
  grade[1] = infile.ReadItem( grade[1], "GRAINDIAM2" );

  taudim= RHO*GRAV;
  refs = (RHOSED-RHO)*GRAV*grade[0];
  refg = (RHOSED-RHO)*GRAV*grade[1];
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
 **
 **  tSedTransWilcock::TransCapacity
 **
 **  This function uses the sediment transport model developed by
 **  P. Wilcock to calculate sediment transport rates of the sand and
 **  gravel fraction individually.  This function should only be used with
 **  two grain sizes and it is assumed that the grain size one is in the
 **  sand range and grain size 2 is in the gravel range.  The sediment
 **  transport rate of both grain sizes is calculated, and the sum of
 **  these two rates is returned. (rate here is in m^3/yr)
 **
 **  Modifications:
 **   - Reverted to earlier computation of tau using roughness and
 **     width
\***********************************************************************/
double tSedTransWilcock::TransCapacity( tLNode *nd )
{
  if( nd->calcSlope() < 0 ){
    nd->setQs(0, 0.);
    nd->setQs(1, 0.);
    nd->setQs(0.);
    return 0.0;
  }

  double tau;
  double taucrit;
  double persand=nd->getLayerDgrade(0,0)/(nd->getLayerDepth(0));
  double factor=nd->getLayerDepth(0)/nd->getMaxregdep();

  // units of Q are m^3/yr; convert to m^3/sec
  //NIC you are doing a test here to see what is causing the
  //downstream coarsening.
  tau = taudim*pow(nd->getHydrRough()*nd->getQ()/SECPERYEAR/nd->getHydrWidth(), 0.6)*pow( nd->calcSlope(), 0.7);

  if(0) { //DEBUG
    std::cout << "hydrrough is " << nd->getChanRough() << std::endl;
    std::cout << "q is " << nd->getQ() << std::endl;
    std::cout << "slope is " << nd->calcSlope() << std::endl;
    std::cout << "taudim is " << taudim << std::endl;
  }

  //Calculate Sand transport rates first
  if(persand<.10)
    taucrit=lowtaucs;
  else if(persand<=.40)
    taucrit=((sands*persand)+sandb);
  else
    taucrit=hightaucs;

  if(tau>taucrit){
    nd->setQs(0, ((0.058/RHOSED)*factor*nd->getHydrWidth()*SECPERYEAR*persand*pow(tau,1.5)*pow((1-sqrt(taucrit/tau)),4.5) ));
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

  if(tau>taucrit){
    nd->setQs(1, (0.058*SECPERYEAR*factor*nd->getHydrWidth()/(RHOSED))*
	      (1-persand)*pow(tau,1.5)*pow((1-(taucrit/tau)),4.5));
  }
  else
    nd->setQs(1,0.0);

  nd->setQs(nd->getQs(0)+nd->getQs(1));
  return nd->getQs();

}


/*********************************************************************\
 **
 **  tSedTransWilcock::TransCapacity
 **
 **  *nd - pointer to the node which you are calculating tranport rates at.
 **  i - the layer which you are basing the transport rate on (for texture)
 **  weight - used in erosion algorithm, a weight based on layer depths.
 **
 **  This function uses the sediment transport model developed by
 **  P. Wilcock to calculate sediment transport rates of the sand and
 **  gravel fraction individually.  This function should only be used with
 **  two grain sizes and it is assumed that grain size one is in the
 **  sand range and grain size 2 is in the gravel range.  The sediment
 **  transport rate of both grain sizes is calculated, and the sum of
 **  these two rates is returned. (rate here is in m^3/yr)
 **  Note that this function assumes that you are looping through layers,
 **  (which is why you need the weight) and so qs total and for each size
 **  was initialized to zero and you just add to it in this function.
 **  It is VERY IMPORTANT that qs is reset to zero before you use this
 **  funtion in a loop.
\***********************************************************************/
double tSedTransWilcock::TransCapacity( tLNode *nd, int i, double weight )
{
  if(0) //DEBUG
    std::cout << "tSedTransWilcock::TransCapacity(tLNode,int,double)\n";

  if( nd->calcSlope() < 0 ){
    nd->setQs(0, 0.);
    if(nd->getNumg()==2)
      nd->setQs(1, 0.);
    nd->setQs(0.);
    return 0.0;
  }

  double tau;
  double taucrit;
  assert( nd->getLayerDepth(i)>0 );
  double persand=nd->getLayerDgrade(i,0)/(nd->getLayerDepth(i));
  double qss, qsg=0.; //gravel and sand transport rate

  // units of Q are m^3/yr; convert to m^3/sec
  //tau = taudim*pow(nd->getHydrRough()*nd->getQ()/nd->getHydrWidth(), 0.6)*pow( nd->calcSlope(), 0.7);
  tau = taudim*pow(0.03, 0.6)*pow(nd->getQ()/SECPERYEAR, 0.3)*pow( nd->calcSlope(), 0.7);

  if(0) { //DEBUG
    std::cout << "channel rough is " << nd->getChanRough() << std::endl;
    std::cout << "channel width is " << nd->getChanWidth() << std::endl;
    std::cout << "q in secs is " << nd->getQ()/SECPERYEAR << std::endl;
    std::cout << "slope is " << nd->calcSlope() << std::endl;
    std::cout << "taudim is " << taudim << std::endl;
  }

  //Calculate Sand transport rates first

  if(persand<.10)
    taucrit=lowtaucs;
  else if(persand<=.40)
    taucrit=((sands*persand)+sandb);
  else
    taucrit=hightaucs;

  if(tau>taucrit){
    qss=(0.058/RHOSED)*weight*nd->getHydrWidth()*SECPERYEAR*persand*pow(tau,1.5)*pow((1-sqrt(taucrit/tau)),4.5) ;
    nd->addQs(0, qss);
  }
  else
    qss=0.;

  //Now calculate Gravel transport rates
  if(nd->getNumg()==2){

    if(persand<.10)
      taucrit=lowtaucg;
    else if(persand<=.40)
      taucrit=((gravs*persand)+gravb);
    else
      taucrit=hightaucg;

    if(tau>taucrit){
      qsg=(0.058*SECPERYEAR*weight*nd->getHydrWidth()/(RHOSED))*
	(1-persand)*pow(tau,1.5)*pow((1-(taucrit/tau)),4.5);
      nd->addQs(1,qsg);
    }
    else
      qsg=0.;
  }

  //NOTE - don't need to update total qs cause this gets updates
  //with update of qs of individual sizes

  return qsg+qss;

}

/*************************************************************************\
 **  FUNCTIONS FOR CLASS tSedTransMineTailings
\**************************************************************************/

/*************************************************************************\
 **
 **  tSedTransMineTailings constructor
 **
 **  note that for now this is the same as tSedTransWilcock
 **  constructor since it is using same taucrit function
\**************************************************************************/
tSedTransMineTailings::tSedTransMineTailings( const tInputFile &infile )
  : grade(2)
{
  int i;
  char add[2], name[20];
  double help;

  if(0) //DEBUG
    std::cout << "tSedTransMineTailings(infile)\n" << std::endl;
  strcpy( add, "1" );  // GT changed from add = '1' to prevent glitch
  for(i=0; i<=1; i++){
    strcpy( name, "GRAINDIAM");
    strcat( name, add );
    help = infile.ReadItem( help, name);
    add[0]++;
    grade[i] = help;
  }

  taudim= RHO*GRAV;
  refs = (RHOSED-RHO)*GRAV*grade[0];
  refg = (RHOSED-RHO)*GRAV*grade[1];
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
 **
 **  tSedTransMineTailings::TransCapacity
 **
 **  This function uses the sediment transport model developed on
 **  mine tailings slopes in Australia by Willgoose and Riley (1998)
 **  to calculate sediment transport rates of the sand and
 **  gravel fraction individually.
 **
 **  For now This function should only be used with
 **  two grain sizes because it uses the Wilcock method to calculate taucrit.
 **  It is assumed that grain size one is in the
 **  sand range and grain size 2 is in the gravel range.  The sediment
 **  transport rate of both grain sizes is calculated, and the sum of
 **  these two rates is returned. (rate here is in m^3/yr)
 **
\***********************************************************************/
double tSedTransMineTailings::TransCapacity( tLNode *nd )
{
  if( nd->calcSlope() < 0 ){
    nd->setQs(0, 0.);
    nd->setQs(1, 0.);
    nd->setQs(0.);
    return 0.0;
  }

  double tau;
  double taucrit;
  double persand=nd->getLayerDgrade(0,0)/(nd->getLayerDepth(0));
  //not sure what the point of factor was
  //double factor=nd->getLayerDepth(0)/nd->getMaxregdep();

  // units of Q are m^3/yr; convert to m^3/sec
  //tau = taudim*pow(nd->getHydrRough()*nd->getQ()/SECPERYEAR/nd->getHydrWidth(), 0.6)*pow( nd->calcSlope(), 0.7);
  tau = taudim*pow(0.03, 0.6)*pow(nd->getQ()/SECPERYEAR, 0.3)*pow( nd->calcSlope(), 0.7);

  if(0) { //DEBUG
    std::cout << "hydrrough is " << nd->getChanRough() << std::endl;
    std::cout << "q is " << nd->getQ() << std::endl;
    std::cout << "slope is " << nd->calcSlope() << std::endl;
    std::cout << "taudim is " << taudim << std::endl;
  }

  //Calculate Sand transport rates first

  if(persand<.10)
    taucrit=lowtaucs;
  else if(persand<=.40)
    taucrit=((sands*persand)+sandb);
  else
    taucrit=hightaucs;

  if(tau>taucrit){
    nd->setQs(0,(0.0541/RHOSED)*SECPERYEAR*persand*pow(nd->getQ()/SECPERYEAR,1.12)*pow(nd->calcSlope(),-0.24)*(tau-taucrit) );
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

  if(tau>taucrit){
    nd->setQs(1, (0.0541/RHOSED)*SECPERYEAR*(1-persand)*pow(nd->getQ()/SECPERYEAR,1.12)*pow(nd->calcSlope(),-0.24)*(tau-taucrit));
  }
  else
    nd->setQs(1,0.0);

  nd->setQs(nd->getQs(0)+nd->getQs(1));
  return nd->getQs();

}

/*********************************************************************\
 **
 **  tSedTransMineTailings::TransCapacity
 **
 **  *nd - pointer to the node which you are calculating tranport rates at.
 **  i - the layer which you are basing the transport rate on (for texture)
 **  weight - used in erosion algorithm, a weight based on layer depths.
 **
 **  This function uses the sediment transport model and parameters in
 **  Willgoose and Riley (1998) to calculate sediment transport rates
 **  of the sand and gravel fraction individually.
 **  For now, this function should only be used with two grain sizes
 **  (it is using the wilcock model to calculate tau crit)
 **  and it is assumed that grain size one is in the
 **  sand range and grain size 2 is in the gravel range.  The sediment
 **  transport rate of both grain sizes is calculated, and the sum of
 **  these two rates is returned. (rate here is in m^3/yr)
 **  Note that this function assumes that you are looping through layers,
 **  (which is why you need the weight) and so qs total and for each size
 **  was initialized to zero and you just add to it in this function.
 **  It is VERY IMPORTANT that qs is reset to zero before you use this
 **  funtion in a loop.
\***********************************************************************/
double tSedTransMineTailings::TransCapacity( tLNode *nd, int i, double weight )
{
  if(0) //DEBUG
    std::cout << "tSedTransMineTailings::TransCapacity(tLNode,int,double)\n";

  if( nd->calcSlope() < 0 ){
    nd->setQs(0, 0.);
    if(nd->getNumg()==2)
      nd->setQs(1, 0.);
    nd->setQs(0.);
    return 0.0;
  }

  double tau;
  double taucrit;
  assert( nd->getLayerDepth(i)>0.0 );
  double persand=nd->getLayerDgrade(i,0)/(nd->getLayerDepth(i));
  double qss, qsg=0.; //gravel and sand transport rate

  // units of Q are m^3/yr; convert to m^3/sec
  //NIC check to see what taudim is -> probably right but U R anal
  tau = taudim*pow(0.03, 0.6)*pow(nd->getQ()/SECPERYEAR, 0.3)*pow( nd->calcSlope(), 0.7);

  if(0) { //DEBUG
    std::cout << "Q is " << nd->getQ() << std::endl;
    std::cout << "slope is " << nd->calcSlope() << std::endl;
    std::cout << "taudim is " << taudim << std::endl;
    std::cout << "persand is " << persand << std::endl;
    std::cout << "weight is " << weight << std::endl;
  }

  //Calculate Sand transport rates first

  // here calculating critical shear stress, same as with wilcock
  if(persand<.10)
    taucrit=lowtaucs;
  else if(persand<=.40)
    taucrit=((sands*persand)+sandb);
  else
    taucrit=hightaucs;

  //remember tau is in units of sec, so compute everything in seconds
  //and transfer back to years in the end.
  //Nic, notice what you did, the commented out qss is what you originally,
  //and still kind of think is right.  Just testing.  Difference is that
  //you aren't multiply by the percentage that is there.
  //now original equation, but multiply by 10 to speed things up.
  if(tau>taucrit){
    //qss=(0.0541/RHOSED)*weight*SECPERYEAR*pow(nd->getQ()/SECPERYEAR,1.12)*pow(nd->calcSlope(),-0.24)*(tau-taucrit);
    qss=(0.0541/RHOSED)*weight*SECPERYEAR*persand*pow(nd->getQ()/SECPERYEAR,1.12)*pow(nd->calcSlope(),-0.24)*(tau-taucrit);
    nd->addQs(0, qss);
  }
  else{
    qss=0.;
  }

  //Now calculate Gravel transport rates
  if(nd->getNumg()==2){
    if(persand<.10)
      taucrit=lowtaucg;
    else if(persand<=.40)
      taucrit=((gravs*persand)+gravb);
    else
      taucrit=hightaucg;

    if(tau>taucrit){
      //qsg=(0.0541/RHOSED)*weight*SECPERYEAR*pow(nd->getQ()/SECPERYEAR,1.12)*pow(nd->calcSlope(),-0.24)*(tau-taucrit);
      qsg=(0.0541/RHOSED)*weight*SECPERYEAR*(1-persand)*pow(nd->getQ()/SECPERYEAR,1.12)*pow(nd->calcSlope(),-0.24)*(tau-taucrit);
      nd->addQs(1,qsg);
    }
    else{
      qsg=0.;
    }

  }

  //NOTE - don't need to update total qs cause this gets updates
  //with update of qs of individual sizes

  return qsg+qss;

}

/**************************************************************************\
 **
 **  Functions for Class tMechanicalWeather; Wr=Wo.exp(-kwc*regdepth)
 **
\**************************************************************************/

//constructor
tMechanicalWeather::tMechanicalWeather()
{
  kwc=1;
  Wo=1;
  K=1;
}

tMechanicalWeather::tMechanicalWeather( const tInputFile &infile)
{
  kwc=infile.ReadItem(kwc,"KWc");
  Wo=infile.ReadItem(Wo,"W0");
  K=infile.ReadItem(K,"KR");
}

tMechanicalWeather::~tMechanicalWeather()
{}

double tMechanicalWeather::getK() const
{
  return K;
}

void tMechanicalWeather::setK(double nk )
{
  K=nk;
}  

double tMechanicalWeather::WeatheringRate(tLNode * n)
{
  double regdepth = n->getLayerDepth(0), Wr;
  double Nlayer = n->getNumLayer();   
  if(Nlayer>1)
    {
      Wr=Wo*exp(-1*kwc*regdepth);
    }
  else
    {
      Wr=Wo;
    }  
  return Wr;
}

/***************************************************************************\
 **  FUNCTIONS FOR CLASS tErosion
\***************************************************************************/

//constructor
tErosion::tErosion( tMesh<tLNode> *mptr, const tInputFile &infile ) :
  meshPtr(mptr),
  bedErode(0), sedTrans(0), weather(infile)
{
  assert( mptr!=0 );

  // Read parameters needed from input file
  kd = infile.ReadItem( kd, "KD" );  // Hillslope diffusivity coefficient
  D_alpha = infile.ReadItem( D_alpha, "D_ALPHA" );  // Hillslope diffusivity exponent
  int optAdaptMesh = infile.ReadItem( optAdaptMesh, "OPTMESHADAPTDZ" );
  if( optAdaptMesh )
    mdMeshAdaptMaxFlux = infile.ReadItem( mdMeshAdaptMaxFlux,
					  "MESHADAPT_MAXNODEFLUX" );

  // Make sure the user wants the detachment and transport options that
  // are compiled in this version
  int optProcessLaw = infile.ReadItem( optProcessLaw,
				       "DETACHMENT_LAW" );
  switch(optProcessLaw){
#define X(a,b) case a: \
               bedErode = new b(infile); \
               break;
    DETACHMENT_LAW_TABLE2
#undef X
      default:
    {
      std::cerr << "\nError: You requested the detachment law '"
	   << optProcessLaw << "' which does not exist."  << std::endl
	   << "Available options:" << std::endl;
      for(int i=0; i!=NUMBER_OF_DETACHMENT_LAWS; ++i ){
	std::cerr << " " << i << ": " << DetachmentLaw[i] << std::endl;
      }
      ReportFatalError( "Requested detachment law not available. "
			"Switch options.\n" );
    }
  }

  std::cout << "DETACHMENT OPTION: "
       << DetachmentLaw[optProcessLaw] << std::endl;

  optProcessLaw = infile.ReadItem( optProcessLaw,
				   "TRANSPORT_LAW" );
  switch(optProcessLaw){
#define X(a,b) case a: \
               sedTrans = new b(infile); \
               break;
    TRANSPORT_LAW_TABLE2
#undef X
      default:
    {
      std::cerr << "\nError: You requested the transport law '"
	   << optProcessLaw << "' which does not exist."  << std::endl
	   << "Available options:" << std::endl;
      for(int i=0; i!=NUMBER_OF_TRANSPORT_LAWS; ++i ){
	std::cerr << " " << i << ": " << TransportLaw[i] << std::endl;
      }
      ReportFatalError( "Requested transport law not available. "
			"Switch options.\n" );
    }
  }
  std::cout << "SEDIMENT TRANSPORT OPTION: "
       << TransportLaw[optProcessLaw] << std::endl;
}

tErosion::~tErosion(){
  meshPtr = 0;
  delete bedErode;
  delete sedTrans;
}

/*****************************************************************************\
 **
 **  tErosion::ErodeDetachLim
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
 **
 **   - added calls to compute channel width (& depth etc) before computing
 **     erosion. This is done because the detachment capacity functions now
 **     require a defined channel width. (GT 2/01)
\*****************************************************************************/
void tErosion::ErodeDetachLim( double dtg, tStreamNet *strmNet,
			       tVegetation * /*pVegetation*/ )
{
  if(1) //DEBUG
    std::cout<<"ErodeDetachLim...";
  double dt,
    dtmax; // time increment
  const double frac = 0.9; //fraction of time to zero slope
  tLNode *cn, *dn;
  //int nActNodes = meshPtr->getNodeList()->getActiveSize();
  tMesh< tLNode >::nodeListIter_t ni( meshPtr->getNodeList() );

  strmNet->FindChanGeom();
  strmNet->FindHydrGeom();

  tArray<double> valgrd(1);
  //TODO: make it work w/ arbitrary # grain sizes

  // Iterate until total time dtg has been consumed
  int debugCount=0;
  do
    {
      //first find erosion rate:
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
	cn->setDzDt( -bedErode->DetachCapacity( cn ) );

      //find max. time step s.t. slope does not reverse:
      dtmax = dtg;
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
	{
	  dn = cn->getDownstrmNbr();
	  const double ratediff = dn->getDzDt() - cn->getDzDt();
	  if( ratediff > 0 )
	    {
	      dt = ( cn->getZ() - dn->getZ() ) / ratediff * frac;
	      if( dt > 0.000005 && dt < dtmax ) dtmax = dt;
	    }
	}
      //assert( dtmax > 0 );

      //apply erosion:
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() ){
	//	if( cn->getX()>7300.0 && cn->getX()<7400.0 && cn->getY()>500.0 && cn->getY()<600.0 )
	//	  cn->TellAll();

	//ng added stuff below to update layering using the other erodep
	valgrd[0]=cn->getDzDt() * dtmax;
	cn->EroDep( 0, valgrd, 0.);
	//cn->EroDep( cn->getQs() * dtmax );
      }

      // Update veg
#if 0
#define NEWVEG 0
      if( pVegetation && NEWVEG ) pVegetation->ErodeVegetation( meshPtr, dtmax );
#undef NEWVEG
#endif

      //update time:
      dtg -= dtmax;
	  
	  if(1) //DEBUG
	  {
	     debugCount++;
		 if( debugCount > 1e6 )
		    ReportFatalError("More than 1e6 iterations in ErodeDetachLim()" );
	  }
	  
    } while( dtg>0.0000001 );

}//end tErosion::ErodeDetachLim( double dtg


/*****************************************************************************\
 **
 **  tErosion::ErodeDetachLim (2 of 2)
 **
 **  This version adds the uplift rate source term to the time-step
 **  estimation. (This could be handled as a default argument, avoiding
 **  the need for two nearly identical versions -- TODO)
 **
 **   - added calls to compute channel width (& depth etc) before computing
 **     erosion. This is done because the detachment capacity functions now
 **     require a defined channel width. (GT 2/01)
\*****************************************************************************/
void tErosion::ErodeDetachLim( double dtg, tStreamNet *strmNet, tUplift const *UPtr )
{
  double dt,
    dtmax; // time increment
  double frac = 0.1; //fraction of time to zero slope
  //Xint i;
  tLNode * cn, *dn;
  //int nActNodes = meshPtr->getNodeList()->getActiveSize();
  tMesh< tLNode >::nodeListIter_t ni( meshPtr->getNodeList() );
  double ratediff;
  //Xdouble dslpdt;
  double dtmin = dtg * 0.0001;
  int debugCount = 0;

  strmNet->FindChanGeom();
  strmNet->FindHydrGeom();

  tArray<double> valgrd(1);
  // Iterate until total time dtg has been consumed
  do
    {
      //first find erosion rate:
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
	cn->setDzDt( -bedErode->DetachCapacity( cn ) );
      dtmax = dtg;
      //find max. time step s.t. slope does not reverse:
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
	{
	  /*slp = cn->calcSlope();
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
		  if(0) //DEBUG
		    std::cout << "time step too small because of node at x,y,z "
			 << cn->getX() << " " << cn->getY() << " " << cn->getZ()
			 << std::endl;
		}
	    }
	}
      //assert( dtmax > 0 );
      //apply erosion:
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() ){
	//ng added stuff below to update layering using the other erodep
	valgrd[0]=cn->getDzDt() * dtmax;
	cn->EroDep( 0, valgrd, 0.);
	//cn->EroDep( cn->getDzDt() * dtmax );
      }
      //update time:
      dtg -= dtmax;
	  
	  if(1) //DEBUG
	  {
	     debugCount++;
		 if( debugCount > 1e6 )
		    ReportFatalError("More than 1e6 iterations in ErodeDetachLim()" );
	  }
	  
    } while( dtg>0 );

}//end tErosion::ErodeDetachLim( double dtg, tUplift *UPtr )


/*****************************************************************************\
 **
 **  tErosion::StreamErode
 **
 **  Solution for general detachment- or transport-limited erosion. Now
 **  replaced by DetachErode.
 **
\*****************************************************************************/
#define kSmallTimeStep 1e-8
void tErosion::StreamErode( double dtg, tStreamNet *strmNet )
{
  double dt,
    dtmax;         // time increment: initialize to arbitrary large val
  double frac = 0.3; // fraction of time to zero slope
  tLNode * cn, *dn;
  tMesh< tLNode >::nodeListIter_t ni( meshPtr->getNodeList() );
  double ratediff,  // Difference in ero/dep rate btwn node & its downstrm nbr
    cap,          // Transport capacity
    pedr,         // Potential erosion/deposition rate
    dcap,         // Bedrock detachment capacity
    dz,           // Depth of deposition/erosion (erosion = negative)
    dzr;          // Potential depth of bedrock erosion
  bool smallflag=false;
  int smallcount=0;

  std::cout << "tErosion::StreamErode\n";

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
	  cap = sedTrans->TransCapacity( cn );
	  pedr = (cn->getQsin() - cap ) / cn->getVArea();
	  //sediment input:
	  if( cn == strmNet->getInletNodePtr() )
	    pedr += strmNet->getInSedLoad() / cn->getVArea();

	  // If we're on bedrock, adjust accordingly
	  if( cn->OnBedrock() && pedr<0 )
	    {
	      // get detachment capacity (this also sets node's drdt)
	      dcap = -bedErode->DetachCapacity( cn );
	      if( dcap > pedr )
                pedr = dcap;
	    }
	  // set the erosion (deposition) rate and send the corresponding
	  // sediment influx downstream
	  cn->setDzDt( pedr );
	  cn->getDownstrmNbr()->addQsin( cn->getQsin() - pedr*cn->getVArea() );
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
		std::cout << "Very small dt " << dt << " at:\n";
                cn->TellAll();
                dn->TellAll();
		}*/

	    }
	}
      dtmax *= frac;  // Take a fraction of time-to-flattening
      if( dtmax < kSmallTimeStep ) dtmax = kSmallTimeStep;
      if( dtmax <= 0.01 && !smallflag )
	{
	  smallflag=true;
	  std::cout << "SMALL STEP: " << dtmax << std::endl;
	}
      if( smallflag )
	{
	  smallcount++;
	  if( smallcount==100 )
	    {
	      std::cout << "TIME REMAINING: " << dtg << std::endl;
	      smallcount=0;
	    }
	}

      //std::cout << "  dt " << dtmax << std::endl;

      // Zero out sed influx again, because depending on bedrock-alluvial
      // interaction it may be modified; if inlet, give it strmNet->inlet.inSedLoad
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
	cn->setQsin( 0.0 );
      //sediment input:
      if(strmNet->getInletNodePtrNC() != NULL)
	{
          strmNet->getInletNodePtrNC()->setQsin( strmNet->getInSedLoad() );
	  std::cout << "Inlet node:\n";
	  strmNet->getInletNodePtr()->TellAll();
	}

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

	  //std::cout << "** THIS node has dz " << dz << std::endl;
	  //cn->TellAll();

	  // Update alluvium thickness and node elevation
	  /*if( cn->getID()==3214 )
	    {
	    cn->TellAll();
	    std::cout << "Applying dz=" << dz << std::endl;
	    }
	    if( cn->getZ() < -1 ) {
            std::cout << "The following node is going below baselevel:\n";
            cn->TellAll();
	    }*/
	  cn->EroDep( dz );
	  dn = cn->getDownstrmNbr();
	  /*if( dn->getID()==3214 )
	    {
            std::cout << "SENDing to 3214: " << cn->getQsin() << " - " << dz*cn->getVArea()/dtmax << std::endl;
            cn->TellAll();
	    }*/

	  // Send sediment downstream: sediment flux is equal to the flux in
	  // plus/minus rate of erosion/deposition times node area
	  assert( dtmax>0 );
	  dn->addQsin( cn->getQsin() - dz*cn->getVArea()/dtmax );
	  //           if( (cn->getQsin() - dz*cn->getVArea()/dtmax) < -0.1 )
	  //           {
	  //              std::cout << "NEG OUTFLUX! (dz=" << dz << ")\n";
	  //              std::cout << "outflux: " << cn->getQsin() - dz*cn->getVArea()/dtmax
	  //                   << std::endl;
	  //              cn->TellAll();
	  //              dn->TellAll();
	  //           }
	}

      // Update time remaining
      dtg -= dtmax;

    } while( dtg>1e-6 ); // Keep going until we've used up the whole time intrvl

  //std::cout << "Leaving StreamErode()\n";
}


/***************************************************************************\
 **  tErosion::StreamErodeMulti
 **
 **  Eventually there should only be one streamerode (not one for multi
 **  and one for single size).  This work is still in progress so both
 **  remain for the time being.  This function calculates the erosion and
 **  deposition rate of each grain size assuming the system is transport
 **  limited.  If there is less than a critical amount of sediment on the
 **  bed, the bedrock is eroded and added to the load.
 **  ...NOW OBSOLETE?
\***************************************************************************/
void tErosion::StreamErodeMulti( double dtg, tStreamNet *strmNet, double time )
{
  double dt,
    dtmax;         // time increment: initialize to arbitrary large val
  double frac = 0.3; //fraction of time to zero slope
  double timegb; //time gone by - for layering time purposes
  tLNode * cn, *dn;
  // int nActNodes = meshPtr->getNodeList()->getActiveSize();
  tMesh< tLNode >::nodeListIter_t ni( meshPtr->getNodeList() );
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
  const tArray <double> sedzero( cn->getNumg() );

  // Sort so that we always work in upstream to downstream order
  strmNet->SortNodesByNetOrder();
  strmNet->FindChanGeom();
  strmNet->FindHydrGeom();

  // Compute erosion and/or deposition until all of the elapsed time (dtg)
  // is used up
  timegb=time;
  do
    {
      //std::cout<<"AT BEGINNING!"<<std::endl;
      // Zero out sed influx of all sizes
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() ){
	//if(cn->getID()==93)
	//  std::cout<<"93 is active"<<std::endl;
	cn->setQsin(0.0); //totals are for ts calculation
	cn->setQs(0.0);
	cn->setQsin( sedzero );
	for( size_t i=0; i<cn->getNumg(); i++ ){
	  cn->setQs(i,0.0);
	}
	//if(cn->getID()==93)
	//    std::cout<<"93 is active"<<std::endl;
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
	      cap = sedTrans->TransCapacity( cn );
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
		  dcap = -bedErode->DetachCapacity(cn)*(1-(cn->getLayerDepth(0)/cn->getMaxregdep()));
		  pedr += dcap;
		  // if(fabs(dcap)>1){
		  //                   std::cout << "huge bedrock erosion of "<< dcap << " at node " << cn->getID()<<std::endl;
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
	      dcap = -bedErode->DetachCapacity( cn );
	      // if(fabs(dcap)>1){
	      //                 std::cout << "huge bedrock erosion of "<< dcap <<" at node " << cn->getID()<< std::endl;
	      //              }

	      //if( dcap > pedr )
	      pedr = dcap;
	      //}
	    }

	  // if(cn->getX()==2.75 && cn->getY()==1){
	  //             std::cout<<cn->getNumLayer()<<" layers" <<std::endl;
	  //             n=0;
	  //             while(n<cn->getNumLayer()){
	  //                std::cout << "layer " << n+1 << std::endl;
	  //                std::cout << "layer creation time is " << cn->getLayerCtime(n) << std::endl;
	  //                std::cout << "layer recent time is " << cn->getLayerRtime(n) << std::endl;
	  //                std::cout << "layer depth is " << cn->getLayerDepth(n) << std::endl;
	  //                std::cout << "layer erodibility is " << cn->getLayerErody(n) << std::endl;
	  //                std::cout << "is layer sediment? " << cn->getLayerSed(n) << std::endl;
	  //                std::cout << "dgrade 1 is " << cn->getLayerDgrade(n,0) << std::endl;
	  //                std::cout << "dgrade 2 is " << cn->getLayerDgrade(n,1) << std::endl;
	  //                n++;
	  //             }
	  //             std::cout<<"qs0 is "<<cn->getQs(0)<<" qs1 is "<<cn->getQs(1)<<std::endl;
	  //             std::cout<<"texture of surface is "<<cn->getLayerSed(0)<<std::endl;
	  //          }

	  // Set the erosion (deposition) rate and send the corresponding
	  // sediment influx downstream
	  cn->setDzDt( pedr );
	  cn->getDownstrmNbr()->addQsin( cn->getQsin() - pedr*cn->getVArea() );
	  // only doing totals here
	  // sediment transport rates for each grn size have been calculated
	  //std::cout << "RATE STEP:\n";
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

		  std::cout << "Very small dt " << dt << " at:\n" << std::endl;
		  //std::cout << "rate dif is " << ratediff << std::endl;
		  //std::cout << "elev dif is " << cn->getZ() - dn->getZ() << std::endl;
		  //std::cout << "dzdt upstream is " << cn->getDzDt() << std::endl;
		  //std::cout << "dzdt downstream is " << dn->getDzDt() << std::endl;
		  // cn->TellAll();
		  //dn->TellAll();
		  //std::cout << "arbitrarily set dt to 0.0015" << std::endl;
		  dtmax=0.005;
		}
	    }
	}
      dtmax *= frac;  // Take a fraction of time-to-flattening

      //std::cout<<"dtmax is "<<dtmax<<std::endl;

      // Zero out sed influx again, because depending on bedrock-alluvial
      // interaction it may be modified; if inlet, give it strmNet->inlet.inSedLoad
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() ){
	cn->setQsin(0.0);
	cn->setQsin( sedzero );
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
	  dzt=0.;
	  //cn->TellAll();
	  for( size_t i=0; i<cn->getNumg(); i++ ){
            dz[i] = ( (cn->getQsin(i)-cn->getQs(i)) / cn->getVArea() ) * dtmax;
            dzt += dz[i];
            retbr[i]=0.;
            retsed[i]=0.;
	  }
	  //  if(cn->getX()==2.75 && cn->getY()==1){
	  //             std::cout<<"qsin0 qs0 "<<cn->getQsin(0)<<" "<<cn->getQs(0)<<std::endl;
	  //             std::cout<<"qsin1 qs1 "<<cn->getQsin(1)<<" "<<cn->getQs(1)<<std::endl;
	  //             std::cout<<"area "<<cn->getVArea()<<" time "<<dtmax<<std::endl;
	  //          }


	  // If we're on bedrock, scour the bedrock
	  dzrt = 0.;
	  if(cn->getLayerSed(0)<1){
            // Bedrock at surface
            for( size_t i=0; i<cn->getNumg(); i++ ){
	      dzr[i] = cn->getDrDt()*cn->getLayerDgrade(0,i)/cn->getLayerDepth(0)*dtmax;
	      dzrt += dzr[i];
	      // if(cn->getX()==2.75 && cn->getY()==1){
	      //                   std::cout<<"no sed drdt is "<<dzrt<<" of node "<<cn->getID()<<std::endl;
	      //                   std::cout<<"should erode some bedrock"<<std::endl;
	      //                }
            }

            if(dzrt<0)
	      retbr=cn->EroDep(0,dzr,timegb);
            //if(cn->getID()==93){
            // std::cout<<"done with erosion nic"<< std::endl;
            //}

	  }
	  else if( fabs(cn->getLayerDepth(0)-cn->getMaxregdep())>0.001 && dzt<0.0 && cn->getLayerSed(1)<1)
	    {
	      // Bedrock not at surface, but not enough sediment either
	      // nic this should work if only regolith and bedrock
	      // and also IF layering is working correctly
	      for( size_t i=0; i<cn->getNumg(); i++ ){
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
	      //                std::cout<<"sed on top drdt is "<<dzrt<<" of node "<<cn->getID()<<std::endl;
	      //                std::cout<<"should erode some bedrock"<<std::endl;
	      //             }

	      if(dzrt<0)
                retbr=cn->EroDep(1,dzr,timegb);
	    }

	  //std::cout << "** THIS node has dzs " << dzs << " & dzr " << dzr
	  //     << " & net change " << dzr+dzs << std::endl;
	  //cn->TellAll();

	  // Update alluvium thickness and node elevation
	  if(fabs(dzt)>0){
            // if(cn->getX()==2.75 && cn->getY()==1){
	    //                 std::cout<<"eroding sediment at specified node, total of "<<dzt<<std::endl;
	    //                 std::cout<<"Total number of layers is "<<cn->getNumLayer()<<std::endl;
	    //                 std::cout<<"dz0 is "<<dz[0]<<" dz1 is "<<dz[1]<<std::endl;

	    //             }

            //std::cout << "dzt is " << dzt << std::endl;
            retsed=cn->EroDep(0, dz, timegb );
	  }

	  dn = cn->getDownstrmNbr();

	  // Send sediment downstream: sediment flux is equal to the flux in
	  // plus/minus rate of erosion/deposition times node area
	  for(size_t i=0; i<cn->getNumg(); i++){
            dn->addQsin( i, cn->getQsin(i) - (retbr[i]+retsed[i])*cn->getVArea()/dtmax );
	  }
	  //cn->TellAll();
	}



      // Update time remaining
      dtg -= dtmax;
    } while( dtg>1e-6 ); // Keep going until we've used up the whole time intrvl

}


/***********************************************************************\
 **
 **  tErosion::DetachErode
 **
 **  Algorithm for eroding sediment and bedrock.  Material is only detached
 **  if the stream has the capacity to carry it. Handles multiple grain
 **  sizes. Replaces StreamErode and StreamErodeMulti.
 **
\************************************************************************/

void tErosion::DetachErode(double dtg, tStreamNet *strmNet, double time,
			   tVegetation * /*pVegetation*/ )
{

  //Added 4/00, if there is no runoff, this would crash, so check
  if(strmNet->getRainRate()-strmNet->getInfilt()>0){

    double dtmax;       // time increment: initialize to arbitrary large val
    double frac = 0.3;  //fraction of time to zero slope
    double timegb=time; //time gone by - for layering time purposes
    bool flag;
    tLNode * cn, *dn;
    // int nActNodes = meshPtr->getNodeList()->getActiveSize();
    tMesh< tLNode >::nodeListIter_t ni( meshPtr->getNodeList() );
    double ratediff,  // Difference in ero/dep rate btwn node & its downstrm nbr
      drdt,
      dz,
      depck,
      qs,
      excap;
    tLNode * inletNode = strmNet->getInletNodePtrNC();
    double insedloadtotal = strmNet->getInSedLoad();
	int debugCount = 0;

    cn = ni.FirstP();

    tArray <double> ret( cn->getNumg() ); //amt actually ero'd/dep'd
    tArray <double> erolist( cn->getNumg() );
    const tArray <double> sedzero( cn->getNumg() );
    tArray <double> insed( strmNet->getInSedLoadm() );

    // Sort so that we always work in upstream to downstream order
    strmNet->SortNodesByNetOrder();
    strmNet->FindChanGeom();
    strmNet->FindHydrGeom();


    // Compute erosion and/or deposition until all of the elapsed time (dtg)
    // is used up
    do
      {
	// Zero out sed influx of all sizes
	for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
	  {
            cn->setQs(0.0);
            if( cn!=inletNode )
	      {
		cn->setQsin(0.0); //totals are for ts calculation
		cn->setQsin( sedzero );
		for( size_t i=0; i<cn->getNumg(); i++ ){
                  cn->setQs(i,0.0);
		}
	      }
            else
	      {
		cn->setQsin(insedloadtotal); //totals are for ts calculation
		cn->setQsin( insed );
		for( size_t i=0; i<cn->getNumg(); i++ ){
                  cn->setQs(i,0.0);
		  //std::cout << "inlet qsin size " << i << "=" << insed[i] << std::endl;
		}
	      }
	  }

	// Estimate erosion rates and time-step size
	// NOTE - in this first loop we are only dealing with
	// totals for time-step calculations, however transport
	// rates for each size are also set within the function call.
	for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
	  {
            depck=0.;
            int i=0;
            qs=0.;

            assert(cn->getChanDepth()<1000);

            while((cn->getChanDepth()-depck)>0.0001)
	      {
		// Total transport capacity is a weighted average
		// of the transport capacity calculated from each
		// layer within the channel depth.
		// sediment and bedrock treated the same
		// units on qs are l^3/t
		if((depck+cn->getLayerDepth(i))<=cn->getChanDepth()){
                  //TransportCapacity function should keep running
                  //sum of qs of each grain size.
                  //qs returned is in m^3/yr; qs stored in tLNode has same units
                  qs+=sedTrans->TransCapacity(cn,i,cn->getLayerDepth(i)/cn->getChanDepth());
		}
		else{
                  qs+=sedTrans->TransCapacity(cn,i,1-(depck/cn->getChanDepth()));
		}
		depck+=cn->getLayerDepth(i); //need to keep this here for qs calc
		i++;
	      }

            //NIC this detachcapacity returns the correct thing, but
            //it also sets within the layer the drdt of each size.
            //You don't want to use detach capacity this way, so
            //I don't think that will affect anything, just be careful of
            //using those values!!!

            if(depck>cn->getChanDepth()) //which layer are you basing detach on?
	      drdt=-bedErode->DetachCapacity( cn, i-1 );
            else
	      drdt=-bedErode->DetachCapacity( cn, i );//[m^3/yr]

            cn->setDrDt(drdt);
            cn->setDzDt(drdt);

	    /*std::cout << "*** EROSION ***\n";
	      if( cn->getID()==342 ) {
	      std::cout << "**** Trans Cap 342 = " << qs << std::endl;
	      cn->TellAll();
	      }*/

            excap=(qs - cn->getQsin())/cn->getVArea();//[m/yr]
            //excap negative = deposition; positive = erosion
            //Note that signs are opposite to what one
            //might expect.  This works out for Qsin addition.
            //Limit erosion to capacity of flow or deposition
            if( -drdt > excap ){
	      cn->setDzDt(-excap);
            }
            cn->getDownstrmNbr()->addQsin(cn->getQsin()-cn->getDzDt()*cn->getVArea());
	  }//ends for( cn = ni.FirstP...

	//Find local time-step based on dzdt
	dtmax = dtg/frac;
	for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
	  {
            //Not for time step calculations, just utilizing loop
            if( cn!=inletNode )
	      {
		//Note-seting qsini for each size should automatically
		//properly set qstotal to be the sum of all qsini
		cn->setQsin( sedzero );
	      }
            else
	      {
		cn->setQsin( insed );
	      }

            dn = cn->getDownstrmNbr();
            ratediff = dn->getDzDt() - cn->getDzDt(); //Are the pts converging?
            if( ratediff > 0. && (cn->calcSlope()) > 1e-7 )  // if yes, get time
            {                                              //  to zero slope
               if(1) {
		 double dt;
		 dt = ( cn->getZ() - dn->getZ() ) / ratediff;
		 if( dt < dtmax ) dtmax = dt;
	       }
	       if( ratediff*dtmax > (cn->getZ() - dn->getZ() ) )
		 {
		 dtmax = ( cn->getZ() - dn->getZ() ) / ratediff;
		 assert( dtmax > 0.0 );
		 if( dtmax < 0.0001 && dtmax < dtg )
		   {
		     if(1) { // debug
		       std::cout << "Very small dtmax " << dtmax <<  std::endl;
		       std::cout << "rate dif is " << ratediff << std::endl;
		       std::cout << "elev dif is " << cn->getZ()-dn->getZ() << std::endl;
		       std::cout << "dzdt upstream is " << cn->getDzDt() << std::endl;
		       std::cout << "dzdt downstream is " << dn->getDzDt() << std::endl;
		       cn->TellAll();
		       dn->TellAll();
		     }
		     //dtmax=0.0001; //GREG I added this just because things
		     // were taking forever.  I kept it for now just for
		     // testing stuff.  Maybe we should discuss this.
		     //this was set to 0.0001 by NIC, I'm making it smaller
		     //because I'm having instability problems.. my storms are
		     // 0.0001! 
		   }
		 }
	    }
	  }// End for( cn = ni.FirstP()..
	dtmax *= frac;  // Take a fraction of time-to-flattening
	timegb+=dtmax;
         
         //At this point: we have drdt and qs for each node, plus dtmax
         
         // Do erosion/deposition
	for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
	  {
            //need to recalculate cause qsin may change due to time step calc
            excap=(cn->getQs() - cn->getQsin())/cn->getVArea();

            //std::cout<<"actual erosion excap = "<<excap<<std::endl;
            //std::cout<<"drdt is "<<cn->getDrDt()<<std::endl;
            //again, excap pos if eroding, neg if depositing
            //nic here is where drdt comes in again
            //flag is used to determine the texture of what should be eroded.
            //If detach limited, just erode what is there, but always limit
            //it by what flow has capacity to transport.  If transport limited,
            //the texture of what erode is determined by the calculated values
            //of qs.
            if( -cn->getDrDt() < excap ){
	      dz = cn->getDrDt()*dtmax; // detach-lim
	      flag = false;
            }
            else{
	      dz = -excap*dtmax; // trans-lim
	      flag = true;
            }

            for(size_t i=0; i<cn->getNumg(); i++)
	      cn->getDownstrmNbr()->addQsin(i,cn->getQsin(i));
            //What goes downstream will be what comes in + what gets ero'd/dep'd
            //This should always be negative or zero since max amt
            //to deposit is what goes in.
            //i.e. send (qsin[i]-ret[i]*varea/dtmax) downstream
            //Note: I think need to do the add in here and possibly take out later
            //because of looping through layers for the same erosion pass.

            /*DEBUG double l0, l1;
              if( cn->getX()>50.0 && cn->getX()<51.0
              && cn->getY()>29.0 && cn->getY()<30.0 )
              {
              std::cout << "f (" << cn->getID() << " ld = " << cn->getLayerDepth(0) << std::endl;
              l0 = cn->getLayerDgrade(0,0);
              l1 = cn->getLayerDgrade(0,1);
              }*/

            if( dz<0 ) //total erosion
	      {
		if(!flag)
		  { // detach-lim
		    int i=0;
		    depck=0.;
		    while(dz<-0.000000001&&depck<cn->getChanDepth()&&i<cn->getNumLayer())
		      {
			depck+=cn->getLayerDepth(i);
			if(-dz<=cn->getLayerDepth(i))
			  {//top layer can supply total depth
			    for(size_t j=0;j<cn->getNumg();j++)
			      {
				erolist[j]=dz*cn->getLayerDgrade(i,j)/cn->getLayerDepth(i);
				if(erolist[j]<(cn->getQsin(j)-cn->getQs(j))*dtmax/cn->getVArea())
				  {
				    //decrease total dz because of capacity limitations
				    erolist[j]=(cn->getQsin(j)-cn->getQs(j))*dtmax/cn->getVArea();
				    cn->setQsin(j,0.0);
				    cn->setQs(j,0.0);
				  }
			      }
			    ret=cn->EroDep(i,erolist,timegb);
			    for(size_t j=0;j<cn->getNumg();j++)
			      {
				cn->getDownstrmNbr()->addQsin(j,-ret[j]*cn->getVArea()/dtmax);
			      }
			    dz=0.;
			  }
			else
			  {//top layer is not deep enough, need to erode more layers
			    flag=false;
			    for(size_t j=0;j<cn->getNumg();j++)
			      {
				erolist[j]=-cn->getLayerDgrade(i,j);
				if(erolist[j]<(cn->getQsin(j)-cn->getQs(j))*dtmax/cn->getVArea())
				  {
				    //decrease total dz because of capacity limitations
				    erolist[j]=(cn->getQsin(j)-cn->getQs(j))*dtmax/cn->getVArea();
				    cn->setQsin(j,0.0);
				    cn->setQs(j,0.0);
				    //need to set these to zero since the capacity has
				    //now been filled by the stuff in this layer
				    flag=true;
				    //Since not taking all of the material from the
				    //surface, surface layer won't be removed-must inc i
				  }
				dz-=erolist[j];
			      }
			    ret=cn->EroDep(i,erolist,timegb);
			    for(size_t j=0;j<cn->getNumg();j++)
			      {
				//if * operator was overloaded for arrays, no loop necessary
				cn->getDownstrmNbr()->addQsin(j,-ret[j]*cn->getVArea()/dtmax);
			      }
			    if(flag)
			      {
				i++;
			      }
			  }
		      }
		  }
		else
		  {//trans-lim
		    //std::cout<<"X "<<cn->getX()<<" Y "<<cn->getY();
		    for(size_t j=0;j<cn->getNumg();j++){
		      erolist[j]=(cn->getQsin(j)-cn->getQs(j))*dtmax/cn->getVArea();
		      // std::cout<<" j "<<j<<" "<<erolist[j];
		    }
		    //std::cout<<std::endl;
		    
		    int i=0;
		    depck=0.;
		    while(depck<cn->getChanDepth())
		      {
			depck+=cn->getLayerDepth(i);
			int flag=cn->getNumLayer();
			ret=cn->EroDep(i,erolist,timegb);
			double sum=0.;
			for(size_t j=0;j<cn->getNumg();j++)
			  {
			    cn->getDownstrmNbr()->addQsin(j,-ret[j]*cn->getVArea()/dtmax);
			    erolist[j]-=ret[j];
			    sum+=erolist[j];
			  }
			if(sum>-0.0000001)
			  depck=cn->getChanDepth();
			if(flag==cn->getNumLayer())
			  i++;
		      }
		  }//end if( trans-limited )
	      }//ends(if dz<0)
            else if(dz>0) //total deposition -> need if cause erodep chokes with 0
	      {
		//Get texture of stuff to be deposited
		for(size_t j=0;j<cn->getNumg();j++)
		  erolist[j]=(cn->getQsin(j)-cn->getQs(j))*dtmax/cn->getVArea();
		ret=cn->EroDep(0,erolist,timegb);
		for(size_t j=0;j<cn->getNumg();j++)
		  {
		    cn->getDownstrmNbr()->addQsin(j,-ret[j]*cn->getVArea()/dtmax);
		  }
		
	      }

	  } // Ends for( cn = ni.FirstP()...

	 // Erode vegetation
#if 0
#define NEWVEG 0
	if( pVegetation && NEWVEG ) pVegetation->ErodeVegetation( meshPtr, dtmax );
#undef NEWVEG
#endif

	// Update time remainig
	dtg -= dtmax;
	
	if(1) //DEBUG
	  {
	    debugCount++;
	    if( debugCount > 1e6 )
	      ReportFatalError("More than 1e6 iterations in ErodeDetachLim()" );
	  }
	
	//std::cout<<"Time remaining now "<<dtg<<std::endl;
      } while( dtg>1e-6 );  //Keep going until we've used up the whole time intrvl
  }//end if rainrate-infilt>0


   //std::cout<<"ending detach erode"<<std::endl;
}// End erosion algorithm

/*****************************************************************************\
**
**  tErosion::WeatherNodes  
**
\*****************************************************************************/

void tErosion::WeatherNodes(double rt, double tt )
{
  tLNode * cn; 
  tMesh< tLNode >::nodeListIter_t ni( meshPtr->getNodeList() );
  double dz;

  for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )                
    {
      dz= weather.WeatheringRate(cn)*rt;  

      tArray<double> depth;
      depth.setSize(1);
      depth[0] =-dz;
      double newK;

      tArray<double> depthe;
      depthe.setSize(1);
      depthe[0] = 1.5*dz;

      if (cn->getNumLayer()>1)     
	{
	  // descending the weathering front to bedrock
	  cn->EroDep(cn->getNumLayer()-1,depth,tt);  
	  //producing soil to the soil store
	  cn->EroDep(cn->getNumLayer()-2,depthe,tt);   
	}
      else
	{    
	  cn->EroDep(cn->getNumLayer()-1,depth,tt); 
	  newK = weather.getK(); 
	  cn->makeNewLayerBelow(-1,tLayer::kSed,newK,depthe,tt);
	}
    } 
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
 **  Inputs:  rt -- time duration over which to compute diffusion
 **           noDepoFlag -- if true, material is only eroded, never
 **                             deposited
 **  Modifies:  node elevations (z); node Qsin (sed influx) is reset and
 **             used in the computations
 **  Notes:  as of 3/98, does not differentiate between rock and sediment
 **          as of 5/06, does differentiate for 1 regolith layer (hf)
\*****************************************************************************/
//#define kVerySmall 1e-6
#define kEpsOver2 0.1
void tErosion::Diffuse( double rt, bool noDepoFlag , double time)
{
  tLNode * cn;
  tEdge * ce;
  double volout,  // Sediment volume output from a node (neg=input)
    delt,       // Max local step size
    dtmax;      // Max global step size (initially equal to total time rt)
  tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );
  tMesh< tLNode >::edgeListIter_t edgIter( meshPtr->getEdgeList() );

  double Kd; //new diffusivity coeff. funct. of vegetation cover.

#ifdef TRACKFNS
  std::cout << "tErosion::Diffuse()" << std::endl;
#endif

  if( kd==0 ) return;

  // Compute maximum stable time-step size based on Courant condition
  // for FTCS (here used as an approximation).
  // (Note: for a fixed mesh, this calculation only needs to be done once;
  // performance could be improved by having this block only called if
  // mesh has changed since last time through)
  dtmax = rt;  // Initialize dtmax to total time rt
  for( ce=edgIter.FirstP(); edgIter.IsActive(); ce=edgIter.NextP() )
    {
      cn = static_cast<tLNode *>(ce->getOriginPtrNC());  //origin node
      Kd = kd*exp(-D_alpha*(cn->getVegCover().getVeg())); //new Kd(veg) [hf]
      if( 0 ) //DEBUG
	{
	  std::cout << "In Diffuse(), large vedglen detected: " << ce->getVEdgLen() << std::endl;
	  ce->TellCoords();
	}
      //Xif( (denom=Kd*ce->getVEdgLen() ) > kVerySmall )
      // Evaluate DT <= DX^2 / Kd
      assert( Kd > 0.0 );
      delt = kEpsOver2 * ce->getLength()*ce->getLength() / Kd;
      if( delt < dtmax )
	{
	  dtmax = delt;
	  if(0) { //DEBUG
	    std::cout << "TIME STEP CONSTRAINED TO " << dtmax << " AT EDGE:\n";
	    ce->TellCoords(); }
	}
    }


  // Loop until we've used up the entire time interval rt
  do
    {
      // Reset sed input for each node for the new iteration
      for( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
	cn->setQsin( 0. );

      // Compute sediment volume transfer along each edge
      for( ce=edgIter.FirstP(); edgIter.IsActive(); ce=edgIter.NextP() )
	{
	  //quick fix below (if) to prevent bedrock from diffusing -HF may 2005:still need to figure out what to do if run out of regolith...
	if (((tLNode *)ce->getOriginPtr())->getNumLayer() < 2)
	  {
	    volout = 0;
	  }
	else
	  {
	    cn = static_cast<tLNode *>(ce->getOriginPtrNC());  //origin node
	    Kd = kd*exp(-D_alpha*(cn->getVegCover().getVeg())); //new Kd(veg) [hf]
	  volout = Kd*ce->CalcSlope()*ce->getVEdgLen()*dtmax;
	  }
	  // Record outgoing flux from origin
	  cn = static_cast<tLNode *>(ce->getOriginPtrNC());
	  cn->addQsin( -volout );
	  // Record incoming flux to dest'n
	  cn = static_cast<tLNode *>(ce->getDestinationPtrNC());
	  cn->addQsin( volout );

	  if( 0 ) { //DEBUG
	    std::cout << volout << " mass exch. from " << ce->getOriginPtr()->getID()
		 << " to "
		 << ce->getDestinationPtr()->getID()
		 << " on slp " << ce->getSlope() << " ve " << ce->getVEdgLen()
		 << "\nvp " << ce->getRVtx().at(0)
		 << " " << ce->getRVtx().at(1) << std::endl;
	    static_cast<tLNode *>(ce->getOriginPtrNC())->TellAll();
	    static_cast<tLNode *>(ce->getDestinationPtrNC())->TellAll();
	    std::cout << std::endl;
	  }

	  /*ce =*/ edgIter.NextP();  // Skip complementary edge
	}

      // Compute erosion/deposition for each node
      for( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
	{
	  if( 0 ) //DEBUG
	    std::cout << "Node " << cn->getID() << " Qsin: " << cn->getQsin()
		 << " dz: " << cn->getQsin() / cn->getVArea() << std::endl;
	  if( noDepoFlag && cn->getQsin() > 0.0 )
	      cn->setQsin( 0.0 );
	  tArray<double> netflux;
	  netflux.setSize(1);
	  netflux[0] = cn->getQsin() / cn->getVArea();
	  if(-netflux[0] < cn->getLayerDepth(0)) //works only for single regolith layer
	    {//simple case where there is enough material to be eroded
	      cn->EroDep(0, netflux, time );  // add or subtract net flux/area
	    }
	  else
	    {//case where there's no enough material to be eroded
	      //netflux has to be negative, and larger in mag. than layer
	      //depth, so in the command line below we are substracting
	      //the difference from qsin of the downstream cell
	      //and the net flux is the layer depth.
	      cn->getDownstrmNbr()->addQsin( (netflux[0] + cn->getLayerDepth(0)));
	      netflux[0] = cn->getLayerDepth(0);
	      cn->EroDep(0, netflux, time ); 
	    }
	  if( 0 ) //DEBUG
	    std::cout<<cn->getZ()<<" Q: "<<cn->getQ()
		<<" dz "<<cn->getQsin() / cn->getVArea()
		<<" dt "<<dtmax<<std::endl;
	  /*if( cn->id==700 ) {
	    cn->TellAll();
	    }*/
	}

      rt -= dtmax;
      if( dtmax>rt ) dtmax=rt;

    } while( rt>0.0 );


}


/***********************************************************************\
 ** tErosion::UpdateExposureTime( double dtg)                         **
 **                                                                   **
 ** This function increments the exposure time of the top layer at    **
 ** every node by the amount dtg.                                     **
 ** Called from main loop.                                            **
 **                                                                   **
 ** created 3/1999 ng                                                 **
 **********************************************************************/
void tErosion::UpdateExposureTime( double dtg)
{
  tLNode * cn;
  tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );

#ifdef TRACKFNS
  std::cout << "tErosion::UpdateExposureTime()" << std::endl;
#endif

  for( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() ){
    cn->addLayerEtime(0, dtg);
  }
}


/***********************************************************************\
 **
 ** tErosion::DensifyMesh
 **
 ** Called only when the option for adaptive remeshing is invoked,
 ** this routine is used to increase mesh resolution in areas where
 ** erosion (or deposition) is especially rapid. For each point,
 ** a check is made to see whether the current erosion rate (ie, the
 ** most recently recorded value) times the Voronoi area is greater
 ** than a user-specified threshold, mdMeshAdaptMaxFlux. Note that
 ** the dimensions are L3/T, ie sediment flux being eroded from
 ** (or deposited into) the node; thus the threshold is a maximum
 ** allowable sediment flux resulting from local erosion. If the
 ** flux exceeds this threshold at a given node, new nodes are
 ** added at each of the node's Voronoi vertices.
 **
 **   Created: 2/2000 gt for gully erosion study
 **   Assumptions: assumes node dzdt value is correct
 **
 **********************************************************************/
void tErosion::DensifyMesh( double time )
{
  tMesh< tLNode >::nodeListIter_t niter( meshPtr->getNodeList() );  // node list iter.
  tLNode *cn;              // Current node being checked

  double dbgnf, dbgmax=0.;

  std::cout << "Checking nodes...\n";

  // Check all active nodes
  for( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
    {
      dbgnf = fabs(cn->getVArea()*cn->getDzDt());
      if( dbgnf>dbgmax ) dbgmax = dbgnf;

      // If local flux (ero rate * varea) exceeds threshold, add new nodes
      if( fabs(cn->getVArea()*cn->getDzDt()) > mdMeshAdaptMaxFlux )
	{
	  meshPtr->AddNodesAround( cn, time );
	  //std::cout << "*** Adding points here:\n";
	  //cn->TellAll();
	}
    }

  std::cout << "Max node flux: " << dbgmax << std::endl;

}


//#undef channel
