/***************************************************************************/
/**
 **  @file erosion.cpp
 **  @brief Functions for equilibrium checking, sediment transport,
 **          bed erosion (detachment), physical and chemical 
 **          weathering, and debris flow objects.
 **
 **    tEquilibCheck
 **  Transport objects:
 **    tSedTransPwrLaw
 **    tSedTransPwrLaw2
 **    tSedTransBridgeDom
 **    tSedTransWilcock
 **    tSedTransMineTailings //added 4/00 ng
 **    tSedTransPwrLawSimp
 **  Detachment objects:
 **    tBedErodePwrLaw
 **    tBedErodePwrLaw2
 **    tBedErodeAParabolic1
 **    tBedErodeGeneralFQS
 **  Physical weathering objects (added 6/2010 sl):
 **    tPhysicalWeatheringNone
 **    tPhysicalWeatheringExpLaw
 **    tPhysicalWeatheringDensityDependent
 **  Chemical weathering objects (added 6/2010 sl):
 **    tChemicalWeatheringNone
 **    tChemicalWeatheringDissolution
 **  Debris flow objects (added 9/2010 sl):
 **    tDebrisFlow
 **  Debris flow runout objects (added 9/2010 sl):
 **    tDF_RunOut
 **    tDF_RunOutNone
 **    tDF_RunOutNoStop
 **  Debris flow scour objects (added 9/2010 sl):
 **    tDF_Scour
 **    tDF_ScourNone
 **    tDF_ScourAllSediment
 **  Debris flow deposition objects (added 9/2010 sl):
 **    tDF_Deposit
 **    tDF_DepositNone
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
 **     - added nonlinear creep function DiffuseNonlinear (GT 8/07)
 **     - Weathering: Added chemical and physical weathering; choice of 
 **       laws, determined at run time as with tBedErode and tSedTrans
 **       (SL, 7/10)
 **     - Added "dummy" tBedErode and tSedTrans versions to provide a 
 **       way to turn off fluvial erosion and transport completely (SL 7/10)
 **     - Nonlinear depth-dependent supply-limited diffusion (SL 8/10)
 **     - Landsliding of node clusters and debris flows; choice of rules
 **       governing debris flow runout, scour, and deposition are chosen
 **       at run time as with tBedErode and tSedTrans, etc. (SL 9/10)
 **
 **    Known bugs:
 **     - ErodeDetachLim assumes 1 grain size. If multiple grain sizes
 **       are specified in the input file and the detachment limited
 **       option is used, a crash will result when tLNode::EroDep
 **       attempts to access array indices above 1. TODO (GT 3/00)
 **
 **  $Id: erosion.cpp,v 1.144 2007-08-21 00:13:46 childcvs Exp $
 */
/***************************************************************************/

#include <math.h>
#include <assert.h>
# include <iomanip>
#include <vector>  // first added for DiffuseNonlinear()
#include <queue> // first added for Landslides()
using namespace std;   // also added for DiffuseNonlinear() to use vector class from STL
//#include <string>
#include "erosion.h"

// Here follows a table for transport, detachment, and physical and chemical
// weathering laws, which are chosen at run time via "X()" trick in 
// tErosion constructor.
//
// ("X()" trick exposed in:
// The New C: X Macros, Randy Meyers, C/C++ Users Journal,
// 19(5), May 2001)

// Transport laws:
#define TRANSPORT_LAW_TABLE \
X(PowerLaw1,"Power-law transport formula"), \
X(PowerLaw2,"Power-law transport formula, form 2"), \
X(BridgeDominic,"Bridge-Dominic form of Bagnold bedload formula"), \
X(Wilcock,"Wilcock sand-gravel formula"), \
X(PowerLawMulti,"Multi-size power-law formula"), \
X(MineTailings,"Willgoose/Riley mine tailings formula"), \
X(PowerLaw3, "Ultra-Simplified power-law transport formula"), \
X(NoSedTrans, "Dummy law for no fluvial transport")

#define TRANSPORT_LAW_TABLE2 \
X(PowerLaw1,tSedTransPwrLaw) \
X(PowerLaw2,tSedTransPwrLaw2) \
X(BridgeDominic,tSedTransBridgeDom) \
X(Wilcock,tSedTransWilcock) \
X(PowerLawMulti,tSedTransPwrLawMulti) \
X(MineTailings,tSedTransMineTailings) \
X(PowerLaw3,tSedTransPwrLawSimp) \
X(NoSedTrans,tSedTransNone)

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

// Detachment laws:
#define DETACHMENT_LAW_TABLE \
X(DetachPwrLaw1,"Power law, form 1"), \
X(DetachPwrLaw2,"Power law, form 2"), \
X(DetachAParabolic1, "Almost Parabolic Law"), \
X(DetachGeneralFQS,"Generalized f(Qs) Detachment-rule"), \
X(DetachNone, "Dummy law for no fluvial erosion")

#define DETACHMENT_LAW_TABLE2 \
X(DetachPwrLaw1,tBedErodePwrLaw) \
X(DetachPwrLaw2,tBedErodePwrLaw2) \
X(DetachAParabolic1, tBedErodeAParabolic1) \
X(DetachGeneralFQS,tBedErodeGeneralFQS) \
X(DetachNone,tBedErodeNone)

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

// Physical weathering (production) laws:
#define PRODUCTION_LAW_TABLE \
X(NoPhysWeath,"No soil production"), \
X(ExpLaw,"Simple exponential production law"), \
X(DensityDependent,"Density-dependent exponential production law")

#define PRODUCTION_LAW_TABLE2 \
X(NoPhysWeath,tPhysicalWeatheringNone) \
X(ExpLaw,tPhysicalWeatheringExpLaw) \
X(DensityDependent,tPhysicalWeatheringDensityDependent)

#define X(a,b) a
enum {
  PRODUCTION_LAW_TABLE
};
#undef X

#define X(a,b) b
char const * const ProductionLaw[] =
{
  PRODUCTION_LAW_TABLE
};
#undef X

const int NUMBER_OF_PRODUCTION_LAWS =
sizeof(ProductionLaw)/sizeof(ProductionLaw[0]);

// Chemical weathering:
#define CHEM_WEATHERING_TABLE \
X(NoChemWeath,"No chemical weathering"), \
X(Dissolution,"Simple dissolution law")

#define CHEM_WEATHERING_TABLE2 \
X(NoChemWeath,tChemicalWeatheringNone) \
X(Dissolution,tChemicalWeatheringDissolution)

#define X(a,b) a
enum {
  CHEM_WEATHERING_TABLE
};
#undef X

#define X(a,b) b
char const * const ChemWeathering[] =
{
  CHEM_WEATHERING_TABLE
};
#undef X

const int NUMBER_OF_CHEM_WEATHERINGS =
sizeof(ChemWeathering)/sizeof(ChemWeathering[0]);

// Debris flow runout:
#define DF_RUNOUT_TABLE			 \
X(NoDF_RunOut,"No debris flow runout"), \
X(NoDF_Stop,"Runout exits domain")

#define DF_RUNOUT_TABLE2 \
X(NoDF_RunOut,tDF_RunOutNone) \
X(NoDF_Stop,tDF_RunOutNoStop)

#define X(a,b) a
enum {
  DF_RUNOUT_TABLE
};
#undef X

#define X(a,b) b
char const * const DebrisFlowRunOut[] =
{
  DF_RUNOUT_TABLE
};
#undef X

const int NUMBER_OF_DF_RUNOUTS =
sizeof(DebrisFlowRunOut)/sizeof(DebrisFlowRunOut[0]);

// Debris flow scour:
#define DF_SCOUR_TABLE			 \
X(NoDF_Scour,"No debris flow scour"), \
X(AllSediment,"Scour all sediment")

#define DF_SCOUR_TABLE2 \
X(NoDF_Scour,tDF_ScourNone) \
X(AllSediment,tDF_ScourAllSediment)

#define X(a,b) a
enum {
  DF_SCOUR_TABLE
};
#undef X

#define X(a,b) b
char const * const DebrisFlowScour[] =
{
  DF_SCOUR_TABLE
};
#undef X

const int NUMBER_OF_DF_SCOURS =
sizeof(DebrisFlowScour)/sizeof(DebrisFlowScour[0]);

// Debris flow deposition:
#define DF_DEPOSITION_TABLE	      \
X(NoDF_Deposition,"No debris flow deposition")

#define DF_DEPOSITION_TABLE2 \
X(NoDF_Deposition,tDF_DepositNone)

#define X(a,b) a
enum {
  DF_DEPOSITION_TABLE
};
#undef X

#define X(a,b) b
char const * const DebrisFlowDeposit[] =
{
  DF_DEPOSITION_TABLE
};
#undef X

const int NUMBER_OF_DF_DEPOSITIONS =
sizeof(DebrisFlowDeposit)/sizeof(DebrisFlowDeposit[0]);

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
 **   - Changed critical shear stress so that one now sets different 
 **     thresholds for bedrock and regolith. This involved changes here
 **     and in tLNode.h/.cpp. The bed erosion functions read TAUCB (bedrock)
 **     and the sediment transport functions read TAUCR (regolith).
 **     GT, apr 07
 **
 \***************************************************************************/
//constructor: reads and sets the parameters
tBedErodePwrLaw::tBedErodePwrLaw( const tInputFile &infile )
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
  taucd = infile.ReadItem( taucd, "TAUCB" ); // tag changed from TAUCD to TAUCB apr 07 gt
  
  // Add unit conversion factor for kt -- this is required to convert
  // the quantity (Q/W)^mb from units of years to units of seconds.
  kt *= pow( secPerYear, -mb );
}

void tBedErodePwrLaw::Initialize_Copy( tBedErode* oPtr )
{
  tBedErodePwrLaw* ptr = static_cast<tBedErodePwrLaw*>( oPtr );
  kb = ptr->kb;
  kt = ptr->kt;
  mb = ptr->mb; // Specific q exponent
  nb = ptr->nb;
  pb = ptr->pb;
  taucd = ptr->taucd; // tag changed from TAUCD to TAUCB apr 07 gt
}

//copy constructor: 
tBedErodePwrLaw::tBedErodePwrLaw( const tBedErodePwrLaw &orig ) 
  : tBedErode()
{
  kb = orig.kb;
  kt = orig.kt;
  mb = orig.mb; // Specific q exponent
  nb = orig.nb;
  pb = orig.pb;
  taucd = orig.taucd; // tag changed from TAUCD to TAUCB apr 07 gt
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
  const double slp = n->calcSlope();
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw::DetachCapacity(tLNode*,double)");
  const double tau =  kt*pow( n->getQ() / n->getHydrWidth(), mb ) * pow( slp, nb );
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
  const double slp = n->calcSlope();
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw::DetachCapacity(tLNode*)");
  const double tau = kt*pow( n->getQ() / n->getHydrWidth(), mb )*pow( slp, nb );
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
  const double slp = n->calcSlope();
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw::DetachCapacity(tLNode*)");
  const double tau = kt*pow( n->getQ() / n->getHydrWidth(), mb )*pow( slp, nb );
  n->setTau( tau );
  double erorate = tau - n->getTauCrit();
  if(0) //DEBUG
    std::cout << "erorate: " << erorate << std::endl;
  erorate = (erorate>0.0) ? erorate : 0.0;
  erorate = n->getLayerErody(i)*pow( erorate, pb );
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
  taucd = infile.ReadItem( taucd, "TAUCB" );
  
  // Add unit conversion factor for kt -- this is required to convert
  // the quantity (Q/W)^mb from units of years to units of seconds.
  kt *= pow( secPerYear, -mb );
}

void tBedErodePwrLaw2::Initialize_Copy( tBedErode* oPtr )
{
  tBedErodePwrLaw2* ptr = static_cast<tBedErodePwrLaw2*>( oPtr );
  kb = ptr->kb;
  kt = ptr->kt;
  mb = ptr->mb; // Specific q exponent
  nb = ptr->nb;
  pb = ptr->pb;
  taucd = ptr->taucd; // tag changed from TAUCD to TAUCB apr 07 gt
}

//copy constructor: 
tBedErodePwrLaw2::tBedErodePwrLaw2( const tBedErodePwrLaw2 &orig )
  : tBedErode()
{
  kb = orig.kb;
  kt = orig.kt;
  mb = orig.mb; // Specific q exponent
  nb = orig.nb;
  pb = orig.pb;
  taucd = orig.taucd;
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
 **  FUNCTIONS FOR CLASS tBedErodeAParabolic1
 \***************************************************************************/

/***************************************************************************\
 **  tBedErodeAParabolic1 Constructor
 **
 \***************************************************************************/
//constructor: reads and sets the parameters
tBedErodeAParabolic1::tBedErodeAParabolic1( const tInputFile &infile )
{
  
  kb = infile.ReadItem(kb, "KB" );  //used for time-step, but get from
  //layer in calculations
  mb = infile.ReadItem( mb, "MB" ); // Q Exponent
  nb = infile.ReadItem( nb, "NB" ); // Slope Exponent
  beta = infile.ReadItem( beta, "BETA"); // fraction of sediment to bedload
  
}

void tBedErodeAParabolic1::Initialize_Copy( tBedErode* oPtr )
{
  tBedErodeAParabolic1* ptr = static_cast<tBedErodeAParabolic1*>( oPtr );
  kb = ptr->kb;
  mb = ptr->mb; // Specific q exponent
  nb = ptr->nb;
  beta = ptr->beta; // fraction of sediment to bedload
}

//copy constructor: 
tBedErodeAParabolic1::tBedErodeAParabolic1( const tBedErodeAParabolic1 &orig )
  : tBedErode()
{
  kb = orig.kb;  //used for time-step, but get from
                  //layer in calculations
  mb = orig.mb; // Q Exponent
  nb = orig.nb; // Slope Exponent
  beta = orig.beta; // fraction of sediment to bedload
}


/***************************************************************************\
 **  tBedErodeAParabolic1::DetachCapacity (1 of 3)
 **
 **  Computes the depth of erosion over a time interval dt assuming the
 **  erosion rate = kb f(Qs) Q^mb S^nb
 **      f(Qs) = 1 - 4(Qsin/Qc - 0.5)^2 ; Qsin/Qc >0.1
 **      f(Qs) = 2.6(Qsin/Qc + 0.1) ; Qsin/Qc <0.1
 **      f(Qs) = 99999999.9; Qsin/Qc >1 
 **              This last case forces model to go transport limited when
 **              there is more sediment than you can transport
 **
 **  Input: n -- node at which to compute detachment capacity
 **         dt -- time interval
 **  Returns: the detachment depth
 **  Assumptions: n->calcSlope() does not return a negative value (returns neg.
 **               only if infinite loop in calcSlope()); kb, mb,
 **               and nb all >=0.
 **               - Qsin is properly accounted for elsewhere.
 **               - Qc is already calculated elsewhere.
 **  Modifications:
 **   - replaced uniform erodibility coefficient kb with erodibility of
 **     topmost layer at node n (GT 8/98)
 **   - this functions now includes the diffusive load as part of the incoming
 **     load 
 \***************************************************************************/
double tBedErodeAParabolic1::DetachCapacity( tLNode * n, double dt )
{
  if( n->getFloodStatus() != tLNode::kNotFlooded) return 0.0;
  const double slp = n->calcSlope();
  double fqs, ratio, tau;
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodeAParabolic12::DetachCapacity(tLNode*,double)");
  
  if( slp == 0.0 ){
    n->setDrDt(0.0);
    return(0.0);
  }
  
  tau=pow(n->getQ(),mb)*pow(slp,nb);
  ratio = beta * (n->getQsin()+n->getQsdin())/ n->getQs();
  
  if(ratio <0.1){
    fqs=2.6*ratio+0.1;
    return(n->getLayerErody(0)*fqs*tau*dt );
  }
  else if(ratio < 1.0){
    fqs=1-(4*pow(ratio-0.5,2));
    return(n->getLayerErody(0)*fqs*tau*dt );
  }
  else{
    return(99999999.9);  //Force transport-limited case
  }
}

/***************************************************************************\
 **  tBedErodeAParabolic1::DetachCapacity (2 of 3)
 **
 **  Computes the rate of erosion as stated above - ONLY DIFFERENCE
 **  FROM PREVIOUS IS THAT A RATE IS RETURNED
 **  erosion rate = kb f(Qs) Q^mb S^nb
 **      f(Qs) = 1 - 4(Qsin/Qc - 0.5)^2 ; Qsin/Qc >0.1
 **      f(Qs) = 2.6(Qsin/Qc + 0.1) ; Qsin/Qc <0.1
 **      f(Qs) = 999999999.9; Qsin/Qc >1 
 **              This last case forces model to go transport limited when
 **              there is more sediment than you can transport
 **
 **  Input: n -- node at which to compute detachment capacity
 **         dt -- time interval
 **  Returns: the detachment depth
 **  Assumptions: n->calcSlope() does not return a negative value (returns neg.
 **               only if infinite loop in calcSlope()); kb, mb,
 **               and nb all >=0.
 **               - Qsin is properly accounted for elsewhere.
 **               - Qc is already calculated elsewhere.
 **  Modifications:
 **   - replaced uniform erodibility coefficient kb with erodibility of
 **     topmost layer at node n (GT 8/98)
 **
 **  Note:
 **    - This version of the code only accounts for sediment brought in
 **      through fluvial processes.  May want to change this later.
 **   - this functions now includes the diffusive load as part of the incoming
 **     load 
 **
 \***************************************************************************/
double tBedErodeAParabolic1::DetachCapacity( tLNode * n )
{
  if( n->getFloodStatus() != tLNode::kNotFlooded) return 0.0;
  const double slp = n->calcSlope();
  double fqs, ratio, tau, erorate;
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodeAParabolic12::DetachCapacity(tLNode*,double)");
  
  if( slp == 0.0 ){
    n->setDrDt(0.0);
    return(0.0);
  }
  
  tau=pow(n->getQ(),mb)*pow(slp,nb);
  ratio = (beta * (n->getQsin()+n->getQsdin()))/ n->getQs();
  
  if(ratio <0.1){
    //linear beginning on fqs to account for upper boundary condition
    fqs=2.6*ratio+0.1;
    erorate=n->getLayerErody(0)*fqs*tau;
    n->setDrDt( -erorate );
    return(erorate );
  }
  else if(ratio < 1.0){
    fqs=1-(4*pow(ratio-0.5,2));
    erorate=n->getLayerErody(0)*fqs*tau;
    n->setDrDt( -erorate );
    return(erorate );
  }
  else{
    n->setDrDt(-99999999.9);
    return(99999999.9);  //Force transport-limited case
  }
}

/***************************************************************************\
 **  tBedErodeAParabolic1::DetachCapacity (3 of 3)
 **
 **  Computes the rate of erosion as above.  The only difference is that
 **  here erodibility of layer is used as the coefficient for detach capacity
 **
 \***************************************************************************/
double tBedErodeAParabolic1::DetachCapacity( tLNode * n, int i )
{
  //cout<<"Parabolic DetachCapacity #3"<<endl;
  if( n->getFloodStatus() != tLNode::kNotFlooded) return 0.0;
  const double slp = n->calcSlope();
  double fqs, ratio, tau, erorate;
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodeAParabolic12::DetachCapacity(tLNode*,double)");
  
  if( slp == 0.0 ){
    n->setDrDt(0.0);
    return(0.0);
  }
  
  tau=pow(n->getQ(),mb)*pow(slp,nb);
  ratio = (beta * (n->getQsin()+n->getQsdin()))/ n->getQs();
  
  if(ratio <0.1){
    //linear beginning to take care of boundary condition when sed influx is zero
    fqs=2.6*ratio+0.1;
    erorate=n->getLayerErody(i)*fqs*tau;
    n->setDrDt( -erorate );
    //cout<<" fqs "<<fqs<<" erorate "<<erorate<<endl;
    return(erorate );
  }  
  else if(ratio < 1.0){
    fqs=1-(4*pow(ratio-0.5,2));
    erorate=n->getLayerErody(i)*fqs*tau;
    n->setDrDt( -erorate );
    //cout<<" fqs "<<fqs<<" erorate "<<erorate<<endl;
    return(erorate );
  }
  else{
    n->setDrDt(-99999999.9);
    //cout<<" erosion rate =  99999999"<<endl;
    return(99999999.9);  //Force transport-limited case
  }
}

/***************************************************************************\
 **  tBedErodeAParabolic1::SetTimeStep
 **
 **  NOTE - This is based on Streampower model, but don't have another 
 **         method right now.
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
double tBedErodeAParabolic1::SetTimeStep( tLNode * n )
{
  const double slp = n->calcSlope();
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodeAParabolic1::setTimeStep(tLNode*)");
  assert( n->getQ()>=0 );
  double eroterm = kb * pow( n->getQ(), mb ) * pow( slp, nb-1.0 );
  if( eroterm==0 ) return 100000;
  return( 0.2*n->getFlowEdg()->getLength() / eroterm );
  
}

/**************************************************************************\
 ** FUNCTIONS FOR CLASS tBedErodeGeneralFQS
 \*************************************************************************/

/***************************************************************************\
 **  tBedErodeGeneralFQS Constructor
 **
 **  Reads in coefficients and exponents for a general wear model
 **
 **  erosion rate = K B Qs/W (1 - Qs/(Qc)) A^m S^n
 ****************************************************************************/
//constructor
tBedErodeGeneralFQS::tBedErodeGeneralFQS( const tInputFile &infile )
{
  m = infile.ReadItem( m, "MB" ); 
  n = infile.ReadItem( n, "NB" );
  K = infile.ReadItem( K, "KB" );
  beta = infile.ReadItem( beta, "BETA" );
  
}

void tBedErodeGeneralFQS::Initialize_Copy( tBedErode* oPtr )
{
  tBedErodeGeneralFQS* ptr = static_cast<tBedErodeGeneralFQS*>( oPtr );
  m = ptr->m; 
  n = ptr->n;
  K = ptr->K;
  beta = ptr->beta;
}

//copy constructor
tBedErodeGeneralFQS::tBedErodeGeneralFQS( const tBedErodeGeneralFQS &orig )
  : tBedErode()
{
  m = orig.m; 
  n = orig.n;
  K = orig.K;
  beta = orig.beta;
}

/***************************************************************************\
 **  tBedErodeGeneralFQS::DetachCapacity (1 of 3)
 **
 **  Computes the depth of erosion over a time interval dt assuming the
 **  erosion rate =  K Qs/W (1 - Qs/Qc) A^m S^n
 **  erosion rate = infinte if Qs>Qc => this will cause the model to
 **  use the transport-limited solution
 **
 **  Input: nd -- node at which to compute detachment capacity
 **         dt -- time interval
 **  Returns: the detachment depth [L]
 **  Assumptions: assumes that Qsin and Qc and W
 **               calculated elsewhere. 
 **               - Diffusion should be calculated before erosion is called
 **               because Qsdin included in sediment load.
 **               - n->getSlope() does not return a negative value (returns neg.
 **               only if infinite loop in getSlope())
 **
 \***************************************************************************/
double tBedErodeGeneralFQS::DetachCapacity( tLNode * nd, double dt )
{
  double slp = sin(atan(nd->calcSlope())), Qs, Qc, W, erorate;
  // calcSlope() returns tan alpha, use sin alpha instead.
  
  if( nd->getFloodStatus() != tLNode::kNotFlooded) return 0.0;
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodeGeneralFQS::DetachCapacity(tLNode*,double)");
  Qs= beta * (nd->getQsin()+nd->getQsdin()); //sediment load
  Qc= nd->getQs(); //sediment transport rate
  W= nd->getHydrWidth(); //channel width
  
  if (Qs <= 0){
    nd->setDrDt( 0.0);
    return(0.0);
  }
  else if (Qs<Qc){
    erorate = K*(Qs/W)*(1-(Qs/Qc))*pow(nd->getQ(),m)*pow(slp,n);
    nd->setDrDt( -erorate );
    return( erorate*dt );
  }
  else{ //Qs>Qc
    nd->setDrDt( -999999999);
    return(999999999);
    //Seems strange, but force to go transport-limited when Qc>Qs or
    //Qs = 0 (this is just a boundary condition).
    //Can force transport-limited by making detachment capacity large.
  }
  
}

/***************************************************************************\
 **  tBedErodeGeneralFQS::DetachCapacity (2 of 3)
 **
 **  Computes the rate of erosion= 
 **  erosion rate =  K Qs/W (1 - Qs/Qc) A^m S^n
 **  erosion rate = infinte if Qs>Qc => this will cause the model to
 **  use the transport-limited solution
 **
 **  Input: nd -- node at which to compute detachment capacity
 **
 **  Returns: the detachment rate [L/T]
 **  Assumptions: assumes that Qsin, Qc, and W are 
 **               calculated elsewhere.
 **               n->getSlope() does not return a negative value (returns neg.
 **               only if infinite loop in getSlope()); 
 ** 
 \***************************************************************************/
double tBedErodeGeneralFQS::DetachCapacity( tLNode * nd )
{
  double slp = sin(atan(nd->calcSlope())), Qs, Qc, W, erorate;
  // calcSlope() returns tan alpha, use sin alpha instead.
  
  if( nd->getFloodStatus() != tLNode::kNotFlooded) return 0.0;
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodeGeneralFQS::DetachCapacity(tLNode*,double)");
  Qs= beta * (nd->getQsin()+nd->getQsdin()); //sediment load
  Qc= nd->getQs(); //Transport rate
  W= nd->getHydrWidth();
  
  if (Qs <= 0){
    nd->setDrDt( 0.0);
    return(0.0);
  }
  else if (Qs<Qc){
    erorate = K*(Qs/W)*(1-(Qs/Qc))*pow(nd->getQ(),m)*pow(slp,n);
    nd->setDrDt( -erorate );
    return( erorate );
  }
  else{ //Qs>Qc
    nd->setDrDt( -999999999);
    return(999999999);
    //Seems strange, but force to go transport-limited when Qc>Qs or
    //Qs = 0 (this is just a boundary condition).
    //Can force transport-limited by making detachment capacity large.
  }   
  
}

/***************************************************************************\
 **  tBedErodeGeneralFQS::DetachCapacity (3 of 3)
 **
 **  Computes the rate of erosion= 
 **  erosion rate = K Qs/W (1 - Qs/Qc) A^m S^n
 **  erosion rate = infinte if Qs>Qc => this will cause the model to
 **  use the transport-limited solution
 **
 **  Input: nd -- node at which to compute detachment capacity
 **         i -- layer which you get erodibility of - not working at this point
 **
 **  Returns: the detachment rate [L/T]
 **  Assumptions: assumes that Qsin and Qc calculated elsewhere.
 **               n->getSlope() does not return a negative value (returns neg.
 **               only if infinite loop in getSlope())
 \***************************************************************************/
double tBedErodeGeneralFQS::DetachCapacity( tLNode * nd, int i )
{
  double slp = sin(atan(nd->calcSlope())), Qs, Qc, W, K, erorate;
  // calcSlope() returns tan alpha, use sin alpha instead.
  
  if( nd->getFloodStatus() != tLNode::kNotFlooded) return 0.0;
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodeGeneralFQS::DetachCapacity(tLNode*,double)");
  
  Qs= beta * (nd->getQsin()+nd->getQsdin()); //sediment load
  Qc= nd->getQs(); //transport rate
  W= nd->getHydrWidth();
  K=nd->getLayerErody(i);
  
  if (Qs <= 0){
    nd->setDrDt( 0.0);
    return(0.0);
  }
  else if (Qs<Qc){
    erorate = K*(Qs/W)*(1-(Qs/Qc))*pow(nd->getQ(),m)*pow(slp,n);
    nd->setDrDt( -erorate );
    return( erorate );
  }
  else{ //Qs>Qc
    nd->setDrDt( -999999999);
    return(999999999);
    //Seems strange, but force to go transport-limited when Qc>Qs or
    //Qs = 0 (this is just a boundary condition).
    //Can force transport-limited by making detachment capacity large.
  }   
  
}

/********************************************************************\
 **  THIS FUNCTION IS JUST A DUMMY.
 **  NOT SURE HOW TO SET TIME STEP.
 \********************************************************************/
double tBedErodeGeneralFQS::SetTimeStep( tLNode * nd )
{
  const double slp = nd->calcSlope();
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw::setTimeStep(tLNode*)");
  assert( nd->getQ()>=0 );
  double eroterm = K * pow( nd->getQ(), m ) * pow( slp, n-1.0 );
  if( eroterm==0 ) return 100000;
  return( 0.2*nd->getFlowEdg()->getLength() / eroterm );
  
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
  mf = infile.ReadItem( mf, "MF" );
  nf = infile.ReadItem( nf, "NF" );
  pf = infile.ReadItem( pf, "PF" );
  tauc = infile.ReadItem( tauc, "TAUCR" );  // to TAUCR apr 07 gt
  
  // Add unit conversion factor for kt -- this is required to convert
  // the quantity (Q/W)^mb from units of years to units of seconds.
  kt *= pow( secPerYear, -mf );
}

void tSedTransPwrLaw::Initialize_Copy( tSedTrans* oPtr )
{
  tSedTransPwrLaw *ptr = static_cast<tSedTransPwrLaw*>( oPtr );
  kf = ptr->kf;
  kt = ptr->kt;
  mf = ptr->mf;
  nf = ptr->nf;
  pf = ptr->pf;
  tauc = ptr->tauc;  // to TAUCR apr 07 gt
}

tSedTransPwrLaw::tSedTransPwrLaw( const tSedTransPwrLaw &orig )
  : tSedTrans()
{
  kf = orig.kf;
  kt = orig.kt;
  mf = orig.mf;
  nf = orig.nf;
  pf = orig.pf;
  tauc = orig.tauc;  // to TAUCR apr 07 gt
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
  const double slp = node->calcSlope();
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw::TransCapacity(tLNode*)");
  double tau, tauex, cap = 0.;
  if( node->getFloodStatus() == tLNode::kNotFlooded )
  {
    tau = kt * pow( node->getQ()/node->getHydrWidth(), mf ) * pow( slp, nf );
    node->setTau( tau );
    if(0) //DEBUG
      std::cout << "kt=" << kt << " Q=" << node->getQ() << " W="
      << node->getHydrWidth() << " S=" << node->calcSlope() << std::endl;
    tauex = tau - tauc;
    tauex = (tauex>0.0) ? tauex : 0.0;
    cap = kf * node->getHydrWidth() * pow( tauex, pf );
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
  const double slp = node->calcSlope();
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tSedTransPwrLaw::TransCapacity(tLNode*)");
  double tau, tauex, cap = 0.;
  
  if( node->getFloodStatus() == tLNode::kNotFlooded )
  {
    tau = kt * pow( node->getQ()/node->getHydrWidth(), mf ) * pow( slp, nf );
    node->setTau( tau );
    if(0) //DEBUG
      std::cout << "kt=" << kt << " Q=" << node->getQ() << " W="
      << node->getHydrWidth() << " S=" << node->calcSlope()
      << " tau=" << tau << std::endl;
    tauex = tau - tauc;
    tauex = (tauex>0.0) ? tauex : 0.0;
    cap = weight * kf * node->getHydrWidth() * pow( tauex, pf );
  }
  for(size_t i=0; i<node->getNumg(); i++)
    node->addQs(i, cap*node->getLayerDgrade(lyr,i)/node->getLayerDepth(lyr));
  
  node->setQs( cap );
  return cap;
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
  tauc = infile.ReadItem( tauc, "TAUCR" );  // to TAUCR apr 07 gt
  
  // Add unit conversion factor for kt -- this is required to convert
  // the quantity (Q/W)^mb from units of years to units of seconds.
  kt *= pow( secPerYear, -mf );
}

void tSedTransPwrLaw2::Initialize_Copy( tSedTrans* oPtr )
{
  tSedTransPwrLaw2 *ptr = static_cast<tSedTransPwrLaw2*>( oPtr );
  kf = ptr->kf;
  kt = ptr->kt;
  mf = ptr->mf;
  nf = ptr->nf;
  pf = ptr->pf;
  tauc = ptr->tauc;  // to TAUCR apr 07 gt
}

tSedTransPwrLaw2::tSedTransPwrLaw2( const tSedTransPwrLaw2 &orig )
  : tSedTrans()
{
  kf = orig.kf;
  kt = orig.kt;
  mf = orig.mf;
  nf = orig.nf;
  pf = orig.pf;
  tauc = orig.tauc;  // to TAUCR apr 07 gt
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
 **  FUNCTIONS FOR CLASS tSedTransPwrLawSimp
 \***************************************************************************/

/***************************************************************************\
 **  tSedTransPwrLawSimp Constructor:
 **
 **    Given a tInputFile as an argument, will read relevant parameters from
 **  the input file.
 **
 **
 \***************************************************************************/
tSedTransPwrLawSimp::tSedTransPwrLawSimp( const tInputFile &infile )
{
  
  kf = infile.ReadItem( kf, "KF" );
  mf = infile.ReadItem( mf, "MF" );
  nf = infile.ReadItem( nf, "NF" );
  
}

void tSedTransPwrLawSimp::Initialize_Copy( tSedTrans* oPtr )
{
  tSedTransPwrLawSimp *ptr = static_cast<tSedTransPwrLawSimp*>( oPtr );
  kf = ptr->kf;
  mf = ptr->mf;
  nf = ptr->nf;
}

tSedTransPwrLawSimp::tSedTransPwrLawSimp( const tSedTransPwrLawSimp &orig )
  : tSedTrans()
{
  kf = orig.kf;
  mf = orig.mf;
  nf = orig.nf;
}

/***************************************************************************\
 **  tSedTransPwrLawSimp::TransCapacity
 **
 **  Ultra-simple sediment transport function -> using with parker-wear-rule
 **  Qs = kf Q^mf S^nf
 **
 \***************************************************************************/
double tSedTransPwrLawSimp::TransCapacity( tLNode *node )
{
  const double slp = sin(atan(node->calcSlope()));
  //NMG added this 01/07 to use sin alpha rather than tan alpha,
  //which is returned by caclSlope()
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tBedErodePwrLaw::TransCapacity(tLNode*)");
  double tau=0;
  if( node->getFloodStatus() == tLNode::kNotFlooded )
  {
    tau = kf * pow( node->getQ(), mf ) * pow( slp, nf );
    node->setTau( tau );
    if(0) //DEBUG
      std::cout << "kf=" << kf << " Q=" << node->getQ() << 
      " S=" << node->calcSlope() << std::endl;
  }
  
  node->setQs( tau );
  return tau;
}


/***************************************************************************\
 **  tSedTransPwrLawSimp::TransCapacity
 **
 **  Ultra-simple sediment transport function -> using with parker-wear-rule
 **  Qs = kf Q^mf S^nf
 **  This is a weighted version which is called from DetachErode.
 **  Weight is a weighting by depth of layer.
 **  Here, qsi is set by proportion in layer and threshold is constant
 **  The value returned should be in units of m^3/yr
 **
 \***************************************************************************/
double tSedTransPwrLawSimp::TransCapacity( tLNode *node, int lyr, double weight )
{
  const double slp = sin(atan(node->calcSlope()));
  //NMG added this 01/07 to use sin alpha rather than tan alpha,
  //which is returned by caclSlope()
  
  if( slp < 0.0 )
    ReportFatalError("neg. slope in tSedTransPwrLaw::TransCapacity(tLNode*)");
  double tau=0;
  
  if( node->getFloodStatus() == tLNode::kNotFlooded )
  {
    tau = weight * kf * pow( node->getQ(), mf ) * pow( slp, nf );
    node->setTau( tau );
    // if((node->getVArea() < 2460 && node->getVArea() > 2459)  ||
    //   (node->getDrArea() < 5461 && node->getDrArea() > 5460) ){ //DEBUG
    if(0){ 
      std::cout << "kf=" << kf << " Q=" << node->getQ() << 
      " S=" << node->calcSlope() << " tau=" << tau <<
      " A = "<<node->getDrArea()<<" VA= "<<node->getVArea()<<std::endl;
      std::cout<<"x "<<node->getX()<<" y "<<node->getY()<<std::endl;
      std::cout<<"dwnstrmnbr A "<<node->getDownstrmNbr()->getDrArea()<<std::endl;
    }
  }
  
  for(size_t i=0; i<node->getNumg(); i++)
    node->addQs(i, tau*node->getLayerDgrade(lyr,i)/node->getLayerDepth(lyr));
  
  node->setQs( tau );
  //cout<<"Area " << node->getDrArea()<<" cap "<<cap<<endl;
  return tau;
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
  tauc = infile.ReadItem( tauc, "TAUCR" );  // to TAUCR apr 07 gt
  sqrtTauc = sqrt( tauc );
  
  // Add unit conversion factor for kt -- this is required to convert
  // the quantity (Q/W)^mb from units of years to units of seconds.
  kt *= pow( secPerYear, -mf );
}

void tSedTransBridgeDom::Initialize_Copy( tSedTrans* oPtr )
{
  tSedTransBridgeDom *ptr = static_cast<tSedTransBridgeDom*>( oPtr );
  kf = ptr->kf;
  kt = ptr->kt;
  mf = ptr->mf;
  nf = ptr->nf;
  tauc = ptr->tauc;  // to TAUCR apr 07 gt
  sqrtTauc = ptr->sqrtTauc; 
}

tSedTransBridgeDom::tSedTransBridgeDom( const tSedTransBridgeDom &orig )
  : tSedTrans()
{
  kf = orig.kf;
  kt = orig.kt;
  mf = orig.mf;
  nf = orig.nf;
  tauc = orig.tauc;  // to TAUCR apr 07 gt
  sqrtTauc = orig.sqrtTauc;
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

void tSedTransPwrLawMulti::Initialize_Copy( tSedTrans* oPtr )
{
  tSedTransPwrLawMulti *ptr = static_cast<tSedTransPwrLawMulti*>( oPtr );
  kf = ptr->kf;
  kt = ptr->kt;
  mf = ptr->mf;
  nf = ptr->nf;
  pf = ptr->pf;
  mdGrndiam = ptr->mdGrndiam;
  mdTauc = ptr->mdTauc;
  miNumgrnsizes = ptr->miNumgrnsizes;
  mdHidingexp = ptr->mdHidingexp;
}

tSedTransPwrLawMulti::tSedTransPwrLawMulti( const tSedTransPwrLawMulti &orig )
  : tSedTrans()
{
  kf = orig.kf;
  kt = orig.kt;
  mf = orig.mf;
  nf = orig.nf;
  pf = orig.pf;
  mdGrndiam = orig.mdGrndiam;
  mdTauc = orig.mdTauc;
  miNumgrnsizes = orig.miNumgrnsizes;
  mdHidingexp = orig.mdHidingexp;
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

void tSedTransWilcock::Initialize_Copy( tSedTrans* oPtr )
{
  tSedTransWilcock *ptr = static_cast<tSedTransWilcock*>( oPtr );
  taudim = ptr->taudim;
  refs = ptr->refs;
  refg = ptr->refg;
  lowtaucs = ptr->lowtaucs;
  lowtaucg = ptr->lowtaucg;
  sandb = ptr->sandb;
  hightaucs = ptr->hightaucs;
  hightaucg = ptr->hightaucg;
  sands = ptr->sands;
  gravb = ptr->gravb;
  gravs = ptr->gravs;
  grade = ptr->grade;
}

tSedTransWilcock::tSedTransWilcock( const tSedTransWilcock &orig )
  : tSedTrans()
{
  taudim = orig.taudim;
  refs = orig.refs;
  refg = orig.refg;
  lowtaucs = orig.lowtaucs;
  lowtaucg = orig.lowtaucg;
  sandb = orig.sandb;
  hightaucs = orig.hightaucs;
  hightaucg = orig.hightaucg;
  sands = orig.sands;
  gravb = orig.gravb;
  gravs = orig.gravs;
  grade = orig.grade;
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
  //tau = taudim*pow(nd->getHydrRough()*nd->getQ()/nd->getHydrWidth(), 0.6)*pow( nd->calcSlope(), 0.7);
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

void tSedTransMineTailings::Initialize_Copy( tSedTrans* oPtr )
{
  tSedTransMineTailings *ptr = static_cast<tSedTransMineTailings*>( oPtr );
  taudim = ptr->taudim;
  refs = ptr->refs;
  refg = ptr->refg;
  lowtaucs = ptr->lowtaucs;
  lowtaucg = ptr->lowtaucg;
  sandb = ptr->sandb;
  hightaucs = ptr->hightaucs;
  hightaucg = ptr->hightaucg;
  sands = ptr->sands;
  gravb = ptr->gravb;
  gravs = ptr->gravs;
  grade = ptr->grade;
}

tSedTransMineTailings::tSedTransMineTailings( const tSedTransMineTailings &orig )
  : tSedTrans()
{
  taudim = orig.taudim;
  refs = orig.refs;
  refg = orig.refg;
  lowtaucs = orig.lowtaucs;
  lowtaucg = orig.lowtaucg;
  sandb = orig.sandb;
  hightaucs = orig.hightaucs;
  hightaucg = orig.hightaucg;
  sands = orig.sands;
  gravb = orig.gravb;
  gravs = orig.gravs;
  grade = orig.grade;
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

/***************************************************************************\
 **  FUNCTIONS FOR CLASS tPhysicalWeatheringExpLaw
 \***************************************************************************/

/***************************************************************************\
 **  tPhysicalWeatheringExpLaw Constructor
 **
 **  Reads in parameters for the simple soil production function.
 **  "soilprodK" is the uniform and constant soil production rate for 
 **  zero soil depth; "soilprodH" is the uniform and constant soil depth
 **  scale controlling the rate of decline of soil production rate with
 **  increasing soil depth.
 **     (STL 6/10)
 **
 \***************************************************************************/
//constructor: reads and sets the parameters
tPhysicalWeatheringExpLaw::
tPhysicalWeatheringExpLaw( const tInputFile &infile )
{
  soilprodK = infile.ReadItem( soilprodK, "SOILPRODRATE" );
  soilprodH = infile.ReadItem( soilprodH, "SOILPRODDEPTH" );
}

//copy constructor: 
tPhysicalWeatheringExpLaw::
tPhysicalWeatheringExpLaw( const tPhysicalWeatheringExpLaw &orig )
  : tPhysicalWeathering()
{
  soilprodK = orig.soilprodK;
  soilprodH = orig.soilprodH;
}

// Initialize() for CSDMS IRF interface:
void tPhysicalWeatheringExpLaw::Initialize( const tInputFile &infile )
{
  soilprodK = infile.ReadItem( soilprodK, "SOILPRODRATE" );
  soilprodH = infile.ReadItem( soilprodH, "SOILPRODDEPTH" );
}

void tPhysicalWeatheringExpLaw::Initialize_Copy( tPhysicalWeathering* oPtr )
{
  tPhysicalWeatheringExpLaw* ptr = static_cast<tPhysicalWeatheringExpLaw*>( oPtr );
  soilprodK = ptr->soilprodK;
  soilprodH = ptr->soilprodH;
}

/***************************************************************************\
 **  tPhysicalWeatheringExpLaw::SoilProduction (1 of 3)
 **
 **  Computes the depth of weathering over a time interval dt.
 **  Simply calls 3rd SoilProduction below and multiplies result by dt.
 \***************************************************************************/
double tPhysicalWeatheringExpLaw::
SoilProduction( tLNode * n, double dt, double )
{
  return( SoilProduction( n ) * dt );
}

double tPhysicalWeatheringExpLaw::Run_Step( tLNode * n, double dt, double time ) 
{
  return SoilProduction( n, dt, time );
}

/***************************************************************************\
 **  tPhysicalWeatheringExpLaw::SoilProduction (2 of 3)
 **
 **  Computes the rate of weathering for a layer, i.e., rate of bedrock 
 **  lowering.
 **  Uses Ahnert/Heimsath soil production function.
 \***************************************************************************/
double tPhysicalWeatheringExpLaw::SoilProduction( tLNode * n, int i )
{
  // don't do physical weathering if layer is sediment
  if( n->getLayerSed(i) == tLayer::kSed ) return 0.0;
  // find soil/sediment thickness above bedrock:
  tListIter< tLayer > lI( n->getLayersRefNC() );
  double soilThickness(0.0);
  tLayer *lP=0;
  for( lP=lI.FirstP(); lP->getSed() == tLayer::kSed; lP=lI.NextP() )
    soilThickness += lP->getDepth();
  for( ; lP->getID() < i; lP=lI.NextP() ) ;
  // calculate rate of bedrock lowering (hence negative sign):
  const double slope = n->calcSlope();
  const double costheta = cos( atan( slope ) );
  double rate = -soilprodK * exp( -soilThickness * costheta / soilprodH );
  return rate;
}

double tPhysicalWeatheringExpLaw::Run_Step( tLNode * n, int i ) 
{
  return SoilProduction( n, i );
}

/***************************************************************************\
 **  tPhysicalWeatheringExpLaw::SoilProduction (3 of 3)
 **
 **  Computes the rate of weathering for the topmost rock layer.
 **  Simply calls 3rd SoilProduction below and multiplies result by dt.
 \***************************************************************************/
double tPhysicalWeatheringExpLaw::SoilProduction( tLNode * n )
{
  // find top bedrock layer and depth of soil above it:
  tListIter< tLayer > lI( n->getLayersRefNC() );
  double soilThickness(0.0);
  tLayer *lP=0;
  for( lP=lI.FirstP(); lP->getSed() == tLayer::kSed; lP=lI.NextP() )
    soilThickness += lP->getDepth();
  // calculate rate of bedrock lowering (hence negative sign):
  const double slope = n->calcSlope();
  const double costheta = cos( atan( slope ) );
  double rate = -soilprodK * exp( -soilThickness * costheta / soilprodH );
  return rate;
}

double tPhysicalWeatheringExpLaw::Run_Step( tLNode * n ) 
{
  return SoilProduction(n);
}


void tPhysicalWeatheringExpLaw::Finalize() {}

/***************************************************************************\
 **  FUNCTIONS FOR CLASS tPhysicalWeatheringDensityDependent
 \***************************************************************************/

/***************************************************************************\
 **  tPhysicalWeatheringDensityDependent Constructor
 **
 **  Reads in parameters for the simple soil production function.
 **  "soilprodK" is the uniform and constant soil production rate for 
 **  zero soil depth; "soilprodH" is the uniform and constant soil depth
 **  scale controlling the rate of decline of soil production rate with
 **  increasing soil depth.
 **     (STL 6/10)
 **
 \***************************************************************************/
//constructor: reads and sets the parameters
tPhysicalWeatheringDensityDependent::
tPhysicalWeatheringDensityDependent( const tInputFile &infile )
{
  soilprodK0 = infile.ReadItem( soilprodK0, "SOILPRODRATEINTERCEPT" );
  soilprodK1 = infile.ReadItem( soilprodK1, "SOILPRODRATESLOPE" );
  soilprodH = infile.ReadItem( soilprodH, "SOILPRODDEPTH" );
}

//copy constructor: 
tPhysicalWeatheringDensityDependent::
tPhysicalWeatheringDensityDependent( const tPhysicalWeatheringDensityDependent &orig )
  : tPhysicalWeathering()
{
  soilprodK0 = orig.soilprodK0;
  soilprodK1 = orig.soilprodK1;
  soilprodH = orig.soilprodH;
}

// Initialize() for CSDMS IRF interface:
void tPhysicalWeatheringDensityDependent::Initialize( const tInputFile &infile )
{
  soilprodK0 = infile.ReadItem( soilprodK0, "SOILPRODRATEINTERCEPT" );
  soilprodK1 = infile.ReadItem( soilprodK1, "SOILPRODRATESLOPE" );
  soilprodH = infile.ReadItem( soilprodH, "SOILPRODDEPTH" );
}

void tPhysicalWeatheringDensityDependent::Initialize_Copy( tPhysicalWeathering* oPtr )
{
  tPhysicalWeatheringDensityDependent* ptr = 
    static_cast<tPhysicalWeatheringDensityDependent*>( oPtr );
  soilprodK0 = ptr->soilprodK0;
  soilprodK1 = ptr->soilprodK1;
  soilprodH = ptr->soilprodH;
}


/***************************************************************************\
 **  tPhysicalWeatheringDensityDependent::SoilProduction (1 of 3)
 **
 **  Computes the depth of weathering over a time interval dt.
 **  Simply calls 3rd SoilProduction below and multiplies result by dt.
 \***************************************************************************/
double tPhysicalWeatheringDensityDependent::
SoilProduction( tLNode * n, double dt, double )
{
  return( SoilProduction( n ) * dt );
}

double tPhysicalWeatheringDensityDependent::Run_Step( tLNode * n, double dt, double time ) 
{
  return SoilProduction( n, dt, time );
}


/***************************************************************************\
 **  tPhysicalWeatheringDensityDependent::SoilProduction (2 of 3)
 **
 **  Computes the rate of weathering for a layer.
 **  Uses Ahnert/Heimsath soil production function.
 **  Preferred function only if you need to find rate for a layer other
 **  than the topmost bedrock (and for that it is a bit suspect).
 \***************************************************************************/
double tPhysicalWeatheringDensityDependent::SoilProduction( tLNode * n, int i )
{
  // don't do physical weathering if layer is sediment
  if( n->getLayerSed(i) == tLayer::kSed ) return 0.0;
  // find soil/sediment thickness above bedrock:
  tListIter< tLayer > lI( n->getLayersRefNC() );
  double soilThickness(0.0);
  tLayer *lP=0;
  int j=0;
  for( lP=lI.FirstP(), j=0; lP->getSed() == tLayer::kSed; lP=lI.NextP(), ++j )
    soilThickness += lP->getDepth();
  // keep adding thickness below bedrock surface 
  // (is this the right way to do this? maybe):
  for( ; j < i; lP=lI.NextP(), ++j )
    soilThickness += lP->getDepth();
  // calculate rate of bedrock lowering (hence negative sign):
  const double slope = n->calcSlope();
  const double costheta = cos( atan( slope ) );
  double rate = -( soilprodK0 - soilprodK1 * lP->getBulkDensity() ) 
  * exp( -soilThickness * costheta / soilprodH );
  return rate;
}

double tPhysicalWeatheringDensityDependent::Run_Step( tLNode * n, int i ) 
{
  return SoilProduction( n, i );
}

/***************************************************************************\
 **  tPhysicalWeatheringDensityDependent::SoilProduction (3 of 3)
 **
 **  Computes the rate of weathering for the topmost rock layer.
 **  Uses Ahnert/Heimsath soil production function.
 **  Preferred function because it automatically calculates rate at top of
 **  bedrock.
 \***************************************************************************/
double tPhysicalWeatheringDensityDependent::SoilProduction( tLNode * n )
{
  // find top bedrock layer and depth of soil above it:
  tListIter< tLayer > lI( n->getLayersRefNC() );
  double soilThickness(0.0);
  tLayer *lP=0;
  for( lP=lI.FirstP(); lP->getSed() == tLayer::kSed; lP=lI.NextP() )
    soilThickness += lP->getDepth();
  // get bedrock surface bulk density:
  double rockDensity = lP->getBulkDensity();
  // calculate rate of bedrock lowering (hence negative sign):
  const double slope = n->calcSlope();
  const double costheta = cos( atan( slope ) );
  double rate = -( soilprodK0 - soilprodK1 * rockDensity ) 
  * exp( -soilThickness * costheta / soilprodH );
  return rate;
}

double tPhysicalWeatheringDensityDependent::Run_Step( tLNode * n ) 
{
  return SoilProduction(n);
}


void tPhysicalWeatheringDensityDependent::Finalize() {}

/***************************************************************************\
 **  FUNCTIONS FOR CLASS tChemicalWeatheringDissolution
 \***************************************************************************/

/***************************************************************************\
 **  tChemicalWeatheringDissolution Constructor
 **
 **  Reads in parameters for the simple dissolution chemical weathering.
 **  "maxDissolution" is the uniform and constant dissolution rate at the
 **  bedrock surface; "chemDepth" is the uniform and constant bedrock depth
 **  scale controlling the rate of decline of dissolution rate with
 **  increasing depth below the bedrock surface.
 **     (STL 6/10)
 **
 \***************************************************************************/
//constructor: reads and sets the parameters
tChemicalWeatheringDissolution::
tChemicalWeatheringDissolution( const tInputFile &infile, 
                               tMesh<tLNode> *meshPtr )
{
  Initialize( infile, meshPtr );
}

//copy constructor: 
tChemicalWeatheringDissolution::
tChemicalWeatheringDissolution( const tChemicalWeatheringDissolution &orig )
  : tChemicalWeathering()
{
  maxDissolution = orig.maxDissolution;
  chemDepth = orig.chemDepth;
  rockBulkDensity_0 = orig.rockBulkDensity_0;
  rockLayerDepth = orig.rockLayerDepth;
  numThinLayers = orig.numThinLayers;
}

// Initialize() for CSDMS IRF interface:
void tChemicalWeatheringDissolution::Initialize( const tInputFile &infile, 
                                                tMesh<tLNode> *meshPtr )
{
  maxDissolution = infile.ReadItem( maxDissolution, "MAXDISSOLUTIONRATE" );
  chemDepth = infile.ReadItem( chemDepth, "CHEMDEPTH" );
  rockBulkDensity_0 = infile.ReadItem( rockBulkDensity_0, "ROCKDENSITYINIT" );
  const double weatheringDepth = 10.0 * chemDepth; // >>chemDepth
  rockLayerDepth = chemDepth / 20; // <<chemDepth
  numThinLayers = ROUND( weatheringDepth / rockLayerDepth );
  Initialize( meshPtr );
}

// Initialize() for CSDMS IRF interface:
void tChemicalWeatheringDissolution::Initialize( tMesh<tLNode> *meshPtr )
{
  tMesh< tLNode >::nodeListIter_t ni( meshPtr->getNodeList() ); // node iter.
  // make bedrock layers for each active node:
  for( tLNode* n = ni.FirstP(); ni.IsActive(); n = ni.NextP() )
  {
    // currently one very thick bedrock layer; set its bulk density,
    // and then make copies of that layer; number and thickness are
    // dependent on the chemical weathering depth scale, chemDepth:
    tListIter< tLayer > lI( n->getLayersRefNC() ); // layer iterator
    tLayer *lP=0; // current layer pointer
    // find top bedrock layer and bottom soil layer:
    for( lP=lI.FirstP(); lP->getSed() == tLayer::kSed; lP=lI.NextP() ) ;
    lP->setBulkDensity( rockBulkDensity_0 );
    // insert copies at bottom of layer list:
    for( int i=0; i<numThinLayers; ++i )
      n->getLayersRefNC().insertAtBack( *lP );
    // now have a bunch of thick layers;
    // reset depth for all rock layers but the thick bottom layer:
    {
      int i;
      for( lP=lI.FirstP(), i=0; 
          lP->getSed() == tLayer::kSed; 
          lP=lI.NextP(), ++i ) ;
      const int numLayers = n->getNumLayer();
      for( ; i<numLayers-1; ++i, lP=lI.NextP() )
        lP->setDepth( rockLayerDepth );
    }
  }
}

void tChemicalWeatheringDissolution::Initialize_Copy( tChemicalWeathering* oPtr, 
						      tMesh<tLNode>* mPtr )
{
  tChemicalWeatheringDissolution* ptr = static_cast<tChemicalWeatheringDissolution*>( oPtr );
  maxDissolution = ptr->maxDissolution;
  chemDepth = ptr->chemDepth;
  rockBulkDensity_0 = ptr->rockBulkDensity_0;
  rockLayerDepth = ptr->rockLayerDepth;
  numThinLayers = ptr->numThinLayers;
  Initialize( mPtr );
}

/***************************************************************************\
 **  tChemicalWeatheringDissolution::SoluteFlux (1 of 3)
 **
 **  Computes the chemical weathering for the whole rock column 
 **  over a time interval dt.
 **  Calculates density change (kg/m3) with exponential decay with depth.
 **  Updates density of each rock layer.
 **  Returns total mass flux (kg), typically negative.
 **  Preferred function for doing chemical weathering at a point.
 \***************************************************************************/
double tChemicalWeatheringDissolution::SoluteFlux( tLNode * n, double dt )
{
  // find top of bedrock:
  tListIter< tLayer > lI( n->getLayersRefNC() );
  tLayer *lP=0;
  for( lP=lI.FirstP(); lP->getSed() == tLayer::kSed; lP=lI.NextP() ) ;
  // find bedrock depth and calculate solute flux rate at each layer:
  double bedrockDepth(0.0);
  double flux(0.0);
  for( int i=0; i<numThinLayers; ++i, lP=lI.NextP() )
  {
    // check layer depth; if thick, insert a copy and make it thin 
    // (compare to initial rockLayerDepth, but allow for thickening
    // due to strain associated with weathering):
    if( lP->getDepth() > rockLayerDepth * 20.0 )
    {
      n->getLayersRefNC().insertAtPrev( *lP, lI.NodePtr() );
      lP = lI.PrevP();
      lP->setDepth( rockLayerDepth );
    }
    // calculate bedrock density change (dissolution,  
    // hence negative sign); use depth at top of layer:
    double deltaRho = 
    -maxDissolution * exp( -bedrockDepth / chemDepth ) * dt;
    // update layer bulk density:
    lP->addBulkDensity( deltaRho );
    // and increment total mass flux per unit area:
    flux += deltaRho * lP->getDepth();
    if( lP->getBulkDensity() > 0.0 )
      // increment depth:
      bedrockDepth += lP->getDepth();
    else
    {
      // don't expect this to happen, but if density drops to zero,
      // remove layer and change elevation:
      n->ChangeZ( -lP->getDepth() );
      lP = lI.PrevP();
      tLayer lay;
      n->getLayersRefNC().removeNext( lay, lI.NodePtr() );
      assert( n>0 );
    }
  }
  // multiply by area for total mass flux:
  flux *=  n->getVArea();
  return flux;
}

double tChemicalWeatheringDissolution::Run_Step( tLNode * n, double dt ) 
{
  return SoluteFlux( n, dt );
}

/***************************************************************************\
 **  tChemicalWeatheringDissolution::SoluteFlux (2 of 3)
 **
 **  Computes the rate of weathering for a layer, i.e., rate of bedrock 
 **  density decrease.
 **  Uses exponential decay with depth function.
 **  Returns rate of density change (kg/m3/yr), typically negative.
 **  Preferred function for updating rock density at one layer.
 \***************************************************************************/
double tChemicalWeatheringDissolution::SoluteFlux( tLNode * n, int i )
{
  // don't do chemical weathering if layer is sediment
  if( n->getLayerSed(i) == tLayer::kSed ) return 0.0;
  // find top of bedrock:
  tListIter< tLayer > lI( n->getLayersRefNC() );
  tLayer *lP=0;
  for( lP=lI.FirstP(); lP->getSed() == tLayer::kSed; lP=lI.NextP() ) ;
  // find bedrock depth at layer:
  double bedrockDepth(0.0);
  for( int j=0; j < i; ++j, lP=lI.NextP() )
    bedrockDepth += lP->getDepth();
  // calculate rate of bedrock density change (dissolution, hence negative sign):
  double rate = -maxDissolution * exp( -bedrockDepth / chemDepth );
  return rate;
}

double tChemicalWeatheringDissolution::Run_Step( tLNode * n, int i ) 
{
  return SoluteFlux( n, i );
}

/***************************************************************************\
 **  tChemicalWeatheringDissolution::SoluteFlux (3 of 3)
 **
 **  Computes the rate of weathering for the whole rock column.
 **  Uses exponential decay with depth function.
 **  Does not update rock density, since no dt.
 **  Returns rate of density change (kg/m3/yr), typically negative.
 **  Preferred function only for finding total solute flux at a node.
 \***************************************************************************/
double tChemicalWeatheringDissolution::SoluteFlux( tLNode * n )
{
  // find top of bedrock:
  tListIter< tLayer > lI( n->getLayersRefNC() );
  tLayer *lP=0;
  for( lP=lI.FirstP(); lP->getSed() == tLayer::kSed; lP=lI.NextP() ) ;
  // find bedrock depth and calculate solute flux rate at each layer:
  double bedrockDepth(0.0);
  double rate(0.0);
  for( int i=0; i<numThinLayers; ++i, lP=lI.NextP() )
  {
    // check layer depth; if thick, insert a copy and make it thin 
    // (compare to initial rockLayerDepth, but allow for thickening
    // due to strain associated with weathering):
    if( lP->getDepth() > rockLayerDepth * 20.0 )
    {
      n->getLayersRefNC().insertAtPrev( *lP, lI.NodePtr() );
      lP = lI.PrevP();
      lP->setDepth( rockLayerDepth );
    }
    // calculate and increment flux rate per unit area (dissolution,  
    // hence negative sign); use depth at top of layer:
    rate += -maxDissolution * exp( -bedrockDepth / chemDepth )
    * lP->getDepth();
    // increment depth:
    bedrockDepth += lP->getDepth();
  }
  // multiply rate by area:
  rate *= n->getVArea();
  return rate;
}

double tChemicalWeatheringDissolution::Run_Step( tLNode * n ) 
{
  return SoluteFlux(n);
}

void tChemicalWeatheringDissolution::Finalize() {}

/***************************************************************************\
 **  FUNCTIONS FOR CLASS tDebrisFlow
 \***************************************************************************/
//usual constructor
inline tDebrisFlow::tDebrisFlow( tPtrList<tLNode> &PList, double netForce, 
                                vector<double> &nodeSoilThickness, 
                                vector<double> &nodeWoodDepth,
                                vector<double> &nodeWaterDepth,
                                vector<int> &nodeIndex,
                                tMesh<tLNode> *meshPtr, tErosion *ePtr )
{ Initialize( PList, netForce, nodeSoilThickness, nodeWoodDepth, 
             nodeWaterDepth, nodeIndex, meshPtr, ePtr ); }

//copy constructor
tDebrisFlow::tDebrisFlow( const tDebrisFlow& orig )
{
  areaFailure = orig.areaFailure;
  areaScour = orig.areaScour;
  areaDeposit = orig.areaDeposit;
  volumeFailure = orig.volumeFailure;
  sedimentVolume = orig.sedimentVolume;
  woodVolume = orig.woodVolume;
  waterVolume = orig.waterVolume;
  areaFrontal = orig.areaFrontal;
  widthFrontal = orig.widthFrontal;
  netForce = orig.netForce;
  orgPtr = orig.orgPtr;
  atPtr = orig.atPtr;
  erosion = orig.erosion;
  slideCluster = orig.slideCluster;
  scourCluster = orig.scourCluster;
  scourZoneMesh = orig.scourZoneMesh;
  depositCluster = orig.depositCluster;
  depositZoneMesh = orig.depositZoneMesh;
  wasList = orig.wasList;
  velocityList = orig.velocityList;
}

/***************************************************************************\
 **  void tDebrisFlow::Initialize
 **
 **  Function to initialize debris flow.
 **  Called by: tDebrisFlow constructor 
 **  Takes: 
 **    - PList: a copy of the slideCluster
 **    - netForce: a copy of the initial net driving force 
 **      (not used in this version)
 **    - nodeSoilThickness: a reference to a vector containing soil
 **      thicknesses for each active node
 **    - nodeWoodDepth: similar to above
 **    - nodeWaterDepth: similar to above
 **    - nodeIndex: a reference to a vector containing the indexes to the
 **      nodeSoilThickness vector by node ID
 **    - meshPtr: a pointer to the simulation tMesh<tLNode>
 **    - ePtr: pointer to tErosion object
 **  Calls: 
 **    - tMesh<tLNode>::LocateTriangle
 **    - ScanLineNodesByPublicFlag (global)
 **    - RightUnitVector (global)
 **    - DotProduct2D (global)
 **    - tDF_RunOut::Start: "rules" function governing runout
 **  Changes: 
 **    - this: sets slideCluster, atPtr, volumeFailure, areaFrontal, 
 **      widthFrontal
 **    - nodes in slideCluster: removes soil from and changes elevations
 **
 **  - STL, 9/2010
 \***************************************************************************/
void tDebrisFlow::Initialize( tPtrList<tLNode> &PList, 
                             double initNetForce,
                             vector<double> &nodeSoilThickness,
                             vector<double> &nodeWoodDepth,
                             vector<double> &nodeWaterDepth,
                             vector<int> &nodeIndex,
                             tMesh<tLNode> *meshPtr, 
                             tErosion *ePtr )
{
  // initialize quantities to zero:
  areaFailure = areaScour = areaDeposit = volumeFailure = sedimentVolume 
  = woodVolume = waterVolume = areaFrontal = widthFrontal = 0.0;
  netForce = initNetForce;
  erosion = ePtr;
  slideCluster = new tPtrList<tLNode>( PList );
  tPtrListIter<tLNode> cI( slideCluster );
  scourCluster = 0;
  scourZoneMesh = 0;
  depositCluster = 0;
  depositZoneMesh = 0;
  wasList = 0;
  velocityList = 0;
  
  // find downstream-most node
  double maxArea = 0.0;
  double volWeightSumX = 0.0;
  double volWeightSumY = 0.0;
  // find landslide volume:
  // for each node in failure cluster:
  for( tLNode *cn = cI.FirstP(); !cI.AtEnd(); cn = cI.NextP() )
  {
    // add up landslide area:
    areaFailure += cn->getVArea();
    // add up landslide volume:
    const double nodeVolume = 
    nodeSoilThickness[nodeIndex[cn->getID()]] * cn->getVArea();
    volumeFailure += nodeVolume;
    sedimentVolume += nodeVolume;
    woodVolume += nodeWoodDepth[nodeIndex[cn->getID()]] * cn->getVArea();
    waterVolume += nodeWaterDepth[nodeIndex[cn->getID()]] * cn->getVArea();
    // add up volume-weighted coordinates:
    volWeightSumX += cn->getX() * nodeVolume;
    volWeightSumY += cn->getY() * nodeVolume;
    tLayer cl;
    // remove its sediment (soil) layers and change elevations:
    while( cn->getLayerSed(0) == tLayer::kSed ) 
    {
      cn->getLayersRefNC().removeFromFront( cl );
      cn->ChangeZ( -cl.getDepth() );
    }
    // find downstream-most node:
    if( cn->getDrArea() > maxArea )
    {
      atPtr = cn;
      maxArea = cn->getDrArea();
    }
  }
  orgPtr = atPtr;
  // set "center" of cluster to point with greatest drainage area to start:
  tLNode *gn = orgPtr;
  // then find geometric center if the cluster size if greater than 2 nodes:
  const int numCluster = orgPtr->public1;
  if( slideCluster->getSize() > 2 )
  {
    //       if( slideCluster->getSize() == 2 )
    // 	{ // for two nodes, LocateTriangle can fail 'cos point is on the edge:
    // 	  gn = cI.FirstP();
    // 	  const double nodeVol1 = 
    // 	    nodeSoilThickness[nodeIndex[gn->getID()]] * gn->getVArea();
    // 	  gn = cI.NextP();
    // 	  const double nodeVol2 =
    // 	    nodeSoilThickness[nodeIndex[gn->getID()]] * gn->getVArea();
    // 	  if( nodeVol1 > nodeVol2 ) gn = cI.FirstP();
    // 	}
    //       else
    // 	{
    // find coordinates of center of mass (volume):
    const double wtAvgX = volWeightSumX / volumeFailure;
    const double wtAvgY = volWeightSumY / volumeFailure;
    // find triangle containing center of mass:
    tTriangle *gt = meshPtr->LocateTriangle( wtAvgX, wtAvgY );
    // if geometric center is within bounds:
    if( gt>0 )
    { // find triangle's node with greatest drainage area (and in cluster):
      double gArea = 0.0;
      for( int i=0; i < 3; ++i )
	    { // and if a node near that center is within the cluster:
	      tLNode *cn = static_cast<tLNode*>( gt->pPtr(i) );
	      if( cn->public1 == numCluster && cn->getDrArea() > gArea )
        { // set "center" to node at geometric center 
          gArea = cn->getDrArea();
          gn = cn;
        }
	    }
    }
    // 	}
  }
  // refer to "rules" object to see whether to prepare for runout:
  if( erosion->getDF_RunOutPtr()->Start( this ) )
  { // find nodes along line perpendicular to the flowedge of the 
    // center-of-mass node up to and including the first nodes
    // outside of those with the same value of tLNode::public1:
    tPtrList<tLNode> scanLineNodesList;
    ScanLineNodesByPublicFlag( gn, scanLineNodesList );
    vector<tLNode*> scanLineNode( scanLineNodesList.getSize() );
    {
      int i=0;
      while( tLNode* cn = scanLineNodesList.removeFromFront() )
        scanLineNode[i++] = cn;
    }
    // calculate cross-sectional area and width of landslide at center of mass;
    // note that, although soil already stripped above, vector of soil
    // thicknesses was not changed:
    // find perpendicular distances from left to right with dot product:
    areaFrontal = 0.0;
    widthFrontal = 0.0;
    const int numScanLine = scanLineNode.size();
    tArray<double> rightUnitVec = RightUnitVector( gn->getFlowEdg() );
    for( int i=0; i<numScanLine-1; ++i )
    {
      double perpDist = 
	    fabs( DotProduct2D( scanLineNode[i], scanLineNode[i+1], 
                         rightUnitVec ) );
      double vertDepth=0.0;
      if( i > 0 && i < numScanLine-2 )
        vertDepth = 
	      ( nodeSoilThickness[ nodeIndex[ scanLineNode[i]->getID() ] ] 
         + nodeSoilThickness[ nodeIndex[ scanLineNode[i+1]->getID() ] ] ) 
	      / 2.0;
      else
	    {
	      // if at ends, use soil thickness for node in cluster:
	      if( i == 0 )
          vertDepth = 
          nodeSoilThickness[nodeIndex[scanLineNode[i+1]->getID()]];
	      else
          vertDepth =
          nodeSoilThickness[nodeIndex[scanLineNode[i]->getID()]];
	      // if at ends, use only half the distance:
	      perpDist /= 2.0;
	    }
      areaFrontal += perpDist * vertDepth;
      widthFrontal += perpDist;
    }
  }
}



/***************************************************************************\
 **  void tDebrisFlow::RunScourDeposit
 **
 **  Function to oversee debris flow runout, scour, and deposition.
 **  Called by: tErosion::LandslideClusters 
 **  Takes: 
 **    - Run time
 **  Calls: 
 **    - Global functions to do scanline and zone stuff:
 **      - ScanLineNodesByDistance
 **      - FillCrossSection
 **      - BuildPublicClusterWithMesh
 **    - Runout function: tDF_RunOut::InMotion
 **    - Scour functions:
 **      - tDF_Scour::InScourZone
 **      - tDF_Scour::BedScour
 **      - tDF_Deposit::InDepositionZone
 **      - tDF_Deposit::FormDeposit
 **    - Temporary tMesh constructor: 
 **        tMesh<tLNode>( x-coords, y-coords, z-coords )
 **  Changes: 
 **    - this: sets scourCluster, depositCluster, scourZoneMesh, 
 **      depositZoneMesh, atPtr, volumeFailure, areaFrontal, 
 **      widthFrontal, sedimentVolume, woodVolume, waterVolume;
 **      other objects (wasList, velocityList) will be made in runout,
 **      scour, and deposition objects
 **    - call of scour and deposition objects results in any erosion,
 **      deposition, and elevation changes
 **
 **  - STL, 9/2010
 \***************************************************************************/
void tDebrisFlow::RunScourDeposit()
{
  // initial value of atPtr is orgPtr:
  tEdge* fe = atPtr->getFlowEdg();
  // if at the outlet, don't need to do this:
  if( fe == 0 || !fe->getDestinationPtr()->isNonBoundary() ) 
    {
      erosion->debris_flow_sed_bucket += sedimentVolume;
      erosion->debris_flow_wood_bucket += woodVolume;
//       erosion->landslideAreas.insertAtBack( areaFailure );
      return;
    }
  // move flow downstream node by node; scour soil off nodes in runout
  // track and to each side up to a height such that the area above the
  // surface is equal to areaFrontal:
  tList< tArray<double> > scourZoneCoords;
  tList< tArray<double> > depositZoneCoords;
  // start at landslide:
  while( erosion->getDF_RunOutPtr()->InMotion( this ) )
  {
    if( erosion->getDF_ScourPtr()->InScourZone( this ) )
    {
      tPtrList<tLNode> scanLineNodesList;
      // finds nodes in scan line perpendicular to atPtr->flowedge;
      // use width of flow to limit scan:
      int iStream = ScanLineNodesByDistance( atPtr, widthFrontal, 
                                            scanLineNodesList );
      vector<tLNode*> scanLineNode( scanLineNodesList.getSize() );
      {
        int i=0;
        while( tLNode* cn = scanLineNodesList.removeFromFront() )
          scanLineNode[i++] = cn;
      }
      tArray<double> leftXYZ(3);
      tArray<double> rightXYZ(3);
      // finds endpoints of cross-section filled by areaFrontal:
      FillCrossSection( iStream, areaFrontal, widthFrontal, scanLineNode, 
                       leftXYZ, rightXYZ );
      scourZoneCoords.insertAtBack( leftXYZ );
      scourZoneCoords.insertAtBack( rightXYZ );
    }
    // note that we haven't pre-determined that the scour and deposition
    // zones are mutually exclusive, or that we must be in one or the other
    // AT THIS POINT; they become mutually exclusive below:
    if( erosion->getDF_DepositPtr()->InDepositionZone( this ) )
    {
      tPtrList<tLNode> scanLineNodesList;
      // finds nodes in scan line perpendicular to atPtr->flowedge;
      // use width of flow to limit scan:
      int iStream = ScanLineNodesByDistance( atPtr, widthFrontal, 
                                            scanLineNodesList );
      vector<tLNode*> scanLineNode( scanLineNodesList.getSize() );
      {
        int i=0;
        while( tLNode* cn = scanLineNodesList.removeFromFront() )
          scanLineNode[i++] = cn;
      }
      tArray<double> leftXYZ(3);
      tArray<double> rightXYZ(3);
      // finds endpoints of cross-section filled by areaFrontal:
      FillCrossSection( iStream, areaFrontal, widthFrontal, scanLineNode, 
                       leftXYZ, rightXYZ );
      depositZoneCoords.insertAtBack( leftXYZ );
      depositZoneCoords.insertAtBack( rightXYZ );
    }
    // finally, move one node downstream:
    atPtr = atPtr->getDownstrmNbr();
  }
  { // do scour if there's a scour zone:
    const int numVtcs = scourZoneCoords.getSize();
    if( numVtcs > 0 )
    {
      tArray<double> zoneX( numVtcs );
      tArray<double> zoneY( numVtcs );
      tArray<double> zoneZ( numVtcs );
      { // write out x, y, and z arrays defining scour zone:
        int i=0;
        tArray<double> XYZ(3);
        while( scourZoneCoords.removeFromFront( XYZ ) )
        {
          zoneX[i] = XYZ[0];
          zoneY[i] = XYZ[1];
          zoneZ[i] = XYZ[2];
          ++i;
        }
      }
      // instantiate scourZoneMesh from the scour zone coordinates:
      scourZoneMesh = new tMesh<tLNode>( zoneX, zoneY, zoneZ );
      // instantiate scourCluster:
      scourCluster = new tPtrList<tLNode>();
      // start potential cluster with node downstream of orgPtr 
      // (wouldn't have come this far in the function if it didn't have a 
      // downstream neighbor), 
      // use flag value from slideCluster:
      BuildPublicClusterWithMesh( scourZoneMesh, scourCluster, 
                                 orgPtr->getDownstrmNbr(), 
                                 orgPtr->public1 );
      // remove soil from and change elevation of nodes in scour zone:
      tPtrListIter<tLNode> cI( scourCluster );
      for( tLNode *cn = cI.FirstP(); !cI.AtEnd(); cn = cI.NextP() )
      {
        erosion->getDF_ScourPtr()->BedScour( this, cn );
        areaScour += cn->getVArea();
      }
    }
  } // end scour
  // update buckets:
  erosion->debris_flow_sed_bucket += sedimentVolume;
  erosion->debris_flow_wood_bucket += woodVolume;
//   erosion->landslideAreas.insertAtBack( areaFailure + areaScour );
  { // do deposition if there's a deposition zone:
    const int numVtcs = depositZoneCoords.getSize();
    if( numVtcs > 0 )
    {
      tArray<double> zoneX( numVtcs );
      tArray<double> zoneY( numVtcs );
      tArray<double> zoneZ( numVtcs );
      { // write out x, y, and z arrays defining deposit zone:
        int i=0;
        tArray<double> XYZ(3);
        while( depositZoneCoords.removeFromFront( XYZ ) )
        {
          zoneX[i] = XYZ[0];
          zoneY[i] = XYZ[1];
          zoneZ[i] = XYZ[2];
          ++i;
        }
      }
      // instantiate depositZoneMesh from the deposit zone coordinates:
      depositZoneMesh = new tMesh<tLNode>( zoneX, zoneY, zoneZ );
      // instantiate depositCluster:
      depositCluster = new tPtrList<tLNode>();
      // note that, here, scourCluster and depositCluster ARE mutually
      // exclusive: the following function will not add a point to the 
      // cluster unless its public1 flag is zero:
      // start potential cluster with atPtr, use flag value from slideCluster:
      BuildPublicClusterWithMesh( depositZoneMesh, depositCluster, atPtr, 
                                 orgPtr->public1 );
      // add sediment and change elevation of nodes in deposit zone:
      tPtrListIter<tLNode> cI( depositCluster );
      for( tLNode *cn = cI.FirstP(); !cI.AtEnd(); cn = cI.NextP() )
      {
        erosion->getDF_DepositPtr()->FormDeposit( this, cn );
        areaDeposit += cn->getVArea();
      }
    }
  } // end deposition
}

/***************************************************************************\
 **  overloaded output operator for tDebrisFlow:
 **  - STL, 9/2010
 \***************************************************************************/
std::ostream& operator<<( std::ostream &os, tDebrisFlow const &dflow )
{
  os.setf( ios::fixed, ios::floatfield );
  os.precision(2);
  os << dflow.orgPtr->getX() << " " << dflow.orgPtr->getY() << " " 
  << dflow.orgPtr->getZ() << " "
  << dflow.atPtr->getX() << " " << dflow.atPtr->getY() << " " 
  << dflow.atPtr->getZ() << " "
  << dflow.areaFailure << " " << dflow.areaScour << " "
  << dflow.areaDeposit << " "
  << dflow.volumeFailure << " " 
  << dflow.sedimentVolume << " " << dflow.woodVolume << " "
  << dflow.waterVolume << " "
  << dflow.areaFrontal << " " << dflow.widthFrontal << " "
  << dflow.netForce;
  return os;
}

/***************************************************************************\
 **  FUNCTIONS FOR CLASS tDF_Scour AND DERIVED CLASSES
 \***************************************************************************/
/***************************************************************************\
 **  void tDF_ScourAllSediment::BedScour
 **
 **  Function to do debris flow scour, where all sediment/soil is 
 **    removed from each node; if there are trees, biomass is also 
 **    removed.
 **  Called by: tDebrisFlow::RunScourDeposit 
 **  Takes: 
 **    - DF: pointer to tDebrisFlow object
 **    - cn: pointer to tLNode object
 **  Calls: 
 **    - various "get" and "set" functions
 **    - tLNode::ChangeZ
 **  Changes: 
 **    - DF: adds scoured material to debris flow
 **    - cn: removes soil and wood and changes elevation
 **
 **  - STL, 9/2010
 \***************************************************************************/
void tDF_ScourAllSediment::BedScour( tDebrisFlow* DF, tLNode* cn )
{
  tLayer cl;
  // remove its sediment (soil) layers, add to tally, and change elevations:
  while( cn->getLayerSed(0) == tLayer::kSed ) 
  {
    cn->getLayersRefNC().removeFromFront( cl );
    DF->addSedimentVolume( cl.getDepth() * cn->getVArea() );
    cn->ChangeZ( -cl.getDepth() );
  }
  // remove biomass and add to debris flow tally:
  if( cn->getVegCover().getTrees() > 0 )
  {
    DF->addWoodVolume( ( cn->getVegCover().getTrees()->getBioMassStand()
                        + cn->getVegCover().getTrees()->getBioMassDown() )
                      * cn->getVArea() );
    cn->getVegCover().getTrees()->setBioMassStand(0.0);
    cn->getVegCover().getTrees()->setBioMassDown(0.0);
  }    
}

/***************************************************************************\
 **  void BuildPublicClusterWithMesh( tMesh<tLNode>* mesh, 
 **          			      tPtrList<tLNode>* cluster, 
 **                                   tLNode* seedNode, 
 **                                   const int flagVal )
 **
 **  Global function to build cluster defined by a tMesh and a seed node.
 **  Start with seed node (tLNode* cn) and grow cluster of nodes that
 **    (a) fall within one of the triangles of the mesh, and
 **    (b) have elevations that fall below the plane of that one triangle.
 **  Puts nodes in tPtrList<tLNode> passed by pointer.
 **
 **  Called by: tDebrisFlow::RunScourDeposit.
 **  Takes: 
 **    - mesh: pointer to tMesh<tLNode> that defines cluster;
 **    - cluster: pointer to tPtrList containing cluster nodes;
 **    - seedNode: pointer to tLNode that "seeds" the cluster;
 **    - flagval: value of tNode::public1 to assign nodes in cluster;
 **  Calls: 
 **    - tMesh::LocateTriangle(x,y) to find whether candidate node is
 **      within mesh;
 **    - PlaneFit (global) to find whether candidate node's elevation
 **      is below that of the plane of the triangle in which node is found.
 **  Changes: values in passed arrays leftXYZ and rightXYZ.
 **
 **  - STL, 7/2010
 \***************************************************************************/
void BuildPublicClusterWithMesh( tMesh<tLNode>* mesh, 
                                tPtrList<tLNode>* cluster, 
                                tLNode* seedNode, 
                                const int flagVal )
{
  tPtrList<tLNode> seedList; // list of nodes on edge of cluster
  //   tLNode* cn = seedNode;
  seedList.insertAtBack( seedNode );
  cluster->insertAtBack( seedNode );
  // set generic flag signifying node in cluster:
  seedNode->public1 = flagVal;
  while( !seedList.isEmpty() )
  {
    //       tLNode* cn = seedList.removeFromFront();
    tSpkIter sI( seedList.removeFromFront() );
    for( tEdge* ce = sI.FirstP(); !sI.AtEnd(); ce = sI.NextP() )
    {
      if( ce->getDestinationPtr()->isNonBoundary() )
	    {
	      tLNode* nn = static_cast<tLNode*>( ce->getDestinationPtrNC() );
	      // if node is not already in a cluster AND point coordinates
	      // fall within a triangle in the mesh AND
	      // the node's elevation is at or below the triangle's plane:
	      // NOTE: if we want flow depth (e.g., for determining 
	      // scour), difference between node elevation and 
	      // triangle surface elevation at that point is effectively
	      // the depth of the flow at that point.
	      tTriangle *ct =  mesh->LocateTriangle( nn->getX(), nn->getY() );
	      if( nn->public1 == 0 
           && 
           ct > 0 
           && 
           nn->getZ() <= PlaneFit( nn->getX(), nn->getY(), ct ) )
        {
          // then add the new node to the cluster:
          nn->public1 = flagVal;
          seedList.insertAtBack( nn );
          cluster->insertAtBack( nn );
        }
	    }
    }
  }
}

/***************************************************************************\
 **  FUNCTIONS FOR CLASS tErosion
 \***************************************************************************/

//constructor
tErosion::tErosion( tMesh<tLNode> *mptr, const tInputFile &infile, 
		    bool no_write_mode /* = false */ ) :
meshPtr(mptr),
bedErode(0), sedTrans(0), physWeath(0), chemWeath(0), 
runout(0), scour(0), deposit(0), DF_fsPtr(0), DF_Hyd_fsPtr(0),
track_sed_flux_at_nodes_( false ), water_sed_tracker_ptr_(NULL),
soilBulkDensity(kDefaultSoilBulkDensity),
rockBulkDensity(kDefaultRockBulkDensity),
wetBulkDensity(kDefaultWetBulkDensity), 
woodDensity(450.0), fricSlope(1.0), num_grain_sizes_(1),
debris_flow_sed_bucket(0), debris_flow_wood_bucket(0)
{
  assert( mptr!=0 );
  
  // Read parameters needed from input file
  num_grain_sizes_ = infile.ReadInt( "NUMG" );
  infile.ReadItem( kd_ts, "KD" );  // Hillslope diffusivity coefficient
  kd = kd_ts.calc(0.);
  //kd = infile.ReadItem( kd, "KD" );  // Hillslope diffusivity coefficient
  difThresh = infile.ReadItem( difThresh, "DIFFUSIONTHRESHOLD");
  bool optNonlinearDiffusion = infile.ReadBool( "OPT_NONLINEAR_DIFFUSION", false );
  if( optNonlinearDiffusion )
    mdSc = infile.ReadItem( mdSc, "CRITICAL_SLOPE" );
  beta=infile.ReadItem( beta, "BETA"); //For Sediment-Flux Detach Rules
  bool optDepthDependentDiffusion = 
  infile.ReadBool( "OPT_DEPTH_DEPENDENT_DIFFUSION", false );
  if( optDepthDependentDiffusion)
    diffusionH = infile.ReadDouble( "DIFFDEPTHSCALE", true );
  else
    if( num_grain_sizes_>1 )
	  ReportFatalError( "When using multiple grain sizes, you must select depth-dependent diffusion and give a value for DIFFDEPTHSCALE" );
  
  int optAdaptMesh = infile.ReadItem( optAdaptMesh, "OPTMESHADAPTDZ", false );
  if( optAdaptMesh )
    mdMeshAdaptMaxFlux = infile.ReadItem( mdMeshAdaptMaxFlux,
                                         "MESHADAPT_MAXNODEFLUX" );
  
  // Make sure the user wants the detachment and transport options that
  // are compiled in this version
  
  // set bedrock detachment law:
  optBedErosionLaw = infile.ReadItem( optBedErosionLaw,
                                      "DETACHMENT_LAW" );
  switch(optBedErosionLaw){
#define X(a,b) case a: \
bedErode = new b(infile);			\
break;
      DETACHMENT_LAW_TABLE2
#undef X
    default:
    {
      std::cerr << "\nError: You requested the detachment law '"
      << optBedErosionLaw << "' which does not exist."  << std::endl
      << "Available options:" << std::endl;
      for(int i=0; i!=NUMBER_OF_DETACHMENT_LAWS; ++i ){
        std::cerr << " " << i << ": " << DetachmentLaw[i] << std::endl;
      }
      ReportFatalError( "Requested detachment law not available. "
                       "Switch options.\n" );
    }
  }
  
  std::cout << "DETACHMENT OPTION: "
  << DetachmentLaw[optBedErosionLaw] << std::endl;
  
  // set sediment transport law:
  optSedTransLaw = infile.ReadItem( optSedTransLaw,
                                  "TRANSPORT_LAW" );
  switch(optSedTransLaw){
#define X(a,b) case a: \
sedTrans = new b(infile);			\
break;
      TRANSPORT_LAW_TABLE2
#undef X
    default:
    {
      std::cerr << "\nError: You requested the transport law '"
      << optSedTransLaw << "' which does not exist."  << std::endl
      << "Available options:" << std::endl;
      for(int i=0; i!=NUMBER_OF_TRANSPORT_LAWS; ++i ){
        std::cerr << " " << i << ": " << TransportLaw[i] << std::endl;
      }
      ReportFatalError( "Requested transport law not available. "
                       "Switch options.\n" );
    }
  }
  std::cout << "SEDIMENT TRANSPORT OPTION: "
	    << TransportLaw[optSedTransLaw] << std::endl;

  // set soil production law:
  optPhysWeathLaw = infile.ReadItem( optPhysWeathLaw,
                                  "PRODUCTION_LAW", false );
  
  switch(optPhysWeathLaw)
  {
#define X(a,b) case a:	   \
physWeath = new b(infile);		\
break;
      PRODUCTION_LAW_TABLE2
#undef X
    default:
    {
      std::cerr << "\nError: You requested the soil production law '"
		  << optPhysWeathLaw << "' which does not exist."  << std::endl
		  << "Available options:" << std::endl;
      for(int i=0; i!=NUMBER_OF_PRODUCTION_LAWS; ++i ){
        std::cerr << " " << i << ": " << ProductionLaw[i] << std::endl;
      }
      ReportFatalError( "Requested soil production law not available. "
                       "Switch options.\n" );
    }
  }
  std::cout << "SOIL PRODUCTION OPTION: "
	    << ProductionLaw[optPhysWeathLaw] << std::endl;
  if( optPhysWeathLaw > 0 )
  {
    soilBulkDensity=infile.ReadDouble( "SOILBULKDENSITY", false );
    if( soilBulkDensity<=0.0 ) soilBulkDensity = kDefaultSoilBulkDensity;
  }

  // set chemical weathering law:
  optChemWeathLaw = infile.ReadItem( optChemWeathLaw,
                                  "CHEM_WEATHERING_LAW", false );
  
  switch(optChemWeathLaw)
  {
#define X(a,b) case a:				\
chemWeath = new b(infile, mptr);		\
break;
      CHEM_WEATHERING_TABLE2
#undef X
    default:
    {
      std::cerr << "\nError: You requested the chemical weathering law '"
		  << optChemWeathLaw << "' which does not exist."  << std::endl
		  << "Available options:" << std::endl;
      for(int i=0; i!=NUMBER_OF_CHEM_WEATHERINGS; ++i ){
        std::cerr << " " << i << ": " << ChemWeathering[i] << std::endl;
      }
      ReportFatalError( "Requested chemical weathering law not available. "
                       "Switch options.\n" );
    }
  }
  std::cout << "CHEMICAL WEATHERING OPTION: "
	    << ChemWeathering[optChemWeathLaw] << std::endl;

  // Landsliding:
  bool optLandslides = infile.ReadBool( "OPT_LANDSLIDES", false );
  if( optLandslides )
  {
    int tmp_ = infile.ReadItem( tmp_, "FLOWGEN" );
    tStreamNet::kFlowGen_t flowFlag = static_cast<tStreamNet::kFlowGen_t>(tmp_);
    if( flowFlag != tStreamNet::kSaturatedFlow1 
       && flowFlag != tStreamNet::kSubSurf2DKinematicWave )
    {
      std::cerr << "\nError: You selected the option for landslides, which "
      << "requires you also select FLOWGEN = 1 for 'type 1' "
      << "saturated flow hydrology or FLOWGEN = 6 for " 
      << "subsurface kinematic wave routing" << std::endl;
      ReportFatalError( "Requested landsliding and incompatible hydrology.\n" );
    }
    rockBulkDensity = infile.ReadDouble( "ROCKDENSITYINIT", false );
    if( rockBulkDensity<=0.0 ) rockBulkDensity = kDefaultRockBulkDensity;
    wetBulkDensity = 
      soilBulkDensity + RHO * ( 1.0 - soilBulkDensity / rockBulkDensity );
    double tmp_double = infile.ReadDouble( "WOODDENSITY", false );
    if( tmp_double > 0.0 )
      woodDensity = tmp_double;
    tmp_double = infile.ReadDouble( "FRICSLOPE", false );
    if( tmp_double > 0.0 )
      fricSlope = tmp_double;
    
    // set runout rules:
    optDF_RunOutLaw = infile.ReadItem( optDF_RunOutLaw, "DF_RUNOUT_RULE" );
    switch(optDF_RunOutLaw)
    {
#define X(a,b) case a: \
runout = new b(infile); \
break;
        DF_RUNOUT_TABLE2
#undef X
      default:
      {
        std::cerr << "\nError: You requested the debris flow runout rule '"
        << optDF_RunOutLaw << "' which does not exist."  << std::endl
        << "Available options:" << std::endl;
        for(int i=0; i!=NUMBER_OF_DF_RUNOUTS; ++i ){
          std::cerr << " " << i << ": " << DebrisFlowRunOut[i] << std::endl;
        }
        ReportFatalError( "Requested debris flow runout rule not available. "
                         "Switch options.\n" );
      }
    }
    // set scour rules:
    optDF_ScourLaw = infile.ReadItem( optDF_ScourLaw, "DF_SCOUR_RULE" );
    switch(optDF_ScourLaw)
    {
#define X(a,b) case a: \
scour = new b(infile); \
break;
        DF_SCOUR_TABLE2
#undef X
      default:
      {
        std::cerr << "\nError: You requested the debris flow scour rule '"
        << optDF_ScourLaw << "' which does not exist."  << std::endl
        << "Available options:" << std::endl;
        for(int i=0; i!=NUMBER_OF_DF_SCOURS; ++i ){
          std::cerr << " " << i << ": " << DebrisFlowScour[i] << std::endl;
        }
        ReportFatalError( "Requested debris flow scour rule not available. "
                         "Switch options.\n" );
      }
    }
    // set deposition rules:
    optDF_DepositLaw = infile.ReadItem( optDF_DepositLaw, "DF_DEPOSITION_RULE" );
    switch(optDF_DepositLaw)
    {
#define X(a,b) case a: \
deposit = new b(infile); \
break;
        DF_DEPOSITION_TABLE2
#undef X
	default:
	  {
	    std::cerr << "\nError: You requested the debris flow deposition rule '"
		      << optDF_DepositLaw << "' which does not exist."  << std::endl
		      << "Available options:" << std::endl;
	    for(int i=0; i!=NUMBER_OF_DF_DEPOSITIONS; ++i ){
	      std::cerr << " " << i << ": " << DebrisFlowDeposit[i] << std::endl;
	    }
	    ReportFatalError( "Requested debris flow deposition rule not available. "
			      "Switch options.\n" );
	  }
	}
      // set DF_fsPtr, DF_Hyd_fsPtr by making new std::ofstreams
      // If landslides, create files for writing them
      if( !no_write_mode )
	{
	  DF_fsPtr = new std::ofstream();
	  {
	    char fname[87];
#define THEEXT ".dflow"
	    infile.ReadItem( fname, sizeof(fname)-sizeof(THEEXT), "OUTFILENAME" );
	    strcat( fname, THEEXT );
#undef THEEXT
	    DF_fsPtr->open( fname );
	    if( !DF_fsPtr->good() )
	      std::cerr << "Warning: unable to create debris flow data file '"
			<< fname << "'\n";
	  }
	  DF_Hyd_fsPtr = new std::ofstream();
	  {
	    char fname[87];
#define THEEXT ".dfhyd"
	    infile.ReadItem( fname, sizeof(fname)-sizeof(THEEXT), "OUTFILENAME" );
	    strcat( fname, THEEXT );
#undef THEEXT
	    DF_Hyd_fsPtr->open( fname );
	    if( !DF_Hyd_fsPtr->good() )
	      std::cerr << "Warning: unable to create debris flow hydrograph file '"
			<< fname << "'\n";
	  }
	}
    }
}


// copy constructor for copying to new mesh:
tErosion::tErosion( const tErosion& orig, tMesh<tLNode>* Ptr )
  : meshPtr(Ptr), bedErode(0), sedTrans(0), physWeath(0), chemWeath(0), runout(0),
    scour(0), deposit(0), DF_fsPtr(0), DF_Hyd_fsPtr(0),
    kd(orig.kd),                 // Hillslope transport (diffusion) coef
    kd_ts(orig.kd_ts),
    difThresh(orig.difThresh),   // Diffusion occurs only at areas < difThresh
    mdMeshAdaptMaxFlux(orig.mdMeshAdaptMaxFlux), // For dynamic point addition: max ero flux rate
    mdSc(orig.mdSc),  // Threshold slope for nonlinear diffusion
    diffusionH(orig.diffusionH), // depth scale for depth-dependent diffusion
    beta(orig.beta), // proportion of sediment flux contributing to bedload
    track_sed_flux_at_nodes_(orig.track_sed_flux_at_nodes_), // option for tracking sed flux at nodes
    water_sed_tracker_ptr_(0),
    soilBulkDensity(orig.soilBulkDensity), // dry bulk density of soil, when made from rock (kg/m3)
    rockBulkDensity(orig.rockBulkDensity), // rock density (kg/m3)
    wetBulkDensity(orig.wetBulkDensity), // wet bulk density of soil (kg/m3)
    woodDensity(orig.woodDensity), // density of wood (kg/m3)
    fricSlope(orig.fricSlope), // tangent of angle of repose for soil (unitless)
    num_grain_sizes_(orig.num_grain_sizes_), // # grain size classes
    debris_flow_sed_bucket(0.0), // tally of debris flow sed. volume
    debris_flow_wood_bucket(0.0), // tally of debris flow wood volume
    landslideAreas(),
    optBedErosionLaw(orig.optBedErosionLaw),
    optSedTransLaw(orig.optSedTransLaw),
    optPhysWeathLaw(orig.optPhysWeathLaw),
    optChemWeathLaw(orig.optChemWeathLaw),
    optDF_RunOutLaw(orig.optDF_RunOutLaw),
    optDF_ScourLaw(orig.optDF_ScourLaw),
    optDF_DepositLaw(orig.optDF_DepositLaw)
{
  // pointers to objects governing rules for sediment transport:
  if( orig.bedErode )
    {
      switch(optBedErosionLaw)
	{
#define X(a,b) case a:				\
	  bedErode = new b();			\
	  break;
	  DETACHMENT_LAW_TABLE2
#undef X
	}
      bedErode->Initialize_Copy( orig.bedErode );
    }
  if( orig.sedTrans )
    {
      switch(optSedTransLaw)
	{
#define X(a,b) case a:				\
	  sedTrans = new b();			\
	  break;
	  TRANSPORT_LAW_TABLE2
#undef X
	}
      sedTrans->Initialize_Copy( orig.sedTrans );
    }
  if( orig.physWeath )
    {
      switch(optPhysWeathLaw)
	{
#define X(a,b) case a:				\
	  physWeath = new b();			\
	  break;
	  PRODUCTION_LAW_TABLE2
#undef X
	}
      physWeath->Initialize_Copy( orig.physWeath );
    }
  if( orig.chemWeath )
    {
      switch(optChemWeathLaw)
	{
#define X(a,b) case a:				\
	  chemWeath = new b();			\
	  break;
	  CHEM_WEATHERING_TABLE2
#undef X
	}
      chemWeath->Initialize_Copy( orig.chemWeath, meshPtr ); // chemical weathering object
    }
  if( orig.runout )
    {
      switch(optDF_RunOutLaw)
	{
#define X(a,b) case a:				\
	  runout = new b();			\
	  break;
	  DF_RUNOUT_TABLE2
#undef X
	}
      runout->Initialize_Copy( orig.runout );
    }
  if( orig.scour )
    {
      switch(optDF_ScourLaw)
	{
#define X(a,b) case a:				\
	  scour = new b();			\
	  break;
	  DF_SCOUR_TABLE2
#undef X
	}
      scour->Initialize_Copy( orig.scour );
    }
  if( orig.deposit )
    {
      switch(optDF_DepositLaw)
	{
#define X(a,b) case a:				\
	  deposit = new b();			\
	  break;
	  DF_DEPOSITION_TABLE2
#undef X
	}
      deposit->Initialize_Copy( orig.deposit );
    }
  //   if( orig.DF_fsPtr )
  //     DF_fsPtr = new std::ofstream( *orig.DF_fsPtr ); // pointer to output stream for debris flows
  //   if( orig.DF_Hyd_fsPtr )
  //     DF_Hyd_fsPtr = new std::ofstream( *orig.DF_Hyd_fsPtr ); // pointer to output stream for debris flow tally

  if( orig.water_sed_tracker_ptr_ )
    water_sed_tracker_ptr_ = 
      new tWaterSedTracker( *orig.water_sed_tracker_ptr_, meshPtr );  // -> water&sed tracker object

}

tErosion::~tErosion(){
  meshPtr = 0;
  delete bedErode;
  delete sedTrans;
  delete physWeath;
  delete chemWeath;
  if( runout > 0 ) delete runout;
  if( scour > 0 ) delete scour;
  if( deposit > 0 ) delete deposit;
  if( DF_fsPtr > 0 ) delete DF_fsPtr;
  if( DF_Hyd_fsPtr > 0 ) delete DF_Hyd_fsPtr;
}

/**************************************************************************\
**
**  tErosion::TurnOnOutput, TurnOffOutput
**
**  Open output file so output will be directed to it, 
**  or close output file so there won't be output.
**
**  SL, 10/2010
**
\**************************************************************************/
void tErosion::TurnOnOutput( const tInputFile& infile )
{
  if( !DF_fsPtr )
    {
      DF_fsPtr = new std::ofstream();
      {
	char fname[87];
#define THEEXT ".dflow"
      infile.ReadItem( fname, sizeof(fname)-sizeof(THEEXT), "OUTFILENAME" );
      strcat( fname, THEEXT );
#undef THEEXT
      DF_fsPtr->open( fname );
      if( !DF_fsPtr->good() )
        std::cerr << "Warning: unable to create debris flow data file '"
		    << fname << "'\n";
    }
    DF_Hyd_fsPtr = new std::ofstream();
    {
      char fname[87];
#define THEEXT ".dfhyd"
      infile.ReadItem( fname, sizeof(fname)-sizeof(THEEXT), "OUTFILENAME" );
      strcat( fname, THEEXT );
#undef THEEXT
      DF_Hyd_fsPtr->open( fname );
      if( !DF_Hyd_fsPtr->good() )
        std::cerr << "Warning: unable to create debris flow hydrograph file '"
		    << fname << "'\n";
    }
  }
}

void tErosion::TurnOffOutput()
{
  if( DF_fsPtr )
    {
      DF_fsPtr->close();
      delete DF_fsPtr;
      DF_fsPtr = NULL;
      DF_Hyd_fsPtr->close();
      delete DF_Hyd_fsPtr;
      DF_Hyd_fsPtr = NULL;
    }
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
  if(0) //DEBUG
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
    tArray <double> inletBedSizeFraction( strmNet->getInletSedSizeFraction() );  // TEMP 6/06: stores desired bed sed proportions at inlet
    // fractions must sum to 1 in input file (INSED1, INSED2, etc)
    double inletSlope;	
	    
    //DEBUGGING 
    if(0) {
      std::cout<<"inletSlope = "<< inletSlope <<std::endl;
      for( size_t i=0; i<cn->getNumg(); i++ )
        std::cout<<"sedfrac "<<i<<"="<<inletBedSizeFraction[i]<<std::endl;
    }
    
    // New stuff in progress for dynamic calculation of sed influx at inlet, 5/06
    // Here's what we need: modify tStreamNet to add a function that returns a ref or ptr to the inlet.
    // Use this to access sed influx info for the inlet.
    // Modify code to set erodibility of inlet node to zero, and compute sed influx before loop using call to 
    // TransCapacity. Assign these fluxes to insed ... etc.
    
    // Sort so that we always work in upstream to downstream order
    strmNet->SortNodesByNetOrder();
    strmNet->FindChanGeom();
    strmNet->FindHydrGeom();
    
    // Compute erosion and/or deposition until all of the elapsed time (dtg)
    // is used up
    do
    {
      if(0) std::cout << "DetachErode: top of do loop\n" << std::flush;
      
      // Zero out sed influx of all sizes
      for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      {
        if(0 && cn==inletNode ) std::cout<<"top loop ID="<<cn->getID()<<std::endl;
        cn->setQs(0.0);
        if( cn!=inletNode )
	      {
          cn->setQsin(0.0); //totals are for ts calculation
          cn->setQsin( sedzero );
          for( size_t i=0; i<cn->getNumg(); i++ ){
            cn->setQs(i,0.0);
          }
	      }
        else  // TEMPORARY MODIFICATIONS FOR TEST, 5/06:  
          // AND SET SLOPE TO A FIXED VALUE
	      {
          // We set the inlet node's erodibility values (for each layer) to zero so it can't be eroded.
          // Also set the grain-size distribution in the upper layers to the values specified in the input 		// file by INSED1, 2, etc. Variable inletBedSizeFraction contains INSED1, 2, etc, which are
          // assumed to be fractions that sum to 1.0 (the user can screw this up ... it isn't checked!)	
        double inletSlope = strmNet->getInletSlope( time );
		    if(0) {
                std::cout<<"inletSlope = "<< inletSlope <<std::endl;
                   }
		  size_t numLayersInlet = cn->getNumLayer();
          if(0) std::cout<<numLayersInlet<<" lay inlt\n";
          for( size_t i=0; i<numLayersInlet; i++ ) {
            cn->setLayerErody( i, 0.0 );
            double layThick = cn->getLayerDepth(i);
            for( size_t j=0; j<cn->getNumg(); j++ ) {
              if(0) 
                std::cout<<"set lay "<<i<<", with thickness " << layThick 
                <<", size "<<j
                <<" to "<< layThick*inletBedSizeFraction[j] << std::endl;
              cn->setLayerDgrade(i,j,layThick*inletBedSizeFraction[j] );
            }
          }
          
          // zero out Qs for each size class
          for( size_t i=0; i<cn->getNumg(); i++ ) {
            cn->setQs(i,0.0);  
          }
          
          // Now we adjust the elevation of the inlet node so that 
          // it has the user-defined slope
          double zdown = cn->getDownstrmNbr()->getZ(); //TEMP TEST
          double len = cn->getFlowEdg()->getLength();   // TEMP TEST
          
          //Xdouble temporary_myslope = 0.05;  // Ultimately, read this from input file
          cn->ChangeZ( (zdown+len * strmNet->getInletSlope( time ) )-cn->getZ() );
          
          // Next, we call TransCapacity, which automatically sets Qs in each size class
          insedloadtotal = sedTrans->TransCapacity( cn, 0, 1.0 );
          if(0) std::cout<<"inlet capacity="<<insedloadtotal<<std::endl;
          
          // Store Qs for each size class in the "insed" array so we can 
          // assign these to Qsin
          for( size_t i=0; i<cn->getNumg(); i++ ) {
            insed[i] = cn->getQs(i);   // Capacity for i-th size fraction
            if(0) std::cout<<" insed["<<i<<"]="<<insed[i]<<std::endl;
          }
          
          // Now, we set the influxes at the inlet node, both total and per-size, 
          //to the capacity values we just calculated and stored
          cn->setQsin( insedloadtotal ); // here's the total influx
          cn->setQsin( insed );  // ... and the per-size influx
          
          //double zdown = cn->getDownstrmNbr()->getZ(); //TEMP TEST
          //double len = cn->getFlowEdg()->getLength();   // TEMP TEST
          //double myslope = 0.025; //TEMP TEST
          //tArray <double> testdz(cn->getNumg() );  //TEMP TEST
          //double testdztotal = zdown+myslope*len - cn->getZ(); //TEMp TEST
          //if( 1 ) std::cout<<"adj inlt "<<testdztotal;
          //for( size_t kk=0; kk<cn->getNumg(); kk++ ) //TEMP TEST
          //{
          //  testdz[kk] = testdztotal*(insed[kk]/insedloadtotal ); //TEMP TEST
          //std::cout<<" kk="<<kk<< "testdz="<<testdz[kk];
          //}
          //if( 1 ) std::cout<<std::endl;
          //cn->EroDep( 0, testdz, timegb );  //TEMP TEST
          //if( 1 ) std::cout<<"nl="<<cn->getNumLayer()<<" thick="<<cn->getLayerDepth(0)<<std::endl;
          //cn->setLayerDepth( nl, 100000.0 ); //TEMP TEST
          //cn->setLayerErody( 0, 1e6 ); //TEMP TEST
          //cn->setLayerErody( 1, 1e6 ); //TEMP TEST
          //if(1) std::cout << "inletnode elev " << cn->getZ() << " dsnbr " << zdown << " len " << len << " slp " << (cn->getZ()-zdown)/len << std::endl;
          //for( size_t i=0; i<cn->getNumg(); i++ ){
          //cn->setQs(i,0.0);
          //cn->setLayerDgrade(0,i,cn->getLayerDepth(0)*(insed[i]/insedloadtotal) ); //TEMP TEST
          //if( cn->getNumLayer()>1) cn->setLayerDgrade(1,i,cn->getLayerDepth(1)*(insed[i]/insedloadtotal) ); //TEMP TEST 
          //std::cout << "inlet size " << i << "=" << cn->getLayerDgrade(0,i) << std::endl;
          //}
	      }
      }
      
      // Estimate erosion rates and time-step size
      // NOTE - in this first loop we are only dealing with
      // totals for time-step calculations, however transport
      // rates for each size are also set within the function call.
      if(0) std::cout << "DetachErode: estimating rates\n" << std::flush;
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
            qs += 
            sedTrans->TransCapacity(cn,i,cn->getLayerDepth(i)
                                    /cn->getChanDepth());
            if(0&&cn==inletNode) 
              std::cout<<"1depck="<<depck<<" qs="<<qs
              <<"wt="<<cn->getLayerDepth(i)/cn->getChanDepth()
              <<" qs/wt="<<qs/(cn->getLayerDepth(i)/cn->getChanDepth())
              <<std::endl;
          }
          else{
            qs += sedTrans->TransCapacity(cn,i,1-(depck/cn->getChanDepth()));
            if(0&&cn==inletNode) 
              std::cout<<"2depck="<<depck<<" qs="<<qs
              <<" wt="<< 1-(depck/cn->getChanDepth())
              << " qs/wt="<<qs/(depck/cn->getChanDepth())<<std::endl;
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
        
        //if( cn==inletNode ) drdt = -1e6;  // TEMP TEST
        
        cn->setDrDt(drdt);
        cn->setDzDt(drdt);
        
        excap=(qs - cn->getQsin())/cn->getVArea();//[m/yr]
        //excap negative = deposition; positive = erosion
        //Note that signs are opposite to what one
        //might expect.  This works out for Qsin addition.
        //Limit erosion to capacity of flow or deposition
        if( -drdt > excap ){
          cn->setDzDt(-excap);
        }
        cn->getDownstrmNbr()->addQsin(cn->getQsin()-cn->getDzDt()*cn->getVArea());
        
        //std::cout << "*** EROSION ***\n";
        if( 0 && cn==inletNode ) {
          std::cout << "Trans Cap inlet = " << qs << "excap=" << excap 
          << " drdt=" << drdt<< "DzDt=" << cn->getDzDt() << std::endl;
          //cn->TellAll();
        }
        
      }//ends for( cn = ni.FirstP...
      
      //Find local time-step based on dzdt
      if(0) std::cout << "DetachErode: finding time step size\n" << std::flush;
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
          if(0) {
            std::cout<<"Inlet qsin set to:\n";
            for( size_t i=0; i<cn->getNumg(); i++ )
              std::cout<< " "<<i<<"="<<insed[i]<<std::endl;
          }
	      }
        
        dn = cn->getDownstrmNbr();
        ratediff = dn->getDzDt() - cn->getDzDt(); //Are the pts converging?
        if( ratediff > 0. && (cn->calcSlope()) > 1e-7 )  // if yes, get time
	      {                                              //  to zero slope
          if(0) {
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
              if(0) { // debug
                std::cout << "Very small dtmax " << dtmax <<  std::endl;
                std::cout << "rate dif is " << ratediff << std::endl;
                std::cout << "elev dif is " << cn->getZ()-dn->getZ() << std::endl;
                std::cout << "dzdt upstream is " << cn->getDzDt() << std::endl;
                std::cout << "dzdt downstream is " << dn->getDzDt() << std::endl;
                cn->TellAll();
                dn->TellAll();
              }
              dtmax=0.0001; //GREG I added this just because things
              // were taking forever.  I kept it for now just for
              // testing stuff.  Maybe we should discuss this.
            }
          }
	      }
      }// End for( cn = ni.FirstP()..
      dtmax *= frac;  // Take a fraction of time-to-flattening
      timegb+=dtmax;
      
      //At this point: we have drdt and qs for each node, plus dtmax
      
      // Do erosion/deposition
      if(0) std::cout << "DetachErode: eroding\n" << std::flush;
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
        if( 0 && cn==inletNode ) {
          std::cout << "dz inlet = " << dz << " dz/dt=" << dz/dtmax << std::endl;
          //cn->TellAll();
        }
        
        if( dz<0 ) //total erosion
	      {
          if(!flag){ // detach-lim
            if(0 && cn==inletNode) std::cout << "dlim\n" << std::endl;
            int i=0;
            depck=0.;
            while(dz<-0.000000001&&depck<cn->getChanDepth()&&i<cn->getNumLayer()){
              depck+=cn->getLayerDepth(i);
              if(-dz<=cn->getLayerDepth(i)){//top layer can supply total depth
                for(size_t j=0;j<cn->getNumg();j++){
                  // Figure out how much of size j is liberated by erosion to depth dz
                  erolist[j]=dz*cn->getLayerDgrade(i,j)/cn->getLayerDepth(i);
                  // Check whether there's enough extra capacity to carry this much of size j
                  if(erolist[j]<(cn->getQsin(j)-cn->getQs(j))*dtmax/cn->getVArea()){
                    //decrease total dz because of capacity limitations
                    erolist[j]=(cn->getQsin(j)-cn->getQs(j))*dtmax/cn->getVArea();
                    cn->setQsin(j,0.0); // ??
                    cn->setQs(j,0.0);   // ??
                  }
                }
                if( 0 && cn==inletNode ) std::cout<<"NO ero "<<dz<<" from lyr "<<i<<std::endl;
                if( cn!=inletNode )  //TEMP 6/06
                { 
                  ret=cn->EroDep(i,erolist,timegb); //ORIGINAL
                  for(size_t j=0;j<cn->getNumg();j++){ //ORIGINAL
                    cn->getDownstrmNbr()->addQsin(j,-ret[j]*cn->getVArea()/dtmax); //ORIGINAL
                  } //ORIGINAL
                } //TEMP 6/06
                dz=0.;
              }
              else{//top layer is not deep enough, need to erode more layers
                flag=false;
                for(size_t j=0;j<cn->getNumg();j++){
                  erolist[j]=-cn->getLayerDgrade(i,j);
                  if(erolist[j]<(cn->getQsin(j)-cn->getQs(j))*dtmax/cn->getVArea()){
                    //decrease total dz because of capacity limitations
                    erolist[j]=(cn->getQsin(j)-cn->getQs(j))*dtmax/cn->getVArea();
                    cn->setQsin(j,0.0); // ??
                    cn->setQs(j,0.0);   // ??
                    //need to set these to zero since the capacity has
                    //now been filled by the stuff in this layer
                    flag=true;
                    //Since not taking all of the material from the
                    //surface, surface layer won't be removed-must inc i
                  }
                  dz-=erolist[j];
                }
                if( 0 && cn==inletNode ) std::cout<<"NO Ero "<<erolist[0]<<"+"<<erolist[1]<<"="<<erolist[0]+erolist[1]<<" from lyr "<<i<<std::endl;
                if( cn!=inletNode ) //TEMP 6/06
                {
                  ret=cn->EroDep(i,erolist,timegb);
                  for(size_t j=0;j<cn->getNumg();j++){
                    //if * operator was overloaded for arrays, no loop necessary
                    cn->getDownstrmNbr()->addQsin(j,-ret[j]*cn->getVArea()/dtmax);
                  }
                }
                if(flag){
                  i++;
                }
              }
            }
          }
          else{//trans-lim
            if( 0 && cn==inletNode ) std::cout<<"Inlet X "<<cn->getX()<<" Y "<<cn->getY() <<" tlim\n";
            for(size_t j=0;j<cn->getNumg();j++){
              erolist[j]=(cn->getQsin(j)-cn->getQs(j))*dtmax/cn->getVArea();
              if( 0 && cn==inletNode ) std::cout<<" j "<<j<<" "<<erolist[j];
            }
            if( 0 && cn==inletNode ) std::cout<<"."<<std::endl;
            
            int i=0;
            depck=0.;
            while(depck<cn->getChanDepth()){
              depck+=cn->getLayerDepth(i);
              int flag=cn->getNumLayer();
              if( 0 && cn==inletNode ) std::cout<<"NO depck="<<depck<<" numLayer="<<flag<<" i="<<i<<std::endl;
              if( cn!=inletNode)  // JUNE 06 TEMP HACK: DON"T ERODE INLET!
              {
                ret=cn->EroDep(i,erolist,timegb);
                //if( 1 && cn==inletNode ) std::cout<<"ret0="<<ret[0]<<" ret1="<<ret[1]<<std::endl;
                double sum=0.;
                for(size_t j=0;j<cn->getNumg();j++){
                  cn->getDownstrmNbr()->addQsin(j,-ret[j]*cn->getVArea()/dtmax);
                  erolist[j]-=ret[j];
                  sum+=erolist[j];
                }
                if( 0 && cn==inletNode ) std::cout<<"end for loop"<<std::endl;
                if(sum>-0.0000001)
                  depck=cn->getChanDepth();
                if(flag==cn->getNumLayer())
                  i++;
              } // END TEMP HACK BRACKETS (INTERIOR IS ORIGINAL)
              if( 0 && cn==inletNode ) std::cout<<"end while loop"<<std::endl;
            } //end while
          }//end if( trans-limited )
	      }//ends(if dz<0)
        else if(dz>0) //total deposition -> need if cause erodep chokes with 0
	      {
          //Get texture of stuff to be deposited
          for(size_t j=0;j<cn->getNumg();j++)
            erolist[j]=(cn->getQsin(j)-cn->getQs(j))*dtmax/cn->getVArea();
          if(0 && cn==inletNode ) std::cout<<"NOT about to erodep inlet\n";
          if( cn!=inletNode ) //CLAUSE ADDED TEMP 6/06 (INTERIOR IS ORIGINAL)
          {
            ret=cn->EroDep(0,erolist,timegb);
            for(size_t j=0;j<cn->getNumg();j++){
              cn->getDownstrmNbr()->addQsin(j,-ret[j]*cn->getVArea()/dtmax);
            }
          }
	      }
        
        if( 0 && cn==inletNode ) std::cout<<"end of node FOR loop\n";
        
      } // Ends for( cn = ni.FirstP()...
      
      if( track_sed_flux_at_nodes_ )
      {
        if(0) std::cout << "WE'RE GOIN ALL THE WAY" << endl;
        water_sed_tracker_ptr_->AddSedVolumesAtTrackingNodes( dtmax );
      }
      else
        if(0) std::cout << "NO WAY JOSE!" << endl;
      
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
  
  
  if(0) std::cout<<"ending detach erode\n"<<std::flush;
  
}// End erosion algorithm

/***********************************************************************\
 **
 **  tErosion::DetachErode2
 **
 **  DetachErode2 is an adaption of detachErode that 
 **  includes the diffusive load in the fluvial
 **  load which applies when using tools and cover type rules, or using 
 **  the GeneralFQS erosion law.  Because the diffusive law gets included in the
 **  sediment flux, this could mess up multiple-grain size erosion, because
 **  the diffusive load doesn't properly account for grain size.  So for now
 **  this version of detachErode is seperate.  Better testing and algorithms for
 **  diffusion of multiple grain size should be developed before using with 
 **  multiple grain-sizes.  (NMG June 2005)
 **
 **  Now beta is included in calculating the transport-limited rate of erosion.
 **  That is, beta *(Qsin+Qsdin) is the bedload flux in.  This assumes that
 **  Qs only contains bedload, which is reasonable for all coded transport
 **  algorithms at this point.  (NMG Nov 2006)
 **
 \************************************************************************/
void tErosion::DetachErode2(double dtg, tStreamNet *strmNet, double time,
                            tVegetation * /*pVegetation*/ )
{
  
  //std::cout<<"welcome to detacherode2"<<endl;
  
  //Added 4/00, if there is no runoff, this would crash, so check
  if(strmNet->getRainRate()-strmNet->getInfilt()>0){
    
    double dtmax;       // time increment: initialize to arbitrary large val
    double frac = 0.3;  //fraction of time to zero slope
    double timegb=time; //time gone by - for layering time purposes
    int flag;
    tLNode * cn, *dn;
    // int nActNodes = meshPtr->getNodeList()->getActiveSize();
    tMesh< tLNode >::nodeListIter_t ni( meshPtr->getNodeList() );
    double ratediff,  // Difference in ero/dep rate btwn node & its downstrm nbr
    drdt,
    dz,
    depck,
    qs=0,
    excap;
    tLNode * inletNode = strmNet->getInletNodePtrNC();
    double insedloadtotal = strmNet->getInSedLoad();
    
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
          for(  size_t i=0; i<cn->getNumg(); i++){
            cn->setQs(i,0.0);
          }
        }
        else
        {
          cn->setQsin(insedloadtotal); //totals are for ts calculation
          cn->setQsin( insed );
          for(  size_t i=0; i<cn->getNumg(); i++ ){
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
        depck=0;
        int i=0;
        qs=0;
        
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
        
        if(depck>cn->getChanDepth()) //which layer are you basing detach on?
          drdt=-bedErode->DetachCapacity( cn, i-1 );
        else
          drdt=-bedErode->DetachCapacity( cn, i );//[m/yr]
        
        cn->setDrDt(drdt);
        cn->setDzDt(drdt);
        
        //nmg now including Qsdin
        //excap=(qs - cn->getQsin())/cn->getVArea();//[m/yr]
        excap=(qs - beta*(cn->getQsin()+cn->getQsdin()))/cn->getVArea();//[m/yr]
        
        //excap negative = deposition; positive = erosion
        //Note that signs are opposite to what one
        //might expect.  This works out for Qsin addition.
        //Limit erosion to capacity of flow or deposition
        if( -drdt > excap ){
          cn->setDzDt(-excap);
        }
        //adding Qsdin for the sed-flux rules
        //Don't use beta here because the tracking the total load.
        //cn->getDownstrmNbr()->addQsin(cn->getQsin()-cn->getDzDt()*cn->getVArea());
        cn->getDownstrmNbr()->addQsin(cn->getQsin()+cn->getQsdin()-cn->getDzDt()*cn->getVArea());
        
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
          if(0) {
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
              if(0) { // debug
                std::cout << "Very small dtmax " << dtmax <<  std::endl;
                std::cout << "rate dif is " << ratediff << std::endl;
                std::cout << "elev dif is " << cn->getZ()-dn->getZ() << std::endl;
                std::cout << "dzdt upstream is " << cn->getDzDt() << std::endl;
                std::cout << "dzdt downstream is " << dn->getDzDt() << std::endl;
                cn->TellAll();
                dn->TellAll();
              }
              dtmax=0.0001; //GREG I added this just because things
              // were taking forever.  I kept it for now just for
              // testing stuff.  Maybe we should discuss this.
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
        //excap=(cn->getQs() - cn->getQsin())/cn->getVArea();
        excap=(cn->getQs() - beta*(cn->getQsin()+cn->getQsdin()))/cn->getVArea();//[m/yr]
        //excap pos if eroding, neg if depositing
        drdt=-bedErode->DetachCapacity( cn, 0 );
        cn->setDrDt(drdt);
        cn->setDzDt(drdt);
        
        if( -cn->getDrDt() < excap ){
          dz = cn->getDrDt()*dtmax; // detach-lim
          flag = 0;
          //std::cout<<"D-lim A: "<<cn->getQ()<<" dz/dt: "<<cn->getDzDt()<<endl;
        }
        else{
          dz = -excap*dtmax; // trans-lim
          flag = 1;
          cn->setDzDt(-excap);
          //std::cout<<"T-lim A: "<<cn->getQ()<<" dz/dt: "<<cn->getDzDt()<<endl;
        }
        
        if( dz<0 ) //total erosion
        {
          if(flag==0){ // detach-lim
            //cout<<"detach-lim varea = "<<cn->getVArea();
            //cout<<" X "<<cn->getX()<<" Y "<<cn->getY()<<endl;
            //total erosion so send all material upstream down.
            for(size_t j=0; j<cn->getNumg(); j++)
              cn->getDownstrmNbr()->addQsin(j,cn->getQsin(j));
            int i=0;
            depck=0;
            //while( sig. amt left to erod && eroding < chan depth && have layers left )
            while(dz<-0.000000001 && depck<cn->getChanDepth() && i<cn->getNumLayer()){
              depck+=cn->getLayerDepth(i);
              if(-dz<=cn->getLayerDepth(i)){//top layer can supply total depth
                for(size_t j=0;j<cn->getNumg();j++){
                  //below distributes erosion based on prop of each grain size
                  erolist[j]=dz*cn->getLayerDgrade(i,j)/cn->getLayerDepth(i);
                  if(erolist[j]<(beta*(cn->getQsin(j)+cn->getQsdin(j))-cn->getQs(j))*dtmax/cn->getVArea()){
                    //decrease total dz because of capacity limitations for grain size
                    //SHOULD ONLY ENTER IF MULTIPLE GRAIN SIZES ARE BEING USED
                    erolist[j]=(beta*(cn->getQsin(j)+cn->getQsdin(j))-cn->getQs(j))*dtmax/cn->getVArea();
                    //WHY ARE Qsin and Qs SET TO ZERO?
                    cn->setQsin(j,0.0);
                    cn->setQs(j,0.0);
                  }
                }
                ret=cn->EroDep(i,erolist,timegb);
                for(size_t j=0;j<cn->getNumg();j++){
                  cn->getDownstrmNbr()->addQsin(j,-ret[j]*cn->getVArea()/dtmax);
                }
                dz=0;
              }
              else{//top layer is not deep enough, need to erode more layers
                flag=0;
                for(size_t j=0;j<cn->getNumg();j++){
                  erolist[j]=-cn->getLayerDgrade(i,j);
                  if(erolist[j]<(beta*(cn->getQsin(j)+cn->getQsdin(j))-cn->getQs(j))*dtmax/cn->getVArea()){
                    
                    //decrease total dz because of capacity limitations
                    erolist[j]=(beta*(cn->getQsin(j)+cn->getQsdin(j))-cn->getQs(j))*dtmax/cn->getVArea();
                    cn->setQsin(j,0.0);
                    cn->setQs(j,0.0);
                    //WHY ARE Qsin and Qs SET TO ZERO?
                    //need to set these to zero since the capacity has
                    //now been filled by the stuff in this layer
                    flag=1;
                    //Since not taking all of the material from the
                    //surface, surface layer won't be removed-must inc i
                  }
                  dz-=erolist[j];
                }
                ret=cn->EroDep(i,erolist,timegb);
                for(size_t j=0;j<cn->getNumg();j++){
                  //if * operator was overloaded for arrays, no loop necessary
                  cn->getDownstrmNbr()->addQsin(j,-ret[j]*cn->getVArea()/dtmax);
                }
                if(flag==1){
                  i++;
                }
              }
            }
            
          }
          else{//trans-lim
            //std::cout<<"trans-lim, varea = "<<cn->getVArea();
            //std::cout<<" X "<<cn->getX()<<" Y "<<cn->getY()<<endl;
            for(size_t j=0;j<cn->getNumg();j++){
              erolist[j]=(cn->getQsin(j)-cn->getQs(j))*dtmax/cn->getVArea();
              //send all upstream material downstream, regardless of whether or not
              //there is deposition because if there is deposition, this material
              //will be removed from the downstream load below
              cn->getDownstrmNbr()->addQsin(j,cn->getQsin(j));
            }
            
            int i=0;
            depck=0;
            while(depck<cn->getChanDepth()){
              depck+=cn->getLayerDepth(i);
              flag=cn->getNumLayer();
              ret=cn->EroDep(i,erolist,timegb);
              double sum=0;
              for(size_t j=0;j<cn->getNumg();j++){
                //here you are sending downstream the amount that was eroded
                //or removing from the load the amount that was deposited
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
            erolist[j]=(beta*(cn->getQsin(j))-cn->getQs(j))*dtmax/cn->getVArea();
          ret=cn->EroDep(0,erolist,timegb);
          for(size_t j=0;j<cn->getNumg();j++){
            //send upstream material down, minus the amount that was deposited
            cn->getDownstrmNbr()->addQsin(j,cn->getQsin(j)-(ret[j]*cn->getVArea()/dtmax));
            //cn->getDownstrmNbr()->addQsin(j,-ret[j]*cn->getVArea()/dtmax);
          }   
        }
        
        cn->getDownstrmNbr()->addQsin(0,cn->getQsdin());
        //Above sends upstream diffusional material downstream.  
        //A bit of a cheat since it might be over capacity, but just assume that
        //diffusive amount is relatively small
        //also, no grain size accounting for diffusive material
      } // Ends for( cn = ni.FirstP()...
      
      // Erode vegetation
#if 0
#define NEWVEG 0
      if( pVegetation && NEWVEG ) pVegetation->ErodeVegetation( meshPtr, dtmax );
#undef NEWVEG
#endif
      
      // Update time remainig
      dtg -= dtmax;
      //cout<<"Time remaining now "<<dtg<<endl;
    } while( dtg>1e-6 );  //Keep going until we've used up the whole time intrvl
  }//end if rainrate-infilt>0
  
  
  //std::cout<<"ending detach erode 2"<<endl;
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
 **  Inputs:  rt -- time duration over which to compute diffusion
 **           noDepoFlag -- if true, material is only eroded, never
 **                             deposited
 **  Modifies:  node elevations (z); node Qsin (sed influx) is reset and
 **             used in the computations
 **  Notes:  as of 3/98, does not differentiate between rock and sediment
 **
 \*****************************************************************************/
//#define kVerySmall 1e-6
#define kEpsOver2 0.1
void tErosion::Diffuse( double rt, bool noDepoFlag, double time )
{
  tLNode * cn;
  tEdge * ce;
  double volout,  // Sediment volume output from a node (neg=input)
  delt,       // Max local step size
  dtmax;      // Max global step size (initially equal to total time rt)
  tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );
  tMesh< tLNode >::edgeListIter_t edgIter( meshPtr->getEdgeList() );
  
  // Fix for "diffusion doesn't update layers" bug GT 11/12. We assume
  // that for multi-sizes, we'll call DiffuseMultiSize instead
  static tArray<double> deposition_depth( 1 );
  
#ifdef TRACKFNS
  std::cout << "tErosion::Diffuse()" << std::endl;
#endif
	
  kd = kd_ts.calc( time );
  if(0) std::cout << "kd = " << kd << std::endl;
  
  if( kd==0 ) return;
  //initialize Qsd, which will record the total amount of diffused material
  //fluxing into a node for the entire time-step.
  //if Qsd is negative, then material was deposited in that node.
  for( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
    cn->setQsdin( 0. );
  
  // Compute maximum stable time-step size based on Courant condition
  // for FTCS (here used as an approximation).
  // (Note: for a fixed mesh, this calculation only needs to be done once;
  // performance could be improved by having this block only called if
  // mesh has changed since last time through)
  if(0) std::cout << "About to enter edge loop...\n" << std::flush;
  dtmax = rt;  // Initialize dtmax to total time rt
  for( ce=edgIter.FirstP(); edgIter.IsActive(); ce=edgIter.NextP() )
  {
    assert( ce!=0 );
    if( 0 ) //DEBUG
    {
      std::cout << "In Diffuse(), large vedglen detected: " << ce->getVEdgLen() << std::endl;
      ce->TellCoords();
    }
    //Xif( (denom=kd*ce->getVEdgLen() ) > kVerySmall )
    // Evaluate DT <= DX^2 / Kd
    assert( kd > 0.0 );
    delt = kEpsOver2 * ce->getLength()*ce->getLength() / kd;
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
      volout = kd*ce->CalcSlope()*ce->getVEdgLen()*dtmax;
      // Record outgoing flux from origin
      cn = static_cast<tLNode *>(ce->getOriginPtrNC());
      if( difThresh>0.0 && cn->getDrArea()>difThresh )
        volout=0;
      cn->addQsin( -volout );
      // Record incoming flux to dest'n
      cn = static_cast<tLNode *>(ce->getDestinationPtrNC());
      cn->addQsin( volout );
      
      if( 0 ) { //DEBUG
        std::cout << volout << " mass exch. from " << ce->getOriginPtr()->getID()
        << " to "
        << ce->getDestinationPtr()->getID()
        << " on slp " << ce->getSlope() << " ve " << ce->getVEdgLen()
        << " kd=" << kd << " dtmax=" << dtmax
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
      deposition_depth[0] = cn->getQsin() / cn->getVArea();
      cn->EroDep( 0, deposition_depth, time );  // add or subtract net flux/area    
      //cn->EroDep( cn->getQsin() / cn->getVArea() );  // add or subtract net flux/area   

	  
      cn->getDownstrmNbr()->addQsdin(-1 * cn->getQsin()/dtmax);  
      //this won't work if time steps are varying, because you are adding fluxes
      
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
    
    if(0) std::cout << "bottom of do loop in Diffuse()\n" << std::flush;
    
  } while( rt>0.0 );
  
  
}
#undef kEpsOver2



#define kEpsOver2 0.1
void tErosion::DiffuseMultiSize( double rt, bool noDepoFlag, double time )
{
  tLNode * cn;
  tLNode * dn;
  tLNode * hn;
  tEdge * ce;
  double hst=diffusionH;	//H_star in m
  double volout,  // Sediment volume output from a node (neg=input)
  delt,       // Max local step size
  dtmax;      // Max global step size (initially equal to total time rt)
  tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );
  tMesh< tLNode >::edgeListIter_t edgIter( meshPtr->getEdgeList() );
  int i;
  
  // Here we create arrays to handle flux and deposition in the 
  // various size classes
  static tArray<double> deposition_depth( num_grain_sizes_ );
  static tArray<double> volout_by_size( num_grain_sizes_ );
  
  if (0) std::cout << "tErosion::DiffuseMultiSize()" << std::endl;
	
  kd = kd_ts.calc( time );
  if(0) std::cout << "kd = " << kd << std::endl;
  
  if( kd==0 ) return;
  //initialize Qsd, which will record the total amount of diffused material
  //fluxing into a node for the entire time-step.
  //if Qsd is negative, then material was deposited in that node.
  for( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
    cn->setQsdin( 0. );
  
  // Compute maximum stable time-step size based on Courant condition
  // for FTCS (here used as an approximation).
  // (Note: for a fixed mesh, this calculation only needs to be done once;
  // performance could be improved by having this block only called if
  // mesh has changed since last time through)
  if(0) std::cout << "About to enter edge loop...\n" << std::flush;
  dtmax = rt;  // Initialize dtmax to total time rt
  for( ce=edgIter.FirstP(); edgIter.IsActive(); ce=edgIter.NextP() )
  {
    assert( ce!=0 );
    if( 0 ) //DEBUG
    {
      std::cout << "In Diffuse(), large vedglen detected: " << ce->getVEdgLen() << std::endl;
      ce->TellCoords();
    }
    //Xif( (denom=kd*ce->getVEdgLen() ) > kVerySmall )
    // Evaluate DT <= DX^2 / Kd
    assert( kd > 0.0 );
    delt = kEpsOver2 * ce->getLength()*ce->getLength() / kd;
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
    {
      cn->setQsin( 0. );
      for( i=0; i<num_grain_sizes_; i++ )
        cn->setQsin( i, 0. );
    }
    
    // Compute sediment volume transfer along each edge
    for( ce=edgIter.FirstP(); edgIter.IsActive(); ce=edgIter.NextP() )
    {
	  cn = static_cast<tLNode *>(ce->getOriginPtrNC());  // get node
      dn = static_cast<tLNode *>(ce->getDestinationPtrNC()); 
	  
	  if (cn->getZ() > dn->getZ() )	  
		hn = cn;
	  else
		hn = dn;
		
		//DEBUG
		if (0){
		std::cout << " Height cn: " << cn->getZ() << " Height dn: " << dn->getZ() << " Height hn: " << hn->getZ() << std::endl;
		}
		
			  
      // Record outgoing flux from origin
      volout = kd*ce->CalcSlope()*ce->getVEdgLen()*dtmax * (1.0-exp((-1.0*hn->getRegolithDepth())/hst)); // volume out 
	  
      if( difThresh>0.0 && cn->getDrArea()>difThresh ) 
        volout=0;
     
	
	  for( int g=0; g<num_grain_sizes_; g++ ) // volume+flux by size
		{
			volout_by_size[g] = volout * ( hn->getLayerDgrade( 0, g ) / hn->getLayerDepth(0) );
			cn->addQsin( g, -volout_by_size[g] );
			dn->addQsin( g, volout_by_size[g] );
			
		 }
		 	if (0){
			std::cout  << " Regolith Depth hn: " << hn->getRegolithDepth() << std::endl;
			}
			
		 //DEBUG
		
		if( 0 ) { //DEBUG
        std::cout << volout << " mass exch. from " << ce->getOriginPtr()->getID()
        << " to "
        << ce->getDestinationPtr()->getID()
        << " on slp " << ce->getSlope() << " ve " << ce->getVEdgLen()
        << " kd=" << kd << " dtmax=" << dtmax
        << "\nvp " << ce->getRVtx().at(0)
        << " " << ce->getRVtx().at(1) << std::endl;
        int gg;
        for ( gg=0; gg<num_grain_sizes_; gg++ )
		{
          std::cout << "  size " << gg << " volout is " << volout_by_size[gg] << std::endl;
		  std::cout << "  size " << gg << " percent vol is " << volout_by_size[gg]/volout << std::endl;
		  }
        static_cast<tLNode *>(ce->getOriginPtrNC())->TellAll();
        static_cast<tLNode *>(ce->getDestinationPtrNC())->TellAll();
        std::cout << std::endl;
      }  
		
			 
      /*ce =*/ edgIter.NextP();  // Skip complementary edge
    }
    
    // Compute erosion/deposition for each node
    for( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
    {
    
      if( noDepoFlag && cn->getQsin() > 0.0 )
        cn->setQsin( 0.0 );
      for( i=0; i<num_grain_sizes_; i++ )
	  {
        deposition_depth[i] = cn->getQsin(i) / cn->getVArea();
		}
		
		
      cn->EroDep( 0, deposition_depth, time );  // add or subtract net flux/area    
	  cn->getDownstrmNbr()->addQsdin(-1 * cn->getQsin()/dtmax);  //what does this do? removed or added in, can't see diff
	  
	  
	  if( 0 ) //DEBUG
	   {
        std::cout << "Node " << cn->getID() << " Qsin: " << cn->getQsin()
        << " dz: " << cn->getQsin() / cn->getVArea() << std::endl;
		int gg;
		  for ( gg=0; gg<num_grain_sizes_; gg++ )
		 {
          std::cout << "  size " << gg << " depo depth is " << deposition_depth[gg] << std::endl;
		  std::cout << "  size " << gg << " Qsin " << cn->getQsin(gg) << std::endl;
		  std::cout << "  size " << gg << " percent depo depth " << 
		  deposition_depth[gg] / (cn->getQsin() / cn->getVArea()) << std::endl;
		  }
		}
		
         
    }
    
    rt -= dtmax;
    if( dtmax>rt ) dtmax=rt;
    
    if(0) std::cout << "bottom of do loop in Diffuse()\n" << std::flush;
    
  } while( rt>0.0 );
  
  
}
#undef kEpsOver2



/*****************************************************************************\
 **
 **  tErosion::DiffuseNonlinear
 **
 **  This function implements the nonlinear creep transport function of
 **  Howard (1994), which was also studied by Roering et al. (1999, 2002).
 **  The sediment flux per unit slope width is:
 **
 **                 Dz
 **    qs = Kd -------------
 **            1 - (Dz/Sc)^2
 **
 **  where Dz is the land surface gradient in 2D, Kd is a creep coefficient
 **  (L^2/T), and Sc is a threshold slope.
 **
 **  The numerical implementation on a Voronoi grid uses the same approach
 **  as used by Diffuse(), with the following differences:
 **
 **  1. The estimated Courant condition includes the denominator in the 
 **     above equation, and because this contains local slope, it must be
 **     re-evaluated for every time step, so it is now inside the time loop.
 **     (The additional square root of [1-(Dz/Sc)^2] is based on experiments with
 **     a 1D version of this solver, which shows improved stability when you
 **     factor this in; otherwise, instabilities appear as f grows small).
 **  2. The equation becomes invalid when Dz/Sc>=1.0. Therefore a safety
 **     feature is used: Dz/Sc is not allowed to go above an arbitrary value
 **     that is very close to unity. Because of this, it's possible that the
 **     slope angle at the base of a slope with a very high flux may exceed
 **     Sc. In a 1D steady-state profile, this shows up as a linear increase in
 **     gradient with distance, above Sc, near the base of the slope.
 **  3. To improve performance (or so I hope), both slope and f (the 
 **     denominator in the above equation) are stored in temporary arrays
 **     during the time step calculation loop.
 **  4. For the first time, I depart from previous practice by using
 **     the vector class from the standard template library, instead of our
 **     custom-built tArray class.
 **  5. The time step calculation loop skips complementary edges (that's 
 **     why slope and f only need N/2 elements).
 **
 **  Inputs:  rt -- time duration over which to compute diffusion
 **           noDepoFlag -- if true, material is only eroded, never
 **                             deposited
 **  Created: August 2007, GT
 **  Modifications:
 ** 
 \*****************************************************************************/
//#define kVerySmall 1e-6
#define kEpsOver2 0.1
#define kBeta 0.999    // Dz/Sc isn't allowed to go higher than this
void tErosion::DiffuseNonlinear( double rt, bool noDepoFlag, double time )
{
  tLNode * cn;
  tEdge * ce;
  double volout,  // Sediment volume output from a node (neg=input)
  delt,       // Max local step size
  dtmax;      // Max global step size (initially equal to total time rt)
  tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );
  tMesh< tLNode >::edgeListIter_t edgIter( meshPtr->getEdgeList() );
  int numActiveEdges = meshPtr->getEdgeList()->getActiveSize();
  vector<double> slope;   // Slope of each edge
  vector<double> f;       // = 1 - (slope/Sc)^2
  int k;      // Counter for edges
  double slopeRatio;   // Ratio of slope to critical slope
  
#ifdef TRACKFNS
  std::cout << "tErosion::DiffuseNonlinear()" << std::endl;
#endif
	
  kd = kd_ts.calc( time );
  
  if( kd==0 ) return;
  
  // Set the size of the vectors
  slope.resize( numActiveEdges/2 );
  f.resize( numActiveEdges/2 );
  
  //initialize Qsd, which will record the total amount of diffused material
  //fluxing into a node for the entire time-step.
  //if Qsd is negative, then material was deposited in that node.
  for( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
    cn->setQsdin( 0. );
  
  // Loop until we've used up the entire time interval rt
  do
  {
    // Compute maximum stable time-step size based on modified Courant condition
    // for FTCS (here used as an approximation).
    dtmax = rt;  // Initialize dtmax to total time rt
    k=0;
    for( ce=edgIter.FirstP(); edgIter.IsActive(); ce=edgIter.NextP() )
    {
      if( 0 ) //DEBUG
      {
        std::cout << "In Diffuse(), large vedglen detected: " << ce->getVEdgLen() << std::endl;
        ce->TellCoords();
      }
      
      // Evaluate DT <= DX^2 f^2 / Kd
      slope[k] = ce->CalcSlope();            // compute and store slope for this edge
      slopeRatio = fabs( slope[k] / mdSc );  // calculate slope ratio
      if( slopeRatio > kBeta ) slopeRatio = kBeta;  // don't let it reach 1 or higher
      f[k] = 1.0 - slopeRatio*slopeRatio;           // compute and store the nonlinear factor
      delt = kEpsOver2 * ce->getLength()*ce->getLength()*f[k]*sqrt(f[k]) / kd;  // max. time step this edge
      if( delt < dtmax )
      {
        dtmax = delt;  // remember the smallest delt
        if(0) { //DEBUG
          std::cout << "TIME STEP CONSTRAINED TO " << dtmax << " AT EDGE:\n";
          ce->TellCoords(); }
      }
      edgIter.NextP();  // Skip complementary edge
      k++;   // increment edge counter
    }
    assert( k==numActiveEdges/2 );
    
    // Reset sed input for each node for the new iteration
    for( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
      cn->setQsin( 0. );
    
    // Compute sediment volume transfer along each edge
    k=0;  // reset edge counter
    for( ce=edgIter.FirstP(); edgIter.IsActive(); ce=edgIter.NextP() )
    {
      volout = kd*(slope[k]/f[k])*ce->getVEdgLen()*dtmax;  // specific flux times width times time step
      // Record outgoing flux from origin
      cn = static_cast<tLNode *>(ce->getOriginPtrNC());
      if( difThresh>0. && cn->getDrArea()>difThresh )
        volout=0;
      cn->addQsin( -volout );
      // Record incoming flux to dest'n
      cn = static_cast<tLNode *>(ce->getDestinationPtrNC());
      cn->addQsin( volout );
      
      if( 0 ) { //DEBUG
        std::cout << volout << " mass exch. from " << ce->getOriginPtr()->getID()
        << " to "
        << ce->getDestinationPtr()->getID()
        << " on slp " << ce->getSlope() << " ve " << ce->getVEdgLen()
        << " kd=" << kd << " dtmax=" << dtmax
        << "\nvp " << ce->getRVtx().at(0)
        << " " << ce->getRVtx().at(1) << std::endl;
        static_cast<tLNode *>(ce->getOriginPtrNC())->TellAll();
        static_cast<tLNode *>(ce->getDestinationPtrNC())->TellAll();
        std::cout << std::endl;
      }
      
      edgIter.NextP();  // Skip complementary edge
      k++;  // reset edge counter
    }
    assert( k==numActiveEdges/2 );
    
    // Compute erosion/deposition for each node
    for( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
    {
      if( 0 ) //DEBUG
        std::cout << "Node " << cn->getID() << " Qsin: " << cn->getQsin()
        << " dz: " << cn->getQsin() / cn->getVArea() << std::endl;
      if( noDepoFlag && cn->getQsin() > 0.0 )
        cn->setQsin( 0.0 );
      cn->EroDep( cn->getQsin() / cn->getVArea() );  // add or subtract net flux/area    
      cn->getDownstrmNbr()->addQsdin(-1 * cn->getQsin()/dtmax);
      //this won't work if time steps are varying, because you are adding fluxes
      
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
#undef kEpsOver2
#undef kBeta

/*****************************************************************************\
 **
 **  tErosion::DiffuseNonlinearDepthDep
 **
 **  This function implements the nonlinear creep transport function of
 **  Howard (1994), which was also studied by Roering et al. (1999, 2002)
 **  and modified to be depth-dependent by Roering (2008).
 **  The sediment flux per unit slope width is:
 **
 **                 Dz
 **    qs = Kd -------------
 **            1 - (Dz/Sc)^2
 **
 **  where Dz is the land surface gradient in 2D, Kd is a creep coefficient
 **  (L^2/T), and Sc is a threshold slope. In this version, Kd is dependent
 **  on soil depth: Kd = eta * (1 - exp( h*cos(theta)/H_0 )). And outflux
 **  at each node is limited to supply of soil.
 **
 **  The numerical implementation on a Voronoi grid uses the same approach
 **  as used by Diffuse(), with the following differences:
 **
 **  1. The estimated Courant condition includes the denominator in the 
 **     above equation, and because this contains local slope, it must be
 **     re-evaluated for every time step, so it is now inside the time loop.
 **     (The additional square root of [1-(Dz/Sc)^2] is based on experiments with
 **     a 1D version of this solver, which shows improved stability when you
 **     factor this in; otherwise, instabilities appear as f grows small).
 **  2. The equation becomes invalid when Dz/Sc>=1.0. Therefore a safety
 **     feature is used: Dz/Sc is not allowed to go above an arbitrary value
 **     that is very close to unity. Because of this, it's possible that the
 **     slope angle at the base of a slope with a very high flux may exceed
 **     Sc. In a 1D steady-state profile, this shows up as a linear increase in
 **     gradient with distance, above Sc, near the base of the slope.
 **  3. To improve performance (or so I hope), both slope and f (the 
 **     denominator in the above equation) are stored in temporary arrays
 **     during the time step calculation loop.
 **  4. For the first time, I depart from previous practice by using
 **     the vector class from the standard template library, instead of our
 **     custom-built tArray class.
 **  5. The time step calculation loop skips complementary edges (that's 
 **     why slope and f only need N/2 elements).
 **
 **  Inputs:  rt -- time duration over which to compute diffusion
 **           time -- runtime, for updating layers (needed by EroDep)
 **
 **  Created: July, 2010, SL
 **  Modifications:
 ** 
 \*****************************************************************************/
#define kEpsOver2 0.1
#define kBeta 0.999    // Dz/Sc isn't allowed to go higher than this
void tErosion::DiffuseNonlinearDepthDep( double rt, double time )
{
#ifdef TRACKFNS
  std::cout << "tErosion::DiffuseNonlinear()" << std::endl;
#endif
	
  kd = kd_ts.calc( time );
  
  if( kd==0 ) return;
  
  tLNode * cn;
  tEdge * ce;
  double delt;       // Max local step size
  double dtmax;      // Max global step size (initially equal to total time rt)
  tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );
  tMesh< tLNode >::edgeListIter_t edgIter( meshPtr->getEdgeList() );
  int numActiveEdges = meshPtr->getEdgeList()->getActiveSize();
  int numEdges = meshPtr->getEdgeList()->getSize();
  
  int k;      // Counter for edges
  double slopeRatio;   // Ratio of slope to critical slope
  
  // Set the size of the vectors
  vector<double> slope( numActiveEdges/2 );   // Slope of each edge
  vector<double> f( numActiveEdges/2 );       // = 1 - (slope/Sc)^2
  vector<double> edgeH( numActiveEdges/2 ); // average soil depth for edge
  vector<double> edgeKd( numActiveEdges/2 ); // depth-dependent param.
  vector<double> edgeFlux( numActiveEdges/2 ); // store fluxes along edges
  vector<int> tempArrayIndex( numEdges ); // indexes to above arrays
  
  //initialize Qsd, which will record the total amount of diffused material
  //fluxing into a node for the entire time-step.
  //if Qsd is negative, then material was deposited in that node.
  for( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
    cn->setQsdin( 0. );
  
  // Loop until we've used up the entire time interval rt
  do
  {
    // Compute maximum stable time-step size based on modified Courant condition
    // for FTCS (here used as an approximation).
    dtmax = rt;  // Initialize dtmax to total time rt
    k=0;
    for( ce=edgIter.FirstP(); edgIter.IsActive(); ce=edgIter.NextP() )
    {
      // Evaluate DT <= DX^2 f^2 / Kd
      slope[k] = ce->CalcSlope(); // compute and store slope for this edge
      slopeRatio = fabs( slope[k] / mdSc );  // calculate slope ratio
      if( slopeRatio > kBeta ) slopeRatio = kBeta; // don't let it reach 1 or higher
      f[k] = 1.0 - slopeRatio*slopeRatio; // compute and store the nonlinear factor
      // determine depth-dependent transport coefficient:
      // use regolith depth of upslope endpoint:
      if( ce->getOriginPtr()->getZ() > ce->getDestinationPtr()->getZ() )
        cn = static_cast<tLNode *>(ce->getOriginPtrNC());
      else
        cn = static_cast<tLNode *>(ce->getDestinationPtrNC());
      // soil thickness:
      double nodeSoilThickness(0.0);
      tListIter< tLayer > lI( cn->getLayersRefNC() );
      for( tLayer *lP=lI.FirstP(); lP->getSed() == tLayer::kSed; lP=lI.NextP() )
        nodeSoilThickness += lP->getDepth();
      edgeH[k] = nodeSoilThickness;    
      edgeKd[k] = 
	    kd * ( 1 - exp( -edgeH[k] * cos( atan( slope[k] ) ) / diffusionH ) );
      // max. time step this edge:
      delt = 
	    kEpsOver2 * ce->getLength()*ce->getLength()*f[k]*sqrt(f[k]) / edgeKd[k];  
      if( delt < dtmax ) dtmax = delt;  // remember the smallest delt
      tempArrayIndex[ce->getID()] = k; // store index for this edge
      ce = edgIter.NextP();  // Skip complementary edge
      tempArrayIndex[ce->getID()] = k; // store index for complementary edge
      k++;   // increment edge counter
    }
    assert( k==numActiveEdges/2 );
    
    // Reset sed input and output for each node for the new iteration
    for( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
    {
      cn->setQsin( 0. );
      cn->setQs( 0. );
    }
    
    // Compute and store sediment volume transfer along each edge
    k=0;  // reset edge counter
    for( ce=edgIter.FirstP(); edgIter.IsActive(); ce=edgIter.NextP() )
    {
      // specific flux times width times time step:
      edgeFlux[k] = edgeKd[k] * ( slope[k] / f[k]) * ce->getVEdgLen() * dtmax;
      // Record fluxes at origin and destination:
      tLNode *on = static_cast<tLNode *>(ce->getOriginPtrNC());
      tLNode *dn = static_cast<tLNode *>(ce->getDestinationPtrNC());
      if( difThresh > 0. && on->getDrArea() > difThresh ) edgeFlux[k]=0;
      // account fluxes in and out separately in order to enforce
      // supply limitation:
      if( edgeFlux[k] > 0.0 )
	    {
	      on->addQs( -edgeFlux[k] );
	      dn->addQsin( edgeFlux[k] );
	    }
      else
	    {
	      on->addQsin( -edgeFlux[k] );
	      dn->addQs( edgeFlux[k] );
	    }      
      edgIter.NextP();  // Skip complementary edge
      k++;  // reset edge counter
    }
    assert( k==numActiveEdges/2 );
    
    // enforce supply limitation: outflux no greater than soil depth:
    for( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
    {
      // if any transport out, compare to soil depth:
      if( cn->getQs() < 0.0 )
	    {
	      // for soil thickness at node, find thickness associated with
	      // flowedge; unless node is flooded, then find first edge
	      // pointing downhill (and all nodes within these brackets will
	      // have a downhill neighbor because they have qs<0.0):
	      if( cn->getFloodStatus() == tLNode::kNotFlooded )
          k = tempArrayIndex[ cn->getFlowEdg()->getID() ];
	      else
        {
          tSpkIter sI( cn );
          for( ce = sI.FirstP(); !sI.AtEnd(); ce = sI.NextP() )
            if( ce->getOriginPtr()->getZ() > ce->getDestinationPtr()->getZ() 
               && ce->FlowAllowed() )
              break;
          k = tempArrayIndex[ ce->getID() ];
        }
	      double nodeSoilThickness = edgeH[k];
	      if( -cn->getQs() / cn->getVArea() > nodeSoilThickness )
        {
          // if flux out more than soil depth, 
          // multiply fluxes out by the factor:
          const double reducFactor = 
          -nodeSoilThickness * cn->getVArea() / cn->getQs();
          // go through node's edges:
          tSpkIter sI( cn );
          for( ce = sI.FirstP(); !sI.AtEnd(); ce = sI.NextP() )
          {
            // check that edge is not connected to a closed boundary: 
            if( likely( ce->FlowAllowed() ) )
            {
              double thisEdgeFlux = 
              edgeFlux[ tempArrayIndex[ce->getID()] ];
              if( ce->getID()%2 == 0 ) thisEdgeFlux *= -1.0;
              if( thisEdgeFlux < 0.0 )
              {
                // if flux is out along this edge, reduce influx
                // at destination:
                tLNode *dn = 
                static_cast<tLNode *>(ce->getDestinationPtrNC());
                dn->addQsin( -thisEdgeFlux * ( reducFactor - 1.0 ) );
              }
            }
          }
          // change sediment outflux for node:
          cn->setQs( -nodeSoilThickness * cn->getVArea() );
        }
	    }
    }
    
    // change elevations, etc., in a separate loop, after done adjusting fluxes:
    for( cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
    {
      tArray<double> erolist( cn->getNumg() );
      // elevation change is net flux per area:
      double deltaZ = ( cn->getQs() + cn->getQsin() ) / cn->getVArea();
      // add elevation change, and mind the layers:
      if( deltaZ > 0.0 )
	    {
	      for( size_t j=0; j<cn->getNumg(); ++j )
          erolist[j] = deltaZ * cn->getLayerDgrade(0,j)/cn->getLayerDepth(0);
	      cn->EroDep( 0, erolist, time );  // add or subtract net flux/area    
	    }
      else if( deltaZ < 0.0 )
        while( deltaZ < 0.0 )
	      {
          if( -deltaZ <= cn->getLayerDepth(0) )
          {
            for( size_t j=0; j<cn->getNumg(); ++j )
              erolist[j] = 
              deltaZ * cn->getLayerDgrade(0,j) / cn->getLayerDepth(0);
            deltaZ = 0.0;
          }
          else
          {
            for( size_t j=0; j<cn->getNumg(); ++j )
              erolist[j] = cn->getLayerDgrade(0,j);
            deltaZ += cn->getLayerDepth(0);
          }
          cn->EroDep( 0, erolist, time );
	      }
      cn->getDownstrmNbr()->addQsdin(-1 * cn->getQs()/dtmax);
      //this won't work if time steps are varying, because you are adding fluxes     
    }
    rt -= dtmax;
    if( dtmax>rt ) dtmax=rt;
  } while( rt>0.0 );
}
#undef kEpsOver2
#undef kBeta

/***************************************************************************\
 **  tErosion::ProduceRegolith( double dtg, double time )
 **
 **  Converts rock into soil.
 **  Called by: main routine.
 **  Takes: time step, runtime.
 **  Calls: 
 **    - tPhysicalWeathering::SoilProduction( tLNode* );
 **    - tLNode::EroDep for adding soil to regolith layers;
 **    - lower-level tLayer and tLNode functions for making and changing
 **      layers and changing elevation since many of the contingencies and 
 **      checks in EroDep are unnecessary in the special case of eroding 
 **      bedrock.
 **  Changes: Elevations and layers for active tLNodes.
 **
 **  - STL, 6/2010
 \***************************************************************************/
void tErosion::ProduceRegolith( double dtg, double time )
{
  tMesh< tLNode >::nodeListIter_t ni( meshPtr->getNodeList() ); // node iter.
  // do physical weathering for each active node:
  for( tLNode* n = ni.FirstP(); ni.IsActive(); n = ni.NextP() )
  {
    // find rate of bedrock lowering:
    double rate = physWeath->SoilProduction( n );
    if( rate < 0.0 ) // skip it all if no soil production
    {
      double rockDeltaZ = rate * dtg; // bedrock lowering
      tListIter< tLayer > lI( n->getLayersRefNC() ); // layer iterator
      tLayer *lP=0; // current layer pointer
      tLayer *soilP=0; // pointer to bottom soil layer
      int i=0; // index to layer, ends up set to top bedrock layer
      // find top bedrock layer and bottom soil layer:
      for( lP=lI.FirstP(), i=0; 
          lP->getSed() == tLayer::kSed; 
          lP=lI.NextP(), ++i )
        soilP = lP;
      tLayer *rockP=lP; // pointer to top bedrock layer
      // soil thickening corresponding to rock lowering; note that rock  
      // density is variable, dependent on chemical weathering, but we're
      // using a constant soil bulk density as a parameter in 
      // tPhysicalWeathering,and there may not be a soil layer from which 
      // to get a value:
      double soilDeltaH = -rockDeltaZ * rockP->getBulkDensity() 
	    / soilBulkDensity;
      // change elevation for bedrock lowering:
      n->ChangeZ( rockDeltaZ );
      // do bedrock lowering:
      tArray< double > erolist( n->getNumg() ); // erosion for each grainsize
      // remove rock from top rock layer or layers; decrement rockDeltaZ as
      // we go until none left; most of the time, this amount will be small,
      // and top bedrock layer will accommodate all rock lowering:
      while( rockDeltaZ < 0.0 )
	    {
	      // check to see whether change will deplete top bedrock layer:
	      if( -rockDeltaZ < rockP->getDepth() )
        {
          // simply remove material from top bedrock layer 
          // (use low-level function):
          for( size_t j=0; j<n->getNumg(); ++j ) 
            rockP->addDgrade( j, 
                             rockDeltaZ * rockP->getDgrade(j) 
                             / rockP->getDepth() );
          // all erosion in this layer, so rockDeltaZ goes to zero:
          rockDeltaZ=0.0;
          // should I worry about removing very thin layers here?
        }
	      else
        {
          // erode entire layer; decrement rockDeltaZ by layer depth:
          rockDeltaZ += rockP->getDepth();
          // if it's the last layer, then make a new one here:
          if( unlikely( i == n->getNumLayer()-1 ) )
          {
            // insert copy of rock layer:
            n->getLayersRefNC().insertAtBack( *rockP );
            // reset its depth:
            n->setLayerDepth( n->getNumLayer()-1, n->getMaxregdep() );
          }	  
          n->removeLayer(i);
          // just removed layer pointed to by rockP; 
          // re-find soil and rock layers:
          for( lP=lI.FirstP(), i=0; 
              lP->getSed() == tLayer::kSed; 
              lP=lI.NextP(), ++i )
            soilP = lP;
          rockP = lP;
        }
	    } // rock lowering loop
      // add material to bottom soil layer (elevation changed in EroDep):
      for( size_t j=0; j<n->getNumg(); ++j ) 
        erolist[j] = soilDeltaH * rockP->getDgrade(j)/rockP->getDepth();
      // is there a soil layer?
      if( soilP > 0 )
        // yes; add material to layer above top rock layer:
        n->EroDep( i-1, erolist, time );
      else
	    {
	      // no; call EroDep such that it will make a new layer:
	      n->EroDep( 0, erolist, time );
	      soilP = lI.FirstP(); // set soil layer pointer
	      // need to set soil bulk density (necessary?):
	      soilP->setBulkDensity( soilBulkDensity );
	    }
    } // end if( rate < 0.0 ) 
  } // nodeList loop
} // void tErosion::ProduceRegolith( double dtg, double time )

/***************************************************************************\
 **  tErosion::WeatherBedrock( double dtg )
 **
 **  Changes bedrock (e.g., reduces its density and/or changes its volume.
 **  Called by: main routine.
 **  Takes: time step.
 **  Calls: 
 **    - tChemicalWeathering::SoluteFlux( tLNode*, timestep );
 **    - tChemicalWeathering::StrainRate( tLNode*, timestep );
 **  Changes: Bulk densities and thicknesses of layers.
 **
 **  - STL, 6/2010
 \***************************************************************************/
void tErosion::WeatherBedrock( double dtg )
{
  tMesh< tLNode >::nodeListIter_t ni( meshPtr->getNodeList() ); // node iter.
  // do chemical weathering for each active node:
  double totalFlux=0.0;
  double totalStrain=0.0;
  for( tLNode* n = ni.FirstP(); ni.IsActive(); n = ni.NextP() )
  {
    // find flux; this version updates bulk density of each bedrock layer:
    totalFlux += chemWeath->SoluteFlux( n, dtg );
    // find strain; as of 6/2010, does nothing:
    totalStrain += chemWeath->StrainRate( n, dtg );
  }
}

/***************************************************************************\
 **  tErosion::LandslideClusters(): landslide initiation
 **
 **  IMPORTANT: This function requires the tStreamNet::kFlowSaturated1
 **  option to work correctly, although calling of this function is 
 **  governed by a separate input option.
 **
 **  Adapted for incorporation into CU CHILD from version in OSU CHILD. 
 **  Old version (OSU) failed only individual nodes. This version finds 
 **  clusters of landslide nodes. This is a BIG DEAL!
 **
 **  Uses a priority_queue of nodes sorted according to net downhill force.
 **  Pops the top node of the queue (greatest net downhill force) and builds
 **  a cluster around that node; nodes are added if (a) original node has
 **  factor > 1 and the new node would decrease the factor or (b) original
 **  node has factor < 1 and the new node would not make the factor > 1.
 **  Continues to make another cluster if the current one fails (i.e., 
 **  stops looking for subsequent failure clusters if previously found
 **  clusters don't fail). Each landslide cluster results in the 
 **  construction of a tDebrisFlow object, which is added to a list. 
 **  Runout, etc., including removing soil and changing elevations in
 **  the failure zone is done separately by the tDebrisFlow object itself.
 **
 **  Incorporates lateral cohesive strength due to roots, as in Schmidt 
 **  et al., 2001. Neglects body (acive and passive) constraining stresses.
 **
 **  This is a big (i.e., many lines) function that does a lot of 
 **  calculations, but actual changes to "permanent" data are made 
 **  elsewhere. It would have been nice to use some smaller functions,
 **  but these don't seem feasible. In particular, methods used here
 **  may not be generalizable in the sense that multiple calculations
 **  are required for each candidate node.
 **
 **  Takes: 
 **    - Rain rate, run time, pointer to tStreamNet (but mainly relies on 
 **      hydrologic calc's in tStreamNet::FlowSaturated1
 **  Calls: 
 **    - tDebrisFlow constructor for each debris flow
 **    - tDebrisFlow::RunScourDeposit: main debris flow function for
 **      each debris flow
 **    - Overloaded output operator for each debris flow
 **  Changes: 
 **    - finds failure clusters, but those members of tDebrisFlow are
 **      actually set within that class 
 **    - tNode::public1
 **
 **  - SL, 9/2010
 **
 \***************************************************************************/
void tErosion::LandslideClusters( double rainrate, 
                                 double time )
{
  tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );
  tMesh< tLNode >::edgeListIter_t edgIter( meshPtr->getEdgeList() );
  int numActiveEdges = meshPtr->getEdgeList()->getActiveSize();
  int numActiveNodes = meshPtr->getNodeList()->getActiveSize();
  int numEdges = meshPtr->getEdgeList()->getSize();
  int numNodes = meshPtr->getNodeList()->getSize();
  vector<double> edgeLatCohesion( numActiveEdges/2 ); // lateral cohesion
  vector<double> edgeSlope( numActiveEdges/2 );
  vector<double> nodeSoilThickness( numActiveNodes );
  vector<double> nodeWoodDepth( numActiveNodes );
  vector<double> nodeWaterDepth( numActiveNodes );
  vector<double> nodeRootCohesionLat( numActiveNodes );
  vector<double> nodeBasalStrength( numActiveNodes );
  vector<double> nodeDrivingForce( numActiveNodes );
  vector<double> nodeLatCohesion( numActiveNodes );
  vector<double> nodeNetForce( numActiveNodes );
  vector<int> tempEdgeIndex( numEdges ); // indexes to above arrays
  vector<int> tempNodeIndex( numNodes );  
  
  tPtrList<tDebrisFlow> dfPList; // list of failures
  
  // empty event tallies:
  debris_flow_sed_bucket = 0.0;
  debris_flow_wood_bucket = 0.0;
  
  { // find soil depth and root cohesion and gradient vector at each node:
    tLNode* cn=0;
    int i;
    for( cn=nodIter.FirstP(), i=0; 
	 nodIter.IsActive(); 
	 cn=nodIter.NextP(), ++i )
    {
      nodeSoilThickness[i] = cn->getRegolithDepth();
      if( cn->getVegCover().getTrees() > 0 )
      {
        nodeWoodDepth[i] = 
	  cn->getVegCover().getTrees()->getBioMassStand() 
	  + cn->getVegCover().getTrees()->getBioMassDown();
        nodeRootCohesionLat[i] = 
	  cn->getVegCover().getTrees()->getRootStrengthLat();
      }
      tempNodeIndex[cn->getID()] = i;
      cn->public1 = 0; // initialize flags for cluster membership
    }
  }
  { // calculate lateral cohesion at each edge:
    tEdge* ce=0;
    int k;
    for( ce=edgIter.FirstP(), k=0; 
	 edgIter.IsActive(); 
	 ce=edgIter.NextP(), ++k )
    {
      edgeSlope[k] = ce->CalcSlope();
      const double costheta = cos( atan( edgeSlope[k] ) );
      const int iOrg = tempNodeIndex[ ce->getOriginPtr()->getID() ];
      if( ce->getDestinationPtr()->isNonBoundary() )
      {
	const int iDest = 
	  tempNodeIndex[ ce->getDestinationPtr()->getID() ];
	edgeSlope[k] -= 
	  ( nodeSoilThickness[iOrg] - nodeSoilThickness[iDest] ) 
	  / ce->getLength();
        const double edgeDepth = 
	  ( nodeSoilThickness[iOrg] + nodeSoilThickness[iDest] ) / 2.0;
	// use minimum of root strengths at edge endpoints:
	const double edgeRoots =
	  min( nodeRootCohesionLat[iOrg], nodeRootCohesionLat[iDest] );
// 	    const double edgeRoots =
// 	      ( nodeRootCohesionLat[iOrg] + nodeRootCohesionLat[iDest] ) 
// 	      / 2.0;
        edgeLatCohesion[k] = 
	  edgeRoots * diffusionH * ce->getVEdgLen() / costheta
	  * ( 1.0 - exp( -edgeDepth / diffusionH * costheta ) );
      }
      else
      {
        edgeSlope[k] -= nodeSoilThickness[iOrg] / ce->getLength();
        // at boundaries use only edge origin:
        // (NOTE: above scheme is the same as assuming there is an identical node 
        // on the other side of the boundary; alternative to above is to have
        // zero cohesion at boundaries, but don't want to "favor" having 
        // failure sides at boundaries)
        edgeLatCohesion[k] = 
	  nodeRootCohesionLat[iOrg] * diffusionH * ce->getVEdgLen() 
	  * ( 1.0 - exp( -nodeSoilThickness[iOrg] 
			 / diffusionH * costheta ) );
      }
      // give same index to edge and its complement:
      tempEdgeIndex[ce->getID()] = k;
      ce=edgIter.NextP();
      tempEdgeIndex[ce->getID()] = k;  
    }
  }
  { // calculate basal strength, driving force, and net downhill force for 
    // each node:
    const double porosity = ( wetBulkDensity - soilBulkDensity ) / RHO;
    //     tSpkIter sI;
    tLNode* cn=0;
    int i;
    tSpkIter sI;
    for( cn=nodIter.FirstP(), i=0; 
	 nodIter.IsActive(); 
	 cn=nodIter.NextP(), ++i )
    {
      // use flowedge bedrock slope:
      const int iEdge = cn->getFlowEdg()->getID();
      double slope = edgeSlope[tempEdgeIndex[iEdge]];
      if( iEdge % 2 == 1 ) slope *= -1.0;
      const double slopeangle = atan( slope );
      const double costheta = cos( slopeangle );
      const double sintheta = sin( slopeangle );
      double basalCohesion = 0.0;
      if( cn->getVegCover().getTrees() > 0 )
	basalCohesion = 
	  cn->getVArea() / costheta
	  * cn->getVegCover().getTrees()->getRootStrengthVert()
	  * exp( -nodeSoilThickness[i] * cos( atan( cn->calcSlope() ) )
		 / diffusionH );
      //HYDROLOGY: use subsurface kinematic wave routing in tStreamNet
      // in case sintheta is zero or negative, max depth is soil depth:
      double saturatedDepth = cn->getHydrDepth();
      if( saturatedDepth == 0.0 ) // should only happen if in a pit:
	saturatedDepth = nodeSoilThickness[i];
// 	// make sure not greater than soil depth:
// 	if( saturatedDepth > nodeSoilThickness[i] )
// 	  saturatedDepth = nodeSoilThickness[i];
      nodeWaterDepth[i] = saturatedDepth * porosity; // water fraction
      const double weight = 
	wetBulkDensity * nodeSoilThickness[i] // soil
	+ woodDensity * nodeWoodDepth[i];     // wood
      const double gravTimesArea = GRAV * cn->getVArea();
      const double basalFriction = 
	( weight - RHO * saturatedDepth )              // pore pressure
	* gravTimesArea * costheta * fricSlope;
      nodeBasalStrength[i] = basalCohesion + basalFriction;
      // if slope and, therefore, sintheta are negative, will just 
      // make negative driving force and net force--no problem,
      // since we're using net force rather than factor of safety:
      nodeDrivingForce[i] = weight * gravTimesArea * sintheta;
      sI.Reset(cn);
      // add up lateral cohesion terms from edges:
      for( tEdge* ce=sI.FirstP(); !sI.AtEnd(); ce=sI.NextP() )
	nodeLatCohesion[i] += edgeLatCohesion[ tempEdgeIndex[ce->getID()] ];
      // find net downhill force for each node.
      nodeNetForce[i] = 
	nodeDrivingForce[i] - nodeLatCohesion[i] - nodeBasalStrength[i];
      if(0) //DEBUG
	cout << cn->getSubSurfaceDischarge() << " " << cn->getDrArea() 
	     << " " << slope << "\n ";
      // record forces for nodes:
      cn->setNetDownslopeForce( nodeNetForce[i] );
    }
    //     cout << "\n";
  }
  typedef vector<NodeNetForceIndex> NodeContainer;
  // construct a priority_queue type that uses the custom comparator:
  priority_queue<NodeNetForceIndex, NodeContainer, NodeNetForce_Lesser> 
    nodePQ;
  // put nodes in a priority_queue with greatest net downhill force at the 
  // top.
  for( tLNode* cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
  {
    NodeNetForceIndex curNSFI;
    curNSFI.node = cn;
    curNSFI.netForce = nodeNetForce[tempNodeIndex[cn->getID()]];
    curNSFI.index = tempNodeIndex[cn->getID()];
    nodePQ.push( curNSFI );
  }
  { // pop top node in queue and build cluster...
    bool anyFailuresThisPass;
    int numCluster = 0;
    tSpkIter sI;
    tSpkIter neI;
    do
    {
      anyFailuresThisPass = false;
      // increment number used to flag nodes in a cluster:
      ++numCluster;
      tPtrList<tLNode> slideCluster;
      tPtrList<tLNode> seedList; // list of nodes on edge of cluster
      NodeNetForceIndex curNSFI;
      do
      {
	// start potential cluster with node with greatest net downhill 
	// force:
	curNSFI = nodePQ.top();
	// and remove it from the queue:
	nodePQ.pop();
      } while( curNSFI.node->public1 > 0 ); // if already part of a 
      // cluster, try again
      seedList.insertAtBack( curNSFI.node );
      slideCluster.insertAtBack( curNSFI.node );
      double initNetForce = curNSFI.netForce;
      double initLatCohesion = nodeLatCohesion[curNSFI.index];
      double initBasalStrength = nodeBasalStrength[curNSFI.index];
      double initDrivingForce = nodeDrivingForce[curNSFI.index];
      double finalNetForce = initNetForce;
      double finalLatCohesion = initLatCohesion;
      double finalBasalStrength = initBasalStrength;
      double finalDrivingForce = initDrivingForce;
      // set generic flag signifying node in cluster:
      curNSFI.node->public1 = numCluster;
      while( !seedList.isEmpty() )
      {
	// get seed for cluster growth
	tLNode* cn = seedList.removeFromFront();
	sI.Reset( cn );
	// search seed's neighbors
	for( tEdge* ce = sI.FirstP(); !sI.AtEnd(); ce = sI.NextP() )
	  if( ce->FlowAllowed() && 
	      ce->getDestinationPtr()->isNonBoundary() )
	  { // if not on boundary:
	    tLNode* nn = 
	      static_cast<tLNode*>( ce->getDestinationPtrNC() );
	    // if node is not already in a cluster
	    // AND it's connected to the cluster by a flow edge:
	    if( nn->public1 == 0 
		&& ( cn->flowThrough( ce ) 
		     || nn->flowThrough( ce->getComplementEdge() ) ) )
	    {
	      // increment basal strength and driving force:
	      finalBasalStrength += 
		nodeBasalStrength[tempNodeIndex[nn->getID()]];
	      finalDrivingForce += 
		nodeDrivingForce[tempNodeIndex[nn->getID()]];
	      // go through candidate's neighbors for lateral 
	      // cohesion:
	      neI.Reset( nn );
	      for( tEdge* ne=neI.FirstP(); !neI.AtEnd(); 
		   ne=neI.NextP() )
		if( ne->FlowAllowed() 
		    && ne->getDestinationPtr()->isNonBoundary() )
		{
		  tLNode* nnn = 
		    static_cast<tLNode*>(ne->getDestinationPtrNC());
		  if( nnn->public1 > 0 )
		    // if neighbor of neighbor already in cluster,
		    // subtract edge's lateral cohesion:
		    finalLatCohesion -=
		      edgeLatCohesion[tempEdgeIndex[ne->getID()]];
		  else
		    // if neighbor of neighbor is not in cluster,
		    // add edge's lateral cohesion:
		    finalLatCohesion +=
		      edgeLatCohesion[tempEdgeIndex[ne->getID()]];
		}
	      // determine new net downhill force with the candidate:
	      finalNetForce =
		finalDrivingForce - finalLatCohesion 
		- finalBasalStrength;
	      // if addiiton of the new node doesn't decrease the 
	      // net downhill force below zero OR, if the net downhill
	      // force is already below zero, addition of the new
	      // node increases the net downhill force:
	      if( ( initNetForce > 0.0 && finalNetForce > 0.0 )
		  || ( initNetForce <= 0.0 && 
		       finalNetForce > initNetForce ) )
	      {
		// then add the new node to the cluster:
		nn->public1 = numCluster;
		seedList.insertAtBack( nn );
		slideCluster.insertAtBack( nn );
		// update "initial" terms:
		initNetForce = finalNetForce;
		initLatCohesion = finalLatCohesion;
		initBasalStrength = finalBasalStrength;
		initDrivingForce = finalDrivingForce;
	      }
	      else
	      {
		// if node not added, reset "final" terms:
		finalNetForce = initNetForce;
		finalLatCohesion = initLatCohesion;
		finalBasalStrength = initBasalStrength;
		finalDrivingForce = initDrivingForce;
	      }
	    }
	  }
      }
      // have cluster (which may be only the node popped off the queue);
      // check whether it fails:
      if( initNetForce > 0.0 )
      {
	// if failure, set logical to continue finding failure clusters:
	anyFailuresThisPass = true;
	// initialize new tDebrisFlow object with reference to cluster, 
	// copy of factor of safety, references to soil depth, wood 
	// depth, water depth, and node index vectors, and pointers to 
	// mesh and erosion objects:
	tDebrisFlow *dfPtr = 
	  new tDebrisFlow( slideCluster, initNetForce,
			   nodeSoilThickness, nodeWoodDepth, 
			   nodeWaterDepth, tempNodeIndex,
			   meshPtr, this );
	landslideAreas.insertAtBack( dfPtr->getAreaFailure() );
	// put copy of new debris flow pointer in list:
	dfPList.insertAtBack( dfPtr );
      }
    } while( anyFailuresThisPass ); // search for new cluster if last cluster failed
  }
  if( !dfPList.isEmpty() )
  {  // run out debris flows:
    const int numDFlows = dfPList.getSize();
    // remove debris flows from front of list (start with the
    // presumably least stable nodes):
    while( tDebrisFlow *dfPtr = dfPList.removeFromFront() )
    {
      dfPtr->RunScourDeposit();
      // output using overloaded "friend" operator:
      if( DF_fsPtr )
	*DF_fsPtr << time << " " << *dfPtr << endl;
      // constructed above with "new"; destroy now with "delete":
      delete dfPtr;
    }
    // write time, storm, # failures, and failed volumes to file:
    if( DF_Hyd_fsPtr )
    {
      DF_Hyd_fsPtr->setf( ios::scientific, ios::floatfield );
      DF_Hyd_fsPtr->precision(4);   
      *DF_Hyd_fsPtr 
	<< time << " " << rainrate << " " << numDFlows << " " 
	<< debris_flow_sed_bucket << " " 
	<< debris_flow_wood_bucket << endl;
    }
  }
}
/***************************************************************************\
 **  tErosion::LandslideClusters3D(): landslide initiation
 **
 **  Like LandslideClusters, but incorporates earth pressures in force
 **  balance, and actually does the vector sum of 2D force vectors for each
 **  node. These force vectors are each parallel to node gradient (positive
 **  uphill), and they are also added, weighted by node area, to find an
 **  average gradient for the cluster. At the end, if the dot-product of the
 **  force vector and the average gradient is negative (i.e., force vector
 **  points downhill) then the cluster fails. 
 **   Since buttressing forces may substantially counteract failures, this
 **  function repeats stability calculations as long as at least one
 **  failure is found in the current pass.
 **
 **  - SL, 11/2010
 **
\***************************************************************************/
void tErosion::LandslideClusters3D( double rainrate, 
				  double time )
{
  // initialize nodeList and edgeList iterators:
  tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );
  tMesh< tLNode >::edgeListIter_t edgIter( meshPtr->getEdgeList() );
  // initialize list sizes:
  const int numActiveEdges = meshPtr->getEdgeList()->getActiveSize();
  const int numActiveNodes = meshPtr->getNodeList()->getActiveSize();
  const int numEdges = meshPtr->getEdgeList()->getSize();
  const int numNodes = meshPtr->getNodeList()->getSize();
  // initialize temporary arrays:
  vector<double> edgeLatCohesion( numActiveEdges/2 ); // lateral cohesion
  vector<double> edgeSlope( numActiveEdges/2 );
  // earth pressures are going to be different for edge and complement:
  vector<double> edgeFrictionMag( numActiveEdges ); 
  vector<double> edgeFrictionX( numActiveEdges );
  vector<double> edgeFrictionY( numActiveEdges );
  vector<double> edgeBurden( numActiveEdges );
  vector<double> nodeGradientX( numActiveNodes );
  vector<double> nodeGradientY( numActiveNodes );
  vector<double> nodeSlopeAngle( numActiveNodes );
  vector<double> nodeGradMag( numActiveNodes );
  vector<double> nodeGradientUnitVectorX( numActiveNodes );
  vector<double> nodeGradientUnitVectorY( numActiveNodes );
  vector<double> nodeSoilThickness( numActiveNodes );
  vector<double> nodeSaturatedDepth( numActiveNodes );
  vector<double> nodeWoodDepth( numActiveNodes );
  vector<double> nodeWaterDepth( numActiveNodes );
  vector<double> nodeRootCohesionLat( numActiveNodes );
  vector<double> nodeLateralFriction( numActiveNodes );
  vector<double> nodeLatCohesion( numActiveNodes );
  vector<double> nodeNetForceX( numActiveNodes );
  vector<double> nodeNetForceY( numActiveNodes );
  vector<double> nodeNetForceMag( numActiveNodes );
  vector<int> tempEdgeIndex( numEdges ); // indexes to above arrays
  vector<int> tempFullEdgeIndex( numEdges );// indexes to long edge arrays
  vector<int> tempNodeIndex( numNodes );  

  tPtrList<tDebrisFlow> dfPList; // list of failures
  tPtrListIter<tDebrisFlow> dfI( dfPList ); 

  // find constants:
//   const double Ksat_soil = strmNet->getInfilt();
  const double porosity = ( wetBulkDensity - soilBulkDensity ) / RHO;
  const double frictionAngle = atan( fricSlope );
  const double coeffAtRestEarthPressure = 1.0 - sin( frictionAngle );
  const double coeffActiveEarthPressure = 
    pow( tan( PI/4.0 - frictionAngle / 2.0 ), 2.0 );
  const double coeffPassiveEarthPressure =
    pow( tan( PI/4.0 + frictionAngle / 2.0 ), 2.0 );
  const double cosPhi = cos( frictionAngle );

  // before first pass, initialize flags for cluster membership:
  for( tLNode* cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP() )
    cn->public1 = 0; 
  // BEGIN main loop: look for failures, and keep looking until we stop
  // finding them:
  int numPasses = 0;
  bool anyFailuresThisStorm = true;
  while( anyFailuresThisStorm )
    {
      ++numPasses;
      anyFailuresThisStorm = false;
      // empty event tallies; each pass through this loop will write
      // these tallies:
      debris_flow_sed_bucket = 0.0;
      debris_flow_wood_bucket = 0.0;

      // BEGIN finding force balance for each active node in mesh:
      { // find soil, wood, and water depth and lateral root 
	// cohesion at each node, and initialize public flags:
	tLNode* cn=0;
	int i;
	tSpkIter sI;
	for( cn=nodIter.FirstP(), i=0; 
	     nodIter.IsActive(); 
	     cn=nodIter.NextP(), ++i )
	  { // soil thickness:
	    nodeSoilThickness[i] = cn->getRegolithDepth();
	    //HYDROLOGY: use subsurface kinematic wave routing in tStreamNet
	    // in case magnitude of gradient is zero, max depth is soil depth:
	    nodeSaturatedDepth[i] = cn->getHydrDepth();
	    if( nodeSaturatedDepth[i] == 0.0 )// only happens if in a pit
	      nodeSaturatedDepth[i] = nodeSoilThickness[i];
	    nodeWaterDepth[i] = 
	      nodeSaturatedDepth[i] * porosity; // water fraction
	    if( cn->getVegCover().getTrees() > 0 )
	      { // wood depth and root cohesion
		nodeWoodDepth[i] = 
		  cn->getVegCover().getTrees()->getBioMassStand() 
		  + cn->getVegCover().getTrees()->getBioMassDown();
		nodeRootCohesionLat[i] = 
		  cn->getVegCover().getTrees()->getRootStrengthLat();
	      }
	    tempNodeIndex[cn->getID()] = i;
	  }
      }
      { // calculate bedrock slopes and lateral cohesion along each edge:
	tEdge* ce=0;
	int k;
	for( ce=edgIter.FirstP(), k=0; 
	     edgIter.IsActive(); 
	     ce=edgIter.NextP(), ++k )
	  {
	    edgeSlope[k] = ce->CalcSlope();
	    const double costheta = cos( atan( edgeSlope[k] ) );
	    const int iOrg = tempNodeIndex[ ce->getOriginPtr()->getID() ];
	    if( ce->getDestinationPtr()->isNonBoundary() )
	      {
		const int iDest = 
		  tempNodeIndex[ ce->getDestinationPtr()->getID() ];
		const double edgeDepth = 
		  ( nodeSoilThickness[iOrg] + nodeSoilThickness[iDest] ) 
		  / 2.0;
		// use minimum of root strengths at edge endpoints:
		const double edgeRoots =
		  min( nodeRootCohesionLat[iOrg], nodeRootCohesionLat[iDest] );
		// lateral cohesion uses surface slope of the edge:
		edgeLatCohesion[k] = 
		  edgeRoots * diffusionH * ce->getVEdgLen() / costheta
		  * ( 1.0 - exp( -edgeDepth / diffusionH * costheta ) );
		edgeSlope[k] -= 
		  ( nodeSoilThickness[iOrg] - nodeSoilThickness[iDest] ) 
		  / ce->getLength();
	      }
	    else
	      { // at boundaries use only edge origin:                          
		// (NOTE: this scheme is the same as assuming there is an 
		// identical node on the other side of the boundary;  
		// alternative to above is to have zero cohesion at boundaries, 
		// but don't want to "favor" having failure sides at 
		// boundaries) lateral cohesion uses zero slope at boundary:
		edgeLatCohesion[k] = 
		  nodeRootCohesionLat[iOrg] * diffusionH * ce->getVEdgLen() 
		  * ( 1.0 - exp( -nodeSoilThickness[iOrg] / diffusionH ) );
		edgeSlope[k] -= nodeSoilThickness[iOrg] / ce->getLength();
	      }
	    // give same index to edge and its complement:
	    tempEdgeIndex[ce->getID()] = k;
	    ce=edgIter.NextP();
	    tempEdgeIndex[ce->getID()] = k;  
	  }
      }
      { // find bedrock gradient vector, lateral and vertical root cohesion, 
	// basal friction, and driving force at each node, and add forces
	// to total:
	tLNode* cn=0;
	int i;
	tSpkIter sI;
	for( cn=nodIter.FirstP(), i=0; 
	     nodIter.IsActive(); 
	     cn=nodIter.NextP(), ++i )
	  { // node gradient (positive uphill); this is 
	    // sum( Vor. edge length * edge slope * edge unit vector )
	    // / sum( Vor. Edge length ):
	    sI.Reset( cn );
	    nodeGradientX[i] = 0.0;
	    nodeGradientY[i] = 0.0;
	    double sumWeights=0.0;
	    for( tEdge* ce=sI.FirstP(); !sI.AtEnd(); ce=sI.NextP() )
	      {
		const int iEdge = tempEdgeIndex[ce->getID()];
		double weightedSlopeByLength = 
		  -edgeSlope[iEdge] * ce->getVEdgLen() / ce->getLength();
		sumWeights += ce->getVEdgLen();
		// make slope positive uphill:
		if( ce->getID() % 2 == 0 ) weightedSlopeByLength *= -1.0;
		nodeGradientX[i] += 
		  weightedSlopeByLength * ce->getEVec().at(0);
		nodeGradientY[i] += 
		  weightedSlopeByLength * ce->getEVec().at(1);
		// for lateral cohesion, sum up cohesion along edges:
		if( cn->getVegCover().getTrees() > 0 )
		  nodeLatCohesion[i] += edgeLatCohesion[iEdge];
	      }
	    nodeGradientX[i] /= sumWeights;
	    nodeGradientY[i] /= sumWeights;
	    nodeGradMag[i] = 
	      sqrt( nodeGradientX[i] * nodeGradientX[i] 
		    + nodeGradientY[i] * nodeGradientY[i] );
	    nodeGradientUnitVectorX[i] = nodeGradientX[i] / nodeGradMag[i];
	    nodeGradientUnitVectorY[i] = nodeGradientY[i] / nodeGradMag[i];
	    nodeSlopeAngle[i] = atan( nodeGradMag[i] );
	    const double costheta = cos( nodeSlopeAngle[i] );
	    // add lateral cohesion (resistance) to force balance:
	    nodeNetForceMag[i] = -nodeLatCohesion[i];
	    // add basal root cohesion (resistance) to force balance
	    // (uses flowedge surface slope angle):
	    if( cn->getVegCover().getTrees() > 0 )
	      nodeNetForceMag[i] -= 
		cn->getVArea() / costheta
		* cn->getVegCover().getTrees()->getRootStrengthVert()
		* exp( -nodeSoilThickness[i] * cos( atan( cn->calcSlope() ) )
		       / diffusionH );
	    const double weight = 
	      wetBulkDensity * nodeSoilThickness[i] // soil
	      + woodDensity * nodeWoodDepth[i];     // wood
	    const double gravTimesArea = GRAV * cn->getVArea();
	    // add basal friction (resistance) to force balance:
	    nodeNetForceMag[i] -= 
	      ( weight - RHO * nodeSaturatedDepth[i] ) 
	      * gravTimesArea * costheta * fricSlope;
	    // add basal driving force to force balance:
	    nodeNetForceMag[i] += 
	      weight * gravTimesArea * sin( nodeSlopeAngle[i] );
	  }
      }
      { // calculate "active" and "passive" burden/buttress earth 
	// pressures and "at-rest" lateral friction for each edge; need 
	// to be able to remove friction for internal faces within 
	// clusters, so need to associate it with edge indexes, and need 
	// gradient calculated in above loop through nodes; could 
	// calculate burden/buttress pressures in above node loop, but
	// doing it here saves another spin through the spokes for each
	// node in the above loop and calculating "integWeight" twice for
	// each edge (cost is another temporary vector):         
	tEdge* ce=0;
	int k;
	for( ce=edgIter.FirstP(), k=0; 
	     edgIter.IsActive(); 
	     ce=edgIter.NextP(), ++k )
	  { 
	    const int iOrg = tempNodeIndex[ ce->getOriginPtr()->getID() ];	
	    const int iDest = 
	      tempNodeIndex[ ce->getDestinationPtr()->getID() ];
	    const double integWeight =
	      ( ( 0.5 * wetBulkDensity * nodeSoilThickness[iDest] 
		  + woodDensity * nodeWoodDepth[iDest] )
		* nodeSoilThickness[iDest]
		- 0.5 * RHO * nodeSaturatedDepth[iDest] 
		* nodeSaturatedDepth[iDest] );
	    edgeBurden[k] = 
	      integWeight * ce->getVEdgLen() / ce->getLength();
	    const double dotprod = 
	      ( ce->getEVec().at(0) * nodeGradientUnitVectorX[iOrg] 
		+ ce->getEVec().at(1) * nodeGradientUnitVectorY[iOrg] );
	    if( dotprod > 0 ) 
	      edgeBurden[k] *= dotprod * coeffActiveEarthPressure;
	    else 
	      edgeBurden[k] *= dotprod * coeffPassiveEarthPressure;
	    edgeFrictionMag[k] = 
	      GRAV * coeffAtRestEarthPressure * integWeight * cosPhi 
	      * fricSlope * 
	      fabs( ce->getVVec().at(0) * nodeGradientUnitVectorX[iOrg] 
		    + ce->getVVec().at(1) * nodeGradientUnitVectorY[iOrg] );
	    // lateral friction points uphill; need vectors when making 
	    // clusters:
	    edgeFrictionX[k] = 
	      edgeFrictionMag[k] * nodeGradientUnitVectorX[iOrg];
	    edgeFrictionY[k] = 
	      edgeFrictionMag[k] * nodeGradientUnitVectorY[iOrg];
	    tempFullEdgeIndex[ce->getID()] = k;
	  }
      }
      { // calculate net burden/buttress pressure and lateral friction and 
	// add them to net downhill force for each node:
	tSpkIter sI;
	tLNode* cn=0;
	int i;
	for( cn=nodIter.FirstP(), i=0; 
	     nodIter.IsActive(); 
	     cn=nodIter.NextP(), ++i )
	  { // add up burden/buttress pressures and lateral friction terms 
	    // at edges:                         
	    double nodeBurden=0.0;
	    sI.Reset(cn);
	    for( tEdge* ce=sI.FirstP(); !sI.AtEnd(); ce=sI.NextP() )
	      {
		const int iEdge = tempFullEdgeIndex[ce->getID()];
		nodeBurden += edgeBurden[iEdge];
		nodeLateralFriction[i] += edgeFrictionMag[iEdge];
	      }
	    const double FricSlopeDiff = 
	      fabs( frictionAngle - nodeSlopeAngle[i] );
	    // add sum of burden/buttress pressures to force balance:
	    nodeNetForceMag[i] += 
	      nodeBurden * GRAV * 
	      ( cos( FricSlopeDiff ) - sin( FricSlopeDiff ) * fricSlope );
	    // magnitude of force is positive downhill:
	    // add lateral friction (resistance) to force balance:
	    nodeNetForceMag[i] -= nodeLateralFriction[i];
	    // find vector by multiplying by components of gradient unit 
	    // vector, which is positive uphill; force magnitude is now 
	    // positive downhill; want force vector to point the right 
	    // direction, hence negative sign:
	    nodeNetForceX[i] = 
	      -nodeNetForceMag[i] * nodeGradientUnitVectorX[i];
	    nodeNetForceY[i] = 
	      -nodeNetForceMag[i] * nodeGradientUnitVectorY[i];
	    // record forces for nodes if first pass:
	    if( numPasses == 1 )
	      cn->setNetDownslopeForce( nodeNetForceMag[i] );
	  }
      } 
      // END finding force balance for each active node in mesh.
      // BEGIN finding landslide clusters:
      typedef vector<NodeNetForceIndex> NodeContainer;
      // construct a priority_queue type that uses the custom comparator:
      priority_queue<NodeNetForceIndex, NodeContainer, NodeNetForce_Lesser> 
	nodePQ;
      // put nodes in a priority_queue with greatest net downhill force 
      // at top.
      for(tLNode* cn=nodIter.FirstP(); nodIter.IsActive(); 
	  cn=nodIter.NextP())
	{
	  NodeNetForceIndex curNSFI;
	  const int iNode = tempNodeIndex[cn->getID()];
	  curNSFI.node = cn;
	  curNSFI.netForce = nodeNetForceMag[iNode];
	  curNSFI.index = iNode;
	  nodePQ.push( curNSFI );
	}
      { // pop top node in queue and build cluster...
	bool anyFailuresThisPass = true;
	int numCluster = 0;
	tSpkIter sI;
	tSpkIter neI;
	while( anyFailuresThisPass ) // search again if last cluster failed
	  {
	    anyFailuresThisPass = false;
	    // increment number used to flag nodes in a cluster:
	    ++numCluster;
	    tPtrList<tLNode> slideCluster;
	    tPtrList<tLNode> seedList; // list of nodes on edge of cluster
	    NodeNetForceIndex curNSFI;
	    do
	      {
		// start potential cluster with node with greatest net 
		// downhill force:
		curNSFI = nodePQ.top();
		// and remove it from the queue:
		nodePQ.pop();
	      } while( curNSFI.node->public1 > 0 ); // if already part of a 
	    // cluster, try again
	    seedList.insertAtBack( curNSFI.node );
	    slideCluster.insertAtBack( curNSFI.node );
	    double initWtGradUnitVecX = 
	      nodeGradientUnitVectorX[curNSFI.index] 
	      * curNSFI.node->getVArea();
	    double initWtGradUnitVecY = 
	      nodeGradientUnitVectorY[curNSFI.index] 
	      * curNSFI.node->getVArea();
	    double initSumWt = curNSFI.node->getVArea();
	    double initNetForceX = nodeNetForceX[curNSFI.index];
	    double initNetForceY = nodeNetForceY[curNSFI.index];
	    // dot-product of force vector and gradient unit vector:
	    double initNetForce = -curNSFI.netForce;
	    double initUphillUnitVectorX = 
	      nodeGradientUnitVectorX[curNSFI.index];
	    double initUphillUnitVectorY = 
	      nodeGradientUnitVectorY[curNSFI.index];
	    double initLatCohesion = nodeLatCohesion[curNSFI.index];
	    double initLatFriction = nodeLateralFriction[curNSFI.index];
	    double finalWtGradUnitVecX = initWtGradUnitVecX;
	    double finalWtGradUnitVecY = initWtGradUnitVecY;
	    double finalSumWt = initSumWt;
	    double finalNetForceX = initNetForceX;
	    double finalNetForceY = initNetForceY;
	    double finalNetForce = initNetForce;
	    double finalUphillUnitVectorX = initUphillUnitVectorX;
	    double finalUphillUnitVectorY = initUphillUnitVectorY;
	    double finalLatCohesion = initLatCohesion;
	    double finalLatFriction = initLatFriction;
	    // set generic flag signifying node in cluster:
	    curNSFI.node->public1 = numCluster;
	    while( !seedList.isEmpty() )
	      {
		// get seed for cluster growth
		tLNode* cn = seedList.removeFromFront();
		sI.Reset( cn );
		// search seed's neighbors
		for( tEdge* ce = sI.FirstP(); !sI.AtEnd(); ce = sI.NextP() )
		  if( ce->FlowAllowed() 
		      && ce->getDestinationPtr()->isNonBoundary() )
		    { // if not on boundary:
		      tLNode* nn = 
			static_cast<tLNode*>( ce->getDestinationPtrNC() );
		      // if node is not already in a cluster
		      // AND it's connected to the cluster by a flow edge:
		      if( nn->public1 == 0 && 
			  ( cn->flowThrough( ce ) || 
			    nn->flowThrough( ce->getComplementEdge() ) ) )
			{
			  // decrement force by lateral cohesion 
			  // (note, subtracting something that's also 
			  // pointing in the direction opposing motion, 
			  // so minus * minus means adding back in what 
			  // was subtracted before):
			  finalNetForceX -= 
			    initLatCohesion * initUphillUnitVectorX;
			  finalNetForceY -= 
			    initLatCohesion * initUphillUnitVectorY;
			  const int iNode = tempNodeIndex[nn->getID()];
			  // increment force by that of new node minus its 
			  // lateral cohesion (gradient unit vector points 
			  // uphill, so again minus * minus... ):
			  finalNetForceX += 
			    nodeNetForceX[iNode] 
			    - nodeLatCohesion[iNode] 
			    * nodeGradientUnitVectorX[iNode];
			  finalNetForceY += 
			    nodeNetForceY[iNode]
			    - nodeLatCohesion[iNode] 
			    * nodeGradientUnitVectorY[iNode];
			  // go through candidate's neighbors for lateral 
			  // cohesion and friction:
			  neI.Reset( nn );
			  for( tEdge* ne=neI.FirstP(); !neI.AtEnd(); 
			       ne=neI.NextP() )
			    if( ne->FlowAllowed() 
				&& ne->getDestinationPtr()->isNonBoundary() )
			      {
				tLNode* nnn = 
				  static_cast<tLNode*>(ne->getDestinationPtrNC());
				if( nnn->public1 > 0 )
				  { // if neighbor of neighbor already in 
				    // cluster, subtract edge's lateral 
				    // cohesion:
				    const int iEdge = 
				      tempEdgeIndex[ne->getID()];
				    finalLatCohesion -=
				      edgeLatCohesion[iEdge];
				    // and lateral friction for both edge and 
				    // complement:
				    const int iFricE = 
				      tempFullEdgeIndex[ne->getID()];
				    const int iFricC = 
				      tempFullEdgeIndex[ne->getComplementEdge()
							->getID()];
				    finalNetForceX -= 
				      edgeFrictionX[iFricE] 
				      + edgeFrictionX[iFricC];
				    finalNetForceY -= 
				      edgeFrictionY[iFricE] 
				      + edgeFrictionY[iFricC];
				  }
				else
				  // if neighbor of neighbor is not in cluster,
				  // add edge's lateral cohesion:
				  finalLatCohesion +=
				    edgeLatCohesion[tempEdgeIndex[ne->getID()]];
			      }
			  // find new force's unit vector (despite name, still 
			  // points in direction of net force here):
			  finalNetForce = 
			    sqrt( finalNetForceX * finalNetForceX 
				  + finalNetForceY * finalNetForceY );
			  finalUphillUnitVectorX = 
			    finalNetForceX / finalNetForce;
			  finalUphillUnitVectorY = 
			    finalNetForceY / finalNetForce;
			  // find weighted uphill gradient vector:
			  finalWtGradUnitVecX += 
			    nodeGradientUnitVectorX[iNode] * nn->getVArea();
			  finalWtGradUnitVecY += 
			    nodeGradientUnitVectorY[iNode] * nn->getVArea();
			  finalSumWt += nn->getVArea();
			  const double tmpGradX = 
			    finalWtGradUnitVecX / finalSumWt;
			  const double tmpGradY = 
			    finalWtGradUnitVecY / finalSumWt;
			  // find dot product of gradient and net force:
			  const double dotprod = 
			    finalNetForceX * tmpGradX 
			    + finalNetForceY * tmpGradY;
			  if( dotprod < 0.0 )
			    { // if force vector is opposite weighted gradient 
			      // vector (i.e., points downhill), make force 
			      // unit vector point uphill:
			      finalUphillUnitVectorX *= -1.0;
			      finalUphillUnitVectorY *= -1.0;
			    }
			  // determine new net downhill force with the 
			  // candidate by adding the lateral cohesion back in 
			  // as resistance to net force (i.e., pointing 
			  // uphill):
			  finalNetForceX += 
			    finalLatCohesion * finalUphillUnitVectorX;
			  finalNetForceY += 
			    finalLatCohesion * finalUphillUnitVectorY;
			  // final net force for comparison is dot product of 
			  // force vector with uphill unit vector; failure if 
			  // negative:
			  finalNetForce = 
			    finalNetForceX * finalUphillUnitVectorX 
			    + finalNetForceY * finalUphillUnitVectorY;
			  // if the initial force points downhill and the  
			  // final force points downhill after addiiton of 
			  // the new node OR if the initial force points 
			  // uphill and addition of the new node decreases 
			  // the uphill magnitude or makes it point downhill:
			  if( ( initNetForce < 0.0 && finalNetForce < 0.0 )
			      || ( initNetForce >= 0.0 && 
				   finalNetForce < initNetForce ) )
			    {
			      // then add the new node to the cluster:
			      nn->public1 = numCluster;
			      seedList.insertAtBack( nn );
			      slideCluster.insertAtBack( nn );
			      // update "initial" terms with "final" values:
			      initWtGradUnitVecX = finalWtGradUnitVecX;
			      initWtGradUnitVecY = finalWtGradUnitVecY;
			      initSumWt = finalSumWt;
			      initNetForceX = finalNetForceX;
			      initNetForceY = finalNetForceY;
			      initNetForce = finalNetForce;
			      initUphillUnitVectorX = finalUphillUnitVectorX;
			      initUphillUnitVectorY = finalUphillUnitVectorY;
			      initLatCohesion = finalLatCohesion;
			      initLatFriction = finalLatFriction;
			    }
			  else
			    { // if node not added, reset "final" terms
			      // with "initial" values:
			      finalWtGradUnitVecX = initWtGradUnitVecX;
			      finalWtGradUnitVecY = initWtGradUnitVecY;
			      finalSumWt = initSumWt;
			      finalNetForceX = initNetForceX;
			      finalNetForceY = initNetForceY;
			      finalNetForce = initNetForce;
			      finalUphillUnitVectorX = initUphillUnitVectorX;
			      finalUphillUnitVectorY = initUphillUnitVectorY;
			      finalLatCohesion = initLatCohesion;
			      finalLatFriction = initLatFriction;
			    }
			}
		    }
	      } // END cluster building for current seed node
	    // have cluster (which may be only the node popped off the 
	    // queue); check whether it fails:
	    if( initNetForce < 0.0 )
	      { // if failure, set logical to continue finding failure 
		// clusters:
		anyFailuresThisPass = true;
		// initialize new tDebrisFlow object with reference to 
		// cluster, copy of factor of safety, references to soil 
		// depth, wood depth, water depth, and node index vectors, 
		// and pointers to mesh and erosion objects (send negative 
		// of initNetForce, which is currently positive uphill, 
		// for consistency with other LandslideClusters and original 
		// convention of positive magnitude = net downhill force):
		tDebrisFlow *dfPtr = 
		  new tDebrisFlow( slideCluster, -initNetForce,
				   nodeSoilThickness, nodeWoodDepth, 
				   nodeWaterDepth, tempNodeIndex,
				   meshPtr, this );
		// put copy of new debris flow pointer in list:
		dfPList.insertAtBack( dfPtr );
	      }
	  } // END while( anyFailuresThisPass )
      }
      tDebrisFlow *dfPtr=0;
      // rarely used this way, but sets DF iterator to first pointer if first pass:
      for( dfPtr = dfI.NextP(); !dfI.AtEnd(); dfPtr = dfI.NextP() )
	{ // DFs added to list; run out debris flows:
	  anyFailuresThisStorm = true;
	  dfPtr->RunScourDeposit();
	  // output using overloaded "friend" operator:
	  if( DF_fsPtr )
	    *DF_fsPtr << time << " " << *dfPtr << endl;
	}
      dfI.Last(); // above loop will take it "past" last
      // need to end this loop after runout, since scour might 
      // unbuttress potential failures:
    } // END while( anyFailuresThisStorm )
  // unlike LandslideClusters(), debris flow list is not emptied 
  // during runout so that we can find agglomerated areas of adjacent
  // failures, so that failures triggered by failure and unbuttressing
  // will be added to the area of that initial failure; we want to be 
  // able to compare landslide areas to field observations, which cannot
  // distinguish between failures en masse and sequences of smaller 
  // adjacent failures:
  if( !dfPList.isEmpty() )
    {
      const int numDFlows = dfPList.getSize();
      vector<double> smallAreas( numDFlows );
      vector<tDebrisFlow*> dFlowPtrs( numDFlows );
      vector<int> dfFlags( numDFlows );
      tSpkIter sI;
      {
	int i=0;
	tDebrisFlow* dfPtr = 0;
	for( dfPtr = dfI.FirstP(); !dfI.AtEnd(); dfPtr = dfI.NextP() )
	  {
	    smallAreas[i] = dfPtr->getAreaFailure();
	    dFlowPtrs[i] = dfPtr;
	    dfFlags[i] = dfPtr->getOriginPtr()->public1;
	    ++i;
	  }
      }
      for( int i=0; i<numDFlows; ++i )
	if( smallAreas[i] > 0.0 )
	  { // hasn't already been added to another cluster, so
	    // initialize potentially larger cluster with this one:
	    tPtrList<tLNode> homeCluster( dFlowPtrs[i]->getSlideCluster() );
	    tPtrListIter<tLNode> hcI( homeCluster );
	    // check neighbors of all nodes in cluster
	    for( tLNode* cn = hcI.FirstP(); !hcI.AtEnd(); cn = hcI.NextP() )
	      {
		sI.Reset( cn );
		for( tEdge* ce = sI.FirstP(); !sI.AtEnd(); ce = sI.NextP() )
		  {
		    tLNode* nn = static_cast<tLNode*>( ce->getDestinationPtrNC() );
		    // if neighbor is part of a cluster,
		    // go through other failures to find if it's one of those:
		    if( nn->public1 > 0 )
		      for( int j=i; j<numDFlows; ++j )
			if( smallAreas[j] > 0.0 && nn->public1 == dfFlags[j] )
			  { // hasn't already been added to another cluster, and
			    // has a neighbor in the current cluster, so
			    // add jth cluster to "homeCluster":
			    tPtrListIter<tLNode> lsI( dFlowPtrs[i]->getSlideCluster() );
			    for( tLNode* sn = lsI.FirstP(); !lsI.AtEnd(); sn = lsI.NextP() )
			      homeCluster.insertAtBack( sn );
			    smallAreas[i] += smallAreas[j];
			    smallAreas[j] = 0.0;
			    break; // found a match, so bail on loop through other DFs
			  }
		  }
	      }
	  }
      // add agglomerated areas to landslideAreas and delete DFs:
      for( int i=0; i<numDFlows; ++i )
	{
	  if( smallAreas[i] > 0.0 )
	    landslideAreas.insertAtBack( smallAreas[i] );
	  tDebrisFlow* dP = dFlowPtrs[i];
	  dFlowPtrs[i] = 0;
	  delete dP;
	}
      // write time, storm, # failures, and failed volumes to file:     
      if( DF_Hyd_fsPtr )
	{
	  DF_Hyd_fsPtr->setf( ios::scientific, ios::floatfield );
	  DF_Hyd_fsPtr->precision(4);   
	  *DF_Hyd_fsPtr 
	    << time << " " << rainrate << " " << numDFlows << " " 
	    << debris_flow_sed_bucket << " " 
	    << debris_flow_wood_bucket << endl;
	}
    }
}
// END LandslideClusters3D

// Conceivably useful bit of code to find an average cross-slope node
// width using cross-products between edges and the gradient vector 
// for the cross-slope widths and dot-products between Voronoi edges and 
// the gradient for the slope-parallel weighting. 
// Written for use in LandslideClusters3D, but not needed.
// -SL, 11/10.
//   vector<double> nodeBasalWidth( numActiveNodes );
// 	// do sums to find average basal width perpendicular to gradient;
// 	// this is the sum( abs(edge vectors X gradient unit vector)
// 	// x abs(Vor. edge vectors "dot" gradient unit vector) )
// 	// / sum( abs(Vor. edge vectors "dot" gradient unit vector) ):
// 	double weightedSumEdgeCrossGrad = 0.0;
// 	double sumVEdgeDotGrad = 0.0;
// 	for( tEdge* ce=sI.FirstP(); !sI.AtEnd(); ce=sI.NextP() )
// 	  {
// 	    const double crossp = 
// 	      ce->getEVec().at(0) * nodeGradientUnitVectorY[i] 
// 	      - ce->getEVec.at(1) * nodeGradientUnitVectorX[i];
// 	    const double dotp = 
// 	      ce->getVVec().at(0) * nodeGradientUnitVectorX[i] 
// 	      + ce->getVVec.at(1) * nodeGradientUnitVectorY[i];
// 	    weightedSumEdgeCrossGrad += fabs( crossp ) * fabs( dotp );
// 	    sumVEdgeDotGrad += fabs( dotp );
// 	  }
// 	nodeBasalWidth[i] = weightedSumEdgeCrossGrad / sumVEdgeDotGrad;

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
