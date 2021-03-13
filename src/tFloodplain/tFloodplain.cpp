/**************************************************************************/
/**
**  @file tFloodplain.cpp
**  @brief functions for class tFloodplain and its helper class tFloodNode
**
**  Class tFloodplain simulates overbank deposition in the CHILD model,
**  using a modified version of Howard's (1996) diffusive overbank
**  deposition model.
**
**  The class provides a constructor, which reads in and initializes
**  all necessary parameters, and a DepositOverbank function that
**  implements floodplain deposition.
**
**  The overbank deposition model works as follows. Nodes with a
**  drainage area greater than a given threshold size are considered
**  large enough to initiate an overbank flood with a significant amount
**  of deposition. The minimum drainage area for a "flood node" is
**  read from the input file and recorded in drarea_min. The rate of
**  overbank deposition at any node depends on (1) its distance from the
**  nearest flood node (which is zero if it _is_ a flood node), and
**  (2) its elevation. The distance dependence is meant to simulate the
**  presumed dropoff in average sedimentation rate of suspended sand or
**  silt with distance from the main channel, and is modeled as an
**  exponential as in the original Howard model. The elevation dependence
**  is meant to simulate the increased likelihood of overbank flooding
**  in low spots versus high spots, as well as a presumed increase in
**  the total amount of sediment raining out of a deeper water column in
**  lower spots. In the original Howard model, the elevation dependence
**  was expressed as the difference between local elevation and a
**  maximum floodplain height above the channel. This was later modified
**  to use an exponential elevation weighting that accounts for a
**  distribution of event sizes. In the implementation here, we can take
**  advantage of the fact that we have an explicit distribution of event
**  sizes (i.e., the storm events that drive the model). Thus, the
**  elevation dependence is modeled as the difference between the local
**  elevation and the water surface height at the nearest flood node
**  during the current event. Water surface height is computed as the
**  sum of the flood node elevation plus water depth as computed from a
**  simple hydraulic geometry relation. Only events larger than a
**  specified precipitation rate event_min are assumed to generate
**  overbank deposition. Taken together, the deposition rate at a point
**  for a given event is (assuming P>Pe and WSH>Z):
**
**    F = MU * (WSH-Z) * exp( - Df / LAMDA )
**
**    WSH = Zf + Kdb * Pe^(mqb-mqs) * A^(mqb-mqs) * Q^mqs
**
**  where (with code variable names in parentheses):
**
**     MU (fpmu) = deposition rate per unit depth at distance = 0
**     WSH = water surface height
**     Z = elevation of node
**     Df = distance to nearest flood node
**     LAMDA (fplamda) = e-folding distance for overbank deposition
**     Zf = elevation of nearest flood node
**     Kdb = bankfull depth-discharge coefficient (ie, Db = Kdb Qb ^ mqb)
**     Pe (event_min) = assumed bankfull event precipitation rate
**     mqb = bankfull depth-discharge exponent
**     mqs = at-a-station depth-discharge exponent
**     A = drainage area of flood node
**     Q = flood node discharge
**
**  (Created 1/99 by GT)
**
**  $Id: tFloodplain.cpp,v 1.30 2005-03-15 17:17:29 childcvs Exp $
*/
/**************************************************************************/

#include "tFloodplain.h"
#include "../tListInputData/tListInputData.h"

/**************************************************************************\
**
**  @fn tFloodplain( tInputFile &infile, tMesh<tLNode *mp );
**  @brief Main constructor for tFloodplain
**
**  @param infile Input file from which parameters are read
**  @param mp Pointer to the mesh
**
**  Takes a tInputFile as an argument and reads from it the various
**  necessary parameters. For convenience and speed, mqbmqs = mqb - mqs
**  and kdb = Kdb * Pe^(mqb-mqs).
**
\**************************************************************************/
tFloodplain::tFloodplain( const tInputFile &infile, tMesh<tLNode> *mp )
  :
  chanDriver(0),
  meshPtr(mp)
{
   int numg;

   // Keep a pointer to the mesh in order to access list of nodes
   assert( meshPtr!=0 );

   // Read in parameters
   fpmode = infile.ReadItem( fpmode,"OPTFLOODPLAIN");
   drarea_min = infile.ReadItem( drarea_min, "FP_DRAREAMIN" );
   kdb = infile.ReadItem( kdb, "HYDR_DEP_COEFF_DS" );
   mqs = infile.ReadItem( mqs, "HYDR_DEP_EXP_STN" );
   mqbmqs = infile.ReadItem( mqbmqs, "HYDR_DEP_EXP_DS" ) - mqs;
   event_min = infile.ReadItem( event_min, "FP_BANKFULLEVENT" );
   infile.ReadItem( fpmuVariation, "FP_MU" );
   fplamda = infile.ReadItem( fplamda, "FP_LAMBDA" );

   kdb = kdb*pow( event_min, mqbmqs );

   if (fpmode==2)
     std::cout<<"Floodplain aggradation based on suspension " <<std::endl;

   //std::cout << "kdb: " << kdb << "  mqbmqs " << mqbmqs << std::endl;

   numg = infile.ReadItem( numg, "NUMGRNSIZE" );
   deparr.setSize(numg); // dimension & init to 0 the deparr array
   deparrRect.setSize(numg+1);

   optControlMainChan = infile.ReadBool( "FP_OPTCONTROLCHAN" );
   if( optControlMainChan )
     chanDriver = new tMainChannelDriver( infile );
}

/**************************************************************************\
**
**  @fn ~tFloodplain
**  @brief destructor
**
\**************************************************************************/
tFloodplain::~tFloodplain()
{
  meshPtr = 0;
  delete chanDriver;
}

/**************************************************************************\
**
**  tFloodplain::DepositOverbank
**
**  Implements the overbank deposition model using the following
**  algorithm:
**
**  IF precip > event_min (ie, overbank flood event)
**    Create "flood list" of all flood nodes, computing WSH for each and
**         recording the maximum WSH
**    FOR each landscape node
**      IF node elevation is below maximum WSH
**        Find the closest flood node
**        IF node elevation < WSH at closest flood node
**          Calculate total deposition depth and update node elevation
**
**  Modifications:
**   - hydraulic geometry parameters should be interpreted for discharge
**     units of meters per second, not meters per year. Corrected. GT 3/99
**   - now uses call to EroDep with layer info. deparr is passed to
**     EroDep; all overbank deposits are assumed to be finest fraction,
**     which is assumed to be the 1st size in the array. GT 3/13/99
**   - Deposition rate WITHIN main channel now reduced to 20% of
**     nominal rate, based on the argument that (a) fine material is
**     much less likely to settle out under the swift in-channel current,
**     and (b) the coarse bedload material, which will, is typically
**     ~20% of the total load. This is not ideal, but it allows
**     in-channel deposition even under "detachment limited" conditions,
**     and is a useful approximation for large-scale floodplain sim's.
**     GT 6/99.
**
**    Parameters:
**      precip -- precipitation rate for current storm event
**      delt -- duration of current storm event
**      ctime -- current time
**    Notes:
**
\**************************************************************************/
void tFloodplain::DepositOverbank( double precip, double delt, double ctime )
{
   if( precip < event_min ) return;

   tMesh< tLNode >::nodeListIter_t ni( meshPtr->getNodeList() ); // iterator for nodes
   tList<tFloodNode> floodList;    // list of "flood nodes"
   tFloodNode *fn;       // ptr to current flood node
   tLNode *cn,           // current landscape node
       *closestNode;     // closest flood node
   double maxWSH = 0.0,  // maximum water surface height at any flood node
       minDist,          // minimum distance to a flood node
       dist,             // distance to flood node
       floodDepth,       // local flood depth
       dx,dy,            // x and y distance to flood node
       wsh=0.0,          // water surface height
       drarea;           // drainage area at flood node

   //std::cout << "tFloodplain\n";

   // Make a list of all nodes with a drainage area large enough to be
   // considered "flood generators". Record the highest water surface height.
   // in these channel cells
   for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
   {
      if( ( drarea = cn->getDrArea() ) >= drarea_min )
      {

      	 // Debug:
	     if (0) //DEBUG
	        std::cout<< "Flood Nodes " <<cn->getX()<< ' ' << cn->getY()<<' '<<cn->getZ()
	           <<" Slp= "<<cn->calcSlope()<<" Q= "<<cn->getQ()
	           <<" W= "<<cn->getChanWidth()<< " MStatus="<<cn->Meanders()<<std::endl;

	     tFloodNode floodNode(cn,
			      kdb*pow( drarea, mqbmqs )
			      *pow( cn->getQ()/SECPERYEAR, mqs )
			      + cn->getZ() );
	     if (0) //DEBUG
	        std::cout << "flood depth at " << cn->getX() << ' ' << cn->getY()
		       <<' '<< cn->getZ() << " = " << floodNode.wsh-cn->getZ() << std::endl;
         if( floodNode.wsh > maxWSH )
             maxWSH = floodNode.wsh;
         floodList.insertAtBack( floodNode );
      }
   }

   // Just in case there are no flood nodes, stop here
   if( floodList.isEmpty() ) return;

   // For each node, find the nearest flood node and if it's WSH is above
   // the local elevation, deposit stuff
   
   double fpmu = fpmuVariation.calc(ctime);
   
   for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
   {
#define NOMAINCHANDEP 1 //Jan 05: try dis-allowing deposition within the main channel
#if NOMAINCHANDEP
      if( cn->getDrArea() < drarea_min ) {
#endif
      closestNode = 0;
      minDist = kVeryFar;
      if( cn->getZ() < maxWSH ) // don't bother if node is above max flood ht
      {
         // Sweep through the flood node list to find the closest flood node,
         // recording its distance and wsh
         for( fn=floodList.FirstP(); fn!=0; fn=floodList.NextP() )
         {
            //std::cout << "  check fn " << fn->nodePtr->getID() << std::endl << flush;
            dx = cn->getX() - fn->nodePtr->getX();
            dy = cn->getY() - fn->nodePtr->getY();
            dist = sqrt( dx*dx + dy*dy );
            if( dist<minDist )
            {
               minDist = dist;
               closestNode = fn->nodePtr;
               wsh = fn->wsh;
            }

         }
         assert( closestNode!=0 ); // (should always find one)
	     assert( wsh>0 );   // should always find a closest node & set wsh
         /*std::cout << " got it: " << closestNode->getID() << " dist=" << minDist
           << std::endl << flush;*/

         // If current node is below flood level, do some deposition
         // (note: deposit is assumed to consist of 100% of the finest
         // grain size fraction, which is assumed to be the first fraction.
         // All other entries in deparr are zero)
         //   NOTE: when delt > 1/mu exp( minDist/fplamda ), the solution
         // is invalid because deposition occurs above the water surface!
         // An analytical solution znew = z + D0(1-exp(-fx delt)) could
         // be used in these cases (ie, when delt is large to speed up
         // computation)
         if( ( floodDepth = wsh - cn->getZ() ) > 0.0 )
         {
            // Depth of deposition: if in a main channel, use 20% of
            // dep rate as the bedload contribution, assuming that susp load
            // is not deposited in main channel (obviously an approximation!);
            // otherwise, use Howard formula

#if !NOMAINCHANDEP
            // In the channel..
            if( minDist==0.0  ){

              if (fpmode==2) {	               // Suspension concentration based
            	deparrRect[0] = 0.0;
            	deparrRect[1] = 0.0;
            	deparr[0]     = 0.0;
              } 
			  else {                         // Howard, allows some sedimentation into the channel
                deparrRect[0] = 0.001*floodDepth*fpmu*delt;      // coarse
                deparrRect[1] = 0.0;                             // fine
                deparr[0]=0.001*floodDepth*fpmu*delt;           // coarse
              }
	          if (0) //DEBUG
		         std::cout<<"FP-Channel: wsh= "<< wsh <<" flooddepth= "<<floodDepth
		            <<" dh_channel= "<<deparrRect[0] << '\n';
            }
            // On the floodplain...
            else {
#endif
               if (fpmode ==2) {                // Suspension concentration based, Gross and Small, 1998
		          //std::cout<<"Calling FloodplainDh"<<std::endl;
		          deparrRect[0] = 0.0;
		          deparrRect[1] = FloodplainDh2(minDist,floodDepth,delt,fplamda);
               } 
			   else {                         // Geometrical, Howard, 1996
		          deparrRect[0] = 0.0;                         // coarse
		          deparrRect[1] = floodDepth*fpmu*exp( -minDist/fplamda )*delt;
		          deparr[0]= floodDepth*fpmu*exp( -minDist/fplamda )*delt;
		          //double suspensiontest = FloodplainDh2(minDist,floodDepth,delt,fplamda);
		          //std::cout<<"At dist "<<minDist<<" FP_1= "<<deparrRect[1]<<", FP_2= "<<suspensiontest<<std::endl;
               }

	           //if(minDist <=100.){
	           //  std::cout<<"FP-Overbank at: " <<minDist<<" wsh=" << wsh
	           //      <<" flooddepth= "<<floodDepth <<" dh_overbank= "<<deparrRect[1] << '\n';
	           //}

#if !NOMAINCHANDEP
	        }
#endif

            if( deparr[0]>floodDepth)
	           std::cout << " *WARNING, deposit thicker than flood depth\n";

	        // Modify heights and communicate to stratigraphy tStratGrid
	        cn->IncrementAccummulatedDh(deparrRect);           // new version, recieves 2d arry
	        cn->EroDep( 0, deparr, ctime );                    // this one recieves 1 value
         }
      }
#if NOMAINCHANDEP
      }
#endif
#undef NOMAINCHANDEP
   }

   if (0) //DEBUG
     std::cout << "Floodplain:: done Overbanks...\n";
}

/************************************************************\
 ** FloodplainDh
 **
 **
 \*************************************************************/
double tFloodplain::FloodplainDh(double minDist, double flooddepth,
				 tLNode *fpnode) const
{
  const double C  = getSuspendedConcentration(minDist);
  double dh = ConcentrationToHeight(flooddepth,fpnode,C);

  //DEBUG, what does it produce
  if (0) //DEBUG
    std::cout<<"For node" <<fpnode->getX()<<","<<fpnode->getY()
	<<" at dist "<<minDist
	<<", with flooddepth "
	<<flooddepth<<", dh = "<<dh<<std::endl;

  if(dh > 0.5 * flooddepth){   // never deposit more than 50 % of the local flooddepth!
    dh = 0.5* flooddepth;
  }

  return dh;

}

/*********************************************************\
 **
 ** FloodplainDh2
 **
 ** An alternative way of aggrading the floodplain wings.
 ** Now, dh is dependent on sediment concentration,
 ** which can be linked (hopefully) to emperical data of
 ** water discharge in the main channel.
 **
 ** dCy/dt = Ey*(d^2Cy/dy^2) - (ws/Z)*Cy
 **
 ** follwing Gross and Small, 1998, in WRR, no 34-9/pp. 2365-2376
 **
 ** Quintijn Clevis November 2003 - January 2004
\**********************************************************/
double tFloodplain::FloodplainDh2(double y, double flood_depth,
				  double /*delt*/, double Wy) const
{
  // Arguments:
  // y = closest distance to the meander channel
  // flood-depth = local water depth
  // delt = duration of flood = storm, in yrs
  // Wy = maximum width of the flooded floodplain, e.g lamda
  //
  // Table of usable values:
  // 		 median dPi (mm)	Fall velocity Ws (m/s)		Diff coeff Es (m/s^2)
  // gravel          4.0			0.280				0.31
  // fine sand       0.5			0.070				0.12
  // fine silt       0.016			0.0002				0.10
  //

  double ws=0.0;
  //double Ey=0.0;
  double Co=0.0;
  //double density=0.0;
  //double G=0.0;
  //double n=0.0;
  double Cy=0.0;
  double Dy=0.0;

  // Hard-coded variables (rates etc.)
  ws = 0.002;       // downward settling velocity in m/s for something between fine sand and fine silt
  //Ey = 0.10 ;       // transverse diffusion coefficient in m^2/s for fine silt
  Co =	1.0 ;       // suspended load concentration integrated over the flow depth
  		    // above the levee, in parts per million times meter, ppm-m
  		    // assumtion used in Gross and Small, use 75% of the in-channel concentration
  //density = 1.8e6;   // density of silt in suspension, g/m^3

  //G  = Wy * sqrt(   ws/flood_depth*Ey  );	// gross and small,1998

  //Ey = 0.013;
  //double Ez = 0.0076;
  //double Vs = 0.008;
  //double G2 = (Wy*Vs)*sqrt(Ey*Ez);                     //floodplain concavity Pizzuto,1987

  //n = y/Wy;

  // Calculate the Local concentration
  //Cy = Co * (  ( (sinh(G*n) * exp(-G))  /cosh(G) ) + exp(-G*n)  );
  Co = 0.02;
  //Cy = Co * (  ( (sinh(G*n) * exp(-G))  /cosh(G) ) + exp(-G*n)  );

  Cy= Co * exp(-y/Wy);				// concentration in per m, column water thickness

  // Depth deposited is estimated as the amount of excess sediment that settles
  // within the discrete period (storm)
  // DY = [Cy (ppm-m) * Ws (m/s) * t (s) ] /density (g/m^3) * flooddepth (m)]
  // where ppm-m is part per million, times meter


  //double delt_sec = delt*SECPERYEAR;
  //Dy = (Cy * ws * delt_sec)/ (density * flood_depth);
  //Dy = (pow(Vs,2.0)*Cy*delt_sec)/Ez;
  Dy = Cy * flood_depth * ws;         // (m^-1) * (m) * (s) ?

  return Dy;
}

/**************************************************************\
 ** tFloodplain::getSuspendedConcentration
 **
 ** Calculates the supended sediment concentration in mg/l
 ** at a certain position in the floodplain, with distance
 ** mindist from the closest channel node
 \**************************************************************/
double tFloodplain::getSuspendedConcentration(double minDist) const
{
  //double Co = k*(channelnode->getQ)

  double Co = 100.0;				// hard coded sediment concentration !
  double C = Co * exp (-minDist/fplamda);

  return C;
}

/*********************************************************\
**   ConcentrationToHeight
**   Calculate the dry thickness of the sediment suspended
**   present in the floodwater column
**
\*********************************************************/

double tFloodplain::ConcentrationToHeight(double flooddepth, tLNode *fpnode,
					  double C) const
{
  double Varea = fpnode->getVArea();

  double density = 1.2;                            // density of the silt (g/cm^3)
  double M3toLit = 1000.0;			    // convert M^3 to liters                                             // suspended sediment concentration         (mg/l)
  double nliters = (flooddepth * Varea)* M3toLit;  // nliters of floodwater on this cell       (l)
  double mgsed =  nliters*C;                       // mg suspended sed in that column          (mg)
  double gsed =   mgsed/1000.0;                    // grams                                    (g)
  double cubCmsed =  gsed/density;                  // n cubic cm's sediment                   (cm^3)
  double cubMsed = cubCmsed/1000000.0;		     // cubic m's sediment                      (m^3)
  double sedheight = cubMsed/Varea;		     // height of dry sed in the water column   (m)


  // Or use something as Rd (rate of sedimentation) = C*Vsettling*(1-To/Tc)   (kgm^2s-1

  //DEBUG
  if(sedheight < 0.0 || sedheight > 10.0){
    std::cout<<"ERROR in tFloodplain::CalculateSedHeight, very thick dry sediment column "<<std::endl;
    std::cout<<"sedheight = "<<sedheight<<std::endl;
  }

  return sedheight;
}

/**************************************************************************\
 **
 ** @tMainChannelDriver::Raisebanks
 **
 ** Function that ensures that the elevation of the meandering channel is
 ** always lower than the elevation of the surrounding banks, by raising the
 ** banks artificially.
 **
 ** Called by: tMainchannelDriver::UpdateMainChannelElevation
 ** Calls:     no functions
 **
 ** Created 7/03 (QC)
\**************************************************************************/
void tMainChannelDriver::RaiseBanks(double future_bed, tLNode *cn, tLNode *upstrN, double tm )
{
  assert(cn);
  assert(upstrN);

  tLNode *dstrN     = cn->getDownstrmNbr();
  tArray<double> delz( num_grnsize_fractions );
  tArray<double> delzRect( num_grnsize_fractions + 1);

  // Water depth during this timestep Flood
  double Flooddepth = kdb*pow( cn->getDrArea(), mqbmqs )
    *pow( cn->getQ()/SECPERYEAR, mqs );

  double HFD=0.5*Flooddepth;

  tEdge *ce;
  tSpkIter spokIter( cn );
  //std::cout<<"Channel node is: "<<cn->getID()<<' '<<cn->getX() <<' '<<cn->getY()<<' '<<cn->getZ()<<'\n';

  for( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() )
    {
      if(ce->getDestinationPtr() != NULL){
	if(ce->getDestinationPtr() != dstrN && ce->getDestinationPtr()!=upstrN){

	  tLNode* banknode = static_cast<tLNode *>( ce->getDestinationPtrNC() );
	  if(banknode != NULL){
	    // Make sure w're not blocking an active meander node:
	    if(banknode->Meanders() != 1 && banknode != dstrN && banknode!=upstrN && banknode->getZ() < future_bed + HFD){
	      delzRect[0] = 0.0;
	      delzRect[1] = (future_bed + HFD) - banknode->getZ();  //raise bank with fines
	      delz[0] = (future_bed + HFD) - banknode->getZ();

	      banknode->IncrementAccummulatedDh(delzRect);
	      banknode->EroDep( 0, delz, tm );                      // normal mesh
	    }
	  } // != NULL


	} // != dstrN || ustrN
      } //!= NULL

    } // spoke loop
}

/**************************************************************************\
**
**  @fn OptControlMainChan()
**  @brief Returns TRUE if user has requested control of main channel elev.
**
\**************************************************************************/
bool tFloodplain::OptControlMainChan() const { return optControlMainChan; }


/**************************************************************************\
**
**  @fn UpdateMainChannelHeight()
**  @brief Calls function in tMainChannelDriver to set main channel height.
**
\**************************************************************************/
void tFloodplain::UpdateMainChannelHeight( double tm, tLNode * inletNode )
{
  if (0) //DEBUG
    std::cout << "Floodplain:: start Updating Main Channel..."<<std::endl;

  chanDriver->UpdateMainChannelElevation( tm, inletNode );

  if (0) //DEBUG
    std::cout << "Floodplain:: Updated Main Channel..."<<std::endl;
}



/**************************************************************************\
**
**  @fn tMainChannelDriver( tInputFile &infile )
**  @brief tMainChannelDriver constructor
**
**  @param infile Reference to the input file for the run
**
\**************************************************************************/
tMainChannelDriver::tMainChannelDriver( const tInputFile &infile )
{
  drop = infile.ReadItem( drop, "FP_VALDROP" );
  num_grnsize_fractions = infile.ReadItem( num_grnsize_fractions,
					   "NUMGRNSIZE" );
  if( num_grnsize_fractions <= 0 )
    ReportFatalError( "NUMGRNSIZE must be >= 1" );

  infile.WarnObsoleteKeyword("FP_OPTCHANALITUDEVARIATION", "FP_INLET_ELEVATION");
  infile.WarnObsoleteKeyword("FP_PERIOD", "FP_INLET_ELEVATION");
  infile.WarnObsoleteKeyword("FP_AMPLITUDE", "FP_INLET_ELEVATION");

  infile.ReadItem( InletElevationVariation, "FP_INLET_ELEVATION");

   kdb = infile.ReadItem( kdb, "HYDR_DEP_COEFF_DS" );
   mqs = infile.ReadItem( mqs, "HYDR_DEP_EXP_STN" );
   mqbmqs = infile.ReadItem( mqbmqs, "HYDR_DEP_EXP_DS" ) - mqs;

   // Create the file used EVERY timestep for writing channel sinuosity, channel belt width etc.

   char fname[87];
#define THEEXT ".meander"
   infile.ReadItem( fname, sizeof(fname)-sizeof(THEEXT), "OUTFILENAME" );
   strcat( fname, THEEXT );
#undef THEEXT
   meanderfile.open( fname );
   if( !meanderfile.good() )
     std::cerr << "Warning: unable to create meander geometry data file '"
	       << fname << "'\n";
}


/**************************************************************************\
**
**  @fn UpdateMainChannelElevation( double tm )
**  @brief Updates elevations along the main stream channel.
**
**  @param tm Current time in simulation
**
**  This function sets elevations along the main channel as a function of
**  (a) time (via sinusoidal, step-wave, or user-specified variation), and
**  (b) a user-specified valley gradient (contained in the parameter
**      "drop", which is the height drop between the top and bottom of
**      the main valley.
**
**  Created: GT, May 2003
**
**  Assumes: (1) delz array elements initialized to zero;
**           (2) highest grain size fraction is the coarsest
**
\**************************************************************************/
void tMainChannelDriver::UpdateMainChannelElevation( double tm,
						     tLNode * inletNode )
{
  int counter;
  double newInletElevation,
    elev,
    chanslp;
  tEdge * fe;    // Ptr to flow edge of current node
  tLNode * cn /* ,*pn */ ;   // Current node along main channel

  assert( inletNode != 0 );



  // Update elevation at head of channel reach (inlet)
  newInletElevation = InletElevationVariation.calc( tm );

  // Compute length and slope of channel
  cn = inletNode;

  if(cn->getZ() > newInletElevation){
    std::cout<<"Channel goes down " <<cn->getZ()-newInletElevation<< " m"<<std::endl;
  } else if(cn->getZ()< newInletElevation){
    std::cout<<"Channel goes up " << newInletElevation-cn->getZ()<< " m"<<std::endl;
  }

  counter=0;
  double totlen = 0.0;  // Total length of channel found so far
  double maxY = inletNode->getY();
  double minY = inletNode->getY();
  do
    {
      fe = cn->getFlowEdg();
      assert( fe );
      totlen += fe->getLength();

      if(cn->getY() > maxY){ maxY = cn->getY(); }
      if(cn->getY() < minY){ minY = cn->getY(); }

      //DEBUG
      if( cn->getZ() < cn->getDownstrmNbr()->getZ()  ){
      	std::cout<<"  "<<std::endl;
      	std::cout<<"Bump In Floodplain loop "<<std::endl;
      	std::cout<<"for ID "<<cn->getID()<<" z= "<<cn->getZ()<<std::endl;
      	std::cout<<"Its downstream " << cn->getDownstrmNbr()->getID() << " z = "<<cn->getDownstrmNbr()->getZ()<<std::endl;
      	std::cout<<"  "<<std::endl;
      	//exit(1);
      }
      //END DEBUG

      counter++;

      cn = cn->getDownstrmNbr();
      assert( cn );
      if (0) //DEBUG
	std::cout<< "Meander Nodes " << cn->getID() << ' ' << cn->getX()<< ' ' << cn->getY()<<' '<<cn->getZ()
	    <<" Q= "<<cn->getQ()<<" W= "<<cn->getChanWidth()
	    << " MStatus="<<cn->Meanders()<<std::endl;


    }
  while( cn->getBoundaryFlag()==kNonBoundary && counter < 1000  );

  if(counter == 1000){                                                    //DEBUG
    std::cout<<"   \n";
    std::cout<<"--------------------------------------------------------------------\n";
    std::cout<<" WARNING:Channel does not find a boundary in function               \n";
    std::cout<<" tMainChannelDriver::UpdateChannelElevation in floodplain.cpp       \n";
    std::cout<<" 								       \n";
    std::cout<<" Recursion along meander path ?				       \n";
    std::cout<<" When this happens might also have problem in tLNode->getSlope()    \n";
    std::cout<<"                                                                    \n";
    std::cout<<" Showing the first 100 nodes below:                                 \n";
    std::cout<<"--------------------------------------------------------------------\n";

    int counter2=0;
    cn=inletNode;
    do
      {
	fe = cn->getFlowEdg();
	assert( fe );
	totlen += fe->getLength();
	counter++;

	cn = cn->getDownstrmNbr();
	assert( cn );
	std::cout<< "M-Nodes " <<cn->getID()<<' '<<cn->getX()<< ' ' << cn->getY()<<' '<<cn->getZ()<<" Q= "<<cn->getQ()<<" W= "<<cn->getChanWidth()<< " MStatus="<<cn->Meanders()<<std::endl;
	counter2++;
      }
    while( cn->getBoundaryFlag()==kNonBoundary && counter2 < 100  );
    //exit(1);

  } //end debug

  // FIXME BIG HACK
  if(counter==1000) 
  {
	totlen=7000.;
	std::cout << "INVOKING BIG HACK IN tMainChannelDriver::UpdateMainChannelElevation()" << std::endl;
  }
  chanslp = drop / totlen;


  std::cout  << "Channel Belt Geometry:       "<< std::endl;
  std::cout  << "1) Channel length= "<<totlen << " m"<< std::endl;
  std::cout  << "2) Sinuosity = "<<totlen/5000.0 << " - "<< std::endl;
  std::cout  << "3) Channel belt width = "<< maxY - minY<< " m"<< std::endl;
  std::cout  << "4) Chanslope ="<<chanslp<<" - "<< std::endl;
  std::cout  << "5) InletElev ="<<inletNode->getZ()<<std::endl;
  std::cout  << "       "<< std::endl;


  meanderfile<<tm<<" "<<totlen<<" "<<totlen/5000.<<" "<<maxY-minY
	     <<" "<<chanslp<<inletNode->getZ()<<std::endl;

  // Update elevations of all channel points below inlet
  tArray<double> delz( num_grnsize_fractions );           // for the triangular nodes
  tArray<double> delzRect( num_grnsize_fractions + 1);    // CAREFUL, temp Hack by Q., this way you have 2 classes in StratGrid, but only one in the mesh
  elev = newInletElevation;
  cn = inletNode;
  //pn = cn;
  //std::cout<< "Inlet Node = " <<cn->getX()<< ' ' << cn->getY()<<std::endl;


  do
    {
      //RaiseBanks(elev,cn,pn,tm);

#define BUGFIX 1   // undo unwanted side effects from hack below (2 places)
#if BUGFIX
      if(elev - cn->getZ() >= 0.0){
#else
      if(elev - cn->getZ() > 0.0){
#endif
      	delzRect[0] = elev - cn->getZ();   // coarse
      	delzRect[1] = 0.0;                 // fine

      	delz[0] = elev - cn->getZ();

      }
      // 'Erosion' because channel goes down
      // Erode only the fines, you have enough of these
      else if( elev - cn->getZ() < 0.0){
      	delzRect[0] = 0.0;
      	delzRect[1] = elev - cn->getZ();

      	delz[0] = elev - cn->getZ();
      }

      //std::cout <<" Channeldriver "<< cn->getX() << ' '<< cn->getY()<<"dh0= "<< delzRect[0] << " dh1= " <<delzRect[1] << '\n';

      if(0) //DEBUG
	     if( cn->getID()==8121 || cn->getID()==8122 )
		    std::cout<<"UpdateMCE: "<<cn->getID()<<" z before "<<cn->getZ()
			      <<" flowedg len "<<cn->getFlowEdg()->getLength()<<" elev "<<elev<<" dz "<<elev-cn->getZ()
				  <<" delz "<<delz[0]<<std::endl;

      cn->IncrementAccummulatedDh(delzRect);
      cn->EroDep( 0, delz, tm );

      if(0) //DEBUG
	     if( cn->getID()==8121 || cn->getID()==8122 )
		    std::cout<<"UpdateMCE: "<<cn->getID()<<" z after "<<cn->getZ()
			      <<" flowedg len "<<cn->getFlowEdg()->getLength()<<" elev "<<elev<<" dz "<<elev-cn->getZ()
				  <<" delz "<<delz[0]<<std::endl;

      // Calculate new elevation of downstream node. 
	  // THE FIX DESCRIBED BELOW HAS UNWANTED SIDE EFFECTS - DISABLING! -GT 1/05
	  // However
      // First, find the local slope along this segment. If the local slope is much large than the
      // average slope, we most likely had an avulsion, and the river is flowing down
      // laterally over one of the floodplain/levee wings. In this case DO NOT raise
      // the channel bed.

#if BUGFIX
      fe = cn->getFlowEdg();
      assert( fe );
      elev = elev - chanslp * fe->getLength();          // this will raise the slope of downstream nodes
#else
      double local_slope = (cn->getZ()) - (cn->getDownstrmNbr()->getZ())/ cn->getFlowEdg()->getLength();

      if(local_slope < 3.0*chanslp){
        fe = cn->getFlowEdg();
        assert( fe );
        elev = elev - chanslp * fe->getLength();          // this will raise the slope of downstream nodes
      }
      else if(local_slope >= 3.0*chanslp){
        elev= cn->getDownstrmNbr()->getZ();              // this will keep them, locally,  unchanged
      }
#endif
#undef BUGFIX

      // Move to downstream node
      //pn = cn;
      cn = cn->getDownstrmNbr();
      assert( cn );
    }
  while( cn->getBoundaryFlag()==kNonBoundary
	 && counter < 1000 );              // until we hit the boundary

  // FIXME BIG HACK
  if(counter==1000) 
  {
	totlen=7000.;
	std::cout << "INVOKING __SECOND__ BIG HACK IN tMainChannelDriver::UpdateMainChannelElevation()" << std::endl;
  }
  chanslp = drop / totlen;


}
