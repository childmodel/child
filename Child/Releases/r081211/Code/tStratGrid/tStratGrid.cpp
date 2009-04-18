/**************************************************************************/
/**
 **  @file tStratGrid.cpp
 **  @brief functions for class tStratgrid
 **
 **  Class tStratGrid creates a second mesh type, an equidistant rectangular grid
 **  underneath the standard Child triangular mesh. The idea behind StratGrid
 **  is to have a simple static rectangular grid available for connection with
 **  the alluvial Stratigraphy linked lists and their layer bookkeeping during
 **  the evolution of the meandering river. In contrast to the stratigraphic
 **  interpolation previously used while relocating meander bends,this simple
 **  data structure should avoid layer smearing and incorrectly shuffeling the
 **  deck of stratigraphic layers, or at least I hope so....
 **
 **  (Created 5/2003 by QC, AD and GT)
 **
 **  $Id: tStratGrid.cpp,v 1.18 2005/03/15 17:17:30 childcvs Exp $
 */
/**************************************************************************/
#include <assert.h>
#include <math.h>
#include "../errors/errors.h"
#include "tStratGrid.h"
#include "../tLNode/tLNode.h"
#include "../tMesh/tMesh.h"

#include <iostream>


/**************************************************************************\
 **								 tSTRATGRID
 **  @fn tStratGrid( tInputFile &infile, tMesh<tLNode *mp)
 **  @brief Main constructor for tStratGrid
 **
 **  @param infile Input file from which parameters are read
 **  @param mp Pointer to the mesh
 **
 **  Takes the tInputFile as an argument and reads from it the various
 **  necessary parameters and spans the grid. The necessary input values are:
 **
 **  xllc = x(m) lower left corner
 **  yllc = y(m) lower left corner
 **
 **  griddx   = distance(m) between the tStratGrid nodes, ideally these should
 **  		   be close or slightly smaller than the width of the main channel
 **  grwidth  = width (m) of the floodplain valley
 **  grlength = length(m) of the floodplain valley, parallel to the axis of the river
 **
 **
\**************************************************************************/
tStratGrid::tStratGrid( tInputFile const &infile, tMesh<tLNode> *mp_)
  : mp(mp_), StratNodeMatrix(0),
    StratConnect(0),
    imax(0), jmax(0)    //,tArray<double>surface(0),tArray<double>subsurface(0)
{
  // Keep a pointer to the mesh in order to access list of nodes
  assert( mp!=0 );

  int i,j,k;					  // x and y position indices (-)
  double x,y;					  // x and y coordinates (m)
  int stepx,stepy;
  tMesh<tLNode>::nodeListIter_t ni( mp->getNodeList() ); // iterator for tLNodes

  i=0;
  x=y=0.;

  // Read in values related to dimensions and resolution of the mesh
  // and desired output format of the stratigraphic sections
  xcorner  = infile.ReadItem(xcorner,"XCORNER");
  ycorner  = infile.ReadItem(ycorner,"YCORNER");
  griddx   = infile.ReadItem(griddx,"GRIDDX");
  grwidth  = infile.ReadItem(grwidth,"GR_WIDTH");
  grlength = infile.ReadItem(grlength,"GR_LENGTH");

  int endtime  = infile.ReadItem(endtime,"RUNTIME");
  nWrite   = infile.ReadItem(nWrite,"OPINTRVL");

  int arraysize = endtime/nWrite;
  surface.setSize(arraysize+1);           // array with timeslice specific surface area
  subsurface.setSize(arraysize+1);        // array with timeslice specific subsurface cummulative height
  subsurface_mbelt.setSize(arraysize+1);  // same, but for zone close to the channel
  outputTime.setSize(arraysize+1);        // array with the storm dependent output times

  for(i=0;i<arraysize;i++){	          // preservation potential arrays
    setSurface(i,0.0);
    setSubsurface(i,0.0);
    setSubsurface_mbelt(i,0.0);
    setOutputTime(i,0.0);
  }

  imax=int(grwidth/griddx);
  jmax=int(grlength/griddx);

  std::cout<<"     "<<std::endl;
  std::cout <<"StratGrid: number of nodes in x-direction = " << imax <<'\n';
  std::cout <<"StratGrid: number of nodes in y-direction = " << jmax <<'\n';


  // Call the constructor for the matrix of stratNodes
  StratNodeMatrix = new tMatrix<tStratNode>(imax,jmax);
  StratConnect = new tMatrix<tTriangle*>(imax,jmax);

  // Fill one stratnode with the initial layerlist properties
  const tStratNode a_StratNode(infile);

  // Construct the Grid by assigning coordinates and the initial layerlist
  // to the nodes present in the tStratNode Matrix
  for(i=0;i<imax;i++){
    for(j=0;j<jmax;j++){

      x = xcorner + i*griddx;
      y = ycorner + j*griddx;

      (*StratNodeMatrix)(i,j) = a_StratNode;    //to copy the layer list of a_stratnode in every node ?

      (*StratNodeMatrix)(i,j).setX(x);
      (*StratNodeMatrix)(i,j).setY(y);
      (*StratNodeMatrix)(i,j).setI(i);
      (*StratNodeMatrix)(i,j).setJ(j);
      (*StratNodeMatrix)(i,j).setZ(0.);
    }
  } // Loop over nodes....

  // Output is done in the form of 10 cross-section. 5 in the x-direction
  // 5 in the y-direction. Write the locations to an array which is data
  // member of Stradgrid and readible in tOutput
  stepx=int(imax/5.0);
  stepy=int(jmax/5.0);

  std::cout<< "Locations of the 10 cross-sections are \n";
  section.setSize(10);
  for(k=0;k<5;k++){
    section[k] = int(0.5*stepx + k*stepx);
    std::cout <<"X-Section "<<k+1<<" at i= " <<section[k] <<" or X= "<<(*StratNodeMatrix)(section[k],0).getX()<< '\n';
  }
  for(k=5;k<10;k++){
    section[k] = int(0.5*stepy + (k-5)*stepy);
    std::cout <<"Y-Section "<<k+1<<" at j= " <<section[k] <<" or Y= " <<(*StratNodeMatrix)(0,section[k]).getY()<< '\n';
  }


  // Build StratConnect
  if (1) { //DEBUG
    std::cout
      <<"   \n"
      <<"Building the StratConnect table of triangles in constructor...."
      <<std::endl;
  }
  updateConnect();
  if (1) { //DEBUG
    std::cout
      <<"Building the StratConnect table of triangles in constructor finished"
      <<"\n    "<<std::endl;
  }

  // Initialize the elevations of the StratNodes by interpolating between the
  // Triangles of the tMesh
  if (1) //DEBUG
    std::cout<<" Initializing tStratGrid elevations by interpolation, for the first Time "<<std::endl;
  InterpolateElevations();

  setSectionBase();   // DEBUG FUNCTION, all stratnodes have to know their initial, stratigraphy basis
  if (1) {//DEBUG
    std::cout
      <<" Finished Initializing tStratGrid elevations by interpolation, for the first Time "
      <<"\n    "<<std::endl;
  }
} // tStratGrid constructor

/**************************************************************************\
 **
 **  @fn tStratGrid
 **  @brief destructor
 **
\**************************************************************************/
tStratGrid::~tStratGrid()
{
  mp = 0;
  delete StratNodeMatrix;
  delete StratConnect;
}

// SET SECTION BASE
// This is a debug function, records the elevations of the stratgrid nodes
// at the beginning of a simulation.
void tStratGrid::setSectionBase()
{
  for(int i=0; i<imax; ++i){
    for(int j=0; j<jmax; ++j){
      tStratNode &sn = (*StratNodeMatrix)(i,j);
      sn.setSectionBase(sn.getZ());
    }
  }
}

/**************************************************************************\
 **
 **  @class
 **  @brief Find a rectangular box in the stratigraphy grid contained within
 **  a given triangle.
 **
 ** AD - 26 March 2004
\**************************************************************************/
class TriBox{
  TriBox();
public:
  int imin;
  int imax;
  int jmin;
  int jmax;
  TriBox(tTriangle const *, double, double, double);
  bool containsNone() const {
    return (imin > imax) || (jmin > jmax);
  }
};

TriBox::TriBox(tTriangle const *ct,
	       double xcorner, double ycorner, double griddx){
  tNode const * const node1Ptr = ct->pPtr(0);
  double
    maxx = node1Ptr->getX(),
    maxy = node1Ptr->getY();
  double
    minx = maxx,
    miny = maxy;
#define COMPUTE_MINMAX(NODEID) \
    do { \
      tNode const * const nodePtr = ct->pPtr(NODEID); \
      const double xx = nodePtr->getX(); \
      const double yy = nodePtr->getY(); \
      maxx = max( maxx, xx ); \
      maxy = max( maxy, yy ); \
      minx = min( minx, xx ); \
      miny = min( miny, yy ); \
    } while(0)

  COMPUTE_MINMAX(1);
  COMPUTE_MINMAX(2);
#undef COMPUTE_MINMAX
  imin = int( ceil((minx-xcorner)/griddx));
  imax = int(floor((maxx-xcorner)/griddx));
  jmin = int( ceil((miny-ycorner)/griddx));
  jmax = int(floor((maxy-ycorner)/griddx));
}

/**************************************************************************\
 **
 **  tStratGrid::updateConnect
 **  @brief update connectivity table StratConnect
 **
\**************************************************************************/
void tStratGrid::updateConnect()
{
  // nullify table
  {
    for(int i=0; i<imax; ++i)
      for(int j=0; j<jmax; ++j)
	(*StratConnect)(i,j) = NULL;
  }
  //
  tTriangle *ct;
  tMesh< tLNode >::triListIter_t triIter( mp->getTriList() );
  for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() ) {
    // Find box containing the current triangle
    const TriBox thisBox( ct, xcorner, ycorner, griddx );
    // triangle does not contain any stratNode
    if (thisBox.containsNone())
      continue;
    // Set StratConnect for the stratNode within the current triangle
    {
      // clip bounds within the actual bounds of the StratGrid
      const int bimin = max(0,thisBox.imin);
      const int bimax = min(getImax()-1,thisBox.imax);
      const int bjmin = max(0,thisBox.jmin);
      const int bjmax = min(getJmax()-1,thisBox.jmax);
      // Find which StratNode is contained within the current triangle.
      for(int i=bimin; i<=bimax; ++i) {
	for(int j=bjmin; j<=bjmax; ++j) {
	  if ((*StratConnect)(i,j) != NULL)
	    continue;
	  if (ct->containsPoint( (*StratNodeMatrix)(i,j).getX(),
				 (*StratNodeMatrix)(i,j).getY() )
	      ){
	    (*StratConnect)(i,j) = ct;
	  }
	}
      }
    }
  }
}

/***********************************************************************\
 ** tStratGrid::UpdateStratGrid(int mode, double time)
 **
 **
 ** 27-10-2003 (QC)
\***********************************************************************/
void tStratGrid::UpdateStratGrid(tUpdate_t mode, double time)
{
  // 0 = Initialisation at the beginning of a time step
  switch(mode){
  case k0:
    {
      std::cout<<"+++ 0 -UPDATESTRATGRID, START INITIALIZING...."<<std::endl;
      updateConnect();
      ResetAccummulatedDh();
      InterpolateErodepFromElevations(time);          // does the same as  the sweep function
      CheckSectionBase(mode);
      std::cout<<"   "<<std::endl;
      std::cout<<"+++ 0 -UPDATESTRATGRID, FINISHED INITIALIZING...."<<std::endl;
      std::cout<<"   "<<std::endl;
    }
    break;
    // 1 = Update after streampower-type erosion and deposition
  case k1:
    {
      std::cout<<"+++ 1 -UPDATESTRATGRID, START STREAMPOWER...."<<std::endl;
      updateConnect();
      InterpolateErodep(time);
      ResetAccummulatedDh();
      CheckSectionBase(mode);
      std::cout<<"   "<<std::endl;
      std::cout<<"+++ 1 -UPDATESTRATGRID, FINISHED STREAMPOWER...."<<std::endl;
      std::cout<<"   "<<std::endl;
    }
    break;
    // 2 = Update after meander migration (may decapitate stratNodes)
  case k2:
    {
      std::cout<<"+++ 2 -UPDATESTRATGRID, START MIGRATION...."<<std::endl;
      updateConnect();
      SweepChannelThroughRectGrid(time);  	// does NOT use accumDh
      CheckSectionBase(mode);
      ResetAccummulatedDh();
      std::cout<<"   "<<std::endl;
      std::cout<<"+++ 2 -UPDATESTRATGRID, FINISHED  MIGRATION...."<<std::endl;
      std::cout<<"   "<<std::endl;
    }
    break;
    // 3 = Update after geometrical meander-related erosion and deposition
  case k3:
    {
      std::cout<<"+++ 3-UPDATESTRATGRID, START CHANNEL DRIVER...."<<std::endl;
      updateConnect();
      InterpolateErodep(time);
      CheckSectionBase(mode);
      ResetAccummulatedDh();
      std::cout<<"    "<<std::endl;
      std::cout<<"+++ 3-UPDATESTRATGRID, FINISHED CHANNEL DRIVER...."<<std::endl;
      std::cout<<"   "<<std::endl;
    }
    break;
  case k4:
    {
      std::cout<<"+++ 4 -UPDATESTRATGRID, START  FLOODPLAIN...."<<std::endl;
      updateConnect();
      InterpolateErodep(time);
      CheckSectionBase(mode);
      ResetAccummulatedDh();
      std::cout<<"    "<<std::endl;
      std::cout<<"+++ 4 -UPDATESTRATGRID, FINISHED  FLOODPLAIN...."<<std::endl;
      std::cout<<"   "<<std::endl;
    }
    break;
  default:
    ReportFatalError("tStratGrid::UpdateStratGrid(): unknown mode.");
  }
}

/************************************************************************\
 **	ResetAccummulatedDh
 **
 **     Loop through the entire tMesh and reset all dh values to 0.0
 \*************************************************************************/
void tStratGrid::ResetAccummulatedDh()
{
  tMesh<tLNode>::nodeListIter_t ni( mp->getNodeList() );         // iterator for tMesh nodes
  tLNode *cn;

  for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
    {
      cn->ResetAccummulatedDh();
    }
}

/**************************************************************************\
 **
 ** tStratGrid::InterpolateElevations( )
 **
 ** Loop over the stratGrid and find for every node the corresponding
 ** elevation by interpolation between the neighbouring nodes of the
 ** tMesh Triange it is associated to.
 **
 ** 28-10-2003 (QC)
\**************************************************************************/
void tStratGrid::InterpolateElevations()
{
  int i,j;

  for(i=0; i<imax; i++){
    for(j=0; j<jmax; j++){

      tTriangle *ct =   (*StratConnect)(i,j);				// fetch  Triangle
      const double sx = (*StratNodeMatrix)(i,j).getX();                // i,j's  X-value
      const double sy = (*StratNodeMatrix)(i,j).getY();	        // i,j's  Y-vlaue

      if(ct != NULL){				                    // If Triangle present
       	int p;
       	tLNode *lnds[3];					// put the nodes in an array
       	int numnodes=0;
       	for(p=0; p<=2; p++){
	  lnds[numnodes] = static_cast<tLNode *>(ct->pPtr(p));
	  numnodes++;
        }
	// and interpolate a new z value

    	const tArray<double> tri_elevs(lnds[0]->getZ(),lnds[1]->getZ(), lnds[2]->getZ() );
        const double newz = PlaneFit(sx, sy, lnds[0]->get2DCoords(), lnds[1]->get2DCoords(),
				     lnds[2]->get2DCoords(), tri_elevs );

	// Set the elevation to the interpolated value, e.g force the StratNode columns up or down
	(*StratNodeMatrix)(i,j).setZ( newz );

	if (0) //DEBUG
	  std::cout<<"i= "<<i<<" j= "<<j<<" elev= "
	      <<newz<<" ,surr 3 nodes have Z "<<lnds[0]->getZ()
	      <<" "<<lnds[1]->getZ()<<" "<<lnds[2]->getZ()<<std::endl;

      } // node is connected to a triangle

      //if(ct==NULL){     //DEBUG:
      //std::cout <<"For i,j= "<<i<<" "<<j<<" No triangle - DEBUG"<<std::endl;
      //}

    } //j
  } //i

} // end of InterpolateElevations

/*****************************************************************\
 **  DEBUG FUNCTION    CheckSectionBase
 **
 **
\*****************************************************************/
void tStratGrid::CheckSectionBase(int mode)
{
  int i;
  int j;

  for(i=1; i<imax; i++){
    for(j=1; j<jmax; j++){

      tStratNode &sn   = (*StratNodeMatrix)(i,j);        // fetch the corresponding stratmode
      double currentZ  = sn.FindCurrentSectionBase();    // the current one
      double startZ    = sn.getSectionBase();            // the original one

      if(currentZ < startZ && startZ-currentZ > 0.001){   // More than 1 cm difference detected ?
	std::cout<<"        "<<std::endl;
	if(mode==0)     { std::cout<<"After Initialisation.. "<<std::endl;}
	else if(mode==1){ std::cout<<"After Streampower erodep.. "<<std::endl;}
	else if(mode==2){ std::cout<<"After Migrate.."<<std::endl;}
	else if(mode==3){ std::cout<<"After ChannelDriver.."<<std::endl;}
	else if(mode==4){ std::cout<<"After Floodplain wings"<<std::endl;}

	else { std::cout<<" UpdateStratGrid Mode not recognized in CheckSectionBase..."<<std::endl; }
	std::cout<<"At i,j "<<i<<", "<<j<<" x,y,z "<<sn.getX()<<" "<<sn.getY()<<" "<<sn.getZ()<<std::endl;
	std::cout<<"WARNING: Stratigraphic column pushed down"<<std::endl;
	std::cout<<"Original Base Z = "<<startZ<<", New Base Z = "<<currentZ<<", DZ= "<<currentZ-startZ<<std::endl;
	std::cout<<"Printing layerlist: "<<std::endl;
	int numlayers = sn.getNumLayer();
	int l=1;
	double totalthickness = 0.0;
	while(l<numlayers){
	  const double thickness = sn.getLayerDepth(l);
	  if(thickness < 100){
	    totalthickness += thickness;
	    //std::cout<<"Layer "<<l<< " has thickness "<<thickness<<std::endl;
	  }
	  std::cout<<"Layer "<<l<< " has thickness "<<thickness<<std::endl;
	  l++;
	}
	std::cout<<"Number of layers is "<<numlayers<<", Total thickness is: "<<totalthickness<<std::endl;
	if(totalthickness > 0.0){
	  exit(1);
	}
      }


    }
  }
}

/*****************************************************************\
 ** tStratGrid::InterpolateErodep(double time)
 **
 **
\******************************************************************/
void tStratGrid::InterpolateErodep(double time)
{
  int i,j;

  for(i=0; i<imax; i++){
    for(j=0; j<jmax; j++){

      // Step 1: get the stratNode location and fill the Tringle point array
      tTriangle *ct =   (*StratConnect)(i,j);				// fetch  Triangle
      double sx = (*StratNodeMatrix)(i,j).getX();                  // i,j's  X-value
      double sy = (*StratNodeMatrix)(i,j).getY();	                // i,j's  Y-vlaue

      if(ct != NULL){				                        // If Triangle present
	tLNode *lnds[3];
	int p;					                // put the Triangle nodes in an array
	int numnodes=0;
	for(p=0; p<=2; p++){					        // Euh,...what about boundary nodes ?
	  lnds[numnodes] = static_cast<tLNode *>(ct->pPtr(p));
	  numnodes++;
	}

	// Step 2: Interpolate the values of the incremented Coarse-grained material at the StratNode's location
	tArray<double>Dh(2);
	const tArray<double> tri_CoarseDh(lnds[0]->getAccCoarse(),
					  lnds[1]->getAccCoarse(),
					  lnds[2]->getAccCoarse() );
	Dh[0] = PlaneFit(sx, sy, lnds[0]->get2DCoords(),
			 lnds[1]->get2DCoords(),
			 lnds[2]->get2DCoords(), tri_CoarseDh );

	// Step 3: Interpolate the values of the incremented Fine-grained material at the StratNode's location
	const tArray<double> tri_FineDh(lnds[0]->getAccFine(),
					lnds[1]->getAccFine(),
					lnds[2]->getAccFine() );
	Dh[1] = PlaneFit(sx, sy, lnds[0]->get2DCoords(),
			 lnds[1]->get2DCoords(),
			 lnds[2]->get2DCoords(), tri_FineDh );

	// Step 4: Pass the interpolated values for dh[0](gravel) and dh[1](sand,fine) to a tStratNode
	//         specific erosion and deposition function, ErodepSimple
	//         This deals with z and the stratigraphy

	double dhtotal = Dh[0] + Dh[1];

	//DEBUG
	//if(dhtotal < 0.0){
	//   std::cout<<" ERODING "<<i<<" "<<j<<" ,with Dh[0,Dh[1] "<<Dh[0]<< " "<<Dh[1]<<", and total "<<dhtotal<<std::endl;
	//}

	if(  (Dh[0] > 0.0 && Dh[1] < 0.0) || (Dh[0] < 0.0 && Dh[1] > 0.0) ){
	  //std::cout<<"    "<<std::endl;
	  //std::cout<<"Doing both erosion and deposition "<<std::endl;
	  //std::cout<<"For node "<<i<<" "<<j<<" at "<<sx<<" "<<sy<<" "<<(*StratNodeMatrix)(i,j).getZ()<<std::endl;
	  //std::cout<<"Dh[0] = "<<Dh[0]<<"  "<<" Dh[1]= "<<Dh[1]<<std::endl;
	  //std::cout<<"Contributing tMesh nodes are: "<<std::endl;
	  //int c;
	  //for(c=0; c<=2; c++){
	  // std::cout<<" x= "<<lnds[c]->getX()<<" y= "<<lnds[c]->getY()<<" z= "<<lnds[c]->getZ()<<" Coarse= "<<lnds[c]->getAccCoarse()<<" Fine= "<<lnds[c]->getAccFine()<<std::endl;
	  //}
	  //std::cout<<"    "<<std::endl;

	  if( dhtotal > 0.0){
	    if( Dh[0] > Dh[1]) {	      //Gravel is more important, make all deposition gravel
              Dh[0] = dhtotal;
              Dh[1] = 0.0;
	    }
	    else if( Dh[0]<=Dh[1]){     //sand is more important, make all deposition sand
              Dh[0] = 0.0;
              Dh[1] = dhtotal;
	    }
	  }
	  else if( dhtotal < 0.0){
	    Dh[0]=0.0;
	    Dh[1]=dhtotal;
	  }

	}

	double current = CalculateMeanderCurrent(ct,sx,sy);
	tStratNode &sn = (*StratNodeMatrix)(i,j);
	sn.EroDepSimple( 0, Dh, time,current );
      }
    } //j
  } // i

} //end of InterpolateErodep

//--------------------old stuff below


void tStratGrid::setSurface( int t, double val)
{
  surface[t]=val;
}
void tStratGrid::setSubsurface( int t, double val)
{
  subsurface[t]=val;
}

void tStratGrid::setSubsurface_mbelt(int t, double val)
{
  subsurface_mbelt[t]=val;
}

double tStratGrid::getSurface( int t ) const
{
  return surface[t];
}

double tStratGrid::getSubsurface (int t ) const
{
  return subsurface[t];
}

double tStratGrid::getSubsurface_mbelt(int t) const
{
  return subsurface_mbelt[t];
}

int tStratGrid::getnWrite() const
{
  return nWrite;
}

double tStratGrid::getOutputTime(int t) const
{
  return outputTime[t];
}

void tStratGrid::setOutputTime( int t, double val)
{
  outputTime[t]=val;
}

int tStratGrid::getSectionLocation(int i) const
{
  return section[i];	  // nb first 5 section are in x, the next 5 are in y
}
/*************************************************************************\
 void tStratGrid::SweepChannelThroughStratGrid

 In the function Movenodes the tLNodes are relocated due to the meandering
 but the elevation of the (new) channel node is identical to that of its old
 position. So as a result of meandering the elevations along the main channel
 will be changed. In order to find the amount of erosion (and depostion)
 needed to transform the elevation change to the startgrid, the elevation  This routine is called after the Meander->Migrate routine to evo

 Created 06/2003 (QC)

\*************************************************************************/

void tStratGrid::SweepChannelThroughRectGrid(double time)
{
  int i,j;
  
  if(1) std::cout<<"SCTRG";

  for(i=0; i<imax; i++){
    for(j=0; j<jmax; j++){

      // make the triangle nodes array

      tTriangle *ct =   (*StratConnect)(i,j);				// fetch  Triangle
      const double sx = (*StratNodeMatrix)(i,j).getX();                  // i,j's  X-value
      const double sy = (*StratNodeMatrix)(i,j).getY();	                // i,j's  Y-vlaue

      if(ct != NULL){				                        // If Triangle present
	tLNode *lnds[3];
	int p;					                // put the Triangle nodes in an array
	int numnodes=0;
	for(p=0; p<=2; p++){					        // Euh,...what about boundary nodes ?
	  lnds[numnodes] = static_cast<tLNode *>(ct->pPtr(p));
	  numnodes++;
	}

	// find the new elevation, after migrate
	// and interpolate a new z value

	const tArray<double> tri_elevs(lnds[0]->getZ(),lnds[1]->getZ(), lnds[2]->getZ() );
	const double newz = PlaneFit(sx, sy, lnds[0]->get2DCoords(), lnds[1]->get2DCoords(),
				     lnds[2]->get2DCoords(), tri_elevs );

	// compare the new, interpolated elevation with the existing elevation for the StratNode
	// What is the elevation difference ?
	tStratNode &sn = (*StratNodeMatrix)(i,j);
	const double oldz = sn.getZ();
	const double dh   = newz-oldz;

	// ---------DEPOSITION----------------------------
	if(dh > 0.0){                       // Depositing, = du to channel, so drop dh as coarse material only
	  const tArray< double > erode(dh,0.);
	  const double current = CalculateMeanderCurrent(ct,sx,sy);
	  sn.EroDepSimple(0,erode,time,current);
	  //std::cout<< "We changed layerlist through depostion \n";
	}
	else if( dh == 0.0){

	}
	// --------EROSION--------------------------------
	else if( dh < 0.0){                 // Erode, most likely situation as channel sweeps away, decapitates banknodes
	  const tArray< double > erode(0.5*dh, 0.5*dh);
	  const double current = CalculateMeanderCurrent(ct,sx,sy);
	  sn.EroDepSimple(0,erode,time,current);
	}

      } // Triangle is not NULL

    } // i
  }   // j

  if(1) std::cout<<".End" << std::endl;

} // end of function tStratGrid SweepChannelThroughRectGrid(double time)



void tStratGrid::InterpolateErodepFromElevations(double time)
{
  int i,j;

  for(i=0; i<imax; i++){
    for(j=0; j<jmax; j++){

      // make the triangle nodes array

      tTriangle *ct =   (*StratConnect)(i,j);        // fetch  Triangle
      tStratNode const &sn_ = (*StratNodeMatrix)(i,j);
      const double sx = sn_.getX();                  // i,j's  X-value
      const double sy = sn_.getY();                  // i,j's  Y-vlaue

      if(ct != NULL){				                        // If Triangle present
	tLNode *lnds[3];
	int p;					                // put the Triangle nodes in an array
	int numnodes=0;
	for(p=0; p<=2; p++){					// Euh,...what about boundary nodes ?
	  lnds[numnodes] = static_cast<tLNode *>(ct->pPtr(p));
	  numnodes++;
	}

	// find the new elevation, after migrate
	// and interpolate a new z value

	const tArray<double> tri_elevs(lnds[0]->getZ(),lnds[1]->getZ(), lnds[2]->getZ() );
	const double newz = PlaneFit(sx, sy, lnds[0]->get2DCoords(), lnds[1]->get2DCoords(),
				     lnds[2]->get2DCoords(), tri_elevs );

	// compare the new, interpolated elevation with the existing elevation for the StratNode
	// What is the elevation difference ?
	tStratNode &sn = (*StratNodeMatrix)(i,j);
	const double oldz = sn.getZ();
	const double dh   = newz-oldz;

	// ---------DEPOSITION----------------------------
	if(dh > 0.0){                       // Depositing, = du to channel, so drop dh as coarse material only
	  const tArray< double > erode(dh, 0.);
	  const double current = CalculateMeanderCurrent(ct,sx,sy);
	  sn.EroDepSimple(0,erode,time,current);
	  //std::cout<< "We changed layerlist through depostion \n";
	}
	else if( dh == 0.0){

	}
	// --------EROSION--------------------------------
	else if( dh < 0.0){                 // Erode, most likely situation as channel sweeps away, decapitates banknodes
	  const tArray< double > erode(0.5*dh, 0.5*dh);
	  const double current = CalculateMeanderCurrent(ct,sx,sy);
	  sn.EroDepSimple(0,erode,time,current);
	}

      } // Triangle is not NULL

    } // i
  }   // j


} // end of function

/*****************************************************************************\
 ** double CalculateMeanderCurrent()
 **
 ** Finds the averaged flow direction for a stratNode, based on interpolation
 ** of surrounding meander nodes
 **
 **
\******************************************************************************/

double tStratGrid::CalculateMeanderCurrent( tTriangle *tri, double sx, double sy) const
{

  double CompassAv = -1;

  if(tri == NULL)  // If Triangle not present
    return CompassAv;

  tLNode *lnds[3];
  int numnodes=0;
  for(size_t p=0; p<=2; ++p){
    if(tri->pPtr(p)->isMobile()){		// only fill the array with meandering nodes
      lnds[numnodes] = static_cast<tLNode *>(tri->pPtr(p));
      numnodes++;
    }
  }

  switch(numnodes){
  case 3:     // There are three surrounding meander nodes
    {
      const tArray<double>
	CompassAngles(
		      CompassAngle(lnds[0],lnds[0]->getDownstrmNbr()),
		      CompassAngle(lnds[1],lnds[1]->getDownstrmNbr()),
		      CompassAngle(lnds[2],lnds[2]->getDownstrmNbr())
		      );
      CompassAv = PlaneFit(sx, sy, lnds[0]->get2DCoords(), lnds[1]->get2DCoords(),
			   lnds[2]->get2DCoords(), CompassAngles);
    }
    break;
  case 2:    // There are two surrounding meander nodes
    {
      double dist[2];    //distance b/w new point and the three triangle points
      dist[0] = DistanceBW2Points(sx, sy, lnds[0]->getX(), lnds[0]->getY());
      dist[1] = DistanceBW2Points(sx, sy, lnds[1]->getX(), lnds[1]->getY());

      double CompassAngles[2];
      CompassAngles[0] = CompassAngle(lnds[0],lnds[0]->getDownstrmNbr() );
      CompassAngles[1] = CompassAngle(lnds[1],lnds[1]->getDownstrmNbr() );

      CompassAv = LineFit(0., CompassAngles[0], dist[0]+dist[1],
			  CompassAngles[1], dist[0]);

    }
    break;
  case 1:    // There is one surrounding meander node
    {
      CompassAv = CompassAngle(lnds[0],lnds[0]->getDownstrmNbr() );
    }
    break;
  default:    // There is no surrounding meander node
    CompassAv = -1;
    break;
  }

  return CompassAv;
}

/*********************************************************\
 **  Returns compass angle
 **  Takes a tLNode and its downstream neighbour and
 **  returns a compass angle. Assumes Y+ = North and X+ = East
 **
\*********************************************************/
double tStratGrid::CompassAngle(tLNode *cn, tLNode *dn) const
{
  double cnx = cn->getX();              // node itself
  double cny = cn->getY();

  double dnx = dn->getX();              // and its downstream neighbour
  double dny = dn->getY();

  double delx = dnx - cnx;
  double dely = dny - cny;
  double compass = -1;

  double RtoD = (180/PI);
  if (delx == 0.0 && dely == 0.0){
    compass = -1;
  }
  else{
    compass = 270 - (  atan2(dely,delx) *RtoD );
  }

  return compass;
}




//------------------END OF CLASS STRATGRID---------------------------------










//--------------------CLASS STRATNODE-----------------------------------------

/**************************************************************************\
 **
 **  fn tStratNode( tInputFile &infile, tMesh<tLNode *mp );
 **  @brief Main constructor for tStratNode
 **
 **  Modifications:
 **    - changed input parameter MAXREGDEPTH to SG_MAXREGDEPTH to 
 **      distinguish from the Voronoi-based layering procedure. GT 1/05
 **
\**************************************************************************/

// initialize static members numg and grade
int tStratNode::numg = 0;
tArray<double> tStratNode::grade = 1;
double tStratNode::maxregdep = 1.;
double tStratNode::KRnew = 1.0;



// 1) tStratNode default constructor
tStratNode::tStratNode() :
  layerlist(),ClosestNode(0),
  x(0.),y(0.), z(0.),sectionZ(0.),newz(0.),i(0),j(0)
{}

tStratNode::tStratNode( int ) :
  layerlist(),ClosestNode(0),
  x(0.),y(0.), z(0.),sectionZ(0.),newz(0.),i(0),j(0)
{}

// 2) tStratNode constructor for indices & coordinates
tStratNode::tStratNode( double x_, double y_):
  layerlist(),ClosestNode(0),
  x(x_), y(y_), z(0.),sectionZ(0.),newz(0.),i(0),j(0)
{}

// 3) tStratNode constructor for indices, coordinates and assigning the initial
//    layerlist
tStratNode::tStratNode(double x_, double y_,
		       const tStratNode &orig):
  layerlist(),ClosestNode(0),
  x(x_), y(y_), z(0.),sectionZ(0.),newz(0.),i(0),j(0)
{
  layerlist = orig.layerlist;
}


// 4) tStratNode constructor for layerlist initialisation by inputfile
tStratNode::tStratNode( tInputFile const &infile ) :
  layerlist(),ClosestNode(0),
  x(0.),y(0.), z(0.),sectionZ(0.),newz(0.),i(0),j(0)
{
  int i;
  char add[2], name[20];
  double help, extra, sum, sumbr;
  tLayer layhelp, niclay;
  tArray<double> dgradehelp;
  tArray<double> dgradebrhelp;

  //std::cout << "=>STRATNODE( infile )" << std::endl;
  numg = infile.ReadItem( numg, "NUMGRNSIZE" );

  // This is a --HACK-- by QUINTIJN !!! Do not use this value
  // while trying to sort grainsizes properly on the tMesh, it only applies to tStratGRid
  // with is more geomettrical deposition functions in the tFloodplain class
  numg++;
  //-Carefull, this is done to give a StratNode 2 grainsizes classes,
  //while the Mesh is still running with 1 grain size as demanded
  //by detachment limited conditions during meandering

  grade.setSize( numg );
  maxregdep = infile.ReadItem( maxregdep, "SG_MAXREGDEPTH" );
  KRnew = infile.ReadItem( KRnew, "KR" );
  if( KRnew<0.0 )
    ReportFatalError( "Erodibility factor KR must be positive." );

  i=0;
  add[0]='1';
  add[1]='\0';

  while ( i<numg ){
    // Reading in grain size diameter info
    strcpy( name, "GRAINDIAM");
    strcat( name, add );
    help = infile.ReadItem( help, name);
    grade[i] = help;
    i++;
    add[0]++;
  }



  // Initializing the layering below:

  bool lay = infile.ReadBool( "OPTREADLAYER" );

  if(!lay){

    dgradehelp.setSize( numg );
    dgradebrhelp.setSize( numg );
    sum = 0;
    sumbr = 0;
    i=0;
    add[0]='1';

    while ( i<numg ){
      // Reading in proportions for initial regolith and bedrock
      strcpy( name, "REGPROPORTION");
      strcat( name, add );
      help = infile.ReadItem( help, name);
      dgradehelp[i]=help;
      sum += help;
      strcpy( name, "BRPROPORTION");
      strcat( name, add );
 	  std::cout<<"In tStratGrid, reading '"<<name<<std::endl;
     help = infile.ReadItem( help, name);
      dgradebrhelp[i]=help;
      sumbr += help;
      i++;
      add[0]++;
    }

    if(fabs(sum-1.0)>0.001)
      ReportFatalError("Problem with the regolith proportion of grain sizes in input file");
    if(fabs(sumbr-1.0)>0.001)
      ReportFatalError("Problem with the bedrock proportion of grain sizes in input file");

    help = infile.ReadItem( help, "REGINIT" );
    layhelp.setCtime( 0.0 );
    layhelp.setEtime( 0.0 );
    layhelp.setRtime( 0.0 );
    layhelp.setPaleoCurrent( -1 );
    if( help > 0){
      // Make a bedrock and regolith layer, possibly two depending
      // on the depth of the regolith layer.  The code will decide
      // the total number of layers needed.  By default the regolith
      // layer(s) is/are put on top of the bedrock.
      // Bedrock layer items read in and set
      help = infile.ReadItem( help, "KB");
      if( help<0.0 )
	ReportFatalError( "Erodibility factor KB must be positive." );
      layhelp.setErody(help);
      layhelp.setSed(tLayer::kBedRock);

      // in the regolith layer dgrade is saving
      // the proportion of grain size available from the bedrock

      layhelp.setDgradesize(numg);
      i=0;
      help = infile.ReadItem( help, "BEDROCKDEPTH");
      while(i<numg){
	layhelp.setDgrade(i, help*dgradebrhelp[i]);
	i++;
      }

      layerlist.insertAtBack( layhelp );

      // Regolith layer items read in and set
      layhelp.setSed(tLayer::kSed);
      help = infile.ReadItem( help, "KR" );
      if( help<0.0 )
	ReportFatalError( "Erodibility factor KR must be positive." );
      layhelp.setErody(help);
      help = infile.ReadItem( help, "REGINIT");
      extra = 0;
      if(help > maxregdep){
	// too much regolith, create two layers the bottom layer is made here
	extra = help - maxregdep;
	//layhelp.setDepth(extra);
	//layhelp.setDgradesize(numg);
	i=0;
	while(i<numg){
	  layhelp.setDgrade(i, extra*dgradehelp[i]);
	  i++;
	}
	if(extra>10000) //TODO, bettter criteria for setting deep sed. flag
	  layhelp.setRtime(-1.);
	layerlist.insertAtFront( layhelp );
	// the top regolith layer is now made
	layhelp.setRtime(0.);
	i=0;
	while(i<numg){
	  layhelp.setDgrade(i, maxregdep*dgradehelp[i]);
	  i++;
	}
	layerlist.insertAtFront( layhelp );
      }
      else{
	// create only one regolith layer
	i=0;
	while(i<numg){
	  layhelp.setDgrade(i, help*dgradehelp[i]);
	  i++;
	}

	layerlist.insertAtFront( layhelp );
      }
    }

    else{
      // no regolith, so by default everything is bedrock
      // Bedrock layer items set
      help = infile.ReadItem( help, "KB");
      if( help<0.0 )
	ReportFatalError( "Erodibility factor KB must be positive." );
      layhelp.setErody(help);
      layhelp.setSed(tLayer::kBedRock);
      layhelp.setDgradesize(numg);
      i=0;
      help = infile.ReadItem( help, "BEDROCKDEPTH");
      while(i<numg){
	layhelp.setDgrade(i, help*dgradebrhelp[i]);
	i++;
      }
      layerlist.insertAtBack( layhelp );
    }

    //std::cout << layerlist.getSize() << "  A Layerlistscreated in STRATNODE constructor " << std::endl;
  }

}

// 5) So, when do I use this one ?
tStratNode::tStratNode( const tStratNode &orig )
  :
  layerlist(),
  ClosestNode(orig.ClosestNode),
  x(orig.x),
  y(orig.y),
  z(orig.z),
  sectionZ(orig.sectionZ),
  newz(orig.newz),
  i(orig.i),
  j(orig.j)
{
  layerlist = orig.layerlist;
}

tStratNode &tStratNode::operator=( const tStratNode &right )
{
  if( &right != this )
    {
      layerlist = right.layerlist;
      ClosestNode = right.ClosestNode;
      x = right.x;
      y = right.y;
      z = right.z;
      sectionZ = right.sectionZ;
      newz = right.newz;
      i = right.i;
      j = right.j;
    }
  return *this;
}

tStratNode::~tStratNode()
{
  if (0) //DEBUG
    std::cout << "    ~STRATNODE()" << std::endl;
}



//"set" and "get" functions for the coordinates:
void tStratNode::setX( double val ) {x = val;}
void tStratNode::setY( double val ) {y = val;}
void tStratNode::setZ( double val ) {z = val;}
void tStratNode::setNewZ( double val ) { newz = val;}
void tStratNode::setSectionBase( double val ) { sectionZ = val; }

void tStratNode::setI( int val ) {i = val;}
void tStratNode::setJ( int val ) {j = val;}

/**************************************************************************\
 **
 **  TellAll
 **
 **  ...is a debugging routine that prints out a bunch of info about a stratNode
 **
\**************************************************************************/

tLNode * tStratNode::getClosestNode()
{
  return ClosestNode;
}


void tStratNode::InsertLayerBack( tLayer const & lyr )
{
  layerlist.insertAtBack( lyr );
}

double tStratNode::FindCurrentSectionBase()
{
  return z - AlluvialColumnThickness();
}

double tStratNode::AlluvialColumnThickness()
{
  double totalthickness = 0.;
  tListIter<tLayer> ly ( layerlist );
  ly.First();
  // skip first layer (hence initialisation of 'layer' with NextP())
  for (tLayer *layer = ly.NextP(); !( ly.AtEnd() ); layer = ly.NextP() ) {
    const double thickness = layer->getDepth();
    if(thickness < 100){
      totalthickness += thickness;
    }
  }
  // the depth of the base is the current elevation - thickness

  return totalthickness;
}
/**********************************************\
 **  getAgeAtDepth
 **
 **  returns the average age of the sediment
 **  at a given depth below the current surface
\**********************************************/

double tStratNode::getAgeAtDepth( double val) const
{
  const int sz = layerlist.getSize();
  int l=1;
  double thickness = 0;
  double totalthickness =0.;
  double age = 0.0;
  double Rsed = 0.0;
  while(l<sz && totalthickness < val){
    thickness = getLayerDepth(l);
    if(thickness < 100){
      totalthickness += thickness;
    }
    l++;
  }

  // Out of the loop
  // 1) in bedrock
  // 2) out of alluvial layers/ at alluvial layer boundary
  // 3) within an alluvial layer

  if(l==sz && totalthickness == 0.0){
    age = 0.0;
  }
  else if (totalthickness == val ){
    age = getLayerCtime(l-1);
  }
  else if (totalthickness > val && thickness < 100){
    if(getLayerRtime(l-1) > getLayerCtime(l-1)){
      // DEBUG
      //std::cout<<"Depth= "<<getLayerDepth(l-1)<<std::endl;
      //std::cout<<"Ctime= "<<getLayerCtime(l-1)<<std::endl;
      //std::cout<<"Rtime= "<<getLayerRtime(l-1)<<std::endl;

      Rsed = getLayerDepth(l-1)/( getLayerRtime(l-1) - getLayerCtime(l-1) );
      //std::cout<<"Rsed= "<<Rsed<<std::endl;

      age = getLayerCtime(l-1) + (totalthickness - val)*Rsed;

      //std::cout<<"age= "<<age<<std::endl;

      // DEBUG
      if(age < 0.0 || age > 1000000){
	std::cout<<"in getAgeAtDepth"<<std::endl;
	std::cout<<"l-1= "<<l-1<<" thickness= "<<thickness
	    <<" totalthick = "<<totalthickness<<std::endl;
	std::cout<<"val = "<<val<<std::endl;
	std::cout<<"l-1= "<<l-1<<" Rtime = "<<getLayerRtime(l-1)
	    <<"Ctime= "<<getLayerCtime(l-1)<<" "<<std::endl;
	std::cout<<"Rsed = "<<Rsed<<std::endl;
	exit(1);
      }

    }
    else{
      age =0.0;
    }
  }

  return age;
}



void tStratNode::setGrade( int i, double size ) const
{
  if(i>=numg)
    ReportFatalError("Trying to set a grain size for an index which is too large");
  grade[i] = size;
}

double tStratNode::getMaxregdep() const
{
  return maxregdep;
}

int tStratNode::getNumg() const
{
  return numg;
}

void tStratNode::setNumg( int size ) const
{
  numg = size;
}

double tStratNode::getGrade( int i ) const
{
  return grade[i];
}

const tArray< double >& tStratNode::getGrade( ) const
{
  return grade;
}

double tStratNode::getLayerCtime( int l ) const
{
  return layerlist.getIthDataRef(l).getCtime();
}
double tStratNode::getPaleoCurrent(int l ) const
{
  return layerlist.getIthDataRef(l).getPaleoCurrent();
}
void tStratNode::setLayerCtime( int l, double tt)
{
  tListIter<tLayer> ly ( layerlist );
  tLayer  * hlp;
  hlp=ly.FirstP();

  int n=0;

  while(n<l){
    n++;
    hlp=ly.NextP();
  }

  hlp->setCtime( tt );
}

double tStratNode::getLayerRtime( int l ) const
{
  return layerlist.getIthDataRef(l).getRtime();
}

void tStratNode::setLayerRtime( int l, double tt)
{
  tListIter<tLayer> ly ( layerlist );
  tLayer  * hlp;
  hlp=ly.FirstP();

  int n=0;

  while(n<l){
    n++;
    hlp=ly.NextP();
  }

  hlp->setRtime( tt );
}

double tStratNode::getLayerEtime( int l ) const
{
  return layerlist.getIthDataRef(l).getEtime();
}

void tStratNode::setLayerEtime( int l, double tt)
{
  tListIter<tLayer> ly ( layerlist );
  tLayer  * hlp;
  hlp=ly.FirstP();

  int n=0;

  while(n<l){
    n++;
    hlp=ly.NextP();
  }

  hlp->setEtime( tt );
}

void tStratNode::addLayerEtime( int l, double tt)
{
  tListIter<tLayer> ly ( layerlist );
  tLayer  * hlp;
  hlp=ly.FirstP();

  int n=0;

  while(n<l){
    n++;
    hlp=ly.NextP();
  }

  hlp->addEtime( tt );
}

double tStratNode::getLayerDepth( int l ) const
{
  if( layerlist.isEmpty() )
    {
      std::cout << "** WARNING layer list is empty\n";
      std::cout << "  x=" << x << " y=" << y << " z=" << z;

    }
  return layerlist.getIthDataRef(l).getDepth();
}

void tStratNode::setLayerDepth( int l, double dep)
{
  assert( dep > 0.0 );
  tListIter<tLayer> ly ( layerlist );
  tLayer *hlp;
  hlp=ly.FirstP();

  int n=0;

  while(n<l){
    n++;
    hlp=ly.NextP();
  }

  hlp->setDepth( dep );
}

double tStratNode::getLayerErody( int l ) const
{
  return layerlist.getIthDataRef(l).getErody();
}

void tStratNode::setLayerErody( int l, double ero)
{
  tListIter<tLayer> ly ( layerlist );
  tLayer  * hlp;
  hlp=ly.FirstP();

  int n=0;

  while(n<l){
    n++;
    hlp=ly.NextP();
  }

  hlp->setErody( ero );
}


tLayer::tSed_t tStratNode::getLayerSed( int l ) const
{
  return layerlist.getIthDataRef(l).getSed();
}

void tStratNode::setLayerSed( int i, tLayer::tSed_t s)
{
  tListIter<tLayer> ly ( layerlist );
  tLayer  * hlp;
  hlp=ly.FirstP();

  int n=0;

  while(n<i){
    n++;
    hlp=ly.NextP();
  }

  hlp->setSed( s );
}

double tStratNode::getLayerDgrade( int l, int num ) const
{
  return layerlist.getIthDataRef(l).getDgrade(num);
}

void tStratNode::setLayerDgrade( int i, int g, double val)
{
  tListIter<tLayer> ly ( layerlist );
  tLayer  * hlp;
  hlp=ly.FirstP();


  // assert( val>=0.0 ); For Use in floodplain, quintijn

  int n=0;

  while(n<i){
    n++;
    hlp=ly.NextP();
  }

  hlp->setDgrade(g, val );
}
/*************************************************************\
  tStratNode::ErodeSimple

  Only for use on the rectangular StratGrid !
  
  Modifications:
   - Bug fix: "if(-remainder<thickness)" changed to <=, to 
   catch special case in which -remainder==thickness
   GT&QC 2/05
\**************************************************************/

void tStratNode::EroDepSimple( int l, tArray<double>dh, double tt,double current)
{
  int debugI,debugJ;
  int g;
  double dhtotal,remainder,thickness;
  double tofill;
  tArray<double> fill(numg);
  tArray<double> ratio(numg);
  tArray<double> newl(numg);

  g=0;
  dhtotal=0.0;
  remainder=0.0;
  thickness=0.0;
  tofill=0.0;

  debugI=19;
  debugJ=19;

  // SimpleErodep does not use an active top layer.
  // real layering starts with index 1, the rest is zero
  g=0;
  while(g<numg){
    setLayerDgrade(l,g,0.0);
    g++;
  }

  //------------Space for debug prints---------------------------------
  // CHECK the elevation balance before erosion and depostion, it should be that
  // base
  double sectionBase  = getSectionBase();
  double columnheight = AlluvialColumnThickness();
  double startdifference = (sectionBase)- (z-columnheight);
  if(z-columnheight != sectionBase && columnheight > 0.0){
    if(startdifference < -0.001 || startdifference > 0.001){
      std::cout<<"    "<<std::endl;
      std::cout<<"Entering ErodepSimple, no elevation balance"<<std::endl;
      std::cout<<"at "<<getI()<<" "<<getJ()<<std::endl;
      std::cout<<"z = "<<z<<" alluv= "<<columnheight<<" base= "<<sectionBase<<std::endl;
      std::cout<<"z-columnheight= "<<z-columnheight<<std::endl;
    }
  }

  if( getI()== debugI && getJ()== debugJ && 0 ){
    std::cout<<" "<<std::endl;
    std::cout<<"Start erodep for node "<<getI()<<" "<<getJ()<<std::endl;
    std::cout<<"dh[0]= "<<dh[0]<<" dh[1]= "<<dh[1]<<" dhtotaal= "<<dh[0] + dh[1]<<std::endl;
    std::cout<<"Start elev= "<<z<<" startbase= "<< getSectionBase()<<std::endl;
    std::cout<<"Start allu-thickness= "<<AlluvialColumnThickness()<<std::endl;
  }



  //-----------end of space for debug prints---------------------------

  l++;                                   // Go to the real stratigraphic layer (start with 1)
  dhtotal = dh[0] + dh[1];               // Total amount of erosion/deposition
  z+=dhtotal;                            // Modify elevation of 'this' node

  /************************************************************\
                        DEPOSITION
 \************************************************************/
  if( (dh[0] > 0.0 || dh[1] > 0.0) && dhtotal > 0.0){
    //1)We have an top sediment layer
    if(getLayerSed(l) != 0 && getLayerDepth(l) > 0.0){


      //std::cout << "Deposition at  "<< x <<' '<< y <<" dh[0]= " << dh[0] << " dh[1]= " << dh[1] << " dhtotal= " << dhtotal << " time= " <<tt<< '\n';

      if(getLayerDepth(l) >=maxregdep){        // The top layer is already too thick
        if(getLayerDepth(l) >=1000.){         // its the initialisation layer
	  makeNewLayerBelow(l-1, getLayerSed(l), getLayerErody(l), dh, tt,current);
        }
        else if(getLayerDepth(l) < 1000.){    // its a true sedimentary layer, just add a new one
          makeNewLayerBelow(l-1, getLayerSed(l), getLayerErody(l), dh, tt,current);
        }
      }
      else if(getLayerDepth(l) < maxregdep){                   //potential for adding material in top layer
     	if(getLayerDepth(l) + dhtotal <= maxregdep){          // enough space in top layer, put it here
	  g=0;
	  while(g<numg){
	    addtoLayer(l,g,dh[g],tt,current);                       // just put everything here
	    g++;
	  }
     	}
     	else if(getLayerDepth(l) + dhtotal > maxregdep){     // you would fill more than one layer
	  tofill = maxregdep-getLayerDepth(l);              // remaining space available in top layer
	  tArray<double>ratio(2);
	  ratio[0] = dh[0]/dhtotal;                         // fraction coarse in added dh
	  ratio[1] = dh[1]/dhtotal;                         // fraction fine

	  g=0;
	  while(g<numg){
     	    fill[g]= (ratio[g])*tofill;
     	    addtoLayer(l,g,fill[g],tt,current);                     // complete top sediment layer
     	    g++;
	  }

     	  // How much is left ?
     	  remainder=dhtotal-tofill;                // remaining stuff to make a new layer
	  if(remainder <= maxregdep+10.){           // we can only fill one layer extra
	    newl[0]=ratio[0]*remainder;
	    newl[1]=ratio[1]*remainder;
	    makeNewLayerBelow(l-1, getLayerSed(l), getLayerErody(l), newl, tt,current);
	    remainder =0.0;
	  }
	  else if(remainder <= maxregdep+10.){
     	    std::cout<<"Loads of deposition at " << x <<' '<< y <<" time= " <<tt<< '\n';
     	    std::cout<<"Make extra layers ??? \n";
     	    exit(1);
	  }


     	}
      }
    }

    //2)No top sediment layer, deposit everything  on top of bedrock
    else if(getLayerSed(i)==0){
      makeNewLayerBelow(-1,tLayer::kSed,KRnew,dh,tt,current);
    }

  }
  /************************************************************\
                         EROSION
  Just erode the dhtotal, do not consider complicated mass
  balance in the composition
 \************************************************************/
  else if( (dh[0] <= 0.0 && dh[1] <= 0.0) && dhtotal < 0.0){
    // Is there a sediment layer present ?
    //std::cout << "Erosion in "<< x <<' '<< y <<" dhtotal= " << dhtotal << " time= "<< tt << '\n';
    if( getLayerSed(l)!=0){

      if(getLayerDepth(l) >= (-dhtotal) ){            // There is enough to erode in the first layer


	if(getI()==debugI && getJ()==debugJ && 0)
	  {std::cout<<" erode 1"<<std::endl;
	  std::cout<<" thickness = "<<getLayerDepth(l)<<std::endl;
	  }
	// just erode from the top sediment layer
	addtoLayer(l,dhtotal,tt);
	if(getLayerDepth(l) > 1000){  setSectionBase(z);   }
	// Is the current op layer after this operation too small ?, if yes just remove it
	//if(getLayerDepth(l)<1e-7){
	//   removeLayer(l);
	//}
	if(getI()==debugI && getJ()==debugJ && 0){std::cout<<" erode 1-finished"<<std::endl;    }

      }
      else if(getLayerDepth(l) < (-dhtotal) ){            // we want to erode more than the top one layer

	if(getI()==debugI && getJ()==debugJ){std::cout<<" erode 2"<<std::endl;    }

	remainder=dhtotal; 				    // negative value, giving the dh we need to erode
	int debugCount=0;
	while(remainder < 0.0 && getLayerSed(l)!=0){
	  thickness=getLayerDepth(l);
	  if(-remainder <= thickness){                    // remainder is smaller thatn the new layer thickness
	    if(getI()==debugI && getJ()==debugJ && 0){
	      std::cout<<" Want to erode layer thickness= "<<thickness<< " in erode 3"<<std::endl;
	      std::cout<<" With "<<remainder<<std::endl;
	    }

	    addtoLayer(l,remainder,tt);                   // just erode the remainder
	    remainder=0.0;                              // and reset it to zero

	    if(getI()==debugI && getJ()==debugJ && 0){
	      std::cout<<" eroded 3, new thickness= "<<getLayerDepth(l)<<std::endl;
	    }
	    if(thickness > 1000){ setSectionBase(z);}
	  }
	  else if(-remainder > thickness){

	    removeLayer(l);                              // remainder is larger than the layer
	    remainder+=thickness;                        // remove the entire layer, but decrement the remainder
	    if(getI()==debugI && getJ()==debugJ && 0){
	      std::cout<<" erode 4 removing thickness "<< thickness<<std::endl;
	      std::cout<<" getlayerdepth now returns  "<<getLayerDepth(l)<<std::endl;
	      std::cout<<"and the remainder is "<<remainder<<std::endl;
	    }
	  }
      if(1) {
  	     debugCount++;
		 if( debugCount>1e6 ) {
		    std::cout<<"More than 1e6 iterations in tStratNode::EroDepSimple." << std::endl;
			std::cout<<"Node "<<getI()<<","<<getJ()<<" remainder "<<remainder<<" thickness "<<thickness<<std::endl;
			ReportFatalError("Apparent endless loop in tStratNode::EroDepSimple");
		 }
	  }
	}


      }



    } // sed layer present
    // We are on bedrock and still want to erode
    else if(getLayerSed(l) ==0){

      if(getI()==debugI && getJ()==debugJ){std::cout<<" erode 5"<<std::endl;    }
      addtoLayer(l,dhtotal,tt);	                // should also only erode bedrock
      setSectionBase(z);                        // adjust the base of the alluvial stratigraphic section


    }

  } // EROSION
  /************************************************************\
              SIMULTANEOUS DEPOSITION & EROSION
 \************************************************************/
  if( (dh[0] > 0.0 && dh[1] < 0.0) || (dh[0] < 0.0 && dh[1] > 0.0) ){

    std::cout<<"WARNING -- SIMULTANEOUS EROSION AN DEPOSITION "<<std::endl;


  } // **********Deposition & Erosion simultaneously************

  if( getI()== debugI && getJ()== debugJ && 0){
    std::cout<<" "<<std::endl;
    std::cout<<"finish erodep for node "<<getI()<<" "<<getJ()<<std::endl;
    std::cout<<"dh[0]= "<<dh[0]<<" dh[1]= "<<dh[1]<<" dhtotaal= "<<dh[0] + dh[1]<<std::endl;
    std::cout<<"finish elev= "<<z<<" finish-base= "<< getSectionBase()<<std::endl;
    std::cout<<"finish allu-thickness= "<<AlluvialColumnThickness()<<std::endl;
  }





  double sectionBase_after  = getSectionBase();
  double columnheight_after = AlluvialColumnThickness();
  double difference = (sectionBase_after)- (z-columnheight_after);
  if(z-columnheight_after != sectionBase_after && columnheight_after > 0.0){
    if(difference < -0.001 || difference > 0.001){
      std::cout<<"   "<<std::endl;
      std::cout<<"Finished ErodepSimple, no elevation balance"<<std::endl;
      std::cout<<"at "<<getI()<<" "<<getJ()<<" dh was "<<dhtotal<<std::endl;
      std::cout<<"z = "<<z<<" alluv= "<<columnheight_after<<" base= "<<sectionBase_after<<std::endl;
      std::cout<<"z-columnheight= "<<z-columnheight_after<<"diff = "<< (sectionBase_after)- (z-columnheight_after)<<std::endl;
      std::cout<<" "<<std::endl;
      std::cout<<" "<<std::endl;
      exit(1);
    }
  }





  //-----------end of space for debug prints---------------------------



} // End of EroDepSimple

/**************************************************************
 ** tStratNode::addtoLayer(int i, int g, double val, double tt,current)
 **
 ** This function is called from the EroDep which updates layering.
 ** Used when material is added/removed to/from a layer, grain size
 ** by grain size. i.e. texture may change
 ** i = layer to add to
 ** g = grain size class
 ** val = amount of that size to deposit
 ** tt = time of deposition/erosion
 **
 ** created by NG
 **
 **    Modifications:
 **      - 3/30/00: if erosion causes a layer to have zero
 **        thickness, the layer is removed. (GT)
 *****************************************************************/
void tStratNode::addtoLayer(int i, int g, double val, double tt, double current)
{
  // For adding material to a layer, grain size by grain size
  tListIter<tLayer> ly ( layerlist );
  tLayer  * hlp;
  hlp=ly.FirstP();

  int n=0;
  while(n<i){
    n++;
    hlp=ly.NextP();
  }

  // Find the current flow patterns of meandering nodes
  // surrounding this stratnode, compass is  0 < 360 for
  // meanders or -1 for floodplain nodes.

  hlp->addPaleoCurrent(val,current);

  if(tt>0)
    hlp->setRtime( tt );

  hlp->addDgrade(g,val);

  //Although check really should be here, I think there may
  //be an issue because this function is called from a loop
  //and if the layer is removed before the loop is through -> big trouble
  // GT re-added this line, 3/00 ... and then took it away again 8/02.
  // Why? Layer might be temporarily zero while eroding one size but
  // depositing another. Check needs to be outside this routine.
  //if(hlp->getDepth() <= 0)
  //    removeLayer(i);
}

/******************************************************************
 **  tStratNode::addtoLayer(int i, double val,tt)
 **  This addtolayer is only used for removing material from a layer.
 **  This erosion will not change the texture of a layer, only depth.
 **  As always, the depth of each grain size class is also updated.
 **  i = layer to deplete
 **  val = amount to deplete by
 **  Returns an array containing the texture of material removed.
 **
 **  created NG
 **  Since only for erosion, nic modified this so that the time
 **  is not passed, since time will not be reset for erosion.
 ******************************************************************/
tArray<double> tStratNode::addtoLayer(int l, double val, double tt)
{
  assert( val<0.0 ); // Function should only be called for erosion

  tListIter<tLayer> ly ( layerlist );
  tLayer  * hlp;
  hlp=ly.FirstP();
  tArray<double> ret(numg);
  double amt;

  int n=0;
  while(n<l){
    n++;
    hlp=ly.NextP();
  }

  //nb. Erosion. Do not reset or change the paleocurrent values

  if(hlp->getDepth()+val>1e-7)
    {
      // have enough material in this layer to fufill all erosion needs
      amt=hlp->getDepth();			     // Total thickness
      n=0;
      while(n<numg)
	{
	  ret[n]=hlp->getDgrade(n)*val/amt;             // return this, will be passed on to top layer
	  hlp->addDgrade(n,hlp->getDgrade(n)*val/amt);  // modify the stratigraphic layer
	  n++;
	}

      if(tt>0)
	hlp->setRtime( tt );


      return ret;        // Return
    }
  else
    {
      // need to remove entire layer
      n=0;
      while(n<numg)
	{
	  ret[n]=-1*hlp->getDgrade(n);
	  n++;
	}
      removeLayer(l);
      return ret;
    }

}




/*******************************************************
 ** tStratNode::removeLayer(int i)
 **
 ** Removes layer i.
 **
 ** called by EroDep which updates layers.
 **
 **    Modifications:
 **      - 3/30/00: apparently the existing code was
 **        actually removing layer i+1. Bug fixed (GT).
 *******************************************************/
void tStratNode::removeLayer(int i)
{
  //tLayer  * hlp;
  //tLayer niclay;
  //hlp=ly.FirstP();

  // GT added code in bug fix attempt 3/00:
  tLayer lay;

  if( i==0 )
    layerlist.removeFromFront( lay );
  else if( (i+1)==layerlist.getSize() )
    layerlist.removeFromBack( lay );
  else
    {
      tListIter<tLayer> ly ( layerlist );
      int n;
      for( n=1; n<i; n++ )
	ly.Next();
      n = layerlist.removeNext( lay, ly.NodePtr() );
      assert( n>0 );
    }

  // end of GT added code

  /*    int n=0;

  while(n<i){
  n++;
  hlp=ly.NextP();
  }

  if(i+1==layerlist.getSize()){
  layerlist.removeFromBack(*hlp );
  }
  else{
  n=layerlist.removeNext((*hlp), layerlist.getListNode(hlp) );
  }
  if(n==0){
  //       n=0;
  //         while(n<layerlist.getSize()){
  //            std::cout << "layer " << n+1 << " node ID "<< getID()<< std::endl;
  //            niclay = layerlist.getIthData(n);
  //            std::cout << "layer creation time is " << getLayerCtime(n) << std::endl;
  //            std::cout << "layer recent time is " << getLayerRtime(n) << std::endl;
  //            std::cout << "layer depth is " << getLayerDepth(n) << std::endl;
  //            std::cout << "layer erodibility is " << getLayerErody(n) << std::endl;
  //            std::cout << "is layer sediment? " << getLayerSed(n) << std::endl;
  //            std::cout << "dgrade 1 is " << getLayerDgrade(n,0) << std::endl;
  //            std::cout << "dgrade 2 is " << getLayerDgrade(n,1) << std::endl;
  //            n++;
  //         }

  ReportFatalError("couldn't remove next layer");
  }*/

}

/*****************************************************************
 ** tStratNode::makeNewLayerBelow(int i, int sd, double erd, tArray<double> sz, double tt)
 ** Makes a new layer below layer i.
 ** if i<0 then make new top layer.
 ** sd = sediment flag
 ** erd = erodibility
 ** sz = array of depths of each grain size class
 ** tt = time
 ** Creation and recent time set to current time, exposure time updated
 ** in erodep.
 ********************************************************************/
void tStratNode::makeNewLayerBelow(int i, tLayer::tSed_t sd, double erd,
				   tArray<double> const& sz,
				   double tt, double current)
{
  tLayer hlp, niclay;
  int n;

  // Debug, Quintijn
  //std::cout << "In makeNewLayerBelow" << '\n';
  //std::cout << "numg = " << numg << '\n';
  //std::cout << "i= " << i << "size[0]= " << sz[0] << "tt= " << tt << '\n';
  //exit(1);


  hlp.setCtime(tt);
  hlp.setRtime(tt);
  hlp.setEtime(0.);
  n=0;
  hlp.setDgradesize(numg);
  while(n<numg){
    assert( sz[n]>=0.0 );
    hlp.setDgrade(n, sz[n]);
    n++;
  }
  assert( hlp.getDepth()>0.0 );
  hlp.setErody(erd);
  hlp.setSed(sd);


  hlp.setPaleoCurrent(current);

  if(i>=0){

    tListIter<tLayer> ly ( layerlist );
    tLayer  * pt;
    pt=ly.FirstP();

    n=0;

    while(n<i){
      n++;
      pt=ly.NextP();
    }

    layerlist.insertAtNext(hlp, layerlist.getListNode(pt));
  }
  else{
    layerlist.insertAtFront(hlp);
  }
}

int tStratNode::getNumLayer() const
{
  return layerlist.getSize();
}

//Returns depth of all the layers - pretty useless for the current model
double tStratNode::getTotalLayerDepth() const
{
  const double sz = layerlist.getSize();
  int i=0;
  double elev = 0;
  while(i<sz)
    {
      elev += getLayerDepth(i);
      i++;
    }
  return elev;
}


/******************************************************************
 **  tStratNode::CopyLayerList
 **
 **  Assigns the layer list in fromNode to this node. Information
 **  about the pre-existing layer list is lost.
 **
 **  Created in order to assign stratigraphy from a footwall point
 **  to the surface node when the fault plane is exhumed.
 **
 **
 ********************************************************************/
void tStratNode::CopyLayerList( tStratNode const * fromNode )
{
  assert( fromNode!=0 );
  layerlist = fromNode->layerlist;
}


//-----------------END OF CLASS STRATNODE-----------------------------------
