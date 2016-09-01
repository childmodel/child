/**************************************************************************/
/**
 **  @file tLNode.cpp
 **  @brief Functions for derived class tLNode and its member classes
 **
 **  Modifications:
 **    - changes related to addition of Vegetation module (GT 1/2000):
 **       - removed tSurface class & fn
 **       - modified constructors (tLNode) to handle tau, tauc, vegCover
 **    - fixed problem with layer initialization in copy constructor
 **      (gt, 2/2000; see below) -- OOPS! (8/02)
 **    - modifications to EroDep to fix bug related to layers under
 **      simultaneous erosion of one size and deposition of another
 **      (GT, 8/2002)
 **
 **  $Id: tLNode.cpp,v 1.140 2007-08-07 02:23:56 childcvs Exp $
 */
/**************************************************************************/

#include <assert.h>
#include <math.h>
#include "../errors/errors.h"
#include "tLNode.h"
//#define kBugTime 5000000

#include "../tStratGrid/tStratGrid.h"

//Sets the total layer depth.  While updating depth, dgrade info is
//automatically updated to keep the same texture.
//So if texture is going to change too - you MUST update dgrade info.
//(use setDgrade)
//which will automatically update total depth.
void tLayer::setDepth( double dep )
{
  assert( dep>0.0 );
  if(dgrade.getSize()>0 && depth>0){
    double sum=0;
    size_t i=0;
    while(i<dgrade.getSize()){
      const double prop = dgrade[i]/depth;
      dgrade[i]=prop*dep;
      assert( dgrade[i]>=0.0 );
      sum+=prop;
      i++;

    }
    if(fabs(sum-1.)>0.01)
      ReportFatalError("Somewhere grain sizes got messed up");
    depth=dep;
  }
  else
    depth=dep;
}

void tLayer::setDgrade( size_t i, double size )
{
  assert( i<dgrade.getSize() );
  assert( size>=0.0 );
  dgrade[i]=size;
  // Automatically update depth when dgrade is changed
  size_t j=0;
  double sum=0;
  while(j<dgrade.getSize()){
    sum += dgrade[j];
    j++;
  }
  depth=sum;
}

tMeander::tMeander()                                              //tMeander
  :
newx(0.), newy(0.),
deltax(0.), deltay(0.), zoldright(0.), zoldleft(0.), bankrough(0.),
xyzd(4),
reachmember(false),
meander(false)
{
  if (0) //DEBUG
    std::cout << "  tMeander()" << std::endl;
}

tMeander::tMeander( const tMeander &orig )                        //tMeander
  :
newx(orig.newx), newy(orig.newy),
deltax(orig.deltax), deltay(orig.deltay),
zoldright(orig.zoldright), zoldleft(orig.zoldleft), bankrough(orig.bankrough),
xyzd( orig.xyzd ),
reachmember(orig.reachmember),
meander(orig.meander)
{
  if (0) //DEBUG
    std::cout << "  tMeander( orig )" << std::endl;
}

tMeander::tMeander( bool state, double x, double y )                //tMeander
  :
newx(x), newy(y),
deltax(0.), deltay(0.), zoldright(0.), zoldleft(0.), bankrough(0.),
xyzd(4),
reachmember(false),
meander(state)
{
  if (0) //DEBUG
    std::cout << "  tMeander( state, x, y )" << std::endl;
}

tMeander::~tMeander()                                             //tMeander
{
  if (0) //DEBUG
    std::cout << "    ~tMeander()" << std::endl;
}

//assignment
const tMeander &tMeander::operator=( const tMeander &right )     //tMeander
{
  if( &right != this )
    {
      newx = right.newx;
      newy = right.newy;
      deltax = right.deltax;
      deltay = right.deltay;
      zoldright = right.zoldright;
      zoldleft = right.zoldleft;
      bankrough = right.bankrough;
      xyzd = right.xyzd;
      reachmember = right.reachmember;
      meander = right.meander;
    }
  return *this;
}

tBedrock::tBedrock()                                             //tBedrock
  :
erodibility(0.)
{
  if (0) //DEBUG
    std::cout << "  tBedrock()" << std::endl;
}

tBedrock::tBedrock( const tBedrock &orig )                       //tBedrock
  :
erodibility(orig.erodibility)
{
  if (0) //DEBUG
    std::cout << "  tBedrock( orig )" << std::endl;
}

tBedrock::~tBedrock()                                            //tBedrock
{
  if (0) //DEBUG
    std::cout << "    ~tBedrock()" << std::endl;
}

//assignment
const tBedrock &tBedrock::operator=( const tBedrock &right )     //tBedrock
{
  if( &right != this )
    {
      erodibility = right.erodibility;
    }
  return *this;
}


/***** Functions for class tRegolith **************************************/

/**************************************************************************\
 **
 **  Constructors:
 **   - Default: initializes variables to zero.
 **   - Input file: takes a tInputFile reference as an argument, and
 **              reads the default initial thickness from the file.
 **   - Copy: copies all fields of tLNode and calls copy constructors for
 **              embedded objects.
 **
\**************************************************************************/

tRegolith::tRegolith()                                          //tRegolith
  :
thickness(0.),
dgrade()
{
  if (0) //DEBUG
    std::cout << "  tRegolith()" << std::endl;
}

tRegolith::tRegolith( const tRegolith &orig )                   //tRegolith
  :
thickness(orig.thickness),
dgrade( orig.dgrade )
{
  if (0) //DEBUG
    std::cout << "  tRegolith( orig ) " << thickness << std::endl;
}


tRegolith::~tRegolith()                                         //tRegolith
{
  if (0) //DEBUG
    std::cout << "    ~tRegolith()" << std::endl;
}

//assignment
const tRegolith &tRegolith::operator=( const tRegolith &right )     //tRegolith
{
  if( &right != this )
    {
      thickness = right.thickness;
      dgrade = right.dgrade;
    }
  return *this;
}

tChannel::tChannel()                                             //tChannel
  :
drarea(0.),q(0.),
mdFlowPathLength(0.),
chanwidth(0.),hydrwidth(0.),
channrough(0.),hydrnrough(0.),
chandepth(0.),hydrdepth(0.),
chanslope(0.),hydrslope(0.),
diam(kVeryHigh),
migration()
{
  if (0) //DEBUG
    std::cout << "  tChannel()" << std::endl;
}

tChannel::tChannel( const tChannel &orig )                       //tChannel
  :
drarea(orig.drarea),q(orig.q),
mdFlowPathLength(orig.mdFlowPathLength),
chanwidth(orig.chanwidth),hydrwidth(orig.hydrwidth),
channrough(orig.channrough),hydrnrough(orig.hydrnrough),
chandepth(orig.chandepth),hydrdepth(orig.hydrdepth),
chanslope(orig.chanslope),hydrslope(orig.hydrslope),
diam(orig.diam),
migration( orig.migration )
{
  if (0) //DEBUG
    std::cout << "  tChannel( orig )" << std::endl;
}

tChannel::~tChannel()                                            //tChannel
{
  if (0) //DEBUG
    std::cout << "    ~tChannel()" << std::endl;
}
//assignment
const tChannel &tChannel::operator=( const tChannel &right )     //tChannel
{
  if( &right != this )
    {
      drarea = right.drarea;
      q = right.q;
      chanwidth = right.chanwidth;
      chandepth = right.chandepth;
      channrough = right.channrough;
      chanslope = right.chanslope;
      hydrwidth = right.hydrwidth;
      hydrdepth = right.hydrdepth;
      hydrnrough = right.hydrnrough;
      hydrslope = right.hydrslope;
      mdFlowPathLength = right.mdFlowPathLength;
      diam = right.diam;
      migration = right.migration;
    }
  return *this;
}

/***** Functions for class tLNode *****************************************/

/**************************************************************************\
 **  Constructors:
 **   - Default: initializes several variables to zero and calls
 **              constructors for embedded objects.
 **   - Input file: takes a tInputFile reference as an argument, and
 **              passes it to the constructor(s) for embedded objects so
 **              that they can read values for default parameters, etc.
 **   - Copy: copies all fields of tLNode and calls copy constructors for
 **              embedded objects.
 **
 **    Modifications:
 **     - fixed an odd error with the layer initialization in the copy
 **       constructor. The use of the embedded object initialization
 **       syntax " : layerlist( orig.layerlist ) ", when compiled with
 **       cxx, appeared to erroneously look for a tList constructor of
 **       type tList( int ), which does not exist. Thus, the layerlist
 **       was not duplicated but instead simply handed by default
 **       memberwise copy, so that the new list pointed to the original
 **       one. Once the problem was discovered, the fix was simple:
 **       use the assignment operator instead. This was probably the
 **       source of the problem noted by nmg below. (gt, 2/2000)
 **     - Changed tLNode(infile) to call new tNode(infile), which reads
 **       static bool tNode::freezeElevations, which is referenced within
 **       tNode::ChangeZ(dz); if freezeElevations = true, ChangeZ will do
 **       nothing. (SL, 11/2010)
 **
\**************************************************************************/

// initialize static members numg and grade
size_t tLNode::numg = 0;
tArray<double> tLNode::grade = 1;
double tLNode::maxregdep = 1;
double tLNode::KRnew = 1.0;
double tLNode::new_sed_bulk_density_ = kDefaultSoilBulkDensity;

tLNode::tLNode()                                                   //tLNode
  :
tNode(), vegCover(), rock(), reg(), chan(),
flood(kNotFlooded), flowedge(0), tracer(0),
dzdt(0.), drdt(0.), tau(0.), taucb(0.), taucr(0.), qs(0.),
qsm(),
qsin(0.),
qsinm(),
qsdin(0.),
qsdinm(),
uplift(0.),
layerlist(),
stratNode(0),
qsubsurf(0.),
is_masked_(false),
netDownslopeForce(0.),
is_moving_(false),
public1(-1)
{
  if (0) //DEBUG
    std::cout << "=>tLNode()" << std::endl;
}

tLNode::tLNode( const tInputFile &infile )
:
tNode( infile ), vegCover(), rock(), reg(), chan(),
flood(kNotFlooded), flowedge(0), tracer(0),
dzdt(0.), drdt(0.), tau(0.), taucb(0.), taucr(0.), qs(0.),
qsm(),
qsin(0.),
qsinm(),
qsdin(0.),
qsdinm(),
uplift(0.),
layerlist(),
stratNode(0),
qsubsurf(0.),
is_masked_(false),
netDownslopeForce(0.),
is_moving_(false),
public1(-1)
{
  char add[2], name[20];
  double help, sum, sumbr;
  tLayer layhelp, niclay;
  tArray<double> dgradehelp;
  tArray<double> dgradebrhelp;
  
  if (0) //DEBUG
    std::cout << "=>tLNode( infile )" << std::endl;
	
  // Modified to read TAUC for both bedrock, TAUCB, and regolith, TAUCR
  taucb = infile.ReadItem( taucb, "TAUCB" );
  taucr = infile.ReadItem( taucr, "TAUCR" );
  
  {
    int tmp_;
    tmp_ = infile.ReadItem( tmp_, "NUMGRNSIZE" );
    assert(tmp_ >= 0);
    numg = tmp_;
  }
  grade.setSize( numg );
  maxregdep = infile.ReadItem( maxregdep, "MAXREGDEPTH" );
  KRnew = infile.ReadItem( KRnew, "KR" );
  if( KRnew<0.0 )
    ReportFatalError( "Erodibility factor KR must be positive." );
  new_sed_bulk_density_ = infile.ReadDouble( "SOILBULKDENSITY", false );
  if( new_sed_bulk_density_<=0.0 ) new_sed_bulk_density_ = kDefaultSoilBulkDensity;
  
  {
    size_t i = 0;
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
  }
  
  qsm.setSize( numg );
  qsinm.setSize( numg );
  
  accumdh.setSize(2);
  accumdh[0]=0.0;
  accumdh[1]=0.0;
  
  bool optReadLayer = infile.ReadBool( "OPTREADLAYER" );
  
  if(!optReadLayer){
    
    dgradehelp.setSize( numg );
    dgradebrhelp.setSize( numg );
    sum = 0;
    sumbr = 0;
    size_t i=0;
    add[0]='1';
    
    while ( i<numg ){
      // Reading in proportions for intital regolith and bedrock
      strcpy( name, "REGPROPORTION");
      strcat( name, add );
      help = infile.ReadItem( help, name);
      dgradehelp[i]=help;
      sum += help;
      strcpy( name, "BRPROPORTION");
      strcat( name, add );
      std::cout<<"In tLNode, reading '"<<name<<std::endl;
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
      help = infile.ReadDouble( "ROCKDENSITYINIT", false );
      if( help<=0.0 ) help = kDefaultRockBulkDensity;
      layhelp.setBulkDensity(help);
      
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
      
      if(0) //DEBUG
      {
        std::cout<<"Just made BR layer thick=" << layhelp.getDepth()
        << " and dgrades:\n";
        for( i=0; i<numg; i++ )
          std::cout << i << "=" << layhelp.getDgrade(i) << std::endl;
      }
      
      // Regolith layer items read in and set
      layhelp.setSed(tLayer::kSed);
      help = infile.ReadItem( help, "KR" );
      if( help<0.0 )
        ReportFatalError( "Erodibility factor KR must be positive." );
      layhelp.setErody(help);
      help=infile.ReadDouble( "SOILBULKDENSITY", false );
      if( help<=0.0 ) help = kDefaultSoilBulkDensity;
      layhelp.setBulkDensity( help );
      help = infile.ReadItem( help, "REGINIT");
      if(help > maxregdep){
        // too much regolith, create two layers the bottom layer is made here
        const double extra = help - maxregdep;
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
        
        if(0) //DEBUG
        {
          std::cout<<"1Just made ALLUV layer thick=" << layhelp.getDepth()
          << " and dgrades:\n";
          for( i=0; i<numg; i++ )
            std::cout << i << "=" << layhelp.getDgrade(i) << std::endl;
        }
        
        // the top regolith layer is now made
        layhelp.setRtime(0.);
        i=0;
        while(i<numg){
          layhelp.setDgrade(i, maxregdep*dgradehelp[i]);
          i++;
        }
        layerlist.insertAtFront( layhelp );
        
        if(0) //DEBUG
        {
          std::cout<<"2Just made ALLUV layer thick=" << layhelp.getDepth()
          << " and dgrades:\n";
          for( i=0; i<numg; i++ )
            std::cout << i << "=" << layhelp.getDgrade(i) << std::endl;
        }
        
      }
      else{
        // create only one regolith layer
        i=0;
        while(i<numg){
          layhelp.setDgrade(i, help*dgradehelp[i]);
          i++;
        }
        
        layerlist.insertAtFront( layhelp );
        
        if(0) //DEBUG
        {
          std::cout<<"3Just made ALLUV layer thick=" << layhelp.getDepth()
          << " and dgrades:\n";
          for( i=0; i<numg; i++ )
            std::cout << i << "=" << layhelp.getDgrade(i) << std::endl;
        }
        
      }
    } // if REGINIT > 0
    
    else{
      // no regolith, so by default everything is bedrock
      // Bedrock layer items set
      help = infile.ReadItem( help, "KB");
      if( help<0.0 )
        ReportFatalError( "Erodibility factor KB must be positive." );
      layhelp.setErody(help);
      layhelp.setSed(tLayer::kBedRock);
      help = infile.ReadDouble( "ROCKDENSITYINIT", false );
      if( help<=0.0 ) help = kDefaultRockBulkDensity;
      layhelp.setBulkDensity(help);
      layhelp.setDgradesize(numg);
      size_t i=0;
      help = infile.ReadItem( help, "BEDROCKDEPTH");
      while(i<numg){
        layhelp.setDgrade(i, help*dgradebrhelp[i]);
        i++;
      }
      layerlist.insertAtBack( layhelp );
    }
    
    if (0) //DEBUG
      std::cout << layerlist.getSize() << " layers created " << std::endl;
  }
  
}

tLNode::tLNode( const tLNode &orig )                               //tLNode
  : tNode( orig ), 
    vegCover( orig.vegCover ), 
    rock( orig.rock ),
    reg( orig.reg ), 
    chan( orig.chan ),
    flood(orig.flood), 
    flowedge(0), 
    tracer(orig.tracer),
    dzdt(orig.dzdt), 
    drdt(orig.drdt), 
    tau(orig.tau), 
    taucb(orig.taucb), 
    taucr(orig.taucr), 
    qs(orig.qs),
    qsm(orig.qsm),
    qsin(orig.qsin),
    qsinm(orig.qsinm ),
    qsdin(orig.qsdin),
    qsdinm(orig.qsdinm ),
    uplift(orig.uplift),
    layerlist(),
    stratNode(0),
    accumdh(orig.accumdh),
    qsubsurf(orig.qsubsurf),
    is_masked_(false),
    netDownslopeForce(orig.netDownslopeForce),
    cumulative_ero_dep_(orig.cumulative_ero_dep_),
    cumulative_sed_xport_volume_(orig.cumulative_sed_xport_volume_),
    is_moving_(orig.is_moving_),
    public1(orig.public1)
{

  //Be aware that the copy constructor should be called using the
  //following syntax
  //tLNode *newnode = new tLNode( *oldtLNode );
  //NOT
  //tLNode newnode(*oldtLNode)
  //The latter seems to call a default constructor which does not
  //properly copy the layerlist

  layerlist = orig.layerlist;
  if (0) //DEBUG
    std::cout << "=>tLNode( orig )" << std::endl;
}

tLNode::~tLNode()                                                  //tLNode
{
  if (0) //DEBUG
    std::cout << "    ~tLNode()" << std::endl;
  flowedge = 0;
  stratNode = 0;
}

//assignment
const tLNode &tLNode::operator=( const tLNode &right )                  //tNode
{
  if( &right != this )
    {
      tNode::operator=( right );
      vegCover = right.vegCover;
      rock = right.rock;
      //surf = right.surf;
      reg = right.reg;
      chan = right.chan;
      flood = right.flood;
      flowedge = right.flowedge;
      tracer = right.tracer;
      dzdt = right.dzdt;
      drdt = right.drdt;
      tau = right.tau;
      taucb = right.taucb;
      taucr = right.taucr;
      qs = right.qs;
      qsm = right.qsm;
      qsin = right.qsin;
      qsinm = right.qsinm;
      qsdin = right.qsdin;
      qsdinm = right.qsdinm;
      uplift = right.uplift;
      layerlist = right.layerlist;
      stratNode = right.stratNode;
      accumdh = right.accumdh;
      qsubsurf = right.qsubsurf;
      is_masked_ = right.is_masked_;
      netDownslopeForce = right.netDownslopeForce;
      cumulative_ero_dep_ = right.cumulative_ero_dep_;
      cumulative_sed_xport_volume_ = right.cumulative_sed_xport_volume_;
	  is_moving_ = right.is_moving_;
      public1 = right.public1;
    }
  return *this;
}

//"get" and "set" functions; most simply return or set a data value, respectively:


//NG changed getDiam 02/1999
double tLNode::getDiam() const {
  double di = 0;
  for(size_t i=0; i<numg; i++){
    di+=grade[i]*getLayerDgrade(0,i)/getLayerDepth(0);
  }
  return di;
}

/************************************************************************\
 **  GetSlopeMeander: Computes and returns the slope of the node's flowedg, or
 **  zero if the slope is less than zero.
 **
 **  The name of this function makes NG cringe.
 **  TODO Change to CALCSLOPE!!!!
 **
 **  Assumptions: edge lengths up to date and nonzero, flowedg's up to
 **    date.
 **  Modified, 2/19/98: to return a reach slope over some distance
 **   independent of the discretization for meandering nodes.
 **   Nominally, we set that distance to 10 channel widths; we limit this
 **   reach slope calculation to meandering nodes because they are the only
 **   ones for which hydraulic geometry is now calculated; if we want to
 **   use a reach slope everywhere, then we need to also calculate the
 **   hydraulic geometry everywhere.
 **   AND make sure that the node 'crn' is not at a lower elevation than
 **   the outlet; in the latter case, it returns slope=0; if it runs into
 **   an infinite loop either searching for the 10-widths-downstream node
 **   or the outlet node, it returns a negative number (-1)
\************************************************************************/
#define kLargeNumber 1000000
double tLNode::CalcSlopeMeander()
{
  int ctr;
  bool breakit;
  double rlen, curlen, slp, delz, downz;
  tLNode *dn, *on;

  assert( flowedge != 0 );
  assert( flowedge->getLength()>0 ); // failure means lengths not init'd

  if( Meanders() )
    {
      rlen = 10.0 * chan.chanwidth;
      curlen = 0.0;
      dn = this;
      ctr = 0;
      //built a couple of 'firewalls' in these loops; first, quit loop if
      //we have a flowedge with zero length; second, if we've somehow gotten
      //an infinite loop, return a negative number as a flag to the calling
      //routine that it needs to update the network (tStreamNet::UpdateNet())
      //and reaches (tStreamMeander::MakeReaches). We've run across this loop
      //bug in a call from tStreamMeander::InterpChannel, which is called by
      //tStreamMeander::MakeReaches.
      while( curlen < rlen && dn->getBoundaryFlag() == kNonBoundary &&
             dn->flowedge->getLength() > 0 )
	{
	  ctr++;
	  if( ctr > kLargeNumber )
	    {
	      ReportFatalError("infinite loop in tLNode::GetSlope(), 1st loop");
	    }
	  curlen += dn->flowedge->getLength();
	  dn = dn->getDownstrmNbr();
	  assert( dn != 0 );
	}
      if(curlen <= 0){
	std::cout<<"Error in getSlope() in 1st loop, curlen is "<<curlen<<std::endl;
	std::cout<<"Do you have an obstruction/ pond along the meander channel ? \n";
	TellAll();
	assert(curlen>0.0);
      }

      assert( curlen > 0 );
      downz = dn->z;
      //if( timetrack >= kBugTime ) std::cout << "GetSlope 1; " << flush;
      delz = z - downz;
      //if( timetrack >= kBugTime ) std::cout << "GS 2; " << flush;
      slp = delz / curlen;
      //if( timetrack >= kBugTime ) std::cout << "GS 3; " << flush;
      on = dn;
      ctr = 0;
      breakit=false;
      while( on->getBoundaryFlag() == kNonBoundary &&
             on->flowedge->getLength() > 0 )
	{
	  ctr++;
	  if( ctr > kLargeNumber )
	    {
	      std::cout<<"At node ID,x,y,z "<< getID()<<' ' <<getX()<<' '<<getY()<<' '<<getZ()<<std::endl;
	      std::cout<<"Not able to find a boundary node of type 1 or 2 (Meander cannot find a gridboundary ?)\n";
	      std::cout<<"breaking the infinite loop, but this is a bug !			\n";
	      //TellAll();
	      //ReportFatalError("an infinite loop in tLNode::GetSlope(), 2nd loop");
	      breakit=true;
	    }
	  if(breakit){
	      tLNode *cn;
	      int counter;
	      counter=0;
	      cn=this;
	      while(cn!=NULL && counter <=50){
	         std::cout<< "M-Nodes "
		     <<cn->getID()<<' '
		     <<cn->getX()<< ' ' << cn->getY()<<' '<<cn->getZ()
		     <<" A= "<<cn->getDrArea()<<" Q= "<<cn->getQ()
		     <<" W= "<<cn->getChanWidth()<< " Mndr="<<cn->Meanders()
		     <<"Flood= "<< FloodName(cn->getFloodStatus()) <<std::endl;
	  	 cn=cn->getDownstrmNbr();
	  	 counter++;
	      }
	      break;
	   //exit(1);
	   }
	                                             // Leave the infinite loop
	  on = on->getDownstrmNbr();
	}
      if( z - on->z < 0.0 ) slp = 0.0;
    }     // for meandering nodes
  else{ // for non-menadering modes
    slp = (z - getDownstrmNbr()->z ) / flowedge->getLength();
  }

  //if( timetrack >= kBugTime ) std::cout << "GS 4; " << std::endl << flush;
  if( slp>=0.0 ) return slp;
  else return 0.0;
}
#undef kLargeNumber

tLNode *tLNode::getDSlopeDtMeander( double &curlen )
{
  double rlen;
  tLNode *dn;
  rlen = 10.0 * chan.chanwidth;
  curlen = 0.0;
  dn = this;
  while( curlen < rlen && dn->getBoundaryFlag() == kNonBoundary )
    {
      curlen += dn->flowedge->getLength();
      dn = dn->getDownstrmNbr();
    }
  assert( curlen > 0 );
  return dn;
}

/****************************************************************************\
**
**  tLNode::FindInitFlowDir
**
**  Initialize flow directions such that each active (non-boundary) node
**  flows to another active node (or open boundary node). This initialization
**  process allows the FlowDirs function to assume that the previous flowedg
**  is always valid in the sense that it doesn't point to a closed boundary
**  (which otherwise could happen when the mesh is first read in).
**
**  Modifies: node flow directions (flowedge)
**  Written 12/1/97 gt.
**  Modified: 1/5/99 SL -- guts moved from tStreamNet::InitFlowDirs()
**  Updated for Oxford CHILD, 9/2003 SL.
**
\****************************************************************************/
#define kMaxSpokes 100
void tLNode::FindInitFlowDir() // overrides tNode; was tStreamNet::
{
   int ctr = 0;
   // Start with the node's default edge
   tEdge* flowedg = getEdg();
   assert( flowedg!=0 );

   // As long as the current edge is a no-flow edge, advance to the next one
   // counter-clockwise
   while( !flowedg->FlowAllowed() )
   {
      flowedg = flowedg->getCCWEdg();
      assert( flowedg != 0 );
      if( ++ctr > kMaxSpokes ) // Make sure to prevent std::endless loops
      {
         std::cerr << "Mesh error: node " << getID()
              << " appears to be surrounded by closed boundary nodes"
              << std::endl;
         ReportFatalError( "Bailing out of InitFlowDirs()" );
      }
   }
   setFlowEdg( flowedg );
   assert( getFlowEdg() != 0 );
}

/****************************************************************************\
**
**  tLNode::FindFlowDir
**
**  Computes flow directions for each node using a
**  steepest-direction routing algorithm. Flow proceeds along the
**  steepest of the directed edges (spokes) around a given node.
**  This edge is assigned as the node's flowedg. Nodes at which
**  the steepest slope is negative (uphill) are flagged as sinks.
**  (If the lake-fill option is active, sink nodes are processed
**  in order to find an outlet; otherwise, all runoff reaching a
**  sink node is assumed to evaporate and/or infiltrate).
**
**      Parameters:     none
**      Called by: tStreamNet::FlowDirs()
**      Modifies: node flow directions and flood status
**      Assumes:
**       - each node points to a valid flow edge (returned by getFlowEdg())
**       - each edge has a valid counter-clockwise edge
**       - edge slopes are up to date
**      Returns: false if flow direction unchanged, true otherwise.
**      Updated: 12/19/97 SL; 12/30/97 GT
**  Modified: 1/5/99 SL -- guts moved from tStreamNet::FlowDirs()
**      Updated for Oxford CHILD, 9/2003 SL.
**
\****************************************************************************/
bool tLNode::FindFlowDir() // was tStreamNet::
{
   setFloodStatus( kNotFlooded );  // Init flood status flag
   tEdge* firstedg = flowedge;
      if( unlikely(firstedg == 0) ) {
          TellAll();
	  assert( 0 );
      }
   double slp = firstedg->getSlope();
   tEdge* nbredg = firstedg;
   tEdge* curedg = firstedg->getCCWEdg();
   int ctr = 0;

   // Check each of the various "spokes", stopping when we've gotten
   // back to the beginning
   while( curedg!=firstedg )
   {
      assert( curedg != 0 );
      if ( curedg->getSlope() > slp && curedg->FlowAllowed() )
      {
         slp = curedg->getSlope();
         nbredg = curedg;
      }
      curedg = curedg->getCCWEdg();
      ++ctr;
      if( unlikely( ctr > kMaxSpokes ) ) // Make sure to prevent std::endless loops
      {
         std::cerr << "Mesh error: node " << getID()
              << " going round and round" << std::endl;
         ReportFatalError( "Bailing out of FlowDirs()" );
      }
   }
   setFlowEdg( nbredg );
   //Nicole changed the line below which is commented out, replaced
   //it with the new if because of the new function WarnSpokeLeaving
   //which can set a node as a boundary node.
   //curnode->setFloodStatus( ( slp>0 ) ? kNotFlooded : kSink );// (NB: opt branch pred?)
   if( slp > 0 && getBoundaryFlag() != kClosedBoundary )
       setFloodStatus( kNotFlooded );
   else
       setFloodStatus( kSink );
   assert( getFlowEdg() != 0 );
   if( nbredg == firstedg ) return false;
   return true;
}

/****************************************************************************\
**
**  tLNode::FindDynamicFlowDir
**
**  Computes flow directions for each node using a
**  steepest-direction routing algorithm. Flow proceeds along the
**  steepest of the directed edges (spokes) around a given node.
**  This edge is assigned as the node's flowedg. Nodes at which
**  the steepest slope is negative (uphill) are flagged as sinks.
**  (If the lake-fill option is active, sink nodes are processed
**  in order to find an outlet; otherwise, all runoff reaching a
**  sink node is assumed to evaporate and/or infiltrate).
**
**      Parameters:     none
**      Called by: tStreamNet::FlowDirs()
**      Modifies: node flow directions and flood status
**      Assumes:
**       - each node points to a valid flow edge (returned by getFlowEdg())
**       - each edge has a valid counter-clockwise edge
**       - edge slopes are up to date
**      Updated: 12/19/97 SL; 12/30/97 GT
**  Modified: 1/5/99 SL -- guts moved from tStreamNet::FlowDirs()
**
**  Only difference from FindFlowDir() is that this calls CalcSlope()
**  for spokes instead of getSlope() the first time.
**  7/2002, SL:
**  9/2003 SL: Also calculates edge length so that UpdateMesh is not needed
**    to redetermine flowedges
**    Also checks that flowedge originates at node and, if not, sets
**    flowedge temporarily to edg. This problem can arise because this
**    function is often called right after node addition or deletion,
**    and the edge that was the flowedge may no longer be connected
**    to this edge!
\****************************************************************************/
bool tLNode::FindDynamicFlowDir() // overrides tNode; was tStreamNet::
{
   bool fechng = false;
   setFloodStatus( kNotFlooded );  // Init flood status flag
   // make sure flowedge starts at "this" node
   if( unlikely( flowedge->getOriginPtr() != this ) ){
      fechng = true;
      FindInitFlowDir();
   }

   tEdge* firstedg = flowedge;
   if( unlikely(firstedg == 0) ) {
      TellAll();
      assert( 0 );
   }
   firstedg->CalcLength();
   double slp = firstedg->CalcSlope();
   tEdge* nbredg = firstedg;
   tEdge* curedg = firstedg->getCCWEdg();
   int ctr = 0;

   // Check each of the various "spokes", stopping when we've gotten
   // back to the beginning
   while( curedg!=firstedg )
   {
      assert( curedg != 0 );
      curedg->CalcLength();
      if ( curedg->CalcSlope() > slp && curedg->FlowAllowed() )
      {
         slp = curedg->getSlope();
         nbredg = curedg;
         fechng = true;
      }
      curedg = curedg->getCCWEdg();
      ++ctr;
      if( unlikely( ctr > kMaxSpokes ) ) // Make sure to prevent std::endless loops
      {
         std::cerr << "Mesh error: node " << getID()
              << " going round and round" << std::endl;
         ReportFatalError( "Bailing out of FlowDirs()" );
      }
   }
   setFlowEdg( nbredg );
   //Nicole changed the line below which is commented out, replaced
   //it with the new if because of the new function WarnSpokeLeaving
   //which can set a node as a boundary node.
   //curnode->setFloodStatus( ( slp>0 ) ? kNotFlooded : kSink );// (NB: opt branch pred?)
   if( slp > 0 && getBoundaryFlag() != kClosedBoundary )
       setFloodStatus( kNotFlooded );
   else
       setFloodStatus( kSink );
   assert( getFlowEdg() != 0 );
   //if( nbredg == firstedg ) return false;
   //return true;
   return fechng;
}

#undef kMaxSpokes

/*****************************************************************************\
 **
 **
 **      DistNew: replaces disnew; the distance of the node's (newx, newy)
 **              from the line
 **              formed by points  p0->(newx, newy) and p1->(newx, newy)
 **      Data members updated:
 **      Called by:
 **      Calls:
 **
 **
\*****************************************************************************/
double tLNode::DistNew(tLNode const * p0,tLNode const * p1 ) const
{
  double a,b,c,res, xp, yp, x0, y0, x1, y1;

  x0 = ( p0->Meanders() ) ? p0->chan.migration.newx : p0->x;
  y0 = ( p0->Meanders() ) ? p0->chan.migration.newy : p0->y;
  x1 = ( p1->Meanders() ) ? p1->chan.migration.newx : p1->x;
  y1 = ( p1->Meanders() ) ? p1->chan.migration.newy : p0->y;
  xp = ( Meanders() ) ? chan.migration.newx : x;
  yp = ( Meanders() ) ? chan.migration.newy : y;

  a = y1 - y0; //(p1->chan.migration.newy) - (p0->chan.migration.newy);
  b = x0 - x1; //-((p1->chan.migration.newx) - (p0->chan.migration.newx));
  c = -( a * x0 + b * y0 );
  //-((a*(p0->chan.migration.newx)) + (b*(p0->chan.migration.newy)));
  res = (a * xp + b * yp + c) / sqrt( a * a + b * b );
  //res = (a*chan.migration.newx + b*chan.migration.newy + c) / sqrt(a*a+b*b);
  if( res < 0 ) res = -res;
  return res;
}



/**************************************************************************\
 **
 **  TellAll
 **
 **  ...is a debugging routine that prints out a bunch of info about a node.
 **
\**************************************************************************/
#ifndef NDEBUG
void tLNode::TellAll() const
{
  tLNode * nbr;

  std::cout << " NODE temp id " << id << " perm id " << permid << ":\n";
  std::cout << "  x=" << x << " y=" << y << " z=" << z;
  if( getEdg() ) {
    std::cout << "  points to edg #" << getEdg()->getID() << std::endl;
    std::cout << "  dr area: " << getDrArea() << "  disch: " << getQ()
	 << "  boundary: " << BoundName(boundary)
	 << "  flood: " << FloodName(flood)
	 << "\n  varea: " << varea << std::endl;

    if( flowedge ) {
      nbr = static_cast<tLNode *>(flowedge->getDestinationPtrNC());
      std::cout << "  Flows along edg " << flowedge->getID() << " to node "
	   << nbr->getID() << " at (" << nbr->getX() << ","
	   << nbr->getY() << "," << nbr->getZ() << ")\n    with vedglen "
	   << flowedge->getVEdgLen() << std::endl;
      flowedge->TellCoords();
      if (0) //DEBUG
	std::cout<<"  ccwedge of flowedge is "<<flowedge->getCCWEdg()->getID()
	    <<" originates at "<<flowedge->getCCWEdg()->getOriginPtrNC()->getID()
	    <<std::endl;
      std::cout << "  qs: " << qs << "  qsin: " << qsin << "  slp: "
	   << flowedge->getSlope() << "  reg: " << reg.thickness << std::endl;
      for(size_t i=0; i<numg; i++)
	std::cout<<"  qsi "<<i<<" "<<qsm[i];
      std::cout<<std::endl;
      //for(i=0; i<numg; i++)
      //std::cout<<"  qsini "<<i<<" "<<qsinm[i];
      //std::cout<<std::endl;
      std::cout<<" creation time top "<<getLayerCtime(0);
      std::cout<<"numlayers is "<<getNumLayer()<<std::endl;
      int j;
      for( j=0; j<getNumLayer(); j++ )
	for( size_t i=0; i<numg; i++)
	  std::cout<<"  dgrade "<<i<<" "<<getLayerDgrade(j,i);
      std::cout << "  dzdt: " << dzdt << "  drdt: " << drdt;
      std::cout<<" meanders "<< Meanders()<<std::endl;
      std::cout << "  width: " << getHydrWidth() << " tau: " << getTau();
      std::cout << " taucrit: " << getTauCrit() << std::endl;
    }
    else std::cout << "  Flowedg is undefined\n";

  }
  else std::cout << "  edg is undefined!\n";

  std::cout << "layerlist addresses are:\n";
  layerlist.DebugTellPtrs();
}
#endif


/**************************************************************************\
 **
 **  EroDep
 **
 **  Erodes or Deposits material of depth dz.
 **
 **  Note: as of now, just changes elev. Should also change alluvial
 **  thickness and set to zero if negative (bedrock erosion).
 **  Alluvial thickness is outdated - just use other EroDep.
 **
 **  Modified, 11/2010: Now calls tLNode::ChangeZ(dz) (which calls 
 **  tNode::ChangeZ(dz)) instead of changing the elevation directly so 
 **  that a new static boolean tNode::freezeElevations is consulted to 
 **  determine whether elevations are meant to change or not.
 **
 **  1/15/98 gt
\**************************************************************************/
void tLNode::EroDep( double dz )
{
//   z += dz;
  ChangeZ( dz );

  cumulative_ero_dep_ += dz;
  if (0) //DEBUG
    std::cout << "  eroding " << id << " by " << dz << std::endl;

  reg.thickness += dz;
  if( reg.thickness < 0. ) reg.thickness = 0.0;
}

void tLNode::InsertLayerBack( tLayer const & lyr )
{
  layerlist.insertAtBack( lyr );
}


/*Xvoid tLNode::setVegErody( double val )
  {surf.vegerody = ( val >= 0.0 ) ? val : 0.0;}

  double tLNode::getVegErody() const {return surf.vegerody;}
*/

void tLNode::setQsinErrorHandler( size_t i ) const
{
  if(i>=numg)
    ReportFatalError( "Trying to index sediment sizes that don't exist ");
  if(i>=qsinm.getSize()){
    std::cout<<"trying to set index "<<i<<" but size of array is "
	<<qsinm.getSize() << std::endl;
    TellAll();
    ReportFatalError( "Index out of bound");
  }
  ReportFatalError("setQsinErrorHandler(): bail out.");
}

void tLNode::addQsin( tArray< double > const &val )
{
  int i;
  const int n = val.getSize();
  for(i=0; i<n; i++){
    qsin +=  val[i];
    qsinm[i] += val[i];
  }

}

void tLNode::setQsdinErrorHandler( size_t i ) const
{
  if(i>=numg)
    ReportFatalError( "Trying to index sediment sizes that don't exist ");
  if(i>=qsdinm.getSize()){
    std::cout<<"trying to set index "<<i<<" but size of array is "
	<<qsdinm.getSize() << std::endl;
    TellAll();
    ReportFatalError( "Index out of bound");
  }
  ReportFatalError("setQsdinErrorHandler(): bail out.");
}

void tLNode::addQsdin( tArray< double > const &val )
{
  int i;
  const int n = val.getSize();
  for(i=0; i<n; i++){
    qsdin +=  val[i];
    qsdinm[i] += val[i];
  }

}

void tLNode::addQs( tArray< double > const & val )
{
  for(size_t i=0; i<val.getSize(); i++){
    qsm[i] += val[i];
    qs += val[i];
  }
}

/************************************************************************\
  tLNode::PrepForAddition

  7/2003 AD after SL
\*************************************************************************/
void tLNode::PrepForAddition( tTriangle const* tri, double time ){
  LayerInterpolation( tri, getX(), getY(), time );
}

/************************************************************************\
  tLNode::PrepForMovement

  7/2003 AD after SL
\*************************************************************************/
void tLNode::PrepForMovement( tTriangle const* tri, double time ){
  const tArray<double> newxy = getNew2DCoords();
  LayerInterpolation( tri, newxy[0], newxy[1], time );
}

/************************************************************************\
  tLNode::LayerInterpolation

  Called by tGrid when new tLNode is added to the grid, or a node is
  moved.
  tri = The old triangle in which the new node will be contained.
        The attributes of its vertices are used for the interpolation.
  tx = The new x location at which the point will be located.
  ty = The new y location at which the point will be located.
  time = current time

  You may wonder why x and y need to be passed.  This is just in case a
  node is moving, in which case you need to reference the newx and newy
  value.  If the node is just being added, only need to reference x and y.
  To make things generic, I just decided to pass in the x and y values.

  This function first creates the new list of layers and then sets the
  layerlist of the node equal to the list which was just created.
  This is done in case a node is moving and its old layers are needed
  for the interpolation.

  LayerInterpolation now checks to make sure that you are interpolating
  only between non-boundary nodes.
  4/1999 Makes sure that top layer depth = maxregdep

  Designed by GT, Coded by NG, last update 4/1999
\*************************************************************************/

#define kAncient -999
void tLNode::LayerInterpolation( tTriangle const* tri, double tx, double ty, double time )
{
  assert(tri!=0);

  if (0) //DEBUG
    std::cout<<std::endl<<"tLNode::LayerInterpolation....";
  if (0) {//DEBUG
    std::cout<<" current x = "<<x<<" current y = "<<y;
    std::cout<<" newx= "<<tx<<" newy= "<<ty<<std::endl;
  }

  tLNode *lnds[3];
  int i;

  int numnodes=0;
  for(i=0; i<=2; i++){
    if(!tri->pPtr(i)->getBoundaryFlag()){
      lnds[numnodes] = static_cast<tLNode *>(tri->pPtr(i));
      numnodes++;
    }
  }
  if (0) //DEBUG
    std::cout<<"numnodes = "<<numnodes<<" newx= "<<tx<<" newy= "<<ty<<std::endl;

  tList< tLayer > helplist; //Make the layer list first.  When
  //the list is made, then set the nodes layerlist equal to helplist.

  if( numnodes==3 ){
    tArray<int> layindex(3);//What layer are you at for each node
    tArray<double> dep(3);//The depth of the layer with the matching age
    tArray<double> age(3);//The age of the current layer at each node
    //The sediment flag of the current layer at each node
    tArray<tLayer::tSed_t> sed(3, tLayer::kBedRock);

    double dist[3]; //distance b/w new point and the three triangle points
    dist[0]=DistanceBW2Points(tx, ty, lnds[0]->getX(), lnds[0]->getY());
    dist[1]=DistanceBW2Points(tx, ty, lnds[1]->getX(), lnds[1]->getY());
    dist[2]=DistanceBW2Points(tx, ty, lnds[2]->getX(), lnds[2]->getY());

    double CA=-1;

    for(i=0; i<=2; i++){
      if(lnds[i]->getLayerRtime(0)>CA){
	CA=lnds[i]->getLayerRtime(0);
      }
      sed[i]=lnds[i]->getLayerSed(0);
      if(sed[i]>0)
	age[i]=lnds[i]->getLayerRtime(0);
      else
	age[i]=kAncient;
      layindex[i]=0;//Initialize layindex
    }

    if (0) //DEBUG
      std::cout<<"Current age is "<<CA<<std::endl;
    //CA now contains the youngest surface layer time of the three nodes.
    //Remember that LayerRtime is the most recent time visited, which
    //implies larger time = younger layer
    double newtex; //texture value of new layer
    double sum;  //helper used to calculate the new texture & erody value
    double newerody; //erodability of new layer
    double newdep; //interpolated depth value
    double newetime; //interpolated exposure time
    tLayer layhelp; //set values in this layer then insert in back of layerlist
    layhelp.setDgradesize(numg);
    layhelp.setCtime(time);

    do
      {
	newtex=0;
	sum=0;
	newerody=0;
	newetime=0;

	for(i=0; i<=2; i++){
	  if(age[i]<=CA+10 && age[i]>=CA-10 && sed[i] != tLayer::kBedRock){
	    //layer is in the window of acceptable ages
	    //set node texture helper to texture of that layer
	    //NOTE - This only works for two sizes right now.
	    //set node depth helper to depth of that layer
	    newtex+=lnds[i]->getLayerDgrade(layindex[i],0)/dist[i];
	    newerody+=lnds[i]->getLayerErody(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
	    newetime+=lnds[i]->getLayerEtime(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
	    sum+=lnds[i]->getLayerDepth(layindex[i])/dist[i];
	    dep[i]=lnds[i]->getLayerDepth(layindex[i]);
	    layhelp.setSed(sed[i]);
	    if( layindex[i] != lnds[i]->getNumLayer()-1 ){
	      //More layers below, continue to search
	      layindex[i]+=1;
	      sed[i]=lnds[i]->getLayerSed(layindex[i]);
	      if(sed[i] != tLayer::kBedRock &&
		 lnds[i]->getLayerRtime(layindex[i])>=0)
		age[i]=lnds[i]->getLayerRtime(layindex[i]);
	      else
		age[i]=kAncient;
	    }
	  }
	  else{
	    dep[i]=0;
	    if (0) //DEBUG
	      std::cout<<"not in correct range, dep set to 0"<<std::endl;
	  }
	}


	//Now find the values of the texture and layer depth in the new
	//location, given what we found out above.  Set everything, and
	//then insert the new layer.

	newdep=PlaneFit(tx, ty, lnds[0]->get2DCoords(), lnds[1]->get2DCoords(),
			lnds[2]->get2DCoords(), dep);
	if(newdep>grade[0]){
	  newtex=newtex/sum;
	  layhelp.setDepth(newdep);
	  layhelp.setDgrade(0,newtex*newdep);
	  layhelp.setErody(newerody/sum);
	  layhelp.setEtime(newetime/sum);
	  layhelp.setRtime(CA);
	  if(numg>1)
	    layhelp.setDgrade(1,(1-newtex)*newdep);
	  helplist.insertAtBack( layhelp );
	}
	CA=age[0];
	for(i=1; i<=2; i++){
	  if(age[i]>CA){
	    CA=age[i];
	  }
	}
	if (0) //DEBUG
	  std::cout<<std::endl<<"after iter, current age is set to "<<CA<<std::endl;

      }while(CA>kAncient);

    //If there is a deep sediment layer, interpolate it
    if(lnds[0]->getLayerRtime(layindex[0])<0){//here assuming if one
      //node has a deep sed layer, they all do
      //TODO fix this!
      newtex=0;
      newerody=0;
      sum=0;
      for(i=0; i<=2; i++){
	//NOTE - This only works for two sizes right now.
	newtex+=lnds[i]->getLayerDgrade(layindex[i],0)/dist[i];
	newerody+=lnds[i]->getLayerErody(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
	sum+=lnds[i]->getLayerDepth(layindex[i])/dist[i];
	dep[i]=lnds[i]->getLayerDepth(layindex[i]);
	layindex[i]++;
      }
      newtex=newtex/sum;
      newdep=PlaneFit(tx, ty, lnds[0]->get2DCoords(),lnds[1]->get2DCoords(),
		      lnds[2]->get2DCoords(), dep);
      layhelp.setDepth(newdep);
      layhelp.setDgrade(0,newtex*newdep);
      if(numg>1)
	layhelp.setDgrade(1,(1-newtex)*newdep);
      layhelp.setSed(tLayer::kSed);
      layhelp.setErody(newerody/sum);
      layhelp.setRtime(-1.);
      layhelp.setEtime(0.);
      helplist.insertAtBack( layhelp );
    }

    newtex=0;
    newerody=0;
    newetime=0;
    sum=0;
    //Now, you need to interpolate bedrock layer here
    for(i=0; i<=2; i++){
      //NOTE - This only works for two sizes right now.
      //debugging routine
      if (0) { //DEBUG
	if(lnds[i]->getNumLayer()<=layindex[i]){
	  lnds[i]->TellAll();
	  for(int j=0; j<lnds[i]->getNumLayer(); j++){
	    std::cout << "layer " << j+1<<std::endl ;
	    std::cout << lnds[i]->getLayerCtime(j);
	    std::cout << " " << lnds[i]->getLayerRtime(j);
	    std::cout << " "<<lnds[i]->getLayerEtime(j)<<std::endl;
	    std::cout << lnds[i]->getLayerDepth(j);
	    std::cout << " " << lnds[i]->getLayerErody(j);
	    std::cout << " " << lnds[i]->getLayerSed(j) << std::endl;
	    std::cout << lnds[i]->getLayerDgrade(j,0);
	  }
	}
      }
      newtex+=lnds[i]->getLayerDgrade(layindex[i],0)/dist[i];
      newerody+=lnds[i]->getLayerErody(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
      newetime+=lnds[i]->getLayerEtime(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
      sum+=lnds[i]->getLayerDepth(layindex[i])/dist[i];
      dep[i]=lnds[i]->getLayerDepth(layindex[i]);
    }
    newtex=newtex/sum;
    newdep=PlaneFit(tx, ty, lnds[0]->get2DCoords(),lnds[1]->get2DCoords(),
		    lnds[2]->get2DCoords(), dep);
    layhelp.setDepth(newdep);
    layhelp.setDgrade(0,newtex*newdep);
    if(numg>1)
      layhelp.setDgrade(1,(1-newtex)*newdep);
    layhelp.setSed(tLayer::kBedRock);
    layhelp.setErody(newerody/sum);
    layhelp.setRtime(0.);
    layhelp.setEtime(newetime/sum);

    helplist.insertAtBack( layhelp );

    //This is the new part that nmg added where the top layers are
    //truncated according to the new elevation and the elevation
    //around the node
    const tArray<double> elevs(
			       lnds[0]->getZ(),
			       lnds[1]->getZ(),
			       lnds[2]->getZ()
			       );
    const double theoryz =
      PlaneFit(tx, ty, lnds[0]->get2DCoords(), lnds[1]->get2DCoords(),
	       lnds[2]->get2DCoords(), elevs );
    double diff=theoryz-getZ();

    /*std::cout<<"new z "<<getZ()<<" theoretical z "<<theoryz<<" diff "<<diff<<std::endl;
      std::cout<<"before removal, layer info"<<std::endl;
      std::cout<<"numlayers = "<<helplist.getSize()<<std::endl;
      std::cout<<"depth of first layer = "<<helplist.FirstP()->getDepth()<<std::endl;*/

    while(diff>0){
      if(helplist.FirstP()->getDepth()<diff){
	diff-=helplist.FirstP()->getDepth();
	helplist.removeFromFront(layhelp);
      }
      else{
	double dpth = helplist.FirstP()->getDepth();
	tLayer * firstlay = helplist.FirstP();
	firstlay->setDepth(dpth-diff);
	diff=0;
      }
    }

  }
  else if( numnodes ==2 ){
    tArray<int> layindex(2);//What layer are you at for each node
    tArray<double> age(2);//The age of the current layer at each node
    //The sediment flag of the current layer at each node
    tArray<tLayer::tSed_t> sed(2, tLayer::kBedRock);

    double dist[2]; //distance b/w new point and the three triangle points
    dist[0]=DistanceBW2Points(tx, ty, lnds[0]->getX(), lnds[0]->getY());
    dist[1]=DistanceBW2Points(tx, ty, lnds[1]->getX(), lnds[1]->getY());
    double distsum=(1/dist[0])+(1/dist[1]);

    double CA=-1;

    for(i=0; i<=1; i++){
      if(lnds[i]->getLayerRtime(0)>CA){
	CA=lnds[i]->getLayerRtime(0);
      }
      sed[i]=lnds[i]->getLayerSed(0);
      if(sed[i] != tLayer::kBedRock)
	age[i]=lnds[i]->getLayerRtime(0);
      else
	age[i]=kAncient;
      layindex[i]=0;//Initialize layindex
    }

    //CA now contains the youngest surface layer time of the three nodes.
    //Remember that LayerRtime is the most recent time visited, which
    //implies larger time = younger layer
    double newtex; //texture value of new layer
    double sum;  //helper used to calculate the new texture & erody value
    double newerody; //erodability of new layer
    double newetime; //exposure time of new layer
    double newdep; //interpolated depth value
    tLayer layhelp; //set values in this layer then insert in back of layerlist
    layhelp.setDgradesize(numg);
    layhelp.setCtime(time);

    do
      {
	newtex=0;
	sum=0;
	newerody=0;
	newetime=0;
	newdep=0;

	for(i=0; i<=1; i++){
	  if(age[i]<=CA+10 && age[i]>=CA-10 && sed[i] != tLayer::kBedRock){
	    //layer is in the window of acceptable ages
	    //set node texture helper to texture of that layer
	    //NOTE - This only works for two sizes right now.
	    newtex+=lnds[i]->getLayerDgrade(layindex[i],0)/dist[i];
	    newerody+=lnds[i]->getLayerErody(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
	    newetime+=lnds[i]->getLayerEtime(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
	    sum+=lnds[i]->getLayerDepth(layindex[i])/dist[i];
	    newdep+=lnds[i]->getLayerDepth(layindex[i])/dist[i]/distsum;

	    layhelp.setSed(sed[i]);
	    if( layindex[i] != lnds[i]->getNumLayer()-1 ){
	      //More layers below, continue to search
	      layindex[i]+=1;
	      sed[i]=lnds[i]->getLayerSed(layindex[i]);
	      if(sed[i]>0 && lnds[i]->getLayerRtime(layindex[i])>=0)
		age[i]=lnds[i]->getLayerRtime(layindex[i]);
	      else
		age[i]=kAncient;
	    }
	  }
	}

	//Now find the values of the texture and layer depth in the new
	//location, given what we found out above.  Set everything, and
	//then insert the new layer.


	if(newdep>grade[0]){
	  newtex=newtex/sum;
	  layhelp.setDepth(newdep);
	  layhelp.setDgrade(0,newtex*newdep);
	  layhelp.setErody(newerody/sum);
	  layhelp.setEtime(newetime/sum);
	  layhelp.setRtime(CA);
	  if(numg>1)
	    layhelp.setDgrade(1,(1-newtex)*newdep);
	  helplist.insertAtBack( layhelp );
	}
	CA=age[0];
	if(age[1]>CA){
	  CA=age[1];
	}

      }while(CA>kAncient);

    //Now, you need to interpolate deep sediment layer (if it exists)
    if(lnds[0]->getLayerRtime(layindex[0])<0){//Just assume that if one
      //node has the deep sediment layer that they all do.
      //TODO fix this!
      newtex=0;
      newerody=0;
      sum=0;
      newdep=0;
      for(i=0; i<=1; i++){
	//NOTE - This only works for two sizes right now.
	newtex+=lnds[i]->getLayerDgrade(layindex[i],0)/dist[i];
	newerody+=lnds[i]->getLayerErody(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
	sum+=lnds[i]->getLayerDepth(layindex[i])/dist[i];
	newdep+=lnds[i]->getLayerDepth(layindex[i])/dist[i]/distsum;
	layindex[i]++;
      }
      newtex=newtex/sum;
      layhelp.setDepth(newdep);
      layhelp.setDgrade(0,newtex*newdep);
      if(numg>1)
	layhelp.setDgrade(1,(1-newtex)*newdep);
      layhelp.setSed(tLayer::kSed);
      layhelp.setErody(newerody/sum);
      layhelp.setRtime(-1.);
      layhelp.setEtime(0.);
      helplist.insertAtBack( layhelp );
    }

    //Finally, interpolate bedrock
    newtex=0;
    newerody=0;
    newetime=0;
    sum=0;
    newdep=0;
    for(i=0; i<=1; i++){
      //NOTE - This only works for two sizes right now.
      newtex+=lnds[i]->getLayerDgrade(layindex[i],0)/dist[i];
      newerody+=lnds[i]->getLayerErody(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
      newetime+=lnds[i]->getLayerEtime(layindex[i])*lnds[i]->getLayerDepth(layindex[i])/dist[i];
      sum+=lnds[i]->getLayerDepth(layindex[i])/dist[i];
      newdep+=lnds[i]->getLayerDepth(layindex[i])/dist[i]/distsum;
    }
    newtex=newtex/sum;
    layhelp.setDepth(newdep);
    layhelp.setDgrade(0,newtex*newdep);
    if(numg>1)
      layhelp.setDgrade(1,(1-newtex)*newdep);
    layhelp.setSed(tLayer::kSed);
    layhelp.setErody(newerody/sum);
    layhelp.setRtime(0.);
    layhelp.setEtime(newetime/sum);
    helplist.insertAtBack( layhelp );


    //This is the new part that nmg added where the top layers are
    //truncated according to the new elevation and the elevation
    //around the node
    const double theoryz=LineFit(0., lnds[0]->getZ(), dist[0]+dist[1],
				 lnds[1]->getZ(), dist[0]);
    double diff=theoryz-getZ();
    while(diff>0){
      if(helplist.FirstP()->getDepth()<diff){
	diff-=helplist.FirstP()->getDepth();
	helplist.removeFromFront(layhelp);
      }
      else{
	double dpth = helplist.FirstP()->getDepth();
	tLayer * firstlay = helplist.FirstP();
	firstlay->setDepth(dpth-diff);
	diff=0;
      }
    }

  }
  else{
    //only one case left - only one non-boundary node
    //check done w/assertion at beginning
    //just set layerlist = to layerlist of non-boundary node
    tLayer layhelp; //set values in this layer then insert in back of layerlist
    layhelp.setDgradesize(numg);
    layhelp.setCtime(time);

    for(i=0; i<lnds[0]->getNumLayer(); i++){
      layhelp.setDepth(lnds[0]->getLayerDepth(i));
      layhelp.setDgrade(0,lnds[0]->getLayerDgrade(i,0));
      if(numg>1)
	layhelp.setDgrade(1,lnds[0]->getLayerDgrade(i,1));
      layhelp.setSed(lnds[0]->getLayerSed(i));
      layhelp.setErody(lnds[0]->getLayerErody(i));
      layhelp.setEtime(lnds[0]->getLayerEtime(i));
      layhelp.setRtime(lnds[0]->getLayerRtime(i));
      helplist.insertAtBack( layhelp );
    }

    //This is the new part that nmg added where the top layers are
    //truncated according to the new elevation and the elevation
    //around the node
    double diff=lnds[0]->getZ()-getZ();
    while(diff>0){
      if(helplist.FirstP()->getDepth()<diff){
	diff-=helplist.FirstP()->getDepth();
	helplist.removeFromFront(layhelp);
      }
      else{
	double dpth = helplist.FirstP()->getDepth();
	tLayer * firstlay = helplist.FirstP();
	firstlay->setDepth(dpth-diff);
	diff=0;
      }
    }
  }

  if(maxregdep-helplist.FirstP()->getDepth()>0.001 &&
     helplist.getIthDataRef(1).getSed()>0){
    //top layer is too small, want to maintain maxregdep for erosion reasons
    //Because NG doesn't really know how to manipulate lists,
    //this is probably done in a cockeyed way.  First remove
    //top layer from list, update it, then add it back to front of list.
    tLayer firstlay;
    helplist.removeFromFront(firstlay);
    double diff=maxregdep-firstlay.getDepth();
    double newerody=firstlay.getErody()*firstlay.getDepth();
    double newetime=firstlay.getEtime()*firstlay.getDepth();
    double newrtime=firstlay.getRtime()*firstlay.getDepth();
    do{
      tLayer * nextlay=helplist.FirstP();
      if(nextlay->getDepth()<diff){
	//add entire contents of layer below to top layer
	for(size_t i=0; i<numg; i++)
	  firstlay.addDgrade(i,nextlay->getDgrade(i));
	newerody+=nextlay->getErody()*nextlay->getDepth();
	newetime+=nextlay->getEtime()*nextlay->getDepth();
	newrtime+=nextlay->getRtime()*nextlay->getDepth();
	diff-=nextlay->getDepth();
	helplist.removeFromFront(*nextlay);
      }
      else{
	//don't remove entire layer below, only take some material
	double prevdep = nextlay->getDepth();
	for(size_t i=0; i<numg; i++){
	  double moving = diff*nextlay->getDgrade(i)/prevdep;
	  firstlay.addDgrade(i,moving);
	  nextlay->addDgrade(i,-1*moving);
	}
	newerody+=nextlay->getErody()*diff;
	newetime+=nextlay->getEtime()*diff;
	//don't want to include negative exposure time flag
	if(nextlay->getRtime()>0)
	  newrtime+=nextlay->getRtime()*diff;
	diff=0;
      }
    }while(diff>0.001);
    firstlay.setErody(newerody/firstlay.getDepth());
    firstlay.setEtime(newetime/firstlay.getDepth());
    firstlay.setRtime(newrtime/firstlay.getDepth());
    helplist.insertAtFront(firstlay);
  }
  else if(helplist.FirstP()->getRtime()<0){
    //top layer is the deep sediment layer.  This messes up things
    //for layering in general.  Form a new top layer with material from
    //the bottom layer.
    tLayer layhelp;
    layhelp.setSed(tLayer::kSed);
    layhelp.setEtime(0.);
    layhelp.setCtime(time);
    layhelp.setRtime(0.);
    layhelp.setDgrade(0,maxregdep*helplist.FirstP()->getDgrade(0)/helplist.FirstP()->getDepth());
    if(numg>1)
      layhelp.setDgrade(1,maxregdep*helplist.FirstP()->getDgrade(1)/helplist.FirstP()->getDepth());
    layhelp.setErody(helplist.FirstP()->getErody());
    helplist.FirstP()->setDepth( helplist.FirstP()->getDepth()-maxregdep);
    helplist.insertAtFront( layhelp );

  }

  layerlist=helplist;
  //     if(getNumLayer()<=2)
  //         std::cout<<"layerinterp made 2 layers at "<<x<<", "<<y<<std::endl;

  //Below  if for debugging purposes
  if(0) { //DEBUG
    /*if(getLayerEtime(0)<0 || (tx>505.0 && tx<506.0 && ty>331.0 && ty<332.0)  ){
     */
    i=0;
    std::cout<<"x "<<tx<<" y "<<ty;
    while(i<layerlist.getSize()){
      std::cout << "layer " << i+1 <<std::endl;
      std::cout << getLayerCtime(i);
      std::cout << " " << getLayerRtime(i);
      std::cout << " " << getLayerEtime(i)<<std::endl;
      std::cout << getLayerDepth(i);
      std::cout << " " << getLayerErody(i);
      std::cout << " " << getLayerSed(i) << std::endl;
      std::cout << getLayerDgrade(i,0) ;
      if( numg>1 ) std::cout << " " << getLayerDgrade(i,1);
      i++;
      std::cout<<std::endl;
    }

    int j;
    for(j=0; j<numnodes; j++){
      std::cout<<std::endl;
      std::cout<<"node "<<j<<" x "<<lnds[j]->getX()<<" y "<<lnds[j]->getY();
      std::cout<<" boundary "<< BoundName(lnds[j]->getBoundaryFlag())<<std::endl;
      tLNode * nicn = lnds[j];
      i=0;
      while(i<nicn->getNumLayer()){
	std::cout << "layer " << i+1<<std::endl ;
	std::cout << nicn->getLayerCtime(i);
	std::cout << " " << nicn->getLayerRtime(i);
	std::cout << " "<<nicn->getLayerEtime(i)<<std::endl;
	std::cout << nicn->getLayerDepth(i);
	std::cout << " " << nicn->getLayerErody(i);
	std::cout << " " << nicn->getLayerSed(i) << std::endl;
	std::cout << nicn->getLayerDgrade(i,0);
	if( numg>1 ) std::cout << " " << nicn->getLayerDgrade(i,1) << std::endl;
	i++;
      }
    }
  }
}
#undef kAncient


/*************************************************************\
  tLNode::WarnSpokeLeaving(tEdge * edglvingptr)

  This function is called when an edge is being removed from the edge list.
  If flowedge is pointing to the edge
  which will be removed, this flowedge must be updated.

  edglvingptr is, as it sais, a pointer to the edge which will be
  removed.

  Called from tGrid::ExtricateEdge

  9/98 NG and GT
\*******************************************************************/
void tLNode::WarnSpokeLeaving(tEdge * edglvingptr)
{
  if (0) //DEBUG
    std::cout<<"tLNode::WarnSpokeLeaving..... node #"<<id<<std::endl;

  //Make sure that edg pointer in tNode won't be affected
  tNode::WarnSpokeLeaving( edglvingptr );

  if( edglvingptr == flowedge ){
    do{
      flowedge = flowedge->getCCWEdg();
    }while( (flowedge->getBoundaryFlag()==kClosedBoundary)
	    && (flowedge != edglvingptr) );
    assert( flowedge->getOriginPtr() == this );

    //After looping around edges, if flow is along a non-flow edge,
    //make this a closedboundary node.  This is potentially a bad
    //thing to do since a node will be marked as a boundary node
    //but not put at the end of the node list.  For now we ignore this
    //TODO fix this - probably in tGrid::ExtricateEdge
    //        if(flowedge->getBoundaryFlag()==kClosedBoundary){
    //           boundary = kClosedBoundary;
    //           std::cout<<"node "<<getID()<<" x "<<getX()<<" y "<<getY()<<" set to boundary in WarnSpokeLeaving"<<std::endl<<flush;
    //        }
  }
}

/**********************************************************************\
 **
 ** tLNode::InitializeNode()
 **
 ** A virtual function.  Used when nodes are creating during the evolution
 ** process, not at the start of a run.
 ** - Assigns one of the spokes as the flow edge.
 ** - sets size of Qsi and Qsini in case not properly set.
 **
 **
 ** Called from tGrid::AddNode
 **
 ** Created 1/1999 NG
 ***********************************************************************/
void tLNode::InitializeNode()
{
  // If we're not a boundary node, make sure we have a valid flowedge
  if( boundary==kNonBoundary )
    {
       // SL 9/2003: call new FindInitFlowDir():
       FindInitFlowDir();
    }

  // Set size of sediment influx and outflux arrays to number of grain sizes
  if( qsm.getSize()!=numg ) {
    qsm.setSize( numg );
    qsinm.setSize( numg );
    qsdinm.setSize( numg );
  }

  accumdh.setSize(2);

  if (0) //DEBUG
    std::cout<<"tLNode::InitializeNode node "<< getID()
	<<" flow edge "<<flowedge->getID()<<std::endl;
}

/*******************************************************************\
  tLNode::FuturePosn() virtual function

  11/98 SL
  5/2003 AD
\*******************************************************************/
tArray< double > tLNode::FuturePosn() {
  //return Meanders() ? getNew2DCoords() : get2DCoords();
  return isMobile() ? getNew2DCoords() : get2DCoords();
}

/***********************************************************************\
  tLNode::splitFlowEdge
  Create a node in the center of the flow edge.
  This node is return dynamically allocated and needs to deleted.
  Its flowedge is set to zero.
  09/2003 AD and QC
\***********************************************************************/
tNode *tLNode::splitFlowEdge() {
  // split only non flippable flow edge
  assert(!getFlowEdg()->isFlippable());

  const tArray< double > zeroArr(4);
  tLNode *nPtr = this->getDownstrmNbr();

  if (0) //DEBUG
    std::cout << "tLNode::splitFlowEdge(): split flowedge between node "
	 << this->getID() << " and node "
	 << nPtr->getID() << "." << std::endl;

  tLNode *nn = new tLNode(*this);
  nn->setXYZD( zeroArr ); //and xyzd ('old' coords)

  // set boundary condition
  const tBoundary_t b0 = this->getBoundaryFlag();
  const tBoundary_t b1 = nPtr->getBoundaryFlag();
  if (b0 != b1){
    if (b0 == kNonBoundary || b1 == kNonBoundary)
      nn->setBoundaryFlag(kNonBoundary);
    else if (b0 == kClosedBoundary || b1 == kClosedBoundary)
      nn->setBoundaryFlag(kClosedBoundary);
    else {
      assert(0);
      abort();
    }
  }

  // set coordinates by averaging
  nn->setX( (this->getX()+nPtr->getX())/2 );
  nn->setY( (this->getY()+nPtr->getY())/2 );
  nn->setZ( (this->getZ()+nPtr->getZ())/2 );

  // set flow edge temporarly to zero, so that it is flippable
  setFlowEdgToZero();

  // return a pointer holding the new node
  return nn;
}

/***********************************************************************\
  tLNode::getUpstrmNbr
  SL, 10/03:
\***********************************************************************/
tLNode* tLNode::getUpstrmNbr()
{
   tSpkIter spI( this );
   double maxarea = 0.0;
   tLNode* mn = NULL;
   for( tEdge* ce = spI.FirstP(); !spI.AtEnd(); ce = spI.NextP() )
   {
      if( ce != flowedge && ce->FlowAllowed() )
      {
         tLNode* un = static_cast<tLNode*>( ce->getDestinationPtrNC() );
         if( un->getDownstrmNbr() == this )
         {
            double curarea = un->getDrArea();
            if( curarea > maxarea )
            {
               mn = un;
               maxarea = curarea;
            }
         }
      }
   }
   return mn;
}

/***********************************************************************\
  tArray<double> tLNode::EroDep

  This function erodes and deposits material and updates
  the layering at the same time.  When eroding, layers are updated
  to maintain the maximum layer depth as long as there is sediment
  to refill them.  Bedrock layers can never receive material.

  Takes : - i is layer to erode/dep from/into
          - valgrd contains the depth of each grain size to erode
          (negative value) or deposit (positive value)
          - tt is the current time - for use if depositing

  Returns : an array containing the actual depth of each grain
            size that was actually ero'd/dep'd.  Only an issue
            if eroding because you might be limited in what you
            can erode (stuff might not be there, only erode from
            one layer at a time).

  Calls : - addtoLayer (both versions)-  for updating the layers below,
          and (if eroding) finding out what to add to the layer you
          are eroding from.
          - removeLayer - sometimes wierd things happens and removing
          everything doesn't really remove everything.  (Why????)
          Don't want really small layers, so remove them.
          - makeNewLayerBelow - have a maximum layer depth.  If you
          are depositing into a layer and it will surpass maxlayerdepth,
          need to move stuff into layer below.  If there is no room
          below, just make a new layer.

  Created 6/98 ng
 
  Modified, 11/2010: Now calls tLNode::ChangeZ(dz) (which calls 
   tNode::ChangeZ(dz)) instead of changing the elevation directly so 
   that a new static boolean tNode::freezeElevations is consulted to 
   determine whether elevations are meant to change or not.

\***********************************************************************/


tArray<double> tLNode::EroDep( int i, tArray<double> valgrd, double tt)
{
  double amt, val, olddep;
  tArray<double> update(numg);
  tArray<double> hupdate(numg);
  
  //NIC these are for testing
  //Xbefore=getLayerDepth(i);
  // int numlay=getNumLayer();
  //Xint h;
  //int nh = 0;
  
  //if(getLayerDgrade(i,0)/getLayerDepth(i)>0.99)
  //{
  //if(x<560.418 && x>560.416){
  // nh=1;
  // std::cout<<"layer "<<i<<" size 0 "<< getLayerDgrade(i,0)<<" size 1 "<<getLayerDgrade(i,1)<<std::endl;
  // if(getNumLayer()==2){
  //    std::cout<<"ERODEP x "<<x<<" y "<<y<<" numlayers "<<getNumLayer();
  //    std::cout<<" to erode b4 "<<valgrd[0]<< " " <<valgrd[1]<<std::endl;
  // }
  //}
    
  size_t g=0;
  val=0;
  double max, min;
  max = -10000.;
  min = 10000.;
  while(g<numg){
    if(-1*valgrd[g]>getLayerDgrade(i,g))
      valgrd[g]=-1*getLayerDgrade(i,g);
    // Checking to see that there is enough stuff
    if(valgrd[g]>max) max = valgrd[g];
    if(valgrd[g]<min) min = valgrd[g];
    val+=valgrd[g];
    g++;
  }
  
  // val is now set to the total amount of erosion or deposition
//   z += val;
  ChangeZ( val );
  // total elevation has now been changed
  cumulative_ero_dep_ += val;
  // elevation change has been added to cumulative (for coupling with other models)
  
  double sume, sumd, olde; //used for calculation of exposure time
  
  // Branch according to whether there is:
  // (1) erosion of all size-classes
  // (2) deposition of all size-classes
  // (3) erosion of some, deposition of others
  if(max <= 0 && min < -0.0000000001)
  {
    // CASE OF EROSION IN ALL SIZE CLASSES
    if(getLayerSed(i) != tLayer::kBedRock &&
	     getLayerSed(i) == getLayerSed(i+1) && getLayerDepth(i)+val<=maxregdep)
	  {
	     // Updating will also be done if entering this statement
	     // No updating of Bedrock layers.
	     // Only update if layer below has the same material
	     sume=(getLayerDepth(i)+val)*getLayerEtime(i); //for averaging
	                                                   //of the exposure time.
	     olde=getLayerEtime(i+1);
	     while(val<-0.000000001 && getLayerSed(i) == getLayerSed(i+1))
	     {
	        // keep eroding until you either get all the material you
	        // need to refill the top layer, or you run out of material
	        hupdate = addtoLayer(i+1, val);//remove stuff from lower layer
	        size_t g=0;
	        sumd=0;
	        while(g<numg)
          {
            sumd-=hupdate[g];
            val-=hupdate[g];//hupdate stores texture of material that will
                            //refil the top layer
              update[g] += hupdate[g];//need update cause might need
                                      //to get material from more than one layer
                g++;
          }
	        sume+=sumd*olde;
       }
       size_t g=0;
       while(g<numg){
         addtoLayer(i, g, valgrd[g], -1.); // Erosion
         addtoLayer(i,g,-1*update[g],-1.);//Updating with material from below
           g++;
       }
       setLayerEtime(i, sume/getLayerDepth(i));
	  }
    else
	  {
	     // No updating, just eroding, don't need to change exposure time
	     //Do this if you have only bedrock below, or if layer you are eroding
	     //from is >maxregdepth-val (could be if lots of deposition)
	     size_t g=0;
	     while(g<numg){
         addtoLayer(i, g, valgrd[g], -1.); // Erosion done on this line
         g++;
	     }
    }
    if(getLayerDepth(i)<1e-7)
	     removeLayer(i);
    assert( getLayerDepth(i)>0.0 );
  }
  else if(min >= 0.0 && max > 0.0000000001)
  {
    // CASE OF DEPOSITION IN ALL SIZE CLASSES
    // method seems to make good sense for surface layers
    // but may not be as appropriate for lower layers.
    // You may want to either make a provision for lower layers
    // or else write a new algorithm for them (will you ever deposit
    // into lower layers?  I don't know)
    // Also, no test done to make sure that you are depositing the right
    // material into the layer, would need to pass the flag for this.
    // For now assume that the test will be done in another place.
    if(getLayerSed(i) != tLayer::kBedRock && val < maxregdep){
      // top layer is sediment, so no issues
      // depositing less than an entire layer of stuff
      olde=getLayerEtime(i);
      if(getLayerDepth(i)+val>maxregdep){
        // Need to move stuff out of top layer to make room for deposited material
        if(getLayerSed(i) == getLayerSed(i+1) && getLayerDepth(i+1)+val<maxregdep)
        {
          // The layer below is of the appropriate material and has space
          amt = getLayerDepth(i)+val-maxregdep;//how much to move out
          setLayerEtime(i+1, (1/(getLayerDepth(i+1)+amt))*
                        (olde*amt + getLayerDepth(i+1)*getLayerEtime(i+1)));//lower layer's etime is now properly set.
            setLayerEtime(i, (1/maxregdep)*(olde*(getLayerDepth(i)-amt)));
            //now etime is set in the layer you are depositing into
            //note deposited material has an etime of 0
            olddep = getLayerDepth(i);
            size_t g=0;
            while(g<numg){
              addtoLayer(i+1,g,amt*getLayerDgrade(i,g)/olddep, -1.);
              // putting material from top layer to layer below
              // nic, at this point you have decided not to change
              // the recent time on the lower layer when you move
              // stuff down into it.  Might think about this.
              addtoLayer(i,g,-1*amt*getLayerDgrade(i,g)/olddep+valgrd[g], tt);
              // changing top layer composition by depositing and removing
              g++;
            }
            assert( getLayerDepth(i)>0.0 );
            assert( getLayerDepth(i+1)>0.0 );
        }
        else
        {
          // Need to create new layer
          amt = getLayerDepth(i)+val-maxregdep;
          assert( amt>=0.0 ); //GT
          setLayerEtime(i, (1/maxregdep)*(olde*(getLayerDepth(i)-amt)));
          //now etime is set in the layer you are depositing into
          //note deposited material has an etime of 0
          olddep = getLayerDepth(i);
          assert( olddep>0.0 ); // if not true, we get div by zero
          size_t g=0;
          while(g<numg){
            update[g]=amt*getLayerDgrade(i,g)/olddep;
            // material which will be moved from top layer
            assert( update[g]>=0.0 ); //GT
            addtoLayer(i,g,-1*amt*getLayerDgrade(i,g)/olddep+valgrd[g], tt);
            // changing top layer composition by depositing
            g++;
          }
          assert( getLayerDepth(i)>0.0 );
          makeNewLayerBelow(i, getLayerSed(i), getLayerErody(i), update, tt,
                            getLayerBulkDensity(i) );
          
          //When new layer is created then you change the recent time.
          // put material into a new layer which is made below
          setLayerEtime(i+1, olde); //now layer below has the proper etime
        }
      }
      
      else
      {
        // No need to move stuff out of top layer, just deposit
        setLayerEtime(i, (1/(getLayerDepth(i)+val))*olde*getLayerDepth(i));
        //set top layers etime
        size_t g=0;
        while(g<numg){
          addtoLayer(i,g,valgrd[g],tt);
          g++;
        }
        assert( getLayerDepth(i)>0.0 );
      }
    }
    else{
      // depositing sediment on top of br - create a new surface layer baby.
      // setting all new sediment layers to have erodibility of KRnew,
      // value read in at begining
      // Also use this if the amount of material deposited is
      // greater than maxregdep
      makeNewLayerBelow(-1, tLayer::kSed, KRnew, valgrd, tt, 
                        new_sed_bulk_density_ );
    }
  }
  else if(max>0.0000000001 && min<-0.0000000001)
  {
    // CASE OFSIMULTANEOUS EROSION AND DEPOSITION:
    // Erosion of one or more sizes, deposition of other(s).
    // First, test whether the net effect is erosion or deposition.
    //Need if for the zero erosion case, in which nothing is done.
    if(val < 0)
    {
      // Case of net erosion
      if(getLayerSed(i) != tLayer::kBedRock)
	    {
	      //layer is sediment
	      if(getLayerSed(i) == getLayerSed(i+1))
        {
          sume=(getLayerDepth(i)+val)*getLayerEtime(i); //for averaging
                                                        //of the exposure time.
          olde=getLayerEtime(i+1);
          //There is material below to update with
          while(val<-0.000000001 && getLayerSed(i) == getLayerSed(i+1))
          {
            // keep getting material from below  until you
            // either get all the material you
            // need to refill the top layer, or you run out of material
            hupdate = addtoLayer(i+1, val);//remove stuff from
                                           //lower layer, hupdate stores texture of material that will
                                           //refil the top layer
            size_t g=0;
            sumd=0;
            while(g<numg)
            {
              sumd-=hupdate[g];
              val-=hupdate[g];
              update[g] += hupdate[g];
              g++;
            }
            sume+=sumd*olde;
          }
          size_t g=0;
          while(g<numg){
            addtoLayer(i, g, valgrd[g], tt); // Erosion and deposition
            addtoLayer(i,g,-1*update[g], tt);//Updating with material from below
                                             //Set layer recent time because some deposition was done
              g++;
          }
          assert( getLayerDepth(i)>0.0 );
          setLayerEtime(i, sume/getLayerDepth(i));
        }
	      else
        {
          //No updating will be done, but can put the
          //deposited material into the surface layer
          size_t g=0;
          sumd=0;
          while(g<numg){
            if(valgrd[g]>0)
              sumd+=valgrd[g];
            addtoLayer(i, g, valgrd[g], tt); // Erosion/Deposition
                                             //Set layer recent time because some deposition was done
            g++;
          }
          assert( getLayerDepth(i)>0.0 );
          setLayerEtime(i, getLayerEtime(i)*(getLayerDepth(i)-sumd)/
                        getLayerDepth(i));
        }
	    }
      else
	    {
	      //Layer is bedrock
	      //First remove material from bedrock, then create a new layer
	      //for the deposited material.
	      for(size_t g=0; g<numg; g++){
          update[g]=valgrd[g];//update stores the composition of new layer
          if(valgrd[g]<0){
            addtoLayer(i, g, valgrd[g], -1.);
            update[g]=0.0;
          }
	      }
	      assert( getLayerDepth(i)>0.0 );
	      makeNewLayerBelow(i-1, tLayer::kSed, KRnew, update, tt, 
                          new_sed_bulk_density_);
	      //New layer made with deposited material
	    }
      assert( getLayerDepth(i) > -1e-7 ); // can't be much < 0
      if(getLayerDepth(i)<1e-7) // If we've eroded through a layer
        removeLayer(i);
    }
    else
    {
      // Case of net deposition (simultaneous ero & dep)
      // First, we loop over all size-fractions. If erosion in a given
      // fraction is occurring, remove the material from the top layer;
      // if deposition, add it to val, which at the end of the loop
      // will record the total depth to be deposited.
      val=0; //val will now contain total amt to deposit
      for(size_t g=0; g<numg; g++) {
        if(valgrd[g]<0) {
          // Here, we erode material in size fraction g
          update[g]=0;//If new layer needs to be made on top of BR
                      //update will be composition.
                      //Not necessarily used unles on bedrock, but I threw it
                      //here anyway since you one should rarely enter this
                      //neck of the woods.
          addtoLayer(i, g, valgrd[g], -1.);
        }
        else {
          // Size fraction g is going to be deposited, so add it to the
          // total amount deposited (val)
          update[g]=valgrd[g];
          val+=valgrd[g];
        }
      }
      // If we've eroded all of layer, remove it (avoid division by zero
      // below) GT 8/02
      if( getLayerDepth(i)<1e-7 ) removeLayer(i);
      assert( getLayerDepth(i)>0.0 );  // replacement should be ok
      if(getLayerSed(i) != tLayer::kBedRock && val<maxregdep) {
        // Case in which top layer is "sediment" and depth to be
        // deposited is less than nominal layer thickness.
        // top layer is sediment, so no issues
        amt = (getLayerDepth(i)+val)-maxregdep; // excess (if pos)
                                                //if((getLayerDepth(i)+val)>maxregdep){
        if( amt>0.0 ){
          // Need to move stuff out of top layer to make room for deposited mat
          //if(getLayerSed(i) == getLayerSed(i+1) && getLayerDepth(i+1)+val<maxregdep)
          if(getLayerSed(i) == getLayerSed(i+1) && (getLayerDepth(i+1)+amt)<maxregdep)
          {
            // The layer below is of the appropriate material and has space
            //amt = getLayerDepth(i)+val-maxregdep;//how much to move out
            olddep = getLayerDepth(i);
            olde=getLayerEtime(i);
            setLayerEtime(i, (1/maxregdep)*olde*(getLayerDepth(i)-amt));
            //exposure time of layer which is getting ero'd/dep'd is set
            setLayerEtime(i+1, (1/(amt+getLayerDepth(i+1)))*
                          (amt*olde+getLayerDepth(i+1)*getLayerEtime(i+1)));//now exposure time of lower layer is set
              size_t g=0;
              while(g<numg){
                addtoLayer(i+1,g,amt*getLayerDgrade(i,g)/olddep, -1.);
                // putting material from top layer to layer below
                // nic, at this point you have decided not to change
                // the recent time on the lower layer when you move
                // stuff down into it.  Might think about this.
                //addtoLayer(i,g,-1*amt*getLayerDgrade(i,g)/olddep+valgrd[g], tt);
                addtoLayer(i,g,-1*amt*getLayerDgrade(i,g)/olddep+update[g], tt);
                // changing top layer composition by depositing and removing
                g++;
              }
              assert( getLayerDepth(i)>0.0 );
          }
          else
          {
            // Need to create new layer
            //amt = getLayerDepth(i)+val-maxregdep;
            olddep = getLayerDepth(i);
            assert( olddep>0.0 ); // otherwise div by zero below
            olde=getLayerEtime(i);
            setLayerEtime(i, (1/maxregdep)*olde*(getLayerDepth(i)-amt));
            size_t g=0;
            while(g<numg){
              update[g]=amt*(getLayerDgrade(i,g)/olddep);
              // update = material which will be moved from top layer
              assert( update[g]>=0.0 );
              // GT added the following Aug 2002 -- previously, erosion
              // was being double-counted
              if( valgrd[g]>0.0 ) // If depositing in this size
                addtoLayer(i,g,valgrd[g]-update[g], tt);
              // changing top layer composition by depositing
              else // eroding in this size: erosion already done above
                addtoLayer(i,g,-update[g],tt);
              g++;
            }
            assert( getLayerDepth(i)>0.0 );
            makeNewLayerBelow(i, getLayerSed(i), getLayerErody(i), update, tt,
                              getLayerBulkDensity(i) );
            setLayerEtime(i+1, olde);
            //When new layer is created then you change the time.
            // put material into a new layer which is made below
          }
        }
        else
	      {
          // No need to move stuff out of top layer, just deposit
          setLayerEtime(i, (1/(getLayerDepth(i)+val))*getLayerDepth(i)*
                        getLayerEtime(i));
          //now exposure time of top layer is properly set
          size_t g=0;
          while(g<numg){
            if(valgrd[g]>0)
              addtoLayer(i,g,valgrd[g],tt);
            g++;
          }
          assert( getLayerDepth(i)>0.0 );
	      }
                                                }
      else
	    {
	      //Layer is bedrock, so make a new layer on top to deposit into
	      //or, depositing more than maxregdep
	      makeNewLayerBelow(i-1, tLayer::kSed, KRnew, update, tt, 
                          new_sed_bulk_density_);
	    }
      }
    }
  
  //Check to make sure that top layer is not too deep
  if(getLayerDepth(i)>1.1*maxregdep && getLayerSed(i) != tLayer::kBedRock ){
    //Make a top layer that is maxregdep deep so that further erosion
    //is not screwed up
    hupdate = addtoLayer(i, -1*maxregdep);
    for(size_t g=0; g<numg; g++) {
      hupdate[g]=-1* hupdate[g];
      assert( hupdate[g]>=0.0 ); //GT
    }
    makeNewLayerBelow(-1, tLayer::kSed, getLayerErody(i), hupdate, tt, 
                      new_sed_bulk_density_);
    setLayerRtime(i,0.);
  }
  
  //   if(getLayerDepth(0)>1.1*maxregdep){
  //   std::cout<<"TOO MUCH SEDIMENT IN TOP LAYER"<<std::endl;
  //   TellAll();
  //}
  
return valgrd;
  }

/**************************************************************
 ** tLNode::addtoLayer(int i, int g, double val, double tt)
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
void tLNode::addtoLayer(int i, int g, double val, double tt)
{
  // For adding material to a layer, grain size by grain size
  tLayer *hlp = layerlist.getIthDataPtrNC( i );

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
 **  tLNode::addtoLayer(int i, double val)
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
tArray<double> tLNode::addtoLayer(int i, double val)
{
  assert( val<0.0 ); // Function should only be called for erosion

  tArray<double> ret(numg);

  tLayer *hlp = layerlist.getIthDataPtrNC( i );

  if(hlp->getDepth()+val>1e-7)
    {
      // have enough material in this layer to fufill all erosion needs
      const double amt=hlp->getDepth();
      size_t n=0;
      while(n<numg)
	{
	  ret[n]=hlp->getDgrade(n)*val/amt;
	  hlp->addDgrade(n,hlp->getDgrade(n)*val/amt);
	  n++;
	}
      return ret;
    }
  else
    {
      // need to remove entire layer
      size_t n=0;
      while(n<numg)
	{
	  ret[n]=-1*hlp->getDgrade(n);
	  n++;
	}
      removeLayer(i);
      return ret;
    }

}

/*******************************************************
 ** tLNode::removeLayer(int i)
 **
 ** Removes layer i.
 **
 ** called by EroDep which updates layers.
 **
 **    Modifications:
 **      - 3/30/00: apparently the existing code was
 **        actually removing layer i+1. Bug fixed (GT).
 *******************************************************/
void tLNode::removeLayer(int i)
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
 ** tLNode::makeNewLayerBelow(int i, int sd, double erd, tArray<double> sz, double tt)
 ** Makes a new layer below layer i.
 ** if i<0 then make new top layer.
 ** sd = sediment flag
 ** erd = erodibility
 ** sz = array of depths of each grain size class
 ** tt = time
 ** Creation and recent time set to current time, exposure time updated
 ** in erodep.
 ********************************************************************/
void tLNode::makeNewLayerBelow(int i, tLayer::tSed_t sd, double erd,
			       tArray<double> const &sz, double tt, 
             double bulk_density )
{
  tLayer hlp, niclay;
  size_t n;

  hlp.setCtime(tt);
  hlp.setRtime(tt);
  hlp.setEtime(0.);
  hlp.setBulkDensity(bulk_density);
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

  if(i>=0){

    tListIter<tLayer> ly ( layerlist );
    tLayer  * pt;
    pt=ly.FirstP();

    int n=0;

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

//Returns depth of all the layers - pretty useless for the current model
double tLNode::getTotalLayerDepth() const
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


double tLNode::DistFromOldXY() const
{
  tArray< double > const &oldpos =  chan.migration.xyzd;
  const double xo = oldpos[0];
  const double yo = oldpos[1];

  return sqrt( (x-xo) * (x-xo) + (y-yo) * (y-yo) );
}

//#undef kBugTime


/******************************************************************
 **  tLNode::CopyLayerList
 **
 **  Assigns the layer list in fromNode to this node. Information
 **  about the pre-existing layer list is lost.
 **
 **  Created in order to assign stratigraphy from a footwall point
 **  to the surface node when the fault plane is exhumed.
 **
 **  Created:  12/99 GT
 **
 ********************************************************************/
void tLNode::CopyLayerList( tLNode const * fromNode )
{
  assert( fromNode!=0 );
  layerlist = fromNode->layerlist;
}

// SL, 8/2010: New function to find total depth of regolith/sediment
// above bedrock with layers (so no need to separately keep track of 
// alluvial thickness).
double tLNode::getRegolithDepth()
{
  tListIter< tLayer > lI( layerlist );
  double soilThickness(0.0);
  for( tLayer *lP=lI.FirstP(); lP->getSed() == tLayer::kSed; lP=lI.NextP() )
    soilThickness += lP->getDepth();
  return soilThickness;
}


