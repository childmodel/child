/**************************************************************************\
**
**  tLNode.cpp
**
**  Functions for derived class tLNode and its member classes
**
**  $Id: tLNode.cpp,v 1.46 1998-06-15 23:10:36 nmgaspar Exp $
\**************************************************************************/

#include <assert.h>
#include <math.h>
#include "../errors/errors.h"
#include "tLNode.h"
#define kBugTime 5000000

/*************************************************************************
**  tLayer::tLayer : Constructor function for tLayer
*************************************************************************/
tLayer::tLayer ()
        : dgrade()
{
   ctime=0;
   rtime=0;
   depth=0;
   erody=0;
   sed=0;
   //cout << "tLayer( num )" << endl;
}

tLayer::tLayer ( int num )
        : dgrade( num )
{
   ctime=0;
   rtime=0;
   depth=0;
   erody=0;
   sed=0;
     //cout << "tLayer( num )" << endl;
}

tLayer::tLayer( const tLayer &orig )                         //tLayer
        :dgrade( orig.dgrade )
{
   ctime=orig.ctime;
   rtime=orig.rtime;
   depth=orig.depth;
   erody=orig.erody;
   sed=orig.sed;
   
  //if( &orig != 0 )
  // {
  //    dpth = orig.dpth;
  // }
}

tLayer::~tLayer()
{
     //cout << "    ~tLayer()" << endl;
}

const tLayer &tLayer::operator=( const tLayer &right )     //tLayer
{
   if( &right != this )
   {
      dgrade = right.dgrade;
      ctime=right.ctime;
      rtime=right.rtime;
      depth=right.depth;
      erody=right.erody;
      sed=right.sed; 
     
   }
   return *this;
}


void tLayer::setCtime( double tt )
{
   ctime = tt; 
}

double tLayer::getCtime() const 
{
   return ctime;
}

void tLayer::setRtime( double tt )
{
   rtime = tt;
}

double tLayer::getRtime() const 
{
   return rtime;
}

void tLayer::setDepth( double dep)
{
   // while updating depth, dgrade info is automatically updated
   // to keep the same texture (as long as size of dgrade > 0)
   // So if texture is going to change too - you MUST update dgrade
   // not depth because depth will be automatically updated in dgrade
   if(dgrade.getSize()>0 && depth>0){      
      double sum=0;
      tArray<double> prop;
      prop.setSize(dgrade.getSize());
      int i=0;
      while(i<dgrade.getSize()){
         prop[i]=dgrade[i]/depth;
         dgrade[i]=prop[i]*dep;
         sum+=prop[i];
         i++;
         
      }
      if(fabs(sum-1)>0.01)
          ReportFatalError("Somewhere grain sizes got messed up");
      depth=dep;
   }
   else
       depth=dep;
}

double tLayer::getDepth() const 
{
   return depth;
}

void tLayer::setErody( double ero)
{
   erody = ero;
}

double tLayer::getErody() const 
{
   return erody;
}

void tLayer::setSed( int rg)
{
   sed = rg;
}

int tLayer::getSed() const 
{
   return sed;
}

void tLayer::setDgradesize( int i )
{
   dgrade.setSize(i);
}

int tLayer::getDgradesize( )
{
   return dgrade.getSize();
}


void tLayer::setDgrade( int i, double size )
{
   if(i>=dgrade.getSize())
      ReportFatalError( "Trying to index sediment sizes that don't exist in layer");
   dgrade[i]=size;
   // Automatically update depth when dgrade is changed
   int j=0;
   double sum=0;
   while(j<dgrade.getSize()){
      sum += dgrade[j];
      j++;
   }
   depth=sum;
}

void tLayer::addDgrade( int i, double size )
{
   if(i>=dgrade.getSize())
      ReportFatalError( "Trying to index sediment sizes that don't exist in layer");
   dgrade[i]+=size;
   depth+=size;
}

double tLayer::getDgrade( int i)
{
   if(i>=dgrade.getSize())
      ReportFatalError( "Trying to index sediment sizes that don't exist in layer");
   return dgrade[i];
}

tArray< double >
tLayer::getDgrade( ) const
{
   return dgrade;
}
   


tErode::tErode()                                                   //tErode
{
     //erodtype = 0;
   sedinput = zp = qs = qsp = qsin = qsinp = tau = dz = 0.0;
   nsmpts = 0;
     //cout << "  tErode()" << endl;
}

tErode::tErode( const tErode &orig )                                //tErode
{
   if( &orig != 0 )
   {
        //erodtype = orig.erodtype;
      sedinput = orig.sedinput;
      zp = orig.zp;
      qs = orig.qs;
      qsp = orig.qsp;
      qsin = orig.qsin;
      qsinp = orig.qsinp;
      tau = orig.tau;
      nsmpts = orig.nsmpts;
      dz = orig.dz;
   }
     //cout << "  tErode( orig )" << endl;
}

tErode::tErode( int numg, int nums )                                //tErode
{
   nsmpts = nums;
   sedinput = zp = qs = qsp = qsin = qsinp = tau = dz = 0.0;
     //cout << "  tErode( numg, nums )" << endl;
}

tErode::~tErode()                                                   //tErode
{
     //cout << "    ~tErode()" << endl;
}

//assignment
const tErode &tErode::operator=( const tErode &right )     //tErode
{
   if( &right != this )
   {
      sedinput = right.sedinput;
      zp = right.zp;
      qs = right.qs;
      qsp = right.qsp;
      qsin = right.qsin;
      qsinp = right.qsinp;
      tau = right.tau;
      nsmpts = right.nsmpts;
      dz = right.dz;
      smooth = right.smooth;
   }
   return *this;
}

tMeander::tMeander()                                              //tMeander
        : xyzd(4)
{
   meander = head = reachmember = 0;
   newx = newy = deltax = deltay = zoldright = zoldleft = bankrough = 0.0;
     //cout << "  tMeander()" << endl;
}

tMeander::tMeander( const tMeander &orig )                        //tMeander
        : xyzd( orig.xyzd )
{
   if( &orig != 0 )
   {
      meander = orig.meander;
      head = orig.head;
      reachmember = orig.reachmember;
      newx = orig.newx;
      newy = orig.newy;
      deltax = orig.deltax;
      deltay = orig.deltay;
      zoldright = orig.zoldright;
      zoldleft = orig.zoldleft;
      bankrough = orig.bankrough;
   }
     //cout << "  tMeander( orig )" << endl;
}

tMeander::tMeander( int state, double x, double y )                //tMeander
        : xyzd(4)
{
     //chanptr = 0;
   meander = state;
   newx = x;
   newy = y;
   head = reachmember = 0;
   deltax = deltay = zoldright = zoldleft = bankrough = 0.0;
     //cout << "  tMeander( state, x, y )" << endl;
}

tMeander::~tMeander()                                             //tMeander
{
     //cout << "    ~tMeander()" << endl;
}

//assignment
const tMeander &tMeander::operator=( const tMeander &right )     //tMeander
{
   if( &right != this )
   {
      meander = right.meander;
      head = right.head;
      reachmember = right.reachmember;
      newx = right.newx;
      newy = right.newy;
      deltax = right.deltax;
      deltay = right.deltay;
      zoldright = right.zoldright;
      zoldleft = right.zoldleft;
      bankrough = right.bankrough;
      xyzd = right.xyzd;
   }
   return *this;
}

tBedrock::tBedrock()                                             //tBedrock
{
   erodibility = 0;
     //cout << "  tBedrock()" << endl;
}

tBedrock::tBedrock( const tBedrock &orig )                       //tBedrock
{
   if( &orig != 0 )
   {
      erodibility = orig.erodibility;
   }
     //cout << "  tBedrock( orig )" << endl;
}

tBedrock::~tBedrock()                                            //tBedrock
{
     //cout << "    ~tBedrock()" << endl;
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

tSurface::tSurface()                                             //tSurface
{
   veg = tauc = vegerody = 0.0;
     //cout << "  tSurface()" << endl;
}

tSurface::tSurface( const tSurface &orig )                       //tSurface
{
   if( &orig != 0 )
   {
      veg = orig.veg;
      tauc = orig.tauc;
      vegerody = orig.vegerody;
   }
     //cout << "  tSurface( orig )" << endl;
}

tSurface::~tSurface()                                            //tSurface
{
     //cout << "    ~tSurface()" << endl;
}

//assignment
const tSurface &tSurface::operator=( const tSurface &right )     //tSurface
{
   if( &right != this )
   {
      veg = right.veg;
      tauc = right.tauc;
      vegerody = right.vegerody;
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
        : dgrade()
{
   thickness = 0;
     //cout << "  tRegolith()" << endl;
}

tRegolith::tRegolith( tInputFile &infile )                     //tRegolith
  : dgrade( )
{
  int i;
  char add, name[20];
  double help, sum, numg;

   cout << "tRegolith(infile)\n";
   thickness = infile.ReadItem( thickness, "REGINIT" );
   numg = infile.ReadItem( numg, "NUMGRNSIZE" );

   if( numg>1 )
   {
     dgrade.setSize( numg );
     sum = 0;
     i=0;
     add='1';
     while ( i<numg ){
       strcpy( name, "REGPROPORTION");
       strcat( name, &add ); 
       help = infile.ReadItem( help, name);
       dgrade[i] = help*thickness;
       sum += dgrade[i];
       i++;
       add++;
     }
     if(fabs(sum-thickness)>0.01)
         ReportFatalError("Problem with the proportion of grain sizes in input file");
     
   }

   cout << "end tRegolith(infile)\n" << flush;
   
}

tRegolith::tRegolith( const tRegolith &orig )                   //tRegolith
        : dgrade( orig.dgrade )
{
   if( &orig != 0 )
   {
      thickness = orig.thickness;
   }
   cout << "  tRegolith( orig ) " << thickness << endl;
}


tRegolith::~tRegolith()                                         //tRegolith
{
     //cout << "    ~tRegolith()" << endl;
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
{
   drarea = q = chanwidth = hydrwidth = channrough = hydrnrough =
       chandepth = hydrdepth = chanslope = hydrslope = diam = 0;
   diam = kVeryHigh;
     //cout << "  tChannel()" << endl;
}

tChannel::tChannel( const tChannel &orig )                       //tChannel
        : /*erosion( orig.erosion ),*/ migration( orig.migration )
{
   if( &orig != 0 )
   {
      drarea = orig.drarea;
      q = orig.q;
      chanwidth = orig.chanwidth;
      chandepth = orig.chandepth;
      channrough = orig.channrough;
      chanslope = orig.chanslope;
      hydrwidth = orig.hydrwidth;
      hydrdepth = orig.hydrdepth;
      hydrnrough = orig.hydrnrough;
      hydrslope = orig.hydrslope;
      diam = orig.diam;
   }
     //cout << "  tChannel( orig )" << endl;
}

tChannel::~tChannel()                                            //tChannel
{
     //cout << "    ~tChannel()" << endl;
}

//assignment
const tChannel &tChannel::operator=( const tChannel &right )     //tChannel
{
   if( &right != this )
   {
      //erosion = right.erosion;
      migration = right.migration;
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
      diam = right.diam;
   }
   return *this;
}


/***** Functions for class tLNode *****************************************/

/**************************************************************************\
**
**  Constructors:
**   - Default: initializes several variables to zero and calls
**              constructors for embedded objects.
**   - Input file: takes a tInputFile reference as an argument, and 
**              passes it to the constructor(s) for embedded objects so
**              that they can read values for default parameters, etc.
**   - Copy: copies all fields of tLNode and calls copy constructors for
**              embedded objects.
**
\**************************************************************************/

// initialize static members numg and grade
int tLNode::numg = 0;
tArray<double> tLNode::grade = 1;
double tLNode::maxregdep = 1;

tLNode::tLNode()                                                   //tLNode
        : tNode(), rock(), surf(), reg(), chan(), qsm(), qsinm(), 
          layerlist()
{
     //cout << "=>tLNode()" << endl;
   flood = 0;
   flowedge = 0;
   tracer = 0;
   dzdt = drdt = qs = qsin = uplift = 0.0;
   //numg=0;
}

tLNode::tLNode( tInputFile &infile )                               //tLNode
        : tNode(), rock(), surf(), reg(), chan(), qsm(), qsinm(), 
          layerlist()
{
   int i;
   char add, name[20];
   double help, helpmax, extra, sum, sumbr;
   tLayer layhelp, niclay;
   tArray<double> dgradehelp;
   tArray<double> dgradebrhelp;
   
   cout << "entered tlnode " << endl;
   
   //cout << "=>tLNode( infile )" << endl;
   flood = 0;
   flowedge = 0;
   tracer = 0;
   dzdt = drdt = qs = qsin = uplift = 0.0;
   numg = infile.ReadItem( numg, "NUMGRNSIZE" );
   cout << "value of numg is " << numg << endl;
   grade.setSize( numg );
   maxregdep = infile.ReadItem( maxregdep, "MAXREGDEPTH" );
   cout << "maximum regolith depth is " << maxregdep << endl;
   
   dgradehelp.setSize( numg );
   dgradebrhelp.setSize( numg );   
   sum = 0;
   i=0;
   add='1';
   
   while ( i<numg ){
      // Reading in grain size info and proportions for intital regolith
      strcpy( name, "GRAINDIAM");
      strcat( name, &add ); 
      help = infile.ReadItem( help, name);
      grade[i] = help;
      strcpy( name, "REGPROPORTION");
      strcat( name, &add ); 
      help = infile.ReadItem( help, name);
      dgradehelp[i]=help;
      sum += help;
      strcpy( name, "BRPROPORTION");
      strcat( name, &add ); 
      help = infile.ReadItem( help, name);
      dgradebrhelp[i]=help;
      sumbr += help;
      i++;
      add++;
   }

   if(fabs(sum-1.0)>0.001)
         ReportFatalError("Problem with the regolith proportion of grain sizes in input file");
   if(fabs(sumbr-1.0)>0.001)
         ReportFatalError("Problem with the bedrock proportion of grain sizes in input file");
//   cout << "nic, you got past fatal errors" << endl;
   
   if( numg > 1 ){
      qsm.setSize( numg );
      qsinm.setSize( numg );
      
   }

   help = infile.ReadItem( help, "REGINIT" );
//    cout << "regolith depth is " << help << endl;

//    cout << "1 grain size is "<<grade[0]<<" 2 grain size is "<<grade[1] << endl;
   
   if( help > 0){
      // Make a bedrock and regolith layer, possibly two depending
      // on the depth of the regolith layer.  The code will decide
      // the total number of layers needed.  By default the regolith
      // layer(s) is/are put on top of the bedrock.
      layhelp.setCtime( 0.0 );
      //nic - need to put real time here (for above and below)-check into this
      layhelp.setRtime( 0.0 );
      // Bedrock layer items read in and set
      help = infile.ReadItem( help, "KB");
      layhelp.setErody(help);
      layhelp.setSed(0);
      
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
      layhelp.setSed(1);
      help = infile.ReadItem( help, "KR" );
      layhelp.setErody(help);
      help = infile.ReadItem( help, "REGINIT");
      extra = 0;
      if(help > maxregdep){
         // too much regolith, create two layers the bottom layer is made here
         extra = help - maxregdep;
         //layhelp.setDepth(extra);
         layhelp.setDgradesize(numg);
         i=0;
         while(i<numg){
            layhelp.setDgrade(i, extra*dgradehelp[i]);
            i++;
         }
         layerlist.insertAtFront( layhelp );

         // the top regolith layer is now made
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
      layhelp.setCtime(0.0);
      //since time isn't global time is just set to zero.
      //assume you are always starting from zero, if not then
      //you use a restart file and layer info will be read from there.
      layhelp.setRtime(0.0);
      // Bedrock layer items set
      help = infile.ReadItem( help, "KR");
      layhelp.setErody(help);
      layhelp.setSed(0);
      layhelp.setDgradesize(numg);
      i=0;
      help = infile.ReadItem( help, "BEDROCKDEPTH");
      while(i<numg){
         layhelp.setDgrade(i, help*dgradebrhelp[i]);
         i++;
      }
      
      layerlist.insertAtBack( layhelp );
      
   }

   cout << layerlist.getSize() << " layers created " << endl;

   //cout << "nic, what is size of grade array ? " << grade.getSize() << endl;
   //cout << "1 is size " << grade[0] << " 2 is size " << grade[1] << endl;   

   i=0;
   while(i<layerlist.getSize()){
      cout << "layer " << i+1 << endl;
//       niclay = layerlist.getIthData(i);
//       cout << "layer creation time is " << niclay.getCtime() << endl;
//       cout << "layer recent time is " << niclay.getRtime() << endl;
//       cout << "layer depth is " << niclay.getDepth() << endl;
//       cout << "layer erodibility is " << niclay.getErody() << endl;
//       cout << "is layer sediment? " << niclay.getSed() << endl;
//       cout << "dgrade 1 is " << niclay.getDgrade(0) << endl;
//       cout << "dgrade 2 is " << niclay.getDgrade(1) << endl;
//       cout << "layer " << i+1 << " now from the getlayer func" << endl;
      cout << "layer creation time is " << getLayerCtime(i) << endl;
      cout << "layer recent time is " << getLayerRtime(i) << endl;
      cout << "layer depth is " << getLayerDepth(i) << endl;
      cout << "layer erodibility is " << getLayerErody(i) << endl;
      cout << "is layer sediment? " << getLayerSed(i) << endl;
      cout << "dgrade 1 is " << getLayerDgrade(i,0) << endl;
      cout << "dgrade 2 is " << getLayerDgrade(i,1) << endl;
      i++;
      
   }

   cout << "does elevation function work?" << endl;
   cout << "elevation is " << getElev() << endl;

   
//      cout << "nic is now removing material from first layer" << endl;
//      dgradehelp[0]=-0.5;
//      dgradehelp[1]=-0.2;
//      dgradebrhelp=updateLayerDepth(0,dgradehelp,0.0);
//      cout << "returned stuff is 1 -> " << dgradebrhelp[0] << endl;
//      cout << "ret stuff 2 is -> " << dgradebrhelp[1] << endl;
   
//      i=0;
//      while(i<layerlist.getSize()){
//         cout << "layer " << i+1 << endl;
//         niclay = layerlist.getIthData(i);
//         cout << "layer creation time is " << getLayerCtime(i) << endl;
//         cout << "layer recent time is " << getLayerRtime(i) << endl;
//         cout << "layer depth is " << getLayerDepth(i) << endl;
//         cout << "layer erodibility is " << getLayerErody(i) << endl;
//         cout << "is layer sediment? " << getLayerSed(i) << endl;
//         cout << "dgrade 1 is " << getLayerDgrade(i,0) << endl;
//         cout << "dgrade 2 is " << getLayerDgrade(i,1) << endl;
//         i++;  
//      }
}

tLNode::tLNode( const tLNode &orig )                               //tLNode
        : tNode( orig ),
          rock( orig.rock ), surf( orig.surf ),
          reg( orig.reg ), chan( orig.chan ), qsm( orig.qsm),
          qsinm( orig.qsinm ), //grade( orig.grade ),
          layerlist( orig.layerlist)
{
   flowedge = orig.flowedge;
   flood = orig.flood;
   tracer = orig.tracer;
   dzdt = orig.dzdt;
   drdt = orig.drdt;
   qs = orig.qs;
   qsin = orig.qsin;
   uplift = orig.uplift;
   //numg = orig.numg;
   
   //cout << "=>tLNode( orig )" << endl;
}

tLNode::~tLNode()                                                  //tLNode
{
     //cout << "    ~tLNode()" << endl;
}

//assignment
const tLNode &tLNode::operator=( const tLNode &right )                  //tNode
{
   if( &right != this )
   {
      tNode::operator=( right );
      rock = right.rock;
      surf = right.surf;
      reg = right.reg;
      chan = right.chan;
      flowedge = right.flowedge;
      flood = right.flood;
      tracer = right.tracer;
      dzdt = right.dzdt;
      drdt = right.drdt;
      qs = right.qs;
      qsin = right.qsin;
      uplift = right.uplift;
      qsm = right.qsm;
      qsinm = right.qsinm;
      layerlist = right.layerlist;
   }
   return *this;
}

//"get" and "set" functions; most simply return or set a data value, respectively:

const tBedrock &tLNode::getRock() const {return rock;}
const tSurface &tLNode::getSurf() const {return surf;}
const tRegolith &tLNode::getReg() const {return reg;}
const tChannel &tLNode::getChan() const {return chan;}

void tLNode::setRock( const tBedrock & val ) {rock = val;}
void tLNode::setSurf( const tSurface & val ) {surf = val;}
void tLNode::setReg( const tRegolith & val ) {reg = val;}
void tLNode::setChan( const tChannel & val ) {chan = val;}

int tLNode::getFloodStatus() {   return flood; }

void tLNode::setFloodStatus( int status )
{
   flood = status;
}

tEdge * tLNode::getFlowEdg() 
{
   return flowedge;
}

void tLNode::setFlowEdg( tEdge * edg )
{
   assert( edg > 0 );  // Fails when passed an invalid edge
   //cout << "Setting flow edge to edge # " << edg->getID() << endl;
   flowedge = edg;
}

void tLNode::setDrArea( double val ) {chan.drarea = val;}
void tLNode::AddDrArea( double val ) {chan.drarea += val;}
void tLNode::AddDischarge( double val ) {chan.q += val;}

tLNode * tLNode::getDownstrmNbr()
{
   //assert( flowedge!=0 );
   if( flowedge == 0 ) return 0;
   return (tLNode *)flowedge->getDestinationPtrNC();     
}

int tLNode::Meanders() const {return chan.migration.meander;}
void tLNode::setMeanderStatus( int val )
{chan.migration.meander = ( val == 0 || val == 1 ) ? val : 0;}


double tLNode::getHydrWidth() const {return chan.hydrwidth;}
double tLNode::getChanWidth() const {return chan.chanwidth;}
double tLNode::getHydrDepth() const {return chan.hydrdepth;}
double tLNode::getChanDepth() const {return chan.chandepth;}
double tLNode::getHydrRough() const {return chan.hydrnrough;}
double tLNode::getChanRough() const {return chan.channrough;}
double tLNode::getHydrSlope() const {return chan.hydrslope;}
double tLNode::getChanSlope() const {return chan.chanslope;}
double tLNode::getDiam() const {return chan.diam;}
double tLNode::getBankRough() const {return chan.migration.bankrough;}

//TODO: suggest doing away with the zero test for performance enhancement
void tLNode::setHydrWidth( double val )  {chan.hydrwidth = ( val > 0 ) ? val : 0;}
void tLNode::setChanWidth( double val )  {chan.chanwidth = ( val > 0 ) ? val : 0;}
void tLNode::setHydrDepth( double val )  {chan.hydrdepth = ( val > 0 ) ? val : 0;}
void tLNode::setChanDepth( double val )  {chan.chandepth = ( val > 0 ) ? val : 0;}
void tLNode::setHydrRough( double val )  {chan.hydrnrough = ( val > 0 ) ? val : 0;}
void tLNode::setChanRough( double val )  {chan.channrough = ( val > 0 ) ? val : 0;}
void tLNode::setHydrSlope( double val )  {chan.hydrslope = ( val > 0 ) ? val : 0;}
void tLNode::setChanSlope( double val )  {chan.chanslope = ( val > 0 ) ? val : 0;}
void tLNode::setBankRough( double val )
{chan.migration.bankrough = ( val > 0 ) ? val : 0;}

double tLNode::getDrArea() const {return chan.drarea;}

tArray< double >
tLNode::getZOld() const
{
   tArray< double > rl(2);
   rl[0] = chan.migration.zoldright;
   rl[1] = chan.migration.zoldleft;
   return rl;
}

void tLNode::setZOld( double right, double left )
{
   chan.migration.zoldright = right;
   chan.migration.zoldleft = left;
}

tArray< double >                                                   //tNode
tLNode::getNew2DCoords() const
{
   tArray< double > xy(2);
   if( Meanders() )
   {
      xy[0] = chan.migration.newx;
      xy[1] = chan.migration.newy;
   }
   else
   {
      xy[0] = x;
      xy[1] = y;
   }
   return xy;
}

tArray< double >                                                   //tNode
tLNode::getNew3DCoords() const
{
   tArray< double > xyz(3);
   if( Meanders() )
   {
      xyz[0] = chan.migration.newx;
      xyz[1] = chan.migration.newy;
   }
   else
   {
      xyz[0] = x;
      xyz[1] = y;
   }
   xyz[2] = z;
   return xyz;
}

void tLNode::setNew2DCoords( double val1, double val2 )                    //tNode
{
   chan.migration.newx = val1;
   chan.migration.newy = val2;
}

tArray< double > tLNode::
getLatDisplace() const
{
   tArray< double > xy(2);
   xy[0] = chan.migration.deltax;
   xy[1] = chan.migration.deltay;
   return xy;
}

void tLNode::setLatDisplace( double dx, double dy )
{
   chan.migration.deltax = dx;
   chan.migration.deltay = dy;
}

void tLNode::addLatDisplace( double dx, double dy )
{
   chan.migration.deltax += dx;
   chan.migration.deltay += dy;
}



void tLNode::setDischarge( double val ) {chan.q = ( val > 0 ) ? val : 0;}

void tLNode::RevertToOldCoords()
{
   chan.migration.newx = x;
   chan.migration.newy = y;
}

void tLNode::UpdateCoords()
{
   x = chan.migration.newx;
   y = chan.migration.newy;
}

// nb: if channel is integrated into node, change this
double tLNode::getQ()
{
   return chan.q;
}

/************************************************************************\
**  GetSlope: Computes and returns the slope of the node's flowedg, or
**  zero if the slope is less than zero.
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
double tLNode::getSlope()
{
   int ctr;
   double rlen, curlen, slp, delz, downz;
   tLNode *dn, *on, *tn;
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
            return (-1);
            //ReportFatalError("infinite loop in tLNode::GetSlope(), 1st loop");
         }
         curlen += dn->flowedge->getLength();
         dn = dn->getDownstrmNbr();
         assert( dn != 0 );
      }
      assert( curlen > 0 );
      downz = dn->z;
      if( timetrack >= kBugTime ) cout << "GetSlope 1; " << flush;
      delz = z - downz;
      if( timetrack >= kBugTime ) cout << "GS 2; " << flush;
      slp = delz / curlen;
      if( timetrack >= kBugTime ) cout << "GS 3; " << flush;
      on = dn;
      ctr = 0;
      while( on->getBoundaryFlag() == kNonBoundary &&
             on->flowedge->getLength() > 0 )
      {
         ctr++;
         if( ctr > kLargeNumber )
         {
            return (-1);
            //ReportFatalError("infinite loop in tLNode::GetSlope(), 2nd loop");
         }
         on = on->getDownstrmNbr();
      }
      if( z - on->z < 0.0 ) slp = 0.0;
   }
   else slp = (z - getDownstrmNbr()->z ) / flowedge->getLength();
   if( timetrack >= kBugTime ) cout << "GS 4; " << endl << flush;
   if( slp>=0.0 ) return slp;
   else return 0.0;
}
#undef kLargeNumber

double tLNode::getDSlopeDt()
{
   double rlen, curlen, slp;
   tLNode *dn;
   assert( flowedge != 0 );
   assert( flowedge->getLength()>0 ); // failure means lengths not init'd
   if( Meanders() )
   {
      rlen = 10.0 * chan.chanwidth;
      curlen = 0.0;
      dn = this;
      while( curlen < rlen && dn->getBoundaryFlag() == kNonBoundary )
      {
         curlen += dn->flowedge->getLength();
         dn = dn->getDownstrmNbr();
      }
      assert( curlen > 0 );
      //slp = (dzdt - dn->dzdt + uplift - dn->uplift ) / curlen;
   }
   
   else
   {
      curlen = flowedge->getLength();
      assert( curlen > 0.0 );
      dn = getDownstrmNbr();
      //slp = (dzdt - dn->dzdt + uplift - dn->uplift ) / curlen;
   }
   slp = ( dzdt - dn->dzdt + uplift - dn->uplift ) / curlen;
   if( slp>=0.0 ) return slp;
   else return 0.0;
}




//void tLNode::insertFrontSpokeList( tEdge *eptr )             //tNode
//{spokeList.insertAtFront( eptr );}

//void tLNode::insertBackSpokeList( tEdge *eptr )              //tNode
//{spokeList.insertAtBack( eptr );}

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
double tLNode::DistNew(tLNode * p0,tLNode * p1 )
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
void tLNode::TellAll()
{
   tLNode * nbr;
   
   cout << " NODE " << id << ":\n";
   cout << "  x=" << x << " y=" << y << " z=" << z;
   if( edg ) {
      cout << "  points to edg #" << edg->getID() << endl;
      cout << "  dr area: " << getDrArea() << "  disch: " << getQ()
           << "  boundary: " << boundary << "  flood: " << flood
           << "\n  varea: " << varea << endl;
      
      if( flowedge ) {
         nbr = (tLNode *)flowedge->getDestinationPtrNC();
         cout << "  Flows along edg " << flowedge->getID() << " to node "
              << nbr->getID() << " at (" << nbr->getX() << ","
              << nbr->getY() << "," << nbr->getZ() << ")\n    with vedglen "
              << flowedge->getVEdgLen() << endl;
         cout << "  qs: " << qs << "  qsin: " << qsin << "  slp: "
              << getSlope() << "  reg: " << reg.thickness << endl;
         cout << "  dzdt: " << dzdt << "  drdt: " << drdt << endl;
      }
      else cout << "  Flowedg is undefined\n";
      
   }
   else cout << "  edg is undefined!\n";
   
}
#endif

/**************************************************************************\
**
**  Tracer-sorting routines:
**
**  These routines are utilities that are used in sorting the nodes
**  according to their position within the drainage network. The main
**  sorting algorithm is implemented in tStreamNet::SortNodesByNetOrder().
**  The sorting method works by introducing a "tracer" at each point,
**  then allowing the tracers to iteratively cascade downstream. At each
**  step any nodes not containing tracers are moved to the back of the
**  list. The result is a list sorted in upstream-to-downstream order.
**
**  These utilities do the following:
**    ActivateSortTracer -- injects a single tracer at a node
**    AddTracer -- adds a tracer to a node (ignored if node is a bdy)
**    MoveSortTracerDownstream -- removes a tracer and sends it to the
**                                downstream neighbor (unless the node is
**                                a sink; then the tracer just vanishes)
**    NoMoreTracers -- reports whether there are any tracers left here
**
**  Created by GT 12/97.
**
\**************************************************************************/

void tLNode::ActivateSortTracer()
{ tracer = 1; }

void tLNode::MoveSortTracerDownstream()
{
   tracer--;
   if( flood!=kSink ) getDownstrmNbr()->AddTracer();
}

void tLNode::AddTracer()
{
   if( !boundary ) tracer++;
}

int tLNode::NoMoreTracers()
{
   assert( tracer>=0 );
   return( tracer==0 );
}


/**************************************************************************\
**
**  EroDep
**
**  Erodes or Deposits material of depth dz.
**
**  Note: as of now, just changes elev. Should also change alluvial
**  thickness and set to zero if negative (bedrock erosion).
**
**  1/15/98 gt
\**************************************************************************/
void tLNode::EroDep( double dz )
{
   z += dz;
   //cout << "  eroding " << id << " by " << dz << endl << flush;
   
   //sed += dz;
   //if( sed<0 ) sed=0;
   reg.thickness += dz;
   if( reg.thickness < 0 ) reg.thickness = 0.0;
}

void tLNode::setAlluvThickness( double val )
{
   //reg.thickness = ( val >= 0.0 ) ? val : 0.0;
   // gt changed for performance speedup (if stmt shouldn't ever be needed)
   assert( val>=0 );
   reg.thickness = val;
}

double tLNode::getAlluvThickness() const {return reg.thickness;}

tArray< double >
tLNode::getAlluvThicknessm( ) const 
{
   return reg.dgrade;
}

void tLNode::setVegErody( double val )
{surf.vegerody = ( val >= 0.0 ) ? val : 0.0;}

double tLNode::getVegErody() const {return surf.vegerody;}

void tLNode::setBedErody( double val )
{rock.erodibility = ( val >= 0.0 ) ? val : 0.0;}

double tLNode::getBedErody() const {return rock.erodibility;}

void tLNode::setReachMember( int val )
{chan.migration.reachmember = ( val == 0 || val == 1 ) ? val : 0;}

int tLNode::getReachMember() const {return chan.migration.reachmember;}

void tLNode::setQs( double val ) {qs = val;}

void tLNode::setQs( int i, double val ) 
{
   if(i>=numg)
       ReportFatalError( "Trying to index sediment sizes that don't exist ");
   qsm[i] = val;
   qs += val;
}

double tLNode::getQs() const {return qs;}

double tLNode::getQs( int i)
{
   if(i>=numg)
       ReportFatalError( "Trying to index sediment sizes that don't exist ");
   return qsm[i];
}

tArray< double >
tLNode::getQsm( ) const
{
   return qsm;
}

void tLNode::setQsin( double val ) {qsin = val;}

void tLNode::setQsin( int i, double val ) 
{
   if(i>=numg)
       ReportFatalError( "Trying to index sediment sizes that don't exist ");
   qsinm[i]=val;
   qsin += val;
}

void tLNode::AddQsin( double val ) 
{
   qsin += val;
}

void tLNode::AddQsin( int i, double val )
{
   if(i>=numg)
       ReportFatalError( "Trying to index sediment sizes that don't exist ");
   qsinm[i] += val;
   qsin += val;
   
}

void tLNode::AddQsinm( tArray< double > val )
{
   int i;
   for(i=0; i<val.getSize(); i++)
       qsinm[i] += val[i];
}

double tLNode::getQsin() const {return qsin;}

double tLNode::getQsin( int i )
{
   if(i>=numg)
       ReportFatalError( "Trying to index sediment sizes that don't exist ");
   return qsinm[i];
}

tArray< double >
tLNode::getQsinm( ) const
{
   return qsinm;
}

void tLNode::setGrade( int i, double size )
{
   if(i>=numg)
       ReportFatalError("Trying to set a grain size for an index which is too large");
   grade[i] = size;
}

int tLNode::getNumg() const 
{
   return numg;
}

void tLNode::setNumg( int size )
{
   numg = size;
}

double tLNode::getGrade( int i )
{
   return grade[i];
}

tArray< double >
tLNode::getGrade( ) const
{
   return grade;
}  

double tLNode::getLayerCtime( int i ) const
{
   tLayer hlp;
   hlp = layerlist.getIthData(i);
   return hlp.getCtime();
}

void tLNode::setLayerCtime( int i, double tt)
 {
    tListIter<tLayer> ly ( layerlist );
    tLayer  * hlp;
    hlp=ly.FirstP();
    
    int n=0;

    while(n<i){
       n++;
       hlp=ly.NextP();
    }

    hlp->setCtime( tt );
}

double tLNode::getLayerRtime( int i ) const
{
   tLayer hlp;
   hlp = layerlist.getIthData(i);
   return hlp.getRtime();
}

void tLNode::setLayerRtime( int i, double tt)
 {
    tListIter<tLayer> ly ( layerlist );
    tLayer  * hlp;
    hlp=ly.FirstP();
    
    int n=0;

    while(n<i){
       n++;
       hlp=ly.NextP();
    }

    hlp->setRtime( tt );
}

double tLNode::getLayerDepth( int i ) const
{
   tLayer hlp;
   hlp = layerlist.getIthData(i);
   return hlp.getDepth();
}

void tLNode::setLayerDepth( int i, double dep)
 {
    tListIter<tLayer> ly ( layerlist );
    tLayer  * hlp;
    hlp=ly.FirstP();
    
    int n=0;

    while(n<i){
       n++;
       hlp=ly.NextP();
    }

    hlp->setDepth( dep );
}

double tLNode::getLayerErody( int i ) const
{
   tLayer hlp;
   hlp = layerlist.getIthData(i);
   return hlp.getErody();
}

void tLNode::setLayerErody( int i, double ero)
 {
    tListIter<tLayer> ly ( layerlist );
    tLayer  * hlp;
    hlp=ly.FirstP();
    
    int n=0;

    while(n<i){
       n++;
       hlp=ly.NextP();
    }

    hlp->setErody( ero );
}


int tLNode::getLayerSed( int i ) const
{
   tLayer hlp;
   hlp = layerlist.getIthData(i);
   return hlp.getSed();
}

void tLNode::setLayerSed( int i, int s)
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

double tLNode::getLayerDgrade( int i, int num ) const
{
   tLayer hlp;
   hlp = layerlist.getIthData(i);
   return hlp.getDgrade(num);
}

void tLNode::setLayerDgrade( int i, int g, double val)
 {
    tListIter<tLayer> ly ( layerlist );
    tLayer  * hlp;
    hlp=ly.FirstP();
    
    int n=0;

    while(n<i){
       n++;
       hlp=ly.NextP();
    }

    hlp->setDgrade(g, val );
}

tArray<double> tLNode::updateLayerDepth( int i, tArray<double> valgrd, double tt)
{
   // Checked on 6/13/98 and everything appears to be working!
   int g;
   double amt, val, olddep;

   g=0;
   val=0;
   while(g<numg){
      if(-1*valgrd[g]>getLayerDgrade(i,g))
          valgrd[g]=-1*getLayerDgrade(i,g);
      // Checking to see that there is enough stuff
      val+=valgrd[g];
      g++;
   }   
   // val is now set to the total amount of erosion or deposition
   
   tArray<double> update;
   update.setSize(numg);
   tArray<double> hupdate;
   hupdate.setSize(numg);
   
   if(val<0)
   {
      // TOTAL EROSION
      if(getLayerSed(i) != 0 && getLayerSed(i) == getLayerSed(i+1))
      {
         // Updating will also be done if entering this statement
         // No updating of Bedrock layers.
         // Only update if layer below has the same material
         while(val<-0.0001 && getLayerSed(i) == getLayerSed(i+1))
         {
            // keep eroding until you either get all the material you
            // need to refill the top layer, or you run out of material
            hupdate = addtoLayer(i+1, val, tt);
            g=0;
            while(g<numg)
            {
               val-=hupdate[g];
               update[g] += hupdate[g];
               g++;
            }
         }
         g=0;
         while(g<numg){
            addtoLayer(i, g, valgrd[g], tt); // Erosion 
            addtoLayer(i,g,-1*update[g],tt);//Updating with material from below
            g++;
         }
      }
      else
      {
         // No updating, just eroding
         g=0;
         while(g<numg){
             addtoLayer(i, g, valgrd[g], tt); // Erosion done on this line
             g++;
          }
      }
   }
   else
   {
      // DEPOSITION
      // NIC, all deposition cases tested and appear to be working
      // method seems to make good sense for surface layers
      // but may not be as appropriate for lower layers.
      // You may want to either make a provision for lower layers
      // or else write a new algorithm for them (will you ever deposit
      // into lower layers?  I don't know)
      // Also, no test done to make sure that you are depositing the right
      // material into the layer, would need to pass the flag for this.
      // For now assume that the test will be done in another place.
      if(getLayerDepth(i)+val>maxregdep){
         // Need to move stuff out of top layer to make room for deposited mat
         if(getLayerSed(i) == getLayerSed(i+1) && getLayerDepth(i+1)+val<maxregdep)
         {
            // The layer below is of the appropriate material and has space
            amt = getLayerDepth(i)+val-maxregdep;
            olddep = getLayerDepth(i);
            g=0;
            while(g<numg){
               addtoLayer(i+1,g,amt*getLayerDgrade(i,g)/olddep, tt);
               // putting material from top layer to layer below
               addtoLayer(i,g,-1*amt*getLayerDgrade(i,g)/olddep+valgrd[g], tt);
               // changing top layer composition by depositing and removing
               g++;
            }
         }
         else
         {
            // Need to create new layer
            amt = getLayerDepth(i)+val-maxregdep;
            olddep = getLayerDepth(i);
            g=0;
            while(g<numg){
               update[g]=amt*getLayerDgrade(i,g)/olddep;
               // material which will be moved from top layer
               addtoLayer(i,g,-1*amt*getLayerDgrade(i,g)/olddep+valgrd[g], tt);
               // changing top layer composition by depositing
               g++;
            }
            makeNewLayerBelow(i, getLayerSed(i), getLayerErody(i), update, tt);
            // put material into a new layer which is made below
         }
      }
      else
      {
         // No need to move stuff out of top layer, just deposit
         g=0;
         while(g<numg){
            addtoLayer(i,g,valgrd[g],tt);
            g++;
         }
      }
   }
   return valgrd;
}


void tLNode::addtoLayer(int i, int g, double val, double tt)
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

    hlp->addDgrade(g,val);
    if(hlp->getDepth() == 0)
        removeLayer(i);
    hlp->setRtime( tt );
    
}

tArray<double> tLNode::addtoLayer(int i, double val, double tt)
{
   // for removing material from a layer, val is amount to remove
   if(val>0)
       ReportFatalError("Using wrong function for depositing into a layer");

   tListIter<tLayer> ly ( layerlist );
   tLayer  * hlp;
   hlp=ly.FirstP();
   tArray<double> ret;
   ret.setSize(numg);
   double amt;
   
   int n=0;
   while(n<i){
      n++;
      hlp=ly.NextP();
   }

   if(-1*val<hlp->getDepth())
   {
      // have enough material in this layer to fufill all erosion needs
      amt=hlp->getDepth();
      n=0;
      while(n<numg)
      {
         ret[n]=hlp->getDgrade(n)*val/amt;
         hlp->addDgrade(n,hlp->getDgrade(n)*val/amt);
         n++;
      }
      hlp->setRtime(tt);
      return ret;
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
      removeLayer(i);
      return ret;
   }
   
}

void tLNode::removeLayer(int i)
{
    tListIter<tLayer> ly ( layerlist );
    tLayer  * hlp;
    hlp=ly.FirstP();
    
    int n=0;

    while(n<i){
       n++;
       hlp=ly.NextP();
    }

    n=layerlist.removeNext((*hlp), layerlist.getListNode(hlp) );
    if(n==0)
        ReportFatalError("couldn't remove next layer");
    
}

void tLNode::makeNewLayerBelow(int i, int sd, double erd, tArray<double> sz, double tt)
{
   tLayer hlp;
   int n;
   
   hlp.setCtime(tt);
   hlp.setRtime(tt);
   n=0;
   hlp.setDgradesize(numg);
   while(n<numg){
      hlp.setDgrade(n, sz[n]);
      n++;
   }
   hlp.setErody(erd);
   hlp.setSed(sd);

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


int tLNode::getNumLayer() const
{
   return layerlist.getSize();
}

double tLNode::getElev() const
{
   double sz = layerlist.getSize();
   int i=0;
   double elev = 0;
   while(i<sz)
   {
      elev += getLayerDepth(i);
      i++;
   }
   return elev;
}


void tLNode::setDzDt( double val ) {dzdt = val;}

double tLNode::getDzDt() {return dzdt;}

void tLNode::setDrDt( double val ) {drdt = val;}

double tLNode::getDrDt() {return drdt;}

void tLNode::setXYZD( tArray< double > arr )
{
   chan.migration.xyzd = ( arr.getSize() == 4 ) ? arr : tArray< double >(4);
     //cout << "setXYZD: " << chan.migration.xyzd[0] 
     //   << " " << chan.migration.xyzd[1] << " " << chan.migration.xyzd[2]
     //   << " " << chan.migration.xyzd[3] << endl;
}

tArray< double >
tLNode::getXYZD() const {return chan.migration.xyzd;}


double tLNode::DistFromOldXY() const
{
   double xo, yo;
   tArray< double > oldpos( chan.migration.xyzd );
   
   xo = oldpos[0];
   yo = oldpos[1];
   
   return sqrt( (x-xo) * (x-xo) + (y-yo) * (y-yo) );
}


// Tests whether bedrock is exposed at a node
int tLNode::OnBedrock()
{
   // For multi-size model, criterion might be active layer thickness less
   // than a nominal thickness; here, it's just an arbitrary alluvial depth
   return ( reg.thickness<0.1 );
}

void tLNode::setUplift( double val ) {uplift = val;}

double tLNode::getUplift() const {return uplift;}
#undef kBugTime
