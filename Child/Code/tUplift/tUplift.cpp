/************************************************************************/
/**
**  @file tUplift.cpp
**  @brief Functions for class tUplift (see tUplift.h).
**
**  $Id: tUplift.cpp,v 1.25 2004-04-27 09:55:05 childcvs Exp $
*/
/************************************************************************/

#include "tUplift.h"
#include "../errors/errors.h"
#include "../Mathutil/mathutil.h"

/************************************************************************\
**
**  tUplift constructor
**
**  Reads in a type code from the input file, then reads parameters
**  needed for the specified type of uplift.
**
**  The current version supports two types of uplift: spatially uniform,
**  and uniform where Y>=Yf (Yf=fault location), zero elsewhere.
**
**  Inputs:  infile -- input file from which to read parameters
**
\************************************************************************/
tUplift::tUplift_t tUplift::DecodeType(int type){
  if( type < 0 || type > 8 )
    {
      cerr << "I don't recognize the uplift type you asked for ("
           << type << ")\n"
	"Valid uplift types are:\n"
	" 0 - none\n"
	" 1 - Spatially and temporally uniform uplift\n"
	" 2 - Uniform uplift at Y >= fault location, zero elsewhere\n"
	" 3 - Block uplift with strike-slip motion along given Y coord\n"
	" 4 - Propagating fold modeled w/ simple error function curve\n"
	" 5 - 2D cosine-based uplift-subsidence pattern\n"
	" 6 - Block, fault, and foreland sinusoidal fold\n"
	" 7 - Two-sided differential uplift\n"
	" 8 - Fault bend fold\n";
      ReportFatalError( "Please specify a valid uplift type and try again." );
    }
  return static_cast<tUplift_t>(type);
}

tUplift::tUplift( const tInputFile &infile ) :
  duration(0.)
{
  int typeCode_;
   // Find out what kind of uplift the user wants
   typeCode_ = infile.ReadItem( typeCode_, "UPTYPE" );
   typeCode = DecodeType(typeCode_);

   if( typeCode==kNoUplift ) return;

   // get the parameters relevant to that type
   duration = infile.ReadItem( duration, "UPDUR" );
   rate = infile.ReadItem( rate, "UPRATE" );
   switch( typeCode )
   {
      case kNoUplift:
      case k1:
          break;
      case k2:
          faultPosition = infile.ReadItem( faultPosition, "FAULTPOS" );
	  rate2 = infile.ReadItem( rate2, "SUBSRATE" );
          break;
      case k3:
          faultPosition = infile.ReadItem( faultPosition, "FAULTPOS" );
          slipRate = infile.ReadItem( slipRate, "SLIPRATE" );
          break;
      case k4:
          faultPosition = infile.ReadItem( faultPosition, "FAULTPOS" );
          slipRate = infile.ReadItem( slipRate, "FOLDPROPRATE" );
          foldParam = infile.ReadItem( foldParam, "FOLDWAVELEN" );
          foldParam = 4.0/foldParam;
          break;
      case k5:
          foldParam = infile.ReadItem( foldParam, "FOLDWAVELEN" );
          slipRate = infile.ReadItem( slipRate, "TIGHTENINGRATE" );
          faultPosition = infile.ReadItem( faultPosition, "ANTICLINEYCOORD" );
          positionParam1 = infile.ReadItem( positionParam1, "ANTICLINEXCOORD");
          deformStartTime1 =
              infile.ReadItem( deformStartTime1, "YFOLDINGSTART" );
          foldParam2 = infile.ReadItem( foldParam2, "UPSUBRATIO" );
	  break;
     case k6:
          foldParam = infile.ReadItem( foldParam, "FOLDWAVELEN" );
          slipRate = infile.ReadItem( slipRate, "FOLDLATRATE" );
          faultPosition = infile.ReadItem( faultPosition, "FAULTPOS" );
	  rate2 = infile.ReadItem( rate2, "FOLDUPRATE" );
	  foldParam2 = infile.ReadItem( foldParam2, "FOLDPOSITION" );
          break;
     case k7:
          rate2 = infile.ReadItem( rate2, "BLFALL_UPPER" );
	  positionParam1 = infile.ReadItem( positionParam1, "BLDIVIDINGLINE" );
	  break;
     case k8:
     	  slipRate = infile.ReadItem( slipRate, "SLIPRATE" );
     	  faultPosition = infile.ReadItem( faultPosition, "FAULTPOS" );
     	  flatDepth = infile.ReadItem( flatDepth, "FLATDEPTH" );
     	  rampDip = infile.ReadItem( rampDip, "RAMPDIP" );
     	  kinkDip = infile.ReadItem( kinkDip, "KINKDIP" );
     	  upperKinkDip = infile.ReadItem( upperKinkDip, "UPPERKINKDIP" );
     	  meanElevation = infile.ReadItem( meanElevation, "MEAN_ELEV" );
     	  break;
   }

}


/************************************************************************\
**
**  tUplift::DoUplift
**
**  Calls the appropriate function to perform the desired type of
**  uplift.
**
**  Inputs:  mp -- pointer to the mesh
**           delt -- duration of uplift
**
\************************************************************************/
void tUplift::DoUplift( tMesh<tLNode> *mp, double delt )
{
   switch( typeCode )
   {
      case kNoUplift:
          break;
      case k1:
          UpliftUniform( mp, delt );
          break;
      case k2:
          BlockUplift( mp, delt );
          break;
      case k3:
          BlockUplift( mp, delt );
          StrikeSlip( mp, delt );
          break;
      case k4:
          FoldPropErf( mp, delt );
          break;
      case k5:
          CosineWarp2D( mp, delt );
          break;
      case k6:
	  BlockUplift( mp, delt );
	  PropagatingFold( mp, delt );
	  break;
      case k7:
	  TwoSideDifferential( mp, delt );
	  break;
      case k8:
          FaultBendFold( mp, delt );
          FaultBendFold2( mp, delt );
          break;
   }

}


/************************************************************************\
**
**  tUplift::UpliftUniform
**
**  Uniform uplift at a constant rate across the entire domain (but not
**  including boundaries).
**
**  Inputs:  mp -- pointer to the mesh
**           delt -- duration of uplift
**
\************************************************************************/
void tUplift::UpliftUniform( tMesh<tLNode> *mp, double delt )
{
   assert( mp!=0 );
   tLNode *cn;
   tMesh<tLNode>::nodeListIter_t ni( mp->getNodeList() );
   const double rise = rate*delt;

   if (0) //DEBUG
     cout << "****UPLIFTUNI: " << rise << endl;

   for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
   {
      cn->ChangeZ( rise );
      cn->setUplift( rate );
   }
}


/************************************************************************\
**
**  tUplift::BlockUplift
**
**  Uplift at a constant rate for all points Y>=Yf, where Yf is the
**  location of a vertical fault (striking parallel to x-axis).
**
**  Inputs:  mp -- pointer to the mesh
**           delt -- duration of uplift
**
**  Modifications:
**   - Added subsidence of area outside of block at rate given by
**     parameter 'rate2' (GT, June 2003)
**
\************************************************************************/
void tUplift::BlockUplift( tMesh<tLNode> *mp, double delt )
{
   assert( mp!=0 );
   tLNode *cn;
   tMesh<tLNode>::nodeListIter_t ni( mp->getNodeList() );
   const double rise = rate*delt,
     sink = rate2*delt;

   for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
   {
     if( cn->getY()>=faultPosition )
     {
        cn->ChangeZ( rise );
        cn->setUplift( rate );
     }
     if( rate2!=0.0 && cn->getY()<faultPosition )
     {
	 cn->ChangeZ( -sink );
	 cn->setUplift( -rate2 );
     }
   }
}


/************************************************************************\
**
**  tUplift::StrikeSlip
**
**  Mesh points at y<faultPosition move to the right relative to other
**  points (ie, left-lateral displacement). Block uplift can also
**  occur (see DoUplift). Note that boundary nodes also move.
**
**  Inputs:  mp -- pointer to the mesh
**           delt -- duration of uplift
**
\************************************************************************/
void tUplift::StrikeSlip( tMesh<tLNode> *mp, double delt ) const
{
   assert( mp!=0 );
   tLNode *cn;
   tMesh<tLNode>::nodeListIter_t ni( mp->getNodeList() );
   double slip = slipRate*delt;

   cout << "StrikeSlip by " << slip << endl;

   for( cn=ni.FirstP(); !(ni.AtEnd()); cn=ni.NextP() )
   {
     if( cn->getY()<faultPosition )
     {
        cn->setMeanderStatus( true );  // redundant: TODO
        cn->setNew2DCoords( cn->getX()+slip, cn->getY() );
     }
   }
   mp->MoveNodes( 0., false );

}

/************************************************************************\
**
**  tUplift::FoldPropErf
**
**  This function provides a very simple vertical kinematic
**  representation of fold propagation. Uplift rate is described by
**  a symmetrical error function curve, with uplift on one side
**  (Y>pivot point) and subsidence on the other (relative to fixed
**  boundary elevation). The faultPosition variable is used to track
**  the position of the pivot point, which migrates in the direction
**  of "decreasing Y" at a rate that is stored in "slipRate".
**  The variable foldParam encodes the fold wavelength, which is
**  defined as the distance between the point where uplift comes close
**  to its maximum value (where U = U0 erf(2)) and the point where
**  uplift comes close to its minimum value (where U = U0 erf(-2)).
**  More precisely, foldParam is defined as 4/wavelength.
**
**  Inputs:  mp -- pointer to the mesh
**           delt -- duration of uplift
**
\************************************************************************/
void tUplift::FoldPropErf( tMesh<tLNode> *mp, double delt )
{
   assert( mp!=0 );
   tLNode *cn;
   tMesh<tLNode>::nodeListIter_t ni( mp->getNodeList() );
   double uprate;

   // For each node, the uplift rate is the maximum rate ("rate") times
   // the error function curve, which depends on the node's Y-coordinate
   // relative to the pivot point, and on the fold wavelength (the
   // reciprocal of which is embedded in "foldParam"). Note that
   // uplift rate is negative (ie, subsidence occurs) where Y < pivot point.
   for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
   {
      uprate = rate * erf( foldParam * ( cn->getY() - faultPosition ) );
      cn->ChangeZ( uprate*delt );
   }

   // Now we "propagate" the fold by moving the pivot point forward
   // (in the direction of decreasing Y coordinate)
   faultPosition -= slipRate * delt;
}

/************************************************************************\
**
**  tUplift::CosineWarp2D
**
**
**  Inputs:  mp -- pointer to the mesh
**           delt -- duration of uplift
**
\************************************************************************/
void tUplift::CosineWarp2D( tMesh<tLNode> *mp, double delt )
{
   assert( mp!=0 );
   tLNode *cn;
   tMesh<tLNode>::nodeListIter_t ni( mp->getNodeList() );
   double uprate;
   static double elapsedTime = 0;

   // For each node, the uplift rate is the uplift rate constant ("rate") times
   // the cosine function in the y- and (if lateral y-directed tightening has
   // begun) x-directions. Note that the variable "faultPosition"
   // is used to store the Y location of the anticline peak (normally, the
   // upper boundary). The ratio of uplift to subsidence rates is controlled
   // by the parameter "foldParam2"; if uplift is positive, the rate is
   // multiplied by this factor. "positionParam1" is used to store the
   // x-location of the anticline.
   if( elapsedTime >= deformStartTime1 )
   {
      for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
      {
         uprate = rate *
             ( cos(PI*(positionParam1 - cn->getX())/foldParam)
               + cos(TWOPI*(faultPosition-cn->getY())/foldParam) );
         if( uprate>0.0 ) uprate *= foldParam2;
         cn->ChangeZ( uprate*delt );
      }
   }
   else
   {
      for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
      {
         uprate = rate *
             cos(TWOPI*(faultPosition-cn->getY())/foldParam);
         if( uprate>0.0 ) uprate *= foldParam2;
         cn->ChangeZ( uprate*delt );
      }
   }
   elapsedTime += delt;

   // The "tightening" of the folds through time is simulated by
   // progressively decreasing the fold wavelength. (Here the variable
   // "slipRate" is used to store the tightening rate in m/yr, and
   // "foldParam" is used to store the fold wavelength in m).
   // As the fold tightens, the anticline axis migrates in the direction
   // of decreasing Y (assumed to be basinward; this means the folding
   // propagates basinward as it tightens).
   foldParam = foldParam - slipRate*delt;
   faultPosition= faultPosition - slipRate*delt;

}



/************************************************************************\
**
**  tUplift::PropagatingFold
**
**  The "propagating fold" has two components: block uplift inboard of
**  a mountain front (which is handled via a separate call to
**  BlockUplift()), and a single laterally propagating fold, which is
**  handled by this function.
**
**  The fold is aligned parallel to the x-axis. The uplift pattern is
**  sinusoidal in the y-direction and symmetrical about the fold axis.
**  The nose of the fold propagates in the direction of increasing x.
**  Parameters are:
**    - foldParam: width of fold in meters
**    - foldParam2: axis location, in meters from y=0
**    - rate2: uplift rate along fold crest
**    - slipRate: lateral propagation rate in m/yr
**
**  Inputs:  mp -- pointer to the mesh
**           delt -- duration of uplift
**
\************************************************************************/
void tUplift::PropagatingFold( tMesh<tLNode> *mp, double delt ) const
{
   assert( mp!=0 );
   tLNode *cn;
   tMesh<tLNode>::nodeListIter_t ni( mp->getNodeList() );
   double uprate;
   const double northEdge = foldParam2 + 0.5*foldParam;
   const double southEdge = northEdge - foldParam;
   static double foldNose = 0.0;
   const double twoPiLam = TWOPI/foldParam;

   // Advance the fold nose
   foldNose += slipRate*delt;

   // For each node, the uplift rate is the uplift rate constant ("rate") times
   // the cosine function in y. The variable "foldParam2" is the location
   // of the fold axis in meters relative to y=0.
   for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
   {
     if( cn->getX()<=foldNose && cn->getY()<=northEdge && cn->getY()>=southEdge )
       {
	 uprate = rate2 * 0.5 * ( cos( twoPiLam*(foldParam2-cn->getY()) )+1.0);
	 cn->ChangeZ( uprate*delt );
       }
   }

}

/************************************************************************\
**
**  tUplift::TwoSideDifferential
**
**  This uplift pattern is designed to simulate a block of land bordered
**  by two active baselevels that are lowering at different rates. It
**  is designed to operate on a rectangular domain with open boundaries
**  along two (opposite) sides. The function was written to study
**  controls on drainage divide migration under conditions of
**  differential relative uplift -- for example in the case range
**  bounded on one side by a rapidly-slipping normal fault, and on the
**  other by a boundary that slowly lowers due to steady erosion.
**
**  The function uses "rate" as the baselevel fall rate along the lower
**  (y=0) boundary. "rate2" is initially read in (in the constructor)
**  as the baselevel fall rate along the upper boundary. However, in order
**  to preserve the lower boundary as a zero and fixed relative elevation
**  boundary, the differential lowering is implemented by having the
**  interior landscape rise at "rate" (lower baselevel fall) and the upper
**  boundary rise or fall at a rate equal to the difference between upper
**  and lower boundary baselevel fall rates (so the relative rates are
**  preserved; so effectively the datum is the lower baselevel). "rate2"
**  is converted in the constructor to the difference between upper and
**  lower baselevel fall rates, and can be positve or negative.
**
**  Since the user might choose any arbitrary boundary condition, the
**  parameter positionParam1 is used to store the y-coordinate of the
**  dividing line between the two base levels. Any open boundary points
**  at y<=positionParam1 are treated as "lower boundary" (ie, fixed at
**  zero elevation), while any open boundary points at y>positionParam1
**  are treated as "upper boundary".
**
**    Created July, 2001 GT
**
**    Parameters:
**      mp -- pointer to mesh object
**      delt -- time interval over which uplift occurs
**
**    Member data accessed:
**      rate, rate2, positionParam1
**
\************************************************************************/
void tUplift::TwoSideDifferential( tMesh<tLNode> *mp, double delt ) const
{
  tMesh<tLNode>::nodeListIter_t ni( mp->getNodeList() );
  tLNode *cn;
  double landraise = rate*delt,
    boundaryupdown = rate2*delt;

  // Raise the landscape
  for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
    cn->ChangeZ( landraise );

  // Raise or lower the "upper" boundary
  for( cn=ni.FirstBoundaryP(); !ni.AtEnd(); cn=ni.NextP() )
    if( cn->getBoundaryFlag()==kOpenBoundary && cn->getY()>positionParam1 )
      cn->ChangeZ( boundaryupdown );
}


/************************************************************************\
**
**  tUplift::FaultBendFold
**
**  Mesh points above ramp and flat move in positive y-direction at
**  sliprate, inclined and horizontally, respectively.  Note that
**  boundary nodes also move.  Right now, get to work if Z at
**  faultPosition (initial break at surface) is same as the mean
**  elevation--and most realistic if mean elevation is zero so we don't
**  have to worry about foreland erosion and deformation of the hanging
**  wall above this foreland topography.  If you make any inclined
**  initial topography, make Z=0 for y-coordinates where y>=faultPosition.
**  Fault strikes parallel to x-axis.
**
**  Created August 2002  (s. miller).  Updated January 2003.
**
**  Inputs:  mp -- pointer to the mesh
**           delt -- duration of uplift
**
\************************************************************************/
void tUplift::FaultBendFold( tMesh<tLNode> *mp, double delt ) const
{
   assert( mp!=0 );
   tLNode *cn;
   tMesh<tLNode>::nodeListIter_t ni( mp->getNodeList() );
   double slip = slipRate*delt;
   static double elapsedTime = delt;


   for( cn=ni.FirstP(); !(ni.AtEnd()); cn=ni.NextP() )
   {

     /* This loop handles the the initial stages of folding, specifically
        that when the active, backlimb-ward dipping upper flat
        axial surface that exists until the forward dipping passive axial
        surface, which migrates up from the lower end of the ramp, reaches
        the upper end of the ramp.  */

     if( elapsedTime < flatDepth/(sin(PI*rampDip/180))/slipRate )
     {

       /* For hinterland (lower ramp) hangingwall */

       if( cn->getY()< (faultPosition-flatDepth/tan(PI*rampDip/180)) &&
       cn->getZ()<=-tan(PI*kinkDip/180)*cn->getY()+(faultPosition-
       (flatDepth/tan(PI*rampDip/180))-(flatDepth*tan(PI*90/180-
       PI*kinkDip/180)))*tan(PI*kinkDip/180) )
       {
          cn->setMeanderStatus( true );  // redundant: TODO
          cn->setNew2DCoords(  cn->getX(), cn->getY()+slip );
       }
       /* For hangingwall above ramp */

       else if( cn->getZ()>=tan(PI*rampDip/180)*cn->getY()-tan(PI*rampDip/180)
       *faultPosition && cn->getZ()>-tan(PI*kinkDip/180)*cn->getY()+
       (faultPosition-(flatDepth/tan(PI*rampDip/180))-(flatDepth*tan(PI*90/180-
       PI*kinkDip/180)))*tan(PI*kinkDip/180) &&
       cn->getZ()>=tan(PI*upperKinkDip/180)*cn->getY()-tan(PI*upperKinkDip/180)*
       faultPosition )
       {
          cn->setMeanderStatus( true );  // redundant: TODO
          cn->setNew2DCoords( cn->getX(), cn->getY()+slip*cos(PI*rampDip/180) );
       }

       /* For foreland (upper ramp) hangingwall.  Note that only nodes
       slightly (here, 3 m) above the mean elevation are moved with the hangingwall. */

       else if( cn->getZ()<tan(PI*upperKinkDip/180)*cn->getY()-
       tan(PI*upperKinkDip/180)*faultPosition && cn->getY()>faultPosition
       && cn->getZ()>(meanElevation + 3) )
       {
          cn->setMeanderStatus( true );  // redundant: TODO
          cn->setNew2DCoords( cn->getX(), cn->getY()+slip );
       }

     }


     /* This loop handles the later stages of deformation, when both
        active axial surface dip towards the forelimb and the passive axial
        surface above the upper flat, if it is not eroded, dips towards
        the backlimb. */

     else
     {

       /* For hinterland (lower ramp) hangingwall */

       if( cn->getY()< (faultPosition-flatDepth/tan(PI*rampDip/180)) &&
       cn->getZ()<=-tan(PI*kinkDip/180)*cn->getY()+(faultPosition-
       (flatDepth/tan(PI*rampDip/180))-(flatDepth*tan(PI*90/180-
       PI*kinkDip/180)))*tan(PI*kinkDip/180) )
       {
          cn->setMeanderStatus( true );  // redundant: TODO
          cn->setNew2DCoords( cn->getX(), cn->getY()+slip );
       }

       /* For hangingwall above ramp */

       else if( cn->getZ()>=tan(PI*rampDip/180)*cn->getY()-tan(PI*rampDip/180)
       *faultPosition && cn->getZ()>-tan(PI*kinkDip/180)*cn->getY()+
       (faultPosition-(flatDepth/tan(PI*rampDip/180))-(flatDepth*tan(PI*90/180-
       PI*kinkDip/180)))*tan(PI*kinkDip/180) &&
       cn->getZ()<=-tan(PI*kinkDip/180)*cn->getY()+tan(PI*kinkDip/180)*
       faultPosition )
       {
          cn->setMeanderStatus( true );  // redundant: TODO
          cn->setNew2DCoords( cn->getX(), cn->getY()+slip*cos(PI*rampDip/180) );
          /* cn->ChangeZ( slip*sin(PI*rampDip/180) );  */
          /* cn->setUplift( slipRate*sin(PI*rampDip/180) );  */
       }

       /* For foreland (upper ramp) hangingwall.  Note only nodes slightly above
       (here, 3 m) the mean elevation are moved. */

       else if( cn->getZ()>-tan(PI*kinkDip/180)*cn->getY()+
       tan(PI*kinkDip/180)*faultPosition && cn->getY()<=faultPosition  &&
       cn->getZ()>(meanElevation + 3) )
       {
          cn->setMeanderStatus( true );  // redundant: TODO
          cn->setNew2DCoords( cn->getX(), cn->getY()+slip );
       }

     }

   }
   mp->MoveNodes( 0., false );


   /*elapsedTime += delt;*/
}


/************************************************************************\
**
**  tUplift::FaultBendFold2
**
**  This is the uplift part of the fault-bend fold uplift function.  It
**  is separate so that the boundary nodes are pinned.
**
**  Created September, 2002  (s. miller)
**
**  Inputs:  mp -- pointer to the mesh
**           delt -- duration of uplift
**
\************************************************************************/
void tUplift::FaultBendFold2( tMesh<tLNode> *mp, double delt ) const
{
   assert( mp!=0 );
   tLNode *cn;
   tMesh<tLNode>::nodeListIter_t ni( mp->getNodeList() );
   double slip = slipRate*delt;
   static double elapsedTime = delt;
   /* Redefinitions so faultPosition and flatDepth are measured with respect
   to where fault intersects z=0 (faultPosition) and depth below z=0. */
   //double faultPosition = faultPosition - meanElevation/tan(rampDip);
   //double flatDepth = meanElevation - flatDepth;

   for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
   {
      if( elapsedTime < flatDepth/(sin(PI*rampDip/180))/slipRate )
      {
          if( cn->getZ()>=tan(PI*rampDip/180)*cn->getY()-tan(PI*rampDip/180)
          *faultPosition && cn->getZ()>-tan(PI*kinkDip/180)*cn->getY()+
          (faultPosition-(flatDepth/tan(PI*rampDip/180))-(flatDepth*tan(PI*90/180-
          PI*kinkDip/180)))*tan(PI*kinkDip/180) &&
          cn->getZ()>=tan(PI*upperKinkDip/180)*cn->getY()-tan(PI*upperKinkDip/180)*
          faultPosition )
          {
             cn->ChangeZ( slip*sin(PI*rampDip/180) );
             cn->setUplift( slipRate*sin(PI*rampDip/180) );
          }
      }

      else
      {
          if( cn->getZ()>=tan(PI*rampDip/180)*cn->getY()-tan(PI*rampDip/180)
          *faultPosition && cn->getZ()>-tan(PI*kinkDip/180)*cn->getY()+
          (faultPosition-(flatDepth/tan(PI*rampDip/180))-(flatDepth*tan(PI*90/180-
          PI*kinkDip/180)))*tan(PI*kinkDip/180) &&
          cn->getZ()<=-tan(PI*kinkDip/180)*cn->getY()+tan(PI*kinkDip/180)*
          faultPosition )
          {
             cn->ChangeZ( slip*sin(PI*rampDip/180) );
             cn->setUplift( slipRate*sin(PI*rampDip/180) );
          }
      }
   }

   elapsedTime += delt;
}



/************************************************************************\
**
**  tUplift "get" functions
**
\************************************************************************/

double tUplift::getDuration() const
{
   return duration;
}

double tUplift::getRate() const
{
   return rate;
}

