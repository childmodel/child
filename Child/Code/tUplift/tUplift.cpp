/************************************************************************\
**
**  tUplift.cpp: Functions for class tUplift (see tUplift.h).
**
**  $Id: tUplift.cpp,v 1.5 1999-04-01 19:09:32 gtucker Exp $
\************************************************************************/

#include "tUplift.h"
#include "../errors/errors.h"

#define kNumUpliftTypes 2
#define kNoUplift 0


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
tUplift::tUplift( tInputFile &infile )
{
   // Find out what kind of uplift the user wants
   typeCode = infile.ReadItem( typeCode, "UPTYPE" );
   if( typeCode < 0 || typeCode > kNumUpliftTypes )
   {
      cerr << "I don't recognize the uplift type you asked for ("
           << typeCode << ")\n";
      cerr << "Valid uplift types are:\n";
      cerr << " 0 - none\n"
           << " 1 - Spatially and temporally uniform uplift\n"
           << " 2 - Uniform uplift at Y >= fault location, zero elsewhere\n";
      ReportFatalError( "Please specify a valid uplift type and try again." );
   }

   if( typeCode==kNoUplift ) return;
   
   // get the parameters relevant to that type
   duration = infile.ReadItem( duration, "UPDUR" );
   rate = infile.ReadItem( rate, "UPRATE" );
   switch( typeCode )
   {
   case 2:
     faultPosition = infile.ReadItem( faultPosition, "FAULTPOS" );
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
**  Inputs:  gp -- pointer to the mesh
**           delt -- duration of uplift
**
\************************************************************************/
void tUplift::DoUplift( tGrid<tLNode> *gp, double delt )
{
   switch( typeCode )
   {
      case 1:
          UpliftUniform( gp, delt );
          break;
      case 2:
          BlockUplift( gp, delt );
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
**  Inputs:  gp -- pointer to the mesh
**           delt -- duration of uplift
**
\************************************************************************/
void tUplift::UpliftUniform( tGrid<tLNode> *gp, double delt )
{
   assert( gp>0 );
   tLNode *cn;
   tGridListIter<tLNode> ni( gp->getNodeList() );
   double rise = rate*delt;

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
**  Inputs:  gp -- pointer to the mesh
**           delt -- duration of uplift
**
\************************************************************************/
void tUplift::BlockUplift( tGrid<tLNode> *gp, double delt )
{
   assert( gp>0 );
   tLNode *cn;
   tGridListIter<tLNode> ni( gp->getNodeList() );
   double rise = rate*delt;

   for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
     if( cn->getY()>=faultPosition )
     {
	cn->ChangeZ( rise );
	cn->setUplift( rate );
     }
}


/************************************************************************\
**
**  tUplift "get" functions
**
\************************************************************************/

double tUplift::getDuration() 
{
   return duration;
}

double tUplift::getRate() const 
{
   return rate;
}

