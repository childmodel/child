/************************************************************************\
**
**  tUplift.cpp: Functions for class tUplift.
**
**  $Id: tUplift.cpp,v 1.1 1998-02-02 20:40:35 gtucker Exp $
\************************************************************************/

#include "tUplift.h"
#include "../errors/errors.h"

#define kNumUpliftTypes 1
#define kNoUplift 0

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
           << " 1 - Spatially and temporally uniform uplift\n";
      ReportFatalError( "Please specify a valid uplift type and try again." );
   }

   if( typeCode==kNoUplift ) return;
   
   // Get the parameters relevant to that type
   duration = infile.ReadItem( duration, "UPDUR" );
   rate = infile.ReadItem( rate, "UPRATE" );
   switch( typeCode )
   {
   }
   
}


void tUplift::DoUplift( tGrid<tLNode> *gp, double delt )
{
   switch( typeCode )
   {
      case 1:
          UpliftUniform( gp, delt );
          break;
   }
   
}


void tUplift::UpliftUniform( tGrid<tLNode> *gp, double delt )
{
   assert( gp>0 );
   tLNode *cn;
   tGridListIter<tLNode> ni( gp->GetNodeList() );
   double rise = rate*delt;

   for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
       cn->ChangeZ( rise );
   
}



double tUplift::GetDuration() 
{
   return duration;
}

