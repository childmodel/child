/***************************************************************************\
**
**  tRunTimer.cpp: functions for tRunTimer objects.
**
**  tRunTimer objects are used to keep track of time in a time-evolving
**  simulation model. Their services include keeping track of when it's
**  time to write output, printing the current time to standard output if
**  desired, and writing the current time to a file every so often.
** 
**  Version 1.0, Greg Tucker, November 1997
**
**  Potential additions/improvements:
**  - add functions to set output interval and time status notification
**    interval
**
**  $Id: tRunTimer.cpp,v 1.3 1998-06-04 21:27:09 gtucker Exp $
\***************************************************************************/

#include <iostream.h>
#include <fstream.h>
#include <assert.h>

#include "../tInputFile/tInputFile.h"

#include "tRunTimer.h"
//****************************************************
// Constructors
//
// - A default constructor
// - A constructor that sets the run duration and
//   output interval (and if desired sets the option
//   to print time steps to stdout; default is 1,
//   see header file)
// - A constructor that reads run duration and
//   output interval from a tInputFile object (and
//   if desired sets the option
//   to print time steps to stdout; default is 1,
//   see header file)
//
// Note that the notifyInterval is always initialized
// at 1000, but of course this could be changed.
//****************************************************

tRunTimer::tRunTimer()
{
	currentTime = 0;
	endTime = 1;
	outputInterval = 1;
	nextOutputTime = 0;
	optPrintEachTime = 1;
	notifyInterval = 1000;
	nextNotify = 0;
}

tRunTimer::tRunTimer( double duration, double opint, int optprint )
{
	currentTime = 0;
	endTime = duration;
	outputInterval = opint;
	nextOutputTime = 0;
	optPrintEachTime = optprint;
	notifyInterval = 1000;
	nextNotify = 0;
}

tRunTimer::tRunTimer( tInputFile &infile, int optprint )
{
	currentTime = 0;
	endTime = infile.ReadItem( endTime, "RUNTIME" );
	outputInterval = infile.ReadItem( outputInterval, "OPINTRVL" );
	nextOutputTime = 0;
	optPrintEachTime = optprint;
	notifyInterval = 1000;
	nextNotify = 0;
}


//****************************************************
// Start
//
// Sets the current time and run duration.
//****************************************************
void tRunTimer::Start( double start, double end )
{
   assert( end > start );
   currentTime = start;
   endTime = end;
}

//****************************************************
// getCurrentTime
//
// Returns the current time.
//****************************************************
double tRunTimer::getCurrentTime()
{
	return currentTime;
}

//****************************************************
// RemainingTime
//
// Returns the remaining time.
//****************************************************
double tRunTimer::RemainingTime()
{
	return endTime - currentTime;
}

//****************************************************
// Advance
//
// Increments the current time by dt.
// Returns 1 if there is still time
// remaining, 0 if the time is up.
//****************************************************
int tRunTimer::Advance( double dt )
{
	currentTime += dt;
	return( currentTime < endTime );
}

//****************************************************
// ReportTimeStatus
//
// Reports the current time to a file and to 
// standard output if desired. Output to file is
// only done at selected intervals.
// Compares currentTime with nextNotify to see
// whether it's time to update the time status file,
// and if so opens the file and writes the current
// time (overwriting any previous contents). 
// Regardless of the status of nextNotify, if
// optPrintEachTime is selected, the current time is
// reported to cout.
//****************************************************
void tRunTimer::ReportTimeStatus()
{
	if( optPrintEachTime ) cout << currentTime << endl;
	if( currentTime >= nextNotify )
	{
		timeStatusFile.open( "run.time" );
		assert( timeStatusFile.good() );
		timeStatusFile << currentTime << endl;
		timeStatusFile.close();
		nextNotify += notifyInterval;
	}
}

//*************************************************
// CheckOutputTime
//
// Checks to see whether it's time to write output
// yet.
//*************************************************
int tRunTimer::CheckOutputTime()
{
	if( currentTime>=nextOutputTime )
	{
		nextOutputTime += outputInterval;
		return 1;
	}
	else return 0;
}

//*************************************************
// IsFinished
//
// Checks to see whether the run is finished.
//*************************************************
int tRunTimer::IsFinished()
{
	return( currentTime >= endTime );
}










