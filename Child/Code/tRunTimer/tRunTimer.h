//-*-c++-*- 

/***************************************************************************/
/**
**  @file tRunTimer.h
**  @brief Header for tRunTimer objects.
**
**  tRunTimer objects are used to keep track of time in a time-evolving
**  simulation model. Their services include keeping track of when it's
**  time to write output, printing the current time to standard output if
**  desired, and writing the current time to a file every so often.
**
**  $Id: tRunTimer.h,v 1.11 2003-05-23 11:47:13 childcvs Exp $
*/
/***************************************************************************/

#ifndef TRUNTIMER_H
#define TRUNTIMER_H


class tRunTimer
{
public:
	tRunTimer( double duration, double opint, bool optprint=true );
	tRunTimer( tInputFile &infile, bool optprint=true );
	tRunTimer();
	double getCurrentTime() const;     // Report the current time
	int Advance( double );             // Advance time by given amount
	int IsFinished() const;            // Are we done yet?
	double RemainingTime() const;      // How much time is left
	void Start( double, double=0.0 );  // Set current and (optionally) end times
	int CheckOutputTime();             // Is it time to write output yet?
	void ReportTimeStatus();           // Report time to file and (opt) screen
	int CheckTSOutputTime();           // Is it time to write time series output yet?
	
private:
	double currentTime;       // current time in simulation
	double endTime;           // time at which simulation ends
	double outputInterval;    // interval between outputs
	double nextOutputTime;    // time of next output
	const bool optPrintEachTime;     // option for reporting time to screen
	ofstream timeStatusFile;  // file "run.time" for tracking current time
	double notifyInterval;    // interval for reporting time to file
	double nextNotify;        // next time for time-reporting
	double nextTSOutputTime;  // time of next time series output
	double TSOutputInterval;  // interval between time series outputs
	int optTSOutput;          // option for time series output
};


#endif
