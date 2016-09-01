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
**  $Id: tRunTimer.h,v 1.15 2004/06/16 13:37:42 childcvs Exp $
*/
/***************************************************************************/

#ifndef TRUNTIMER_H
#define TRUNTIMER_H


class tRunTimer
{
public:
	tRunTimer( double duration, double opint, bool optprint=true );
	tRunTimer( const tInputFile &infile, bool optprint=true );
	tRunTimer();
	double getCurrentTime() const;     // Report the current time
	bool Advance( double );            // Advance time by given amount
	bool IsFinished() const;           // Are we done yet?
	double RemainingTime() const;      // How much time is left
	void Start( double, double=0.0 );  // Set current and (optionally) end times
	bool CheckOutputTime();             // Is it time to write output yet?
	void ReportTimeStatus();           // Report time to file and (opt) screen
	bool CheckTSOutputTime();           // Is it time to write time series output yet?

private:
	std::ofstream timeStatusFile;  // file "run.time" for tracking current time
	double currentTime;       // current time in simulation
	double endTime;           // time at which simulation ends
	double outputInterval;    // interval between outputs
	double nextOutputTime;    // time of next output
	double notifyInterval;    // interval for reporting time to file
	double nextNotify;        // next time for time-reporting
	double nextTSOutputTime;  // time of next time series output
	double TSOutputInterval;  // interval between time series outputs
	bool optTSOutput;          // option for time series output
	const bool optPrintEachTime;     // option for reporting time to screen
};


#endif
