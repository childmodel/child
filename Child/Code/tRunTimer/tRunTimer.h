/***************************************************************************\
**
**  tRunTimer.h: Header for tRunTimer objects.
**
**  tRunTimer objects are used to keep track of time in a time-evolving
**  simulation model. Their services include keeping track of when it's
**  time to write output, printing the current time to standard output if
**  desired, and writing the current time to a file every so often.
**
**  $Id: tRunTimer.h,v 1.3 1998-06-04 21:27:16 gtucker Exp $
\***************************************************************************/

#ifndef TRUNTIMER_H
#define TRUNTIMER_H


class tRunTimer
{
public:
	tRunTimer( double duration, double opint, int optprint=1 );
	tRunTimer( tInputFile &infile, int optprint=1 );
	tRunTimer();
	double getCurrentTime();
	int Advance( double );
	int IsFinished();
	double RemainingTime();
	void Start( double, double );
	int CheckOutputTime();
	void ReportTimeStatus();
	
private:
	double currentTime;
	double endTime;
	double outputInterval;
	double nextOutputTime;
	int optPrintEachTime;
	ofstream timeStatusFile;
	double notifyInterval;
	double nextNotify;
};


#endif
