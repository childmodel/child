/***************************************************************************\
**
**  tRunTimer.h: Header for tRunTimer objects.
**
**  tRunTimer objects are used to keep track of time in a time-evolving
**  simulation model. Their services include keeping track of when it's
**  time to write output, printing the current time to standard output if
**  desired, and writing the current time to a file every so often.
**
**  $Id: tRunTimer.h,v 1.1 1998-01-14 20:44:56 gtucker Exp $
\***************************************************************************/

#ifndef TRUNTIMER_H
#define TRUNTIMER_H

class tRunTimer
{
public:
	tRunTimer( float duration, float opint, int optprint=1 );
	tRunTimer( tInputFile &infile, int optprint=1 );
	tRunTimer();
	float GetCurrentTime();
	int Advance( float );
	int IsFinished();
	float RemainingTime();
	void Start( float, float );
	int CheckOutputTime();
	void ReportTimeStatus();
	
private:
	float currentTime;
	float endTime;
	float outputInterval;
	float nextOutputTime;
	int optPrintEachTime;
	ofstream timeStatusFile;
	float notifyInterval;
	float nextNotify;
};


#endif
