//-*-c++-*-

//--------------------------------------------------------------
// @file tTimeSeries.h
//
// This file implements temporal variable parameters into CHILD.
// See the file TimeSeries.pdf for info how to use it.
//
// Copyright 2000 P. W. Bogaart / VUA
// January 2004: Cleaned up and integrated in Child (AD)
//--------------------------------------------------------------

#ifndef __tTimeSeries_h
#define __tTimeSeries_h

class tTimeSeriesImp;

/** @class
**
**  Hold a time series
**  Use "bridge" design pattern.
*/
class tTimeSeries {
  tTimeSeriesImp *ts;
  // Clunky, but need tag to find out which kind of time series
  // "Imp" to instantiate in copy constructor:
  int tagImp;
public:
  tTimeSeries();
  tTimeSeries(tTimeSeries const &);
  ~tTimeSeries();
  void configure(const char *s);
  void reconfigure(const char *s);  // GT added June 09
  double calc(double time) const;
  void TellAll() const;
private:
  tTimeSeries& operator=(tTimeSeries const &);
};

#endif
//eof
