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
public:
  tTimeSeries();
  ~tTimeSeries();
  void configure(const char *s);
  double calc(double time) const;
  void TellAll() const;
private:
  tTimeSeries(tTimeSeries const &);
  tTimeSeries& operator=(tTimeSeries const &);
};

#endif
//eof
