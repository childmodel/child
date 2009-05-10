/*
 *  tIDGenerator.h
 *  tTracer
 *
 *  Created by D. Nathan Bradley on 9/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef TIDGENERATOR_H
#define TIDGENERATOR_H

class tIDGenerator
{
	public:
	  tIDGenerator();
	  tIDGenerator(int startingValue);
	  int getNextID();
	  
	private:
	   int id;
 };

#endif


