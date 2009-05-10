/*
 *  tIDGenerator.cpp
 *  tTracer
 *
 *  Created by D. Nathan Bradley on 9/12/06.
 *  Modifications:
 *   - GT added default constructor that initializes ID to zero (5/09)
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tIDGenerator.h"

tIDGenerator::tIDGenerator()
{
    id = 0;
}

tIDGenerator::tIDGenerator(int startingValue)
{
	id = startingValue;
}
	
int tIDGenerator::getNextID()
{
	return id++;
}
