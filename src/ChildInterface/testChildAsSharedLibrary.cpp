/*
** testChildAsSharedLibrary.cpp
**
** This is designed to test the ability to use CHILD as a shared library.
**
** GT, Aug '11
*/

#include <iostream>
#include "../ChildInterface/childInterface.h"
#include "../Inclusions.h"


int main( int argc, char **argv )
{
  childInterface myChild;

  myChild.Initialize( argc, argv );

  return 0;
}

