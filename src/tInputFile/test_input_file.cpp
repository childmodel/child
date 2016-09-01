/*
 *  test_input_file.cpp
 *  
 *  This file is just a little driver file designed to test the tInputFile
 *  class. I wrote it to test the new ReadString function, but it can be
 *  modified for other sorts of tests. See notes below on how to use it.
 *  To compile it, use:
 *
 *  g++ -o test_input_file test_input_file.cpp tInputFile.cpp 
 *                 ../tTimeSeries/tTimeSeries.cpp ../errors/errors.cpp
 *
 *  Created by Greg Tucker on 10/18/09.
 *
 */

#include "tInputFile.h"
#include <iostream>
#include <sstream>  // just for testing stringstream functions

using namespace std;

int main( int argc, char** argv )
{
  
  if( argc<2 )
  {
    cout << "Must include name of input file on the command line\n";
    exit(0);
  }
  
  tInputFile inputFile( argv[1] );

  bool mybool;
  mybool = inputFile.ReadBool( "MY_TRUE_BOOL" );
  if( mybool )
    cout << "MY_TRUE_BOOL is true!" << endl;
  else
    cout << "MY_TRUE_BOOL is false!" << endl;
  mybool = inputFile.ReadBool( "MY_FALSE_BOOL" );
  if( mybool )
    cout << "MY_FALSE_BOOL is true!" << endl;
  else
    cout << "MY_FALSE_BOOL is false!" << endl;
  mybool = inputFile.ReadBool( "MY_MISSING_BOOL", false );
  if( mybool )
    cout << "MY_MISSING_BOOL is true!" << endl;
  else
    cout << "MY_MISSING_BOOL is false, which is what it defaults to when it is missing!" << endl;

  int myint;
  myint = inputFile.ReadInt( "MY_INTEGER" );
  cout << "MY_INTEGER=" << myint << endl;
  myint = inputFile.ReadInt( "MY_MISSING_INTEGER", false );
  cout << "MY_MISSING_INTEGER=" << myint << endl;
  
  long mylong;
  mylong = inputFile.ReadLong( "MY_LONG" );
  cout << "MY_LONG=" << mylong << endl;
  mylong = inputFile.ReadInt( "MY_MISSING_LONG", false );
  cout << "MY_MISSING_LONG=" << mylong << endl;
  
  double mydouble;
  mydouble = inputFile.ReadDouble( "MY_DOUBLE" );
  cout << "MY_DOUBLE=" << mydouble << endl;
  mydouble = inputFile.ReadDouble( "MY_MISSING_DOUBLE", false );
  cout << "MY_MISSING_DOUBLE=" << mydouble << endl;
  
  string mystring;
  mystring = inputFile.ReadString( "MY_STRING" );
  cout << "MY_STRING='" << mystring << "'" << endl;
  mystring = inputFile.ReadString( "MY_MISSING_STRING", false );
  cout << "MY_MISSING_STRING='" << mystring << "'" << endl;
  
  mystring = inputFile.ReadString( "TEST_STRINGSTREAM" );
  stringstream ss;
  ss.str( mystring );
  for( int i=0; i<8; i++ )
  {
    if( ss.eof() )
      cout << "Warning: we've hit the end of the string" << endl;
    ss >> mydouble;
    cout << "mydouble=" << mydouble << endl;
  }
  
  exit(0);

}


/*

This test script is designed to be run with an input file which looks like this:

#
# test_input_file.in
#
# This is a template for a generic input file.
#
# Its purpose is to test the tInputFile class. It contains examples of each of the main
# input types, as well as comment lines.
#
# Created by GT, Oct 09
#
OUTFILENAME: need this so the tInputFile knows how to write an output log
my_test_input_file
MY_TRUE_BOOL
1
MY_FALSE_BOOL
0
#MY_MISSING_BOOL
#1
#because "MY_MISSING_BOOL" is commented out, it shouldn't be read
MY_INTEGER
12345
# random comment line here ...
MY_LONG
123456789
MY_DOUBLE
3.1415926
MY_STRING
To be or not to be: is that the question?
# some other missing things:
# MY_MISSING_INTEGER
# 54321
# MY_MISSING_DOUBLE
# 0.54321
# MY_MISSING_STRING
# Is that a dagger before me?


The output should look like this:

MY_TRUE_BOOL is true!
MY_FALSE_BOOL is false!
Cannot find  'MY_MISSING_BOOL' in the input file.
WARNING: Missing parameter in input file, use zero or default value
(Set "CHILD_ABORT_ON_WARNING" to generate a crash.
.e.g. "env CHILD_ABORT_ON_WARNING=1 child ...")
MY_MISSING_BOOL is false, which is what it defaults to when it is missing!
MY_INTEGER=12345
Cannot find  'MY_MISSING_INTEGER' in the input file.
WARNING: Missing parameter in input file, use zero or default value
(Set "CHILD_ABORT_ON_WARNING" to generate a crash.
.e.g. "env CHILD_ABORT_ON_WARNING=1 child ...")
MY_MISSING_INTEGER=0
MY_LONG=123456789
Cannot find  'MY_MISSING_LONG' in the input file.
WARNING: Missing parameter in input file, use zero or default value
(Set "CHILD_ABORT_ON_WARNING" to generate a crash.
.e.g. "env CHILD_ABORT_ON_WARNING=1 child ...")
MY_MISSING_LONG=0
MY_DOUBLE=3.14159
Cannot find  'MY_MISSING_DOUBLE' in the input file.
WARNING: Missing parameter in input file, use zero or default value
(Set "CHILD_ABORT_ON_WARNING" to generate a crash.
.e.g. "env CHILD_ABORT_ON_WARNING=1 child ...")
MY_MISSING_DOUBLE=0
MY_STRING='To be or not to be: is that the question?'
Cannot find  'MY_MISSING_STRING' in the input file.
WARNING: Missing parameter in input file, use zero or default value
(Set "CHILD_ABORT_ON_WARNING" to generate a crash.
.e.g. "env CHILD_ABORT_ON_WARNING=1 child ...")
MY_MISSING_STRING=''

*/