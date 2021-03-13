/***************************************************************************/
/**
**  @File Option.cpp
**  @brief Functions for class tOption
**
**  A. Desitter - March 2004
**
**  $Id: tOption.cpp,v 1.5 2004-06-16 13:37:39 childcvs Exp $
*/
/***************************************************************************/

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <sstream>

#include "tOption.h"
#include "../Definitions.h"

tOption::tOption(int argc, char const * const argv[])
  : exeName(argv[0]),
    silent_mode(false), checkMeshConsistency(true), no_write_mode(false), 
    inputFile(0)
{
  argv++;
  while(argc > 1){
    int i = parseOptions(argv);
    argc -= i;
    argv += i;
  }
  if (inputFile ==NULL){
    usage();
    exit(EXIT_FAILURE);
  }
}

tOption::tOption(string arguments)
: exeName("child"),
silent_mode(false), checkMeshConsistency(true), no_write_mode(false), 
inputFile(0)
{
  ProcessOptionsFromString( arguments );
}

tOption::tOption(const char * args)
: exeName("child"),
silent_mode(false), checkMeshConsistency(true), no_write_mode(false), 
inputFile(0)
{
  string arg_string( args );
  ProcessOptionsFromString( arg_string );
}

void tOption::ProcessOptionsFromString( string arguments )
{
  istringstream ss( arguments );
  string s;
  while( ss>>s )
    parseOptions( s );
  if (inputFile ==NULL){
    usage();
    exit(EXIT_FAILURE);
  }
}


// Parse options one at a time. Returns number of options consumed.
int tOption::parseOptions(char const * const argv[]) {
	
  const char * const thisOption = argv[0];
	
  if (strcmp(thisOption, "--silent-mode") == 0){
    silent_mode = true;
    return 1;
  }
  if (strcmp(thisOption, "--no-write-mode") == 0){
    no_write_mode = true;
    return 1;
  }
  if (strcmp(thisOption, "--no-check") == 0){
    checkMeshConsistency = false;
    return 1;
  }
  if (strcmp(thisOption, "--help") == 0){
    usage();
    exit(EXIT_SUCCESS);
  }
  if (strcmp(thisOption, "--version") == 0){
    version();
    exit(EXIT_SUCCESS);
  }
  if (thisOption[0] == '-') {
    usage();
    exit(EXIT_FAILURE);
  }
  if (inputFile != NULL)  {
    std::cerr << exeName << ": Several input files given." << std::endl;
    exit(EXIT_FAILURE);
  }
	
  inputFile = thisOption;
  return 1;
}

// Parse options one at a time. Returns number of options consumed.
int tOption::parseOptions(string thisOption) {
	
  if (thisOption.compare("--silent-mode") == 0){
    silent_mode = true;
    return 1;
  }
  if (thisOption.compare("--no-write-mode") == 0){
    no_write_mode = true;
    return 1;
  }
  if (thisOption.compare("--no-check") == 0){
    checkMeshConsistency = false;
    return 1;
  }
  if (thisOption.compare("--help") == 0){
    usage();
    exit(EXIT_SUCCESS);
  }
  if (thisOption.compare("--version") == 0){
    version();
    exit(EXIT_SUCCESS);
  }
  if (thisOption[0] == '-') {
    usage();
    exit(EXIT_FAILURE);
  }
  if (inputFile != NULL)  {
    std::cerr << exeName << ": Several input files given." << std::endl;
    exit(EXIT_FAILURE);
  }
	
  inputFileString = thisOption;
  inputFile = inputFileString.c_str();
  return 1;
}

void tOption::usage() const {
  std::cerr
    << "Usage: " << exeName << " [options] <input file>\n"
    << " --help: display this help message.\n"
    << " --no-check: disable CheckMeshConsistency().\n"
    << " --silent-mode: silent mode.\n"
    << " --no-write-mode: no writing to output.\n"
    << " --version: display version.\n"
    << std::endl;
}

void tOption::version() {
  std::cout
    << "\nThis is CHILD, version " << CHILD_VERSION
    << " (compiled " __DATE__ " " __TIME__ ")"
    << '\n' << std::endl;
}
