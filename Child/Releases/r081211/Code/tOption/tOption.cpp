/***************************************************************************/
/**
**  @File Option.cpp
**  @brief Functions for class tOption
**
**  A. Desitter - March 2004
**
**  $Id: tOption.cpp,v 1.5 2004/06/16 13:37:39 childcvs Exp $
*/
/***************************************************************************/

#include <iostream>
#include <stdlib.h>
#include <string.h>

#include "tOption.h"
#include "../Definitions.h"

tOption::tOption(int argc, char const * const argv[])
  : exeName(argv[0]),
    silent_mode(false), checkMeshConsistency(true),
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

// Parse options one at a time. Returns number of options consumed.
int tOption::parseOptions(char const * const argv[]) {

  const char * const thisOption = argv[0];

  if (strcmp(thisOption, "--silent-mode") == 0){
    silent_mode = true;
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

void tOption::usage() const {
  std::cerr
    << "Usage: " << exeName << " [options] <input file>\n"
    << " --help: display this help message.\n"
    << " --no-check: disable CheckMeshConsistency().\n"
    << " --silent-mode: silent mode.\n"
    << " --version: display version.\n"
    << std::endl;
}

void tOption::version() {
  std::cout
    << "\nThis is CHILD, version " << VERSION
    << " (compiled " __DATE__ " " __TIME__ ")"
    << '\n' << std::endl;
}
