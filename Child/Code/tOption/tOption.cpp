/***************************************************************************/
/**
**  @file Option.cpp
**  @brief Functions for class tOption
**
**  A. Desitter - March 2004
**
**  $Id: tOption.cpp,v 1.1 2004-03-31 17:53:52 childcvs Exp $
*/
/***************************************************************************/

#if !defined(HAVE_NO_NAMESPACE)
# include <iostream>
using namespace std;
#else
# include <iostream.h>
#endif
#include <stdlib.h>

#include "tOption.h"

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
  if (thisOption[0] == '-') {
    usage();
    exit(EXIT_FAILURE);
  }
  if (inputFile != NULL)  {
    cerr << exeName << ": Several input files given." << endl;
    exit(EXIT_FAILURE);
  }

  inputFile = thisOption;
  return 1;
}

void tOption::usage(){
  cerr
    << "Usage: " << exeName << " [options] <input file>\n"
    << " --help: display this help message.\n"
    << " --no-check: disable CheckMeshConsistency().\n"
    << " --silent-mode: silent mode.\n"
    << endl;
}
