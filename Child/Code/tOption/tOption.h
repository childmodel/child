//-*-c++-*-

/************************************************************************/
/**
**  @file tOption.h
**  @brief Header file for tOption objects.
**
**  tOption parse the command line.
**
**  A. Desitter - March 2004
**
**  $Id: tOption.h,v 1.1 2004-03-31 17:53:52 childcvs Exp $
*/
/************************************************************************/

#ifndef TOPTION_H
#define TOPTION_H

class tOption {
  char const * const exeName;
public:
  bool silent_mode;      // Option for silent mode (no time output to stdout)
  bool checkMeshConsistency;
  char const *inputFile;

  tOption(int argc, char const * const argv[]);
private:
  int parseOptions(char const * const argv[]);
  void usage();
};

#endif








