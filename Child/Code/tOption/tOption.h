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
**  $Id: tOption.h,v 1.4 2004-04-16 18:33:54 childcvs Exp $
*/
/************************************************************************/

#ifndef TOPTION_H
#define TOPTION_H

#include <vector>

using namespace std;

class tOption {
  char const * const exeName;
public:
  bool silent_mode;      // Option for silent mode (no time output to stdout)
  bool checkMeshConsistency;
  bool no_write_mode; // option to force no writing to files
  char const *inputFile;

  tOption(int argc, char const * const argv[]);
  tOption(string arguments);
  static void version();
private:
  int parseOptions(char const * const argv[]);
  int parseOptions(std::string thisOption);
  void usage() const;
private:
  tOption();
  tOption(tOption const &);
  tOption& operator=(tOption const &);
};

#endif








