#include <iostream.h>

#include "readChildData.h"

int main(void){
  ReadChildData ChildData;
  int nStep = 2;
  const int TypeVariable = 0;
  //const char* basename = "/miami/academic/arnaud/ForGreg/Child/Example/junk";
  const char* basename = "/miami/academic/arnaud/ForGreg/Child/Example/example";

  try {
    if (!ChildData.LoadData(basename, nStep, TypeVariable)) {
      cerr << "Time step \"" << nStep << "\" does not exist." << endl;
      return 0;
    }
  } catch (const BadFile &ex) {
    cerr << "Error when reading file \"" << ex.filename()
	 << "\": " << ex.ExceptionName()
	 << "." << endl;
    return 0;
  } catch (...) {
    cerr << "A unkown exception occured during file processing." << endl;
    return 0;
  }
  
  return 0;
}
