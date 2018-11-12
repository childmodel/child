#ifndef CHILD_H_
#define CHILD_H_

#include "bmi.h"
#include "child.h"


class Model : public bmi::Bmi {
 public:
  Model ();

  // Model control functions.
  void Initialize (const char *);
  void Update ();
  void UpdateUntil (double);
  void Finalize ();

  // Model information functions.
  void GetComponentName(char * const name);
  int GetInputVarNameCount (void);
  int GetOutputVarNameCount (void);
  void GetInputVarNames (char * const * const names);
  void GetOutputVarNames (char * const * const names);

  // Variable information functions
  int GetVarGrid (const char * var_name);
  void GetVarType (const char * var_name, char * const vtype);
  void GetVarUnits (const char * var_name, char * const units);
  int GetVarItemsize(const char * name);
  int GetVarNbytes(const char * name);
  void GetVarLocation(const char * name, char * const location);

  double GetCurrentTime ();
  double GetStartTime ();
  double GetEndTime ();
  double GetTimeStep ();
  void GetTimeUnits (char * const units);

  // Variable getters
  void GetValue (const char * var_name, void *);
  // double * GetValuePtr (const char * var_name) {
  //   throw bmi::FUNCTION_NOT_IMPLEMENTED;
  // }
  // void GetValueStride (const char * var_name, int * const stride) {
  //   throw bmi::FUNCTION_NOT_IMPLEMENTED;
  // }

  // Variable setters
  void SetValue (const char *, void *);

  // Grid information functions
  void GetGridType (const int grid_id, char * const gtype);
  int GetGridRank (const int grid_id);
  int GetGridSize (const int grid_id);

  // void GetGridShape (const int, int *) {
  //   throw bmi::FUNCTION_NOT_IMPLEMENTED;
  // }
  // void GetGridSpacing (const int, double *) {
  //   throw bmi::FUNCTION_NOT_IMPLEMENTED;
  // }
  // void GetGridOrigin (const int, double *) {
  //   throw bmi::FUNCTION_NOT_IMPLEMENTED;
  // }

  void GetGridX (const int, double * const);
  void GetGridY (const int, double * const);
  // void GetGridZ (const int, double * const) {
  //   throw bmi::FUNCTION_NOT_IMPLEMENTED;
  // }

  int GetGridFaceCount(const int);
  int GetGridVertexCount(const int);
  void GetGridConnectivity (const int, int * const );
  void GetGridOffset (const int, int * const);

  int GetGridNumberOfEdges(const int) { return -1; }
  int GetGridNumberOfFaces(const int);

  void GetGridEdgeNodes(const int, int * const) { }
  void GetGridFaceEdges(const int, int * const) { }
  void GetGridFaceNodes(const int, int * const);
  void GetGridNodesPerFace(const int, int * const);

private:
  int input_var_name_count;
  int output_var_name_count;

  char **input_var_names;
  char **output_var_names;

  Child model;
};

#endif
