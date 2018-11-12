#ifndef BMI_H_
#define BMI_H_

namespace bmi {

const int MAX_COMPONENT_NAME = 2048;
const int MAX_VAR_NAME = 2048;
const int MAX_UNITS_NAME = 2048;

typedef enum {
  FAILURE = 1,
  BAD_ARGUMENT,
  BAD_VAR_NAME,
  BAD_FILE,
  CLASS_NOT_INITIALIZED,
} FatalError;

typedef enum {
  FUNCTION_NOT_IMPLEMENTED = -1,
} NonFatalError;


class Bmi {
 public:
  // Bmi();

  // Model control functions.
  virtual void Initialize(const char *) = 0;
  virtual void Update() = 0;
  virtual void UpdateUntil(double) = 0;
  virtual void Finalize() = 0;

  // Model information functions.
  virtual void GetComponentName(char * const name) = 0;
  virtual int GetInputVarNameCount (void) = 0;
  virtual int GetOutputVarNameCount (void) = 0;
  virtual void GetInputVarNames (char * const * const names) = 0;
  virtual void GetOutputVarNames (char * const * const names) = 0;

  // Variable information functions
  virtual int GetVarGrid (const char * var_name) = 0;
  virtual void GetVarType (const char * var_name, char * const vtype) = 0;
  virtual void GetVarUnits (const char * var_name, char * const units) = 0;
  virtual int GetVarItemsize(const char * name) = 0;
  virtual int GetVarNbytes(const char * name) = 0;
  virtual void GetVarLocation(const char * name, char * const location) = 0;

  virtual double GetCurrentTime () = 0;
  virtual double GetStartTime () = 0;
  virtual double GetEndTime () = 0;
  virtual double GetTimeStep () = 0;
  virtual void GetTimeUnits (char * const units) = 0;

  // Variable getters
  virtual void GetValue (const char * var_name, void *) = 0;
  virtual double * GetValuePtr (const char * var_name) {
    throw FUNCTION_NOT_IMPLEMENTED;
    // throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }
  virtual void GetValueStride (const char * var_name, int * const stride) {
    throw FUNCTION_NOT_IMPLEMENTED;
    //throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }

  // Variable setters
  virtual void SetValue (const char *, void *) = 0;

  // Grid information functions
  virtual void GetGridType (const int grid_id, char * const gtype) = 0;
  virtual int GetGridRank (const int grid_id) = 0;
  virtual int GetGridSize (const int grid_id) = 0;

  virtual void GetGridShape (const int, int *) {
    throw FUNCTION_NOT_IMPLEMENTED;
    // throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }
  virtual void GetGridSpacing (const int, double *) {
    throw FUNCTION_NOT_IMPLEMENTED;
    // throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }
  virtual void GetGridOrigin (const int, double *) {
    throw FUNCTION_NOT_IMPLEMENTED;
    // throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }

  virtual void GetGridX (const int, double * const) = 0;
  virtual void GetGridY (const int, double * const) = 0;
  virtual void GetGridZ (const int, double * const) {
    throw FUNCTION_NOT_IMPLEMENTED;
    // throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }

  virtual int GetGridFaceCount(const int) = 0;
  virtual int GetGridVertexCount(const int) = 0;
  virtual void GetGridConnectivity (const int, int * const ) = 0;
  virtual void GetGridOffset (const int, int * const) = 0;

  virtual int GetGridNumberOfEdges(const int) = 0;
  virtual int GetGridNumberOfFaces(const int) = 0;

  virtual void GetGridEdgeNodes(const int, int * const) = 0;
  virtual void GetGridFaceEdges(const int, int * const) = 0;
  virtual void GetGridFaceNodes(const int, int * const) = 0;
  virtual void GetGridNodesPerFace(const int, int * const) = 0;
};

} // namespace bmi

#endif
