#ifndef CHILD_CHILDINTERFACE_BMI_CHILD_H_
#define CHILD_CHILDINTERFACE_BMI_CHILD_H_

#include "child.h"

#define BMI_GET_GRID_CONNECTIVITY

namespace bmi {

const int COMPONENT_NAME_MAX = 2048;
const int VAR_NAME_MAX = 2048;
const int UNITS_NAME_MAX = 2048;

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

typedef enum {
  GRID_TYPE_UNKNOWN = -1,
  GRID_TYPE_UNIFORM,
  GRID_TYPE_RECTILINEAR,
  GRID_TYPE_STRUCTURED,
  GRID_TYPE_UNSTRUCTURED,
} GridType;

typedef enum {
  VAR_TYPE_UNKNOWN = 0,
  VAR_TYPE_CHAR,
  VAR_TYPE_UNSIGNED_CHAR,
  VAR_TYPE_INT,
  VAR_TYPE_LONG,
  VAR_TYPE_UNSIGNED_INT,
  VAR_TYPE_UNSIGNED_LONG,
  VAR_TYPE_FLOAT,
  VAR_TYPE_DOUBLE,
  VAR_TYPE_NUMBER,
} VarType;

// Implementation of the CSDMS Basic Modeling Interface (BMI) for Child.
//    bmi::Model child;
//    child.Initialize (input_file);
//    for (int i=0; i<10; ++i)
//      child.Update ();
//    child.Finalize ();
class Model : public Child {
 public:
  Model () {
    const char *inputs[] = {
      "surface__elevation",
      "sea_floor__elevation",
      "sea_floor_bedrock_surface__elevation",
      "bedrock_surface__elevation",
      "bedrock_surface__elevation_increment",
      NULL,
    };
    const char *outputs[] = {
      "surface__elevation",
      "sea_floor__elevation",
      "surface__elevation_increment",
      "sediment__erosion_rate",
      "channel_water__discharge",
      "bed_load__mass_flow_rate",
      NULL,
    };
    input_var_names = NULL;
    output_var_names = NULL;
    SetInputVarNames (inputs);
    SetOutputVarNames (outputs);
  }
  ~Model () {
    SetInputVarNames (NULL);
    SetOutputVarNames (NULL);
  }

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
  VarType GetVarType (const char * var_name);
  void GetVarUnits (const char * var_name, char * const type);
  int GetVarRank (const char *var_name);

  double GetCurrentTime ();
  double GetStartTime ();
  double GetEndTime ();
  double GetTimeStep ();
  void GetTimeUnits (char * const units);

  // Variable getters
  void GetDouble (const char * var_name, double *);
  void GetDoubleAtIndices (const char *, double *, int *, int);
  double * GetDoublePtr (const char * var_name) {
    throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }
  void GetVarStride (const char * var_name, int * const stride) {
    throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }

  void GetValue (const char *, void *) {
    throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }
  void *GetValuePtr (const char * var_name) {
    throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }
  void GetValueAtIndices (const char *, void *, int *, int) {
    throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }

  // Variable setters
  void SetDouble (const char *, double *);
  void SetDoubleAtIndices (const char *, int *, int, double *);

  // Grid information functions
  int GetVarPointCount (const char *var_name);
  int GetVarVertexCount (const char *var_name);
  int GetVarCellCount (const char *var_name);

  GridType GetGridType (const char *var_name);

  void GetGridShape (const char *, int *) {
    throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }
  void GetGridSpacing (const char *, double *) {
    throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }
  void GetGridOrigin (const char *, double *) {
    throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }

  void GetGridX (const char *, double * const);
  void GetGridY (const char *, double * const);
  void GetGridZ (const char *, double * const) {
    throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }

  void GetGridConnectivity (const char *, int * const );
  void GetGridOffset (const char *, int * const);

 private:
  bool HasInputVar (const char * var_name);
  bool HasOutputVar (const char * var_name);
  void SetInputVarNames (const char **names);
  void SetOutputVarNames (const char **names);

  int input_var_name_count;
  int output_var_name_count;

  char** input_var_names;
  char** output_var_names;
};

} // namespace bmi

#endif /* bmi_child.h is included */

