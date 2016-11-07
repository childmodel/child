#ifndef CHILD_H_
#define CHILD_H_

#define BMI_GET_GRID_CONNECTIVITY

#include "child.h"


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
  int GetVarGrid (const char * var_name);
  void GetVarType (const char * var_name, char * const vtype);
  void GetVarUnits (const char * var_name, char * const units);
  int GetVarItemsize(const char * name);
  int GetVarNbytes(const char * name);

  double GetCurrentTime ();
  double GetStartTime ();
  double GetEndTime ();
  double GetTimeStep ();
  void GetTimeUnits (char * const units);

  // Variable getters
  void GetValue (const char * var_name, void *);
  double * GetValuePtr (const char * var_name) {
    throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }
  void GetValueStride (const char * var_name, int * const stride) {
    throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }

  // Variable setters
  void SetValue (const char *, void *);

  // Grid information functions
  void GetGridType (const int grid_id, char * const gtype);
  int GetGridRank (const int grid_id);
  int GetGridSize (const int grid_id);

  void GetGridShape (const int, int *) {
    throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }
  void GetGridSpacing (const int, double *) {
    throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }
  void GetGridOrigin (const int, double *) {
    throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }

  void GetGridX (const int, double * const);
  void GetGridY (const int, double * const);
  void GetGridZ (const int, double * const) {
    throw bmi::FUNCTION_NOT_IMPLEMENTED;
  }

  void GetGridConnectivity (const int, int * const );
  void GetGridOffset (const int, int * const);

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

#endif
