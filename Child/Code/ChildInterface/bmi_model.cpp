#include <stdio.h>
#include <string.h>
#include <stdexcept>

#include "child.h"
#include "bmi_model.h"


#define VERBOSE (false)
#define DEBUG

#if defined (DEBUG)
# define CHECK_OR_THROW(assertion, err) { if (!(assertion)) throw err; }
#else
# define CHECK_OR_THROW(assertion, err) { }
#endif


Model::Model () {
    input_var_name_count = 5;
    input_var_names = (char**)malloc(sizeof(char*) * input_var_name_count);
    input_var_names[0] = (char*)malloc(sizeof(char) * 2048 * input_var_name_count);
    for (int i=1; i<input_var_name_count; i++)
      input_var_names[i] = input_var_names[i-1] + 2048;

    output_var_name_count = 6;
    output_var_names = (char**)malloc(sizeof(char*) * output_var_name_count);
    output_var_names[0] = (char*)malloc(sizeof(char) * 2048 * output_var_name_count);
    for (int i=1; i<output_var_name_count; i++)
      output_var_names[i] = output_var_names[i-1] + 2048;

    strcpy(input_var_names[0], "land_surface__elevation");
    strcpy(input_var_names[1], "sea_bottom_surface__elevation");
    strcpy(input_var_names[2], "sea_floor_bedrock_surface__elevation");
    strcpy(input_var_names[3], "bedrock_surface__elevation");
    strcpy(input_var_names[4], "bedrock_surface__elevation_increment");

    strcpy(output_var_names[0], "land_surface__elevation");
    strcpy(output_var_names[1], "sea_bottom_surface__elevation");
    strcpy(output_var_names[2], "land_surface__elevation_increment");
    strcpy(output_var_names[3], "sediment__erosion_rate");
    strcpy(output_var_names[4], "channel_water__discharge");
    strcpy(output_var_names[5], "channel_water_sediment~bedload__mass_flow_rate");

    this->model = Child();
}


void Model::Initialize (const char *file) {
  // Child::Initialize(std::string(file));
  this->model.Initialize(file);
}

void Model::Update () {
	if (this->model.initialized==false)
    throw bmi::CLASS_NOT_INITIALIZED;

  this->model.RunOneStorm ();
  this->model.MaskNodesBelowElevation (0);
}

void Model::UpdateUntil (double t) {
  double dt;

  dt = t - this->GetCurrentTime ();
  if (dt > 0)
    this->model.Run (dt);
  this->model.MaskNodesBelowElevation (0);
}

void Model::Finalize () {
  this->model.CleanUp ();
}

void Model::GetComponentName (char * const name) {
  strcpy (name, "child");
}

int Model::GetInputVarNameCount (void) {
  return this->input_var_name_count;
}

int Model::GetOutputVarNameCount (void) {
  return this->output_var_name_count;
}

void Model::GetInputVarNames (char * const * const names) {
  for (int i=0; i<this->input_var_name_count; i++) {
    strncpy (names[i], input_var_names[i], 2048);
  }
}

void Model::GetOutputVarNames (char * const * const names) {
  for (int i=0; i<this->output_var_name_count; i++) {
    strcpy (names[i], output_var_names[i]);
  }
}

void Model::GetVarType (const char * name, char * const vtype) {
  if (strcmp(name, "land_surface__elevation") == 0 ) {
    strncpy(vtype, "double", bmi::MAX_UNITS_NAME);
  } else if (strcmp(name, "sea_bottom_surface__elevation") == 0 ) {
    strncpy(vtype, "double", bmi::MAX_UNITS_NAME);
  } else if (strcmp(name, "sea_floor_bedrock_surface__elevation") == 0 ) {
    strncpy(vtype, "double", bmi::MAX_UNITS_NAME);
  } else if (strcmp(name, "bedrock_surface__elevation") == 0 ) {
    strncpy(vtype, "double", bmi::MAX_UNITS_NAME);
  } else if (strcmp(name, "bedrock_surface__elevation_increment") == 0 ) {
    strncpy(vtype, "double", bmi::MAX_UNITS_NAME);
  } else if (strcmp(name, "land_surface__elevation_increment") == 0 ) {
    strncpy(vtype, "double", bmi::MAX_UNITS_NAME);
  } else if (strcmp(name, "sediment__erosion_rate") == 0 ) {
    strncpy(vtype, "double", bmi::MAX_UNITS_NAME);
  } else if (strcmp(name, "channel_water__discharge") == 0 ) {
    strncpy(vtype, "double", bmi::MAX_UNITS_NAME);
  } else if (strcmp(name, "channel_water_sediment~bedload__mass_flow_rate") == 0 ) {
    strncpy(vtype, "double", bmi::MAX_UNITS_NAME);
  } else {
    throw bmi::BAD_VAR_NAME;
  }
}

void Model::GetVarUnits (const char * name,
                              char * const units) {
  if (strcmp(name, "land_surface__elevation") == 0 ) {
    strncpy(units, "m", bmi::MAX_UNITS_NAME);
  } else if (strcmp(name, "sea_bottom_surface__elevation") == 0 ) {
    strncpy(units, "m", bmi::MAX_UNITS_NAME);
  } else if (strcmp(name, "sea_floor_bedrock_surface__elevation") == 0 ) {
    strncpy(units, "m", bmi::MAX_UNITS_NAME);
  } else if (strcmp(name, "bedrock_surface__elevation") == 0 ) {
    strncpy(units, "m", bmi::MAX_UNITS_NAME);
  } else if (strcmp(name, "bedrock_surface__elevation_increment") == 0 ) {
    strncpy(units, "m", bmi::MAX_UNITS_NAME);
  } else if (strcmp(name, "land_surface__elevation_increment") == 0 ) {
    strncpy(units, "m", bmi::MAX_UNITS_NAME);
  } else if (strcmp(name, "sediment__erosion_rate") == 0 ) {
    strncpy(units, "m / s", bmi::MAX_UNITS_NAME);
  } else if (strcmp(name, "channel_water__discharge") == 0 ) {
    strncpy(units, "m^3 / s", bmi::MAX_UNITS_NAME);
  } else if (strcmp(name, "channel_water_sediment~bedload__mass_flow_rate") == 0 ) {
    strncpy(units, "kg / s", bmi::MAX_UNITS_NAME);
  } else {
    throw bmi::BAD_VAR_NAME;
  }
}

int Model::GetVarGrid (const char * name) {
  int grid_id;

  if (strcmp(name, "land_surface__elevation") == 0 ) {
    grid_id = 0;
  } else if (strcmp(name, "sea_bottom_surface__elevation") == 0 ) {
    grid_id = 0;
  } else if (strcmp(name, "sea_floor_bedrock_surface__elevation") == 0 ) {
    grid_id = 0;
  } else if (strcmp(name, "bedrock_surface__elevation") == 0 ) {
    grid_id = 0;
  } else if (strcmp(name, "bedrock_surface__elevation_increment") == 0 ) {
    grid_id = 0;
  } else if (strcmp(name, "land_surface__elevation_increment") == 0 ) {
    grid_id = 0;
  } else if (strcmp(name, "sediment__erosion_rate") == 0 ) {
    grid_id = 0;
  } else if (strcmp(name, "channel_water__discharge") == 0 ) {
    grid_id = 0;
  } else if (strcmp(name, "channel_water_sediment~bedload__mass_flow_rate") == 0 ) {
    grid_id = 0;
  } else {
    throw std::runtime_error("unknown var name");
    // throw bmi::BAD_VAR_NAME;
  }
  return grid_id;
}

int Model::GetGridRank (const int grid_id) {
  int rank;

  if (grid_id == 0) {
    rank = 2;
  } else {
    throw bmi::BAD_VAR_NAME;
  }
  return rank;
}

int Model::GetGridSize (const int grid_id) {
  int size = 0;
  size = this->model.mesh->getNodeList()->getSize();
  return size;
}

double Model::GetStartTime () {
  return 0.0;
}

double Model::GetCurrentTime () {
  double now = 0.;
  now = this->model.time->getCurrentTime ();
  return now;
}

double Model::GetEndTime () {
  double stop = 0.;
  stop = this->model.time->getCurrentTime () + this->model.time->RemainingTime ();
  return stop;
}

void Model::GetTimeUnits (char * const units) {
  strncpy(units, "year", 2048);
}


double Model::GetTimeStep () {
  return -1.;
}


int Model::GetVarItemsize(const char * name) {
  int itemsize = 0;

  if (strcmp(name, "land_surface__elevation") == 0 ) {
    itemsize = sizeof(double);
  } else if (strcmp(name, "sea_bottom_surface__elevation") == 0 ) {
    itemsize = sizeof(double);
  } else if (strcmp(name, "sea_floor_bedrock_surface__elevation") == 0 ) {
    itemsize = sizeof(double);
  } else if (strcmp(name, "bedrock_surface__elevation") == 0 ) {
    itemsize = sizeof(double);
  } else if (strcmp(name, "bedrock_surface__elevation_increment") == 0 ) {
    itemsize = sizeof(double);
  } else if (strcmp(name, "land_surface__elevation_increment") == 0 ) {
    itemsize = sizeof(double);
  } else if (strcmp(name, "sediment__erosion_rate") == 0 ) {
    itemsize = sizeof(double);
  } else if (strcmp(name, "channel_water__discharge") == 0 ) {
    itemsize = sizeof(double);
  } else if (strcmp(name, "channel_water_sediment~bedload__mass_flow_rate") == 0 ) {
    itemsize = sizeof(double);
  } else {
    throw bmi::BAD_VAR_NAME;
  }
  return itemsize;
}

int Model::GetVarNbytes(const char * name) {
  const int itemsize = GetVarItemsize(name);
  const int id = GetVarGrid(name);
  const int size = GetGridSize(id);

  return itemsize * size;
}

void Model::GetVarLocation(const char * name, char * const location) {
  strcpy(location, "node");
}


void Model::GetValue(const char * name, void * const dest) {
  if (strcmp(name, "land_surface__elevation") == 0 ) {
    this->model.CopyNodeElevations ((double*)dest);
  } else if (strcmp(name, "sea_bottom_surface__elevation") == 0 ) {
    this->model.CopyNodeElevations ((double*)dest);
  } else if (strcmp(name, "land_surface__elevation_increment") == 0 ) {
    this->model.CopyNodeElevations ((double*)dest);
  } else if (strcmp(name, "sediment__erosion_rate") == 0 ) {
    this->model.CopyNodeErosion ((double*)dest);
  } else if (strcmp(name, "channel_water__discharge") == 0 ) {
    this->model.CopyNodeDischarge ((double*)dest);
  } else if (strcmp(name, "channel_water_sediment~bedload__mass_flow_rate") == 0 ) {
    this->model.CopyNodeSedimentFlux ((double*)dest);
  } else {
    throw bmi::BAD_VAR_NAME;
  }
}

void Model::SetValue (const char * name, void *vals) {
  if (strcmp(name, "land_surface__elevation") == 0 ) {
    return this->model.SetNodeElevations ((double*)vals);
  } else if (strcmp(name, "sea_bottom_surface__elevation") == 0 ) {
    return this->model.SetNodeElevations ((double*)vals);
  } else if (strcmp(name, "sea_floor_bedrock_surface__elevation") == 0 ) {
    return this->model.SetNodeElevations ((double*)vals);
  } else if (strcmp(name, "bedrock_surface__elevation") == 0 ) {
    return this->model.SetNodeElevations ((double*)vals);
  } else if (strcmp(name, "bedrock_surface__elevation_increment") == 0 ) {
    return this->model.SetNodeUplift ((double*)vals);
  } else {
    throw bmi::BAD_VAR_NAME;
  }
}

void Model::GetGridX (const int grid_id, double * const x) {
  if (grid_id == 0) {
    tLNode *current_node;
    tMesh<tLNode>::nodeListIter_t ni( this->model.mesh->getNodeList() );

    for (current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP()) {
      x[current_node->getPermID ()] = current_node->getX ();
    }

  } else {
    throw bmi::FAILURE;
  }
}

void Model::GetGridY (const int grid_id, double * const y) {
  if (grid_id == 0) {
    tLNode *current_node;
    tMesh<tLNode>::nodeListIter_t ni( this->model.mesh->getNodeList() );

    for (current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP()) {
      y[current_node->getPermID ()] = current_node->getY ();
    }
  } else {
    throw bmi::FAILURE;
  }
}


int Model::GetGridNumberOfFaces(const int grid_id) {
  if (grid_id == 0) {
    return this->model.mesh->getTriList()->getSize();
  } else {
    throw bmi::FAILURE;
  }
}


int Model::GetGridFaceCount(const int grid_id) {
  if (grid_id == 0) {
    return this->model.mesh->getTriList()->getSize();
  } else {
    throw bmi::FAILURE;
  }
}


int Model::GetGridVertexCount(const int grid_id) {
  if (grid_id == 0) {
    return 3. * GetGridFaceCount(grid_id);
  } else {
    throw bmi::FAILURE;
  }
}


void Model::GetGridNodesPerFace (const int grid_id, int * edges_per_face) {
  if (grid_id == 0) {
    int i = 0;
    const int n_faces = GetGridFaceCount(grid_id);
    for (i=0; i<n_faces; i++)
      edges_per_face[i] = 3;
  } else {
    throw bmi::FAILURE;
  }
}


void Model::GetGridFaceNodes (const int grid_id, int * face_nodes) {
  if (grid_id == 0) {
    // Implement this: connectivity for this grid.
    {
      int i = 0;
      tMesh<tLNode>::triListIter_t ti(this->model.mesh->getTriList());
      for (tTriangle * current_tri=ti.FirstP(); !ti.AtEnd();
           i+=3, current_tri=ti.NextP ()) {
        face_nodes[i] = current_tri->pPtr(0)->getPermID();
        face_nodes[i+1] = current_tri->pPtr(1)->getPermID();
        face_nodes[i+2] = current_tri->pPtr(2)->getPermID();
      }
    }
  } else {
    throw bmi::FAILURE;
  }
}


void Model::GetGridConnectivity (const int grid_id, int * connectivity) {
  if (grid_id == 0) {
    // Implement this: connectivity for this grid.
    {
      int i = 0;
      tMesh<tLNode>::triListIter_t ti(this->model.mesh->getTriList());
      for (tTriangle * current_tri=ti.FirstP(); !ti.AtEnd();
           i+=3, current_tri=ti.NextP ()) {
        connectivity[i] = current_tri->pPtr(0)->getPermID();
        connectivity[i+1] = current_tri->pPtr(1)->getPermID();
        connectivity[i+2] = current_tri->pPtr(2)->getPermID();
      }
    }
  } else {
    throw bmi::FAILURE;
  }
}


void Model::GetGridOffset (const int grid_id, int * const offset) {
  if (grid_id == 0) {
    // Implement this: connectivity for this grid.
    const int n_offsets = this->model.mesh->getTriList()->getSize();
    offset[0] = 3;
    for (int i=1; i<n_offsets; i++)
      offset[i] = offset[i-1] + 3;
  } else {
    throw bmi::FAILURE;
  }
}

void Model::GetGridType (const int grid_id, char * const gtype) {
  if (grid_id == 0) {
    strncpy(gtype, "unstructured_triangular", 2048);
  } else {
    throw bmi::FAILURE;
  }
}
/*
bool Model::HasInputVar (const char * var_name) {
  for (int i=0; i<this->input_var_name_count; i++) {
    if (strcmp (input_var_names[i], var_name) == 0)
      return true;
  }
  return false;
}

bool Model::HasOutputVar (const char * var_name) {
  for (int i=0; i<this->output_var_name_count; i++) {
    if (strcmp (output_var_names[i], var_name) == 0)
      return true;
  }
  return false;
}
*/
