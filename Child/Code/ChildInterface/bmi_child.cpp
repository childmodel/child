#include <stdio.h>
#include <string.h>
#include <stdexcept>
#include <vector>

#include "child.h"
#include "bmi_child.hxx"


#define VERBOSE (false)
#define DEBUG

#if defined (DEBUG)
# define CHECK_OR_THROW(assertion, err) { if (!(assertion)) throw err; }
#else
# define CHECK_OR_THROW(assertion, err) { }
#endif


void Model::Initialize (std::string file) {
  this->model.Initialize(file);
}

void Model::Update () {
	if (this->model.initialized == false)
    throw CLASS_NOT_INITIALIZED;

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

std::string Model::GetComponentName () {
  return "Child";
}

int Model::GetInputItemCount (void) {
  return this->input_var_name_count;
}

int Model::GetOutputItemCount (void) {
  return this->output_var_name_count;
}

std::vector<std::string>
Model::GetInputVarNames () {
  std::vector<std::string> names;

  for (int i=0; i<this->input_var_name_count; i++)
    names.push_back(this->input_var_names[i]);

  return names;
}

std::vector<std::string>
Model::GetOutputVarNames () {
  std::vector<std::string> names;

  for (int i=0; i<this->output_var_name_count; i++)
    names.push_back(this->output_var_names[i]);

  return names;
}


std::string
Model::GetVarType (std::string name) {
  if (name.compare("land_surface__elevation") == 0 ) {
    return "double";
  } else if (name.compare("sea_bottom_surface__elevation") == 0 ) {
    return "double";
  } else if (name.compare("sea_floor_bedrock_surface__elevation") == 0 ) {
    return "double";
  } else if (name.compare("bedrock_surface__elevation") == 0 ) {
    return "double";
  } else if (name.compare("bedrock_surface__elevation_increment") == 0 ) {
    return "double";
  } else if (name.compare("land_surface__elevation_increment") == 0 ) {
    return "double";
  } else if (name.compare("sediment__erosion_rate") == 0 ) {
    return "double";
  } else if (name.compare("channel_water__discharge") == 0 ) {
    return "double";
  } else if (name.compare("channel_water_sediment~bedload__mass_flow_rate") == 0 ) {
    return "double";
  } else {
    throw BAD_VAR_NAME;
  }
}


std::string
Model::GetVarUnits (std::string name) {
  if (name.compare("land_surface__elevation") == 0 ) {
    return "m";
  } else if (name.compare("sea_bottom_surface__elevation") == 0 ) {
    return "m";
  } else if (name.compare("sea_floor_bedrock_surface__elevation") == 0 ) {
    return "m";
  } else if (name.compare("bedrock_surface__elevation") == 0 ) {
    return "m";
  } else if (name.compare("bedrock_surface__elevation_increment") == 0 ) {
    return "m";
  } else if (name.compare("land_surface__elevation_increment") == 0 ) {
    return "m";
  } else if (name.compare("sediment__erosion_rate") == 0 ) {
    return "m / s";
  } else if (name.compare("channel_water__discharge") == 0 ) {
    return "m^3 / s";
  } else if (name.compare("channel_water_sediment~bedload__mass_flow_rate") == 0 ) {
    return "kg / s";
  } else {
    throw BAD_VAR_NAME;
  }
}


int Model::GetVarGrid (std::string name) {
  int grid_id;

  if (name.compare("land_surface__elevation") == 0 ) {
    grid_id = 0;
  } else if (name.compare("sea_bottom_surface__elevation") == 0 ) {
    grid_id = 0;
  } else if (name.compare("sea_floor_bedrock_surface__elevation") == 0 ) {
    grid_id = 0;
  } else if (name.compare("bedrock_surface__elevation") == 0 ) {
    grid_id = 0;
  } else if (name.compare("bedrock_surface__elevation_increment") == 0 ) {
    grid_id = 0;
  } else if (name.compare("land_surface__elevation_increment") == 0 ) {
    grid_id = 0;
  } else if (name.compare("sediment__erosion_rate") == 0 ) {
    grid_id = 0;
  } else if (name.compare("channel_water__discharge") == 0 ) {
    grid_id = 0;
  } else if (name.compare("channel_water_sediment~bedload__mass_flow_rate") == 0 ) {
    grid_id = 0;
  } else {
    throw std::runtime_error("unknown var name");
    // throw BAD_VAR_NAME;
  }
  return grid_id;
}

int Model::GetGridRank (const int grid_id) {
  int rank;

  if (grid_id == 0) {
    rank = 2;
  } else {
    throw BAD_VAR_NAME;
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

std::string
Model::GetTimeUnits () {
  return "year";
}


double Model::GetTimeStep () {
  return -1.;
}


int Model::GetVarItemsize(std::string name) {
  int itemsize = 0;

  if (name.compare("land_surface__elevation") == 0 ) {
    itemsize = sizeof(double);
  } else if (name.compare("sea_bottom_surface__elevation") == 0 ) {
    itemsize = sizeof(double);
  } else if (name.compare("sea_floor_bedrock_surface__elevation") == 0 ) {
    itemsize = sizeof(double);
  } else if (name.compare("bedrock_surface__elevation") == 0 ) {
    itemsize = sizeof(double);
  } else if (name.compare("bedrock_surface__elevation_increment") == 0 ) {
    itemsize = sizeof(double);
  } else if (name.compare("land_surface__elevation_increment") == 0 ) {
    itemsize = sizeof(double);
  } else if (name.compare("sediment__erosion_rate") == 0 ) {
    itemsize = sizeof(double);
  } else if (name.compare("channel_water__discharge") == 0 ) {
    itemsize = sizeof(double);
  } else if (name.compare("channel_water_sediment~bedload__mass_flow_rate") == 0 ) {
    itemsize = sizeof(double);
  } else {
    throw BAD_VAR_NAME;
  }
  return itemsize;
}

int Model::GetVarNbytes(std::string name) {
  const int itemsize = GetVarItemsize(name);
  const int id = GetVarGrid(name);
  const int size = GetGridSize(id);

  return itemsize * size;
}

std::string Model::GetVarLocation(std::string name) {
  return "node";
}


void Model::GetValue(std::string name, void * const dest) {
  if (name.compare("land_surface__elevation") == 0 ) {
    this->model.CopyNodeElevations ((double*)dest);
  } else if (name.compare("sea_bottom_surface__elevation") == 0 ) {
    this->model.CopyNodeElevations ((double*)dest);
  } else if (name.compare("land_surface__elevation_increment") == 0 ) {
    this->model.CopyNodeElevations ((double*)dest);
  } else if (name.compare("sediment__erosion_rate") == 0 ) {
    this->model.CopyNodeErosion ((double*)dest);
  } else if (name.compare("channel_water__discharge") == 0 ) {
    this->model.CopyNodeDischarge ((double*)dest);
  } else if (name.compare("channel_water_sediment~bedload__mass_flow_rate") == 0 ) {
    this->model.CopyNodeSedimentFlux ((double*)dest);
  } else {
    throw BAD_VAR_NAME;
  }
}

void *Model::GetValuePtr(std::string name) {
  throw NotImplemented();
}


void Model::GetValueAtIndices(std::string name, void *dest, int *inds, int count) {
  throw NotImplemented();
}


void Model::SetValue (std::string name, void *vals) {
  if (name.compare("land_surface__elevation") == 0 ) {
    return this->model.SetNodeElevations ((double*)vals);
  } else if (name.compare("sea_bottom_surface__elevation") == 0 ) {
    return this->model.SetNodeElevations ((double*)vals);
  } else if (name.compare("sea_floor_bedrock_surface__elevation") == 0 ) {
    return this->model.SetNodeElevations ((double*)vals);
  } else if (name.compare("bedrock_surface__elevation") == 0 ) {
    return this->model.SetNodeElevations ((double*)vals);
  } else if (name.compare("bedrock_surface__elevation_increment") == 0 ) {
    return this->model.SetNodeUplift ((double*)vals);
  } else {
    throw BAD_VAR_NAME;
  }
}


void Model::SetValueAtIndices(std::string name, int *inds, int len, void *src) {
  throw NotImplemented();
}


void Model::GetGridShape(const int grid, int *shape) {
  throw NotImplemented();
}


void Model::GetGridSpacing(const int grid, double *spacing) {
  throw NotImplemented();
}


void Model::GetGridOrigin(const int grid, double *origin) {
  throw NotImplemented();
}


void Model::GetGridX (const int grid_id, double * const x) {
  if (grid_id == 0) {
    tLNode *current_node;
    tMesh<tLNode>::nodeListIter_t ni( this->model.mesh->getNodeList() );

    for (current_node=ni.FirstP(); !ni.AtEnd(); current_node=ni.NextP()) {
      x[current_node->getPermID ()] = current_node->getX ();
    }

  } else {
    throw bmi::BMI_FAILURE;
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
    throw bmi::BMI_FAILURE;
  }
}


void Model::GetGridZ (const int grid_id, double * const y) {
  throw NotImplemented();
}


int Model::GetGridNodeCount(const int grid_id) {
  return this->model.mesh->getNodeList()->getSize();
}


int Model::GetGridEdgeCount(const int grid_id) {
  return 0;
  // return this->model.mesh->getEdgeList()->getSize();
}


int Model::GetGridFaceCount(const int grid_id) {
  if (grid_id == 0) {
    return 0;
    // return this->model.mesh->getTriList()->getSize();
  } else {
    throw bmi::BMI_FAILURE;
  }
}

void Model::GetGridNodesPerFace (const int grid_id, int * edges_per_face) {
  if (grid_id == 0) {
    int i = 0;
    const int n_faces = GetGridFaceCount(grid_id);
    for (i=0; i<n_faces; i++)
      edges_per_face[i] = 3;
  } else {
    throw bmi::BMI_FAILURE;
  }
}


void Model::GetGridFaceNodes (const int grid_id, int * face_nodes) {
  if (grid_id == 0) {
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
    throw bmi::BMI_FAILURE;
  }
}


void Model::GetGridEdgeNodes(const int grid, int *edge_nodes) {
  throw NotImplemented();
}


void Model::GetGridFaceEdges(const int grid, int *face_edges) {
  throw NotImplemented();
}


std::string Model::GetGridType (const int grid_id) {
  if (grid_id == 0) {
    return "unstructured";
    // return "unstructured_triangular";
  } else {
    throw bmi::BMI_FAILURE;
  }
}
