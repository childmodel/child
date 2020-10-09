#ifndef BMI_CHILD_H_INCLUDED
#define BMI_CHILD_H_INCLUDED

#include <string>
#include <iostream>

#include "bmi.hxx"
#include "child.h"


class NotImplemented : public std::logic_error {
  public:
  NotImplemented() : std::logic_error("Not Implemented") { };
};


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


class Model : public bmi::Bmi {
  public:
    Model() {
      this->input_var_names[0] = "land_surface__elevation";
      this->input_var_names[1] = "sea_bottom_surface__elevation";
      this->input_var_names[2] = "sea_floor_bedrock_surface__elevation";
      this->input_var_names[3] = "bedrock_surface__elevation";
      this->input_var_names[4] = "bedrock_surface__elevation_increment";

      this->output_var_names[0] = "land_surface__elevation";
      this->output_var_names[1] = "sea_bottom_surface__elevation";
      this->output_var_names[2] = "land_surface__elevation_increment";
      this->output_var_names[3] = "sediment__erosion_rate";
      this->output_var_names[4] = "channel_water__discharge";
      this->output_var_names[5] = "channel_water_sediment~bedload__mass_flow_rate";

      this->model = Child();
    };

    void Initialize(std::string config_file);
    void Update();
    void UpdateUntil(double time);
    void Finalize();

    std::string GetComponentName();
    int GetInputItemCount();
    int GetOutputItemCount();
    std::vector<std::string> GetInputVarNames();
    std::vector<std::string> GetOutputVarNames();

    int GetVarGrid(std::string name);
    std::string GetVarType(std::string name);
    int GetVarItemsize(std::string name);
    std::string GetVarUnits(std::string name);
    int GetVarNbytes(std::string name);
    std::string GetVarLocation(std::string name);

    double GetCurrentTime();
    double GetStartTime();
    double GetEndTime();
    std::string GetTimeUnits();
    double GetTimeStep();

    void GetValue(std::string name, void *dest);
    void *GetValuePtr(std::string name);
    void GetValueAtIndices(std::string name, void *dest, int *inds, int count);

    void SetValue(std::string name, void *src);
    void SetValueAtIndices(std::string name, int *inds, int len, void *src);

    int GetGridRank(const int grid);
    int GetGridSize(const int grid);
    std::string GetGridType(const int grid);

    void GetGridShape(const int grid, int *shape);
    void GetGridSpacing(const int grid, double *spacing);
    void GetGridOrigin(const int grid, double *origin);

    void GetGridX(const int grid, double *x);
    void GetGridY(const int grid, double *y);
    void GetGridZ(const int grid, double *z);

    int GetGridNodeCount(const int grid);
    int GetGridEdgeCount(const int grid);
    int GetGridFaceCount(const int grid);

    void GetGridEdgeNodes(const int grid, int *edge_nodes);
    void GetGridFaceEdges(const int grid, int *face_edges);
    void GetGridFaceNodes(const int grid, int *face_nodes);
    void GetGridNodesPerFace(const int grid, int *nodes_per_face);

  private:
    Child model;
    static const int input_var_name_count = 5;
    static const int output_var_name_count = 6;

    std::string input_var_names[6];
    std::string output_var_names[7];
};

#endif
