#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../bmi_model_child.h"

#define CHECK_INTERFACE_FUNC(f_name, call) { \
  fprintf (stdout, "\033[32m%s... \033[39m", f_name); \
  try { \
    call; \
    fprintf (stdout, "\033[32mPASS\033[39m\n"); \
  } \
  catch (bmi::FatalError e) { \
    fprintf (stdout, "\033[31mFAIL\033[39m\n"); \
    exit (EXIT_FAILURE); \
  } \
  catch (bmi::NonFatalError e) { \
    fprintf (stdout, "\033[33mWARN\033[39m\n"); \
  } \
}

int
main (int argc, char *argv[])
{
  bmi::Model child;
  char input_file[2048];

  if (argc>1)
    strncpy (input_file, argv[1], 2048);
  else {
    fprintf (stderr, "ERROR: Incorrect number of arguments (%d).\n", argc);
    exit (EXIT_FAILURE);
  }

  { /* Model control functions */
    fprintf (stdout, "Model control functions\n");
    fprintf (stdout, "=======================\n\n");

    CHECK_INTERFACE_FUNC ("Initialize", child.Initialize (input_file));
    CHECK_INTERFACE_FUNC ("Update", child.Update ());
    CHECK_INTERFACE_FUNC ("Finalize", child.Finalize ());
  }

  child.Initialize (input_file);
  child.Update ();

  { /* Model information functions */
    int i;
    char model_name[bmi::COMPONENT_NAME_MAX];
    int input_var_count;
    int output_var_count;
    char **input_names = NULL;
    char **output_names = NULL;
    char **name;
    double time;

    fprintf (stdout, "\n\n");
    fprintf (stdout, "Model information functions\n");
    fprintf (stdout, "===========================\n");

    CHECK_INTERFACE_FUNC ("Get_component_name", child.GetComponentName (model_name));

    CHECK_INTERFACE_FUNC ("Get_current_time", child.GetCurrentTime ());
    CHECK_INTERFACE_FUNC ("Get_start_time", child.GetStartTime ());
    CHECK_INTERFACE_FUNC ("Get_end_time", child.GetEndTime ());

    CHECK_INTERFACE_FUNC ("Get_input_var_name_count", input_var_count = child.GetInputVarNameCount ());
    CHECK_INTERFACE_FUNC ("Get_output_var_name_count", output_var_count = child.GetOutputVarNameCount ());

    input_names = (char**) malloc (sizeof (char*) * input_var_count);
    for (i=0; i<input_var_count; i++)
      input_names[i] = (char*) malloc (sizeof (char) * bmi::VAR_NAME_MAX);

    CHECK_INTERFACE_FUNC ("Get_input_var_names", child.GetInputVarNames (input_names));

    output_names = (char**) malloc (sizeof (char*) * output_var_count);
    for (i=0; i<output_var_count; i++)
      output_names[i] = (char*) malloc (sizeof (char) * bmi::VAR_NAME_MAX);

    CHECK_INTERFACE_FUNC ("Get_output_var_names", child.GetOutputVarNames (output_names));

    for (i=0, name=input_names; i<input_var_count; i++, name++) {
      bmi::VarType var_type;
      char units[bmi::UNITS_NAME_MAX];
      int rank;
      int number_of_points;
      bmi::GridType grid_type;
      double *buffer = NULL;
      int *shape;
      double *spacing;
      double *origin;

      fprintf (stdout, "\n\n");
      fprintf (stdout, "Input variable: %s\n", *name);
      fprintf (stdout, "==================================\n");

      /* Variable information functions */
      CHECK_INTERFACE_FUNC ("Get_var_type", var_type = child.GetVarType (*name));
      CHECK_INTERFACE_FUNC ("Get_var_rank", rank = child.GetVarRank (*name));

      CHECK_INTERFACE_FUNC ("Get_var_point_count", number_of_points = child.GetVarPointCount (*name));

      buffer = (double*) malloc (sizeof (double) * number_of_points);

      /* Variable getters */
      //CHECK_INTERFACE_FUNC ("bmi_Get_double", bmi_Get_double (child, *name, buffer));

      /* Variable setters */
      CHECK_INTERFACE_FUNC ("Set_double", child.SetDouble (*name, buffer));
      free (buffer);
    }

    for (i=0, name=output_names; i<output_var_count; i++, name++) {
      bmi::VarType var_type;
      char units[bmi::UNITS_NAME_MAX];
      int rank;
      int number_of_points;
      bmi::GridType grid_type;
      double *buffer = NULL;
      int *shape;
      double *spacing;
      double *origin;

      fprintf (stdout, "\n\n");
      fprintf (stdout, "Output variable: %s\n", *name);
      fprintf (stdout, "==================================\n");

      /* Variable information functions */
      var_type = child.GetVarType (*name);
      rank = child.GetVarRank (*name);
      number_of_points = child.GetVarPointCount (*name);

      fprintf (stdout, "            Type: %d\n", var_type);
      fprintf (stdout, "            Rank: %d\n", rank);
      fprintf (stdout, "Number of points: %d\n", number_of_points);
    }

    for (i=0, name=output_names; i<output_var_count; i++, name++) {
      bmi::VarType var_type;
      char units[bmi::UNITS_NAME_MAX];
      int rank;
      int number_of_points;
      bmi::GridType grid_type;
      int *shape;
      double *spacing;
      double *origin;

      fprintf (stdout, "\n\n");
      fprintf (stdout, "Output variable: %s\n", *name);
      fprintf (stdout, "==================================\n");

      /* Variable information functions */
      CHECK_INTERFACE_FUNC ("Get_var_type", var_type = child.GetVarType (*name));
      CHECK_INTERFACE_FUNC ("Get_var_rank", rank = child.GetVarRank (*name));

      CHECK_INTERFACE_FUNC ("Get_var_point_count", number_of_points = child.GetVarPointCount (*name));

      if (number_of_points > 0) {
        double *buffer = (double*) malloc (sizeof (double) * number_of_points);

        /* Variable getters */
        CHECK_INTERFACE_FUNC ("Get_double", child.GetDouble (*name, buffer));

        /* Variable setters */
        //CHECK_INTERFACE_FUNC ("bmi_Set_double", bmi_Set_double (child, *name, buffer));
        free (buffer);

        CHECK_INTERFACE_FUNC ("Get_double_ptr", buffer = child.GetDoublePtr (*name));

        try {
          buffer = child.GetDoublePtr (*name);
          if (rank > 0) {
            int *stride = (int*) malloc (sizeof (int) * rank);
            CHECK_INTERFACE_FUNC ("Get_var_stride", child.GetVarStride (*name, stride));
            free (stride);
          }
        }
        catch (bmi::NonFatalError) { }

        if (rank > 0) {
          CHECK_INTERFACE_FUNC ("Get_grid_type", grid_type = child.GetGridType (*name));

          if (grid_type == bmi::GRID_TYPE_UNIFORM) {
            shape = (int*) malloc (sizeof (int) * rank);
            spacing = (double*) malloc (sizeof (double) * rank);
            origin = (double*) malloc (sizeof (double) * rank);

            CHECK_INTERFACE_FUNC ("Get_grid_shape", child.GetGridShape (*name, shape));
            CHECK_INTERFACE_FUNC ("Get_grid_spacing", child.GetGridSpacing (*name, spacing));
            CHECK_INTERFACE_FUNC ("Get_grid_origin", child.GetGridOrigin (*name, origin));

            free (origin);
            free (spacing);
            free (shape);
          }
          else if (grid_type == bmi::GRID_TYPE_UNSTRUCTURED) {
            CHECK_INTERFACE_FUNC ("GetVarVertexCount", child.GetVarVertexCount (*name));
            CHECK_INTERFACE_FUNC ("GetVarCellCount", child.GetVarCellCount (*name));

            const int number_of_vertices = child.GetVarVertexCount (*name);
            const int number_of_cells = child.GetVarCellCount (*name);
            double * x = new double[number_of_points];
            double * y = new double[number_of_points];
            int * c = new int[number_of_vertices];
            int * o = new int[number_of_cells];

            CHECK_INTERFACE_FUNC ("GetGridX", child.GetGridX (*name, x));
            CHECK_INTERFACE_FUNC ("GetGridY", child.GetGridX (*name, y));
            CHECK_INTERFACE_FUNC ("GetGridConnectivity", child.GetGridConnectivity (*name, c));
            CHECK_INTERFACE_FUNC ("GetGridOffset", child.GetGridOffset (*name, o));

            delete x, y, c, o;
          }
        }
      }
    }
  }

  child.Finalize ();

  fprintf (stdout, "\n\n");

  return EXIT_SUCCESS;
}

