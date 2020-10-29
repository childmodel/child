#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../bmi_child.hxx"

#define CHECK_INTERFACE_FUNC(f_name, call) { \
  fprintf (stdout, "\033[32m%s... \033[39m", f_name); \
  try { \
    call; \
    fprintf (stdout, "\033[32mPASS\033[39m\n"); \
  } \
  catch (FatalError e) { \
    fprintf (stdout, "\033[31mFAIL\033[39m\n"); \
    exit (EXIT_FAILURE); \
  } \
  catch (NonFatalError e) { \
    fprintf (stdout, "\033[33mWARN\033[39m\n"); \
  } \
}

Predicates predicate;


int
main (int argc, char *argv[])
{
  Model child;
  // char input_file[2048];
  std::string input_file;

  if (argc>1) {
    input_file = argv[1];
    // strncpy (input_file, argv[1], 2048);
    /*
    FILE *fp = fopen(argv[1], "r");
    if (fp) {
        input_file = fgets(fp);
    }
    fclose(fp);
    */
  } else {
    fprintf (stderr, "ERROR: Incorrect number of arguments (%d).\n", argc);
    exit (EXIT_FAILURE);
  }

  { /* Model control functions */
    fprintf (stdout, "Model control functions\n");
    fprintf (stdout, "=======================\n\n");

    fprintf(stdout, "%s\n", input_file.c_str());

    CHECK_INTERFACE_FUNC ("Initialize", child.Initialize (input_file));
    CHECK_INTERFACE_FUNC ("Update", child.Update ());
    CHECK_INTERFACE_FUNC ("Finalize", child.Finalize ());
  }

  child.Initialize (input_file);
  child.Update ();

  { /* Model information functions */
    int i;
    int input_var_count;
    int output_var_count;
    vector<std::string> input_names;
    vector<std::string> output_names;
    std::string name;
    double time;

    fprintf (stdout, "\n\n");
    fprintf (stdout, "Model information functions\n");
    fprintf (stdout, "===========================\n");

    CHECK_INTERFACE_FUNC ("Get_component_name", child.GetComponentName());

    CHECK_INTERFACE_FUNC ("Get_current_time", child.GetCurrentTime ());
    CHECK_INTERFACE_FUNC ("Get_start_time", child.GetStartTime ());
    CHECK_INTERFACE_FUNC ("Get_end_time", child.GetEndTime ());

    CHECK_INTERFACE_FUNC ("Get_input_item_count", input_var_count = child.GetInputItemCount ());
    CHECK_INTERFACE_FUNC ("Get_output_item_count", output_var_count = child.GetOutputItemCount ());

    input_names = child.GetInputVarNames();
    output_names = child.GetOutputVarNames();

    for (i=0; i<input_var_count; i++) {
      std::string var_type;
      std::string units;
      std::string grid_type;
      int rank;
      int grid;
      int number_of_points;
      double *buffer = NULL;
      int *shape;
      double *spacing;
      double *origin;

      name = input_names[i];

      fprintf (stdout, "\n\n");
      fprintf (stdout, "Input variable: %s\n", name.c_str());
      fprintf (stdout, "==================================\n");

      /* Variable information functions */
      CHECK_INTERFACE_FUNC ("Get_var_type", var_type = child.GetVarType (name));
      CHECK_INTERFACE_FUNC ("GetVarGrid", grid = child.GetVarGrid (name));
      CHECK_INTERFACE_FUNC ("Get_var_rank", rank = child.GetGridRank (grid));

      CHECK_INTERFACE_FUNC ("GetGridNodeCount", number_of_points = child.GetGridNodeCount (grid));

      buffer = (double*) malloc (sizeof (double) * number_of_points);

      /* Variable setters */
      CHECK_INTERFACE_FUNC ("SetValue", child.SetValue (name, (void*)buffer));
      free (buffer);
    }

    for (i=0; i<output_var_count; i++) {
      std::string var_type;
      std::string units;
      std::string grid_type;
      int rank;
      int grid;
      int number_of_points;
      double *buffer = NULL;
      int *shape;
      double *spacing;
      double *origin;

      name = output_names[i];

      fprintf (stdout, "\n\n");
      fprintf (stdout, "Output variable: %s\n", name.c_str());
      fprintf (stdout, "==================================\n");

      /* Variable information functions */
      var_type = child.GetVarType (name);
      grid = child.GetVarGrid (name);
      rank = child.GetGridRank (grid);
      number_of_points = child.GetGridNodeCount (grid);

      fprintf (stdout, "            Type: %s\n", var_type.c_str());
      fprintf (stdout, "            Rank: %d\n", rank);
      fprintf (stdout, "Number of points: %d\n", number_of_points);
    }

    for (i=0; i<output_var_count; i++) {
      std::string var_type;
      std::string units;
      std::string grid_type;
      int rank;
      int grid;
      int number_of_points;
      int *shape;
      double *spacing;
      double *origin;

      name = output_names[i];

      fprintf (stdout, "\n\n");
      fprintf (stdout, "Output variable: %s\n", name.c_str());
      fprintf (stdout, "==================================\n");

      /* Variable information functions */
      CHECK_INTERFACE_FUNC ("Get_var_type", var_type = child.GetVarType(name));
      CHECK_INTERFACE_FUNC ("Get_var_grid", grid = child.GetVarGrid (name));
      CHECK_INTERFACE_FUNC ("Get_var_rank", rank = child.GetGridRank (grid));

      CHECK_INTERFACE_FUNC ("GetGridNodeCount", number_of_points = child.GetGridNodeCount (grid));

      if (number_of_points > 0) {
        double *buffer = (double*) malloc (sizeof (double) * number_of_points);

        /* Variable getters */
        CHECK_INTERFACE_FUNC ("GetValue", child.GetValue (name, (void*)buffer));

        /* Variable setters */
        free (buffer);
      }
    }
  }

  child.Finalize ();

  fprintf (stdout, "\n\n");

  return EXIT_SUCCESS;
}

