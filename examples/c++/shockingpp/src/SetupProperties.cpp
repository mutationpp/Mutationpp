#include <iostream>
#include <stdlib.h>

#include "SetupProperties.h"

SetupProperties::SetupProperties(std::string& s_file_input)
        : s_problem_type("empty"),
          s_mixture("empty"),
          s_state_model("empty"),
          s_thermo_db("empty"),
          s_mesh("empty"),
          s_free_stream_conditions("empty"){

    // Open File

    // Check Input Line

    // HARDCODE INPUT FILE FOR THE TIME BEING
    s_problem_type = "shocking";
    s_mixture = "air5";
    s_state_model = "ChemNonEqTTv";
    s_thermo_db = "RRHO";
    s_mesh = "NotEmpty!";
    s_free_stream_conditions = "NotEmpty!";

    // Check Members not Empty
    errorInputFileInaccurate();

}

void SetupProperties::errorInputFileInaccurate(){
    if (s_problem_type.compare("empty") == true) {
        std::cerr << "In the input file a problem type should be provided!" << std::endl;
        exit(1);
    }
    if (s_mixture.compare("empty") == true) {
        std::cerr << "In the input file a mixture should be provided!" << std::endl;
        exit(1);
    }
    if (s_state_model.compare("empty") == true) {
        std::cerr << "In the input file a state_model should be provided!" << std::endl;
        exit(1);
    }
    if (s_thermo_db.compare("empty") == true) {
        std::cerr << "In the input file a thermodynamic database should be provided!" << std::endl;
        exit(1);
    }
    if (s_mesh.compare("empty") == true) {
        std::cerr << "In the input file a mesh should be provided!" << std::endl;
        exit(1);
    }
    if (s_free_stream_conditions.compare("empty") == true) {
        std::cerr << "In the input file free stream conditions for should be provided!" << std::endl;
        exit(1);
    }
}

void SetupProperties::cleanMeshInfo(){s_mesh.clear();}

void SetupProperties::cleanFreeStreamConditionsInfo(){s_free_stream_conditions.clear();}

