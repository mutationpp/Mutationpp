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
    readInputFile(s_file_input);

    // Check Input Line

    // HARDCODE INPUT FILE FOR THE TIME BEING
    // s_problem_type = "shocking";
    // s_mixture = "air5";
    // s_state_model = "ChemNonEqTTv";
    // s_thermo_db = "RRHO";
    // s_mesh = "NotEmpty!";
    // s_free_stream_conditions = "NotEmpty!";

    // Check Members not Empty
    errorInputFileInaccurate();

}

void SetupProperties::errorInputFileInaccurate(){
    if (s_problem_type.compare("empty") == 0) {
        std::cerr << "In the input file a problem type should be provided!" << std::endl;
        exit(1);
    }

    if (s_mixture.compare("empty") == 0) {
        std::cerr << "In the input file a mixture should be provided!" << std::endl;
        exit(1);
    }
    if (s_state_model.compare("empty") == 0) {
        std::cerr << "In the input file a state_model should be provided!" << std::endl;
        exit(1);
    }
    if (s_thermo_db.compare("empty") == 0) {
        std::cerr << "In the input file a thermodynamic database should be provided!" << std::endl;
        exit(1);
    }
    
    if (s_problem_type.compare("shocking") == 0) { // SHOCKING also needs other parameters!
        if (s_mesh.compare("empty") == 0) {
            std::cerr << "In the input file for Shocking, a mesh should be provided!" << std::endl;
            exit(1);
        }
        if (s_free_stream_conditions.compare("empty") == 0) {
            std::cerr << "In the input file for Shocking, free stream conditions ";
            std::cerr << "should be provided!" << std::endl;
            exit(1);
        }
    }
}

void SetupProperties::cleanMeshInfo(){s_mesh.clear();}
void SetupProperties::cleanFreeStreamConditionsInfo(){s_free_stream_conditions.clear();}

// ------------------------------------------------------------------------------------

std::string SetupProperties::parseline(std::string line)
{
  // Remove comments starting somewhere in the line
  for(size_t ii = 0; ii < line.length(); ++ii) {
    if(line[ii] == '#') { // then the rest is a comment, trim it!
      line = line.substr(0,ii);
    }
  }
  // Removing whitespaces in the beginning of the line
  while( line[0] == ' ' ) {
    line = line.substr(1, line.length());
  }
  // Removing whitespaces in the end of the line
  while( line[line.length()-1] == ' ' ) {
    line = line.substr(0, line.length()-1);
  }
  return line;
};

// ------------------------------------------------------------------------------------

void SetupProperties::readInputFile(std::string& s_file_input)
{

  std::cout << "Loading Input file '" << s_file_input << "'" << std::endl;

  // Open file
  const char * filename = s_file_input.c_str();

  std::ifstream filein;
  filein.open(filename);
  if (!filein) {
      std::cout << "ERROR: couldn't open input file." << std::endl;
      exit(1);
  }

  // Initializing variables to value "empty"
  s_mixture      = "empty";
  s_state_model  = "empty";
  s_thermo_db    = "empty";

  // Auxiliary variables
  std::string line;
  std::vector<std::string> linesVect;

  // read lines
  while(getline(filein, line)){

      // Parse the line
      line = parseline(line);

      // If the line is not empty, save it in a vector of strings
      if(line.length() > 0) {
         linesVect.push_back(line);
      }

  }

  // Check for known syntax into the vector of saved lines
  std::string lineNow;
  for(size_t ii = 0; ii < linesVect.size(); ++ii) {

    lineNow = linesVect.at(ii);

    if(lineNow.compare("Problem Type:")==0) {
        s_problem_type = linesVect.at(ii+1);
    }
    if(lineNow.compare("Name of the mixture:")==0) {
        s_mixture = linesVect.at(ii+1);
    }
    if(lineNow.compare("State Model:")==0) {
        s_state_model = linesVect.at(ii+1);
    }
    if(lineNow.compare("Thermodynamic Database:")==0) {
        s_thermo_db = linesVect.at(ii+1);
    }
    if(lineNow.compare("Mesh:") == 0) {
        s_mesh = linesVect.at(ii+1);
    }
    if(lineNow.compare("Free Stream Conditions:") == 0) {
        s_free_stream_conditions = linesVect.at(ii+1);
    }
  }

  // Variables to save read strings
  std::cout << "\nInput file was read: "                  << std::endl;
  std::cout << "    Problem Type: "                       << s_problem_type            << std::endl;
  std::cout << "    Mixture: "                            << s_mixture                 << std::endl;
  std::cout << "    State Model: "                        << s_state_model             << std::endl;
  std::cout << "    Thermodynamic Database: "             << s_thermo_db               << std::endl;
  std::cout << "    Mesh (TEMPORARY): "                   << s_mesh                    << std::endl;
  std::cout << "    Free Stream Conditions (TEMPORARY): " << s_free_stream_conditions  << std::endl;
  std::cout << std::endl;

  return;
}

