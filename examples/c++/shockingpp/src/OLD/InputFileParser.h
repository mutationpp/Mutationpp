#ifndef INPUTFILEPARSER_H
#define INPUTFILEPARSER_H

#include <fstream>
#include <string>

#include "mutation++.h"

#include "InputData.h"
#include "MeshOptions.h"

class InputFileParser{

public:
    InputFileParser( const std::string& l_input_file_name );
    ~InputFileParser();

    Mutation::Mixture* getMixture(){ return mp_mixture; }
    InputData* getInputData(){ return mp_input_data; }
    MeshOptions* getMeshOptions(){ return mp_mesh_options; }

private:
    Mutation::MixtureOptions* mp_mixture_options;
    Mutation::Mixture* mp_mixture;

    InputData* mp_input_data;
    MeshOptions* mp_mesh_options;

    std::string m_input_file_name;
    std::ifstream m_input_file;
    std::string line;

    std::string m_mixture;
    std::string m_state_model;
    std::string m_thermo_db;

    inline void locateInputFile();

};

#endif // INPUTFILEPARSER_H
