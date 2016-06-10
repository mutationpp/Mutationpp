#ifndef SETUPPROPERTIES_H
#define SETUPPROPERTIES_H

#include <string>

class SetupProperties {
public:
    SetupProperties(std::string& s_file_input);
    ~SetupProperties(){}

    std::string getProblemType(){ return s_problem_type; }
    std::string getMixture(){ return s_mixture; }
    std::string getStateModel(){ return s_state_model; }
    std::string getThermoDB(){ return s_thermo_db; }
    std::string getMesh(){ return s_mesh; }
    std::string getFreeStreamConditions(){ return s_free_stream_conditions; }

    void cleanMeshInfo();
    void cleanFreeStreamConditionsInfo();

private:
    std::string s_problem_type;
    std::string s_mixture;
    std::string s_state_model;
    std::string s_thermo_db;

    std::string s_mesh;
    std::string s_free_stream_conditions;

    inline void errorInputFileInaccurate();

};

#endif /* SETUPPROPERTIES_H */
