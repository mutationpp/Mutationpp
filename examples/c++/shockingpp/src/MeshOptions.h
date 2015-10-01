#ifndef MESHOPTIONS_H
#define MESHOPTIONS_H

#include <fstream>
#include <string>
#include <vector>

class MeshOptions{
public:
    MeshOptions( std::ifstream& l_input_file );

    double getMeshX0(){ return m_x_init; }
    double getMeshXend(){ return m_x_end; }
    double getMeshdx(){ return m_dx; }

    int getNumberMeshPoints(){ return n_mesh_points; }
    std::vector<double> getMesh(){ return m_mesh; }

private:
    std::string line;

    double m_x_init;
    double m_x_end;
    double m_dx;

    int n_mesh_points;

    std::vector<double> m_mesh;
    
    void createMesh();

};

#endif // MESHOPTIONS_H
