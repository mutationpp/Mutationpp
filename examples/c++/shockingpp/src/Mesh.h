#ifndef MESH_H
#define MESH_H

#include <string>

class Mesh {
public:
    Mesh(const std::string& l_mesh_options);
    ~Mesh(){}

    double getXinitial();
    double getXfinal();
    double getdX();
    double getdXprint();

private:
    double m_Xinitial;
    double m_Xfinal;
    double m_dX;
    double m_dX_print;
};

#endif /* MESH_H */
