#include "Mesh.h"

Mesh::Mesh(const std::string& l_mesh_info)
        : m_Xinitial(0.0),
          m_Xfinal(0.0),
          m_dX(0.0),
          m_dX_print(0.0)
{
    m_Xinitial = 0.0;
    m_Xfinal = 5.e-2;
    m_dX = 1.e-8;
    m_dX_print = 1.e-5;
}

double Mesh::getXinitial(){return m_Xinitial; }
double Mesh::getXfinal(){return m_Xfinal;}
double Mesh::getdX(){return m_dX;}
double Mesh::getdXprint(){return m_dX_print;}
