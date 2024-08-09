#ifndef C3D4_H
#define C3D4_H
#include <vector>
#include <string>
#include <mkl.h>

struct Material{
    std::string name;
    short int where; // start from 0
    double E;
    double mu;
    double rho;
};

/// @brief c3d4 element(1-order)
class c3d4{
public:
    std::vector<double*> xyz0;
    Material *mat;
    std::vector<int> Nodetag;
    double vol;
// public:
    /// @brief initial the element
    /// @param xyz element nodal coordinate
    /// @param e Young's Modules
    /// @param m Poisson ratio
    /// @param nodet the element nodes numbering
    c3d4(std::vector<double*> xyz = {}, std::vector<int> nodet = {}, Material *Mat = NULL){
        xyz0 = xyz;
        mat = Mat;
        Nodetag = nodet;
    };

    /// @brief calculate and save the volume of this element
    /// @return volume
    double getvolume();

    /// @brief get the stiffness matrix and \frac{\partial K_{ij}}{\partial E} of the element
    /// @param K element stifness K^e
    /// @param KdE element \frac{\partial K_{ij}^e}{\partial E}
    void getK_e(double K[][12], double KdE[][12]);

    /// @brief only calculate element stiffness
    /// @param K element stifness K^e
    void getK_e(double K[][12]);

    /// @brief calculate the \frac{\partial \sigma_{ij}}{\partial u_k} and
    /// the variation of strain and the \frac{\partial \sigma_{ij}}{\partial u_k \partial E}
    /// @param deps the variation of strain
    /// @param sigdudE \frac{\partial \sigma_{ij}}{\partial u_k \partial E}
    /// @return Tensor (i,j,k)
    void getdSigdu(double deps[][3][12], double sigdudE[][3][12], double sigdu[][3][12]);

    /// @brief only calculate the variation of strain and \frac{\partial \sigma_{ij}}{\partial u_k}
    /// @param deps the variation of strain
    /// @return Tensor (i,j,k)
    void getdSigdu(double deps[][3][12], double sigdu[][3][12]);

    ~c3d4(){
    };
};

#endif