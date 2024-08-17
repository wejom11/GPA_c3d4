#ifndef P9SF_ELE_H
#define P9SF_ELE_H
#include <vector>
#include <mkl.h>
#include <math.h>

struct Material{
    double beta = 0.0241;
    double alpha = 0.632185/994/4174;
    double T0 = 25;
    double rho = 994;
    double c1 = -1e-5;
    double c2 = 0.0012;
};

// defination of 9 nodes is demonstrate here:
// 4 —— 7 —— 3
// |    |    |
// 8 —— 9 —— 6
// |    |    |
// 1 —— 5 —— 2
class P9SF{
private:
    bool is_fuild;
    int fd_ele; // freedom of this ele
    int Nodetag[9]{0};  // the tag of 9 nodes
    Material* mat;

    std::vector<double*> XY; // Coordinate of 9 nodes
    std::vector<double*> T; // Temperanture of 9 nodes (for simplicity, call T value)
    std::vector<double*> UV; // Velocity of 9 nodes (U value or V value)
    std::vector<double*> P; // Pressure of 1-4 nodes (P value)
public:
    
    P9SF(Material* mater = nullptr, bool is_fuild = true){
        this->is_fuild = is_fuild;
        fd_ele = is_fuild ? 31 : 9;
        if(is_fuild){
            T.resize(9, nullptr);
            UV.resize(9, nullptr);
            P.resize(4, nullptr);
        }
        else{
            T.resize(9, nullptr);
        }
        XY.resize(9, nullptr);
        mat = mater;
    }

    /// @brief Determine whether the element is valid
    /// @return if valid, return true
    bool is_valid();

    /// @brief determinate whether the element if fuild or solid
    /// @return if ture, it's fuild; 
    bool is_f();

    /// @brief set this elements' type to fuild-1 or solid-0
    /// @param is_f 1-fuild and 0-solid
    void set_fs(bool is_f);

    /// @brief add node's coordinate and T, U, V, P value to its node
    /// @param ntag the tag of 9 nodes in this element
    /// @param xyz coordinate of all nodes
    /// @param uvt_p UVTP value of all nodes
    /// @param tn_num specify the i-th_type_node num (by order: 1-st: nodes contain 
    ///               U,V,T,P; 2-nd: nodes contain U,V,T; 3-rd: nodes contain T{solid element})
    void set_nval(int* ntag, double* xyz, double* uvtp, int* tn_num);

    /// @brief calculate the shape function of U,V,T value of one element with 9 Gauss points (N[9 * 9], 2-order)
    /// and the shape function of P value of element with 9 Gauss points (Np[4 * 9], 1-order)
    /// @param N the shape function of U,V,T value
    /// @param Np the shape function of P value
    /// @param weight the weightness in 9 Gauss Points
    void get_N_Np(double N[9 * 9], double Np[9 * 4], double weight[9]);
    
    /// @brief calculate the shape function of U,V,T value of one element with 9 Gauss points (N[9 * 9], 2-order)
    /// @param N the shape function of U,V,T value
    /// @param weight the weightness in 9 Gauss Points
    void get_N(double N[9 * 9], double weight[9]);

    /// @brief calculate the derivative of shape function N against natural coordinate {\xi, \eta} in 9 Gauss Point. 
    /// @param Ndl the derivative of shape function N
    void get_Ndl(double Ndl[9 * 9 * 2]);

    /// @brief calculate the derivative of shape function N against local coordinate {x, y} in 9 Gauss Point.
    /// @param Ndx (output) the derivative of shape function N against local coordinate
    /// @param Det (output) the determinate of Jacobian with N in 9 Gauss Points
    void get_Ndx(double Ndx[9 * 9 * 2], double Det[9]);

    /// @brief get U,V,T,P value field val in 9 Gauss Point {fuild version}
    /// @param UV_fld U,V field value in 9 Gauss Point
    /// @param T_fld T field value in 9 Gauss Point
    /// @param P_fld P field value in 9 Gauss Point
    /// @param N shape function
    /// @param Np P value shape function
    void get_UVTP_field(double UV_fld[9 * 2], double T_fld[9], double P_fld[9],
                        const double N[9 * 9], const double Np[9 * 4]);

    /// @brief get U,V,T,P value field val in 9 Gauss Point {solid version}
    /// @param T_fld T field value in 9 Gauss Point
    /// @param N shape function
    void get_UVTP_field(double T_fld[9], const double N[9 * 9]);

    /// @brief get the derivative of U,V,T,P value field val in 9 Gauss Point {fuild version}
    /// @param UVdx_fld the derivative of U,V field value in 9 Gauss Point
    /// @param Tdx_fld the derivative of T field value in 9 Gauss Point
    /// @param Ndx the derivative of shape function
    void get_UVTdx_field(double UVdx_fld[9 * 2 * 2], double Tdx_fld[9 * 2], const double Ndx[9 * 9 * 2]);

    /// @brief get the derivative of U,V,T,P value field val in 9 Gauss Point {solid version}
    /// @param Tdx_fld the derivative of T field value in 9 Gauss Point
    /// @param Ndx the derivative of shape function
    void get_UVTdx_field(double Tdx_fld[9 * 2], const double Ndx[9 * 9 * 2]);

    /// @brief calculate the element load vector Fint {fluid element}
    /// @param F the element load vector
    /// @param N shape function
    /// @param UV_fld U,V field value in 9 Gauss Point
    /// @param UVdx_fld the derivative of U,V field value in 9 Gauss Point
    /// @param T_fld T field value in 9 Gauss Point
    /// @param Tdx_fld the derivative of T field value in 9 Gauss Point
    /// @param P_fld P field value in 9 Gauss Point
    /// @param Np P value shape function
    /// @param Ndx the derivative of shape function
    /// @param Det the determinate of Jacobian with N in 9 Gauss Points
    /// @param weight the weightness in 9 Gauss Points
    void get_Fint_f(double F[31], const double N[9 * 9], const double UV_fld[9 * 2], const double UVdx_fld[9 * 2 * 2],
                    const double T_fld[9], const double Tdx_fld[9 * 2], const double P_fld[9], const double Np[9 * 4], 
                    const double Ndx[9 * 9 * 2], const double Det[9], const double weight[9]);

    /// @brief calculate the element stiffness matrix K {fluid element}
    /// @param K the element stiffness matrix K
    /// @param N shape function
    /// @param UV_fld U,V field value in 9 Gauss Point
    /// @param UVdx_fld the derivative of U,V field value in 9 Gauss Point
    /// @param T_fld T field value in 9 Gauss Point
    /// @param Tdx_fld the derivative of T field value in 9 Gauss Point
    /// @param P_fld P field value in 9 Gauss Point
    /// @param Pdx_fld the derivative of P field value in 9 Gauss Point
    /// @param Np P value shape function
    /// @param Ndx the derivative of shape function
    /// @param Npdx the derivative of P value shape function
    /// @param Det the determinate of Jacobian with N in 9 Gauss Points
    /// @param weight the weightness in 9 Gauss Points
    void get_K_f(double K[31][31], const double N[9 * 9], const double UV_fld[9 * 2], const double UVdx_fld[9 * 2 * 2],
                    const double T_fld[9], const double Tdx_fld[9 * 2], const double P_fld[9], const double Np[9 * 4], 
                    const double Ndx[9 * 9 * 2], const double Det[9], const double weight[9]);

    /// @brief calculate the element load vector Fint {solid element}
    /// @param F the element load vector
    /// @param N shape function
    /// @param T_fld T field value in 9 Gauss Point
    /// @param Tdx_fld the derivative of T field value in 9 Gauss Point
    /// @param Ndx the derivative of shape function
    /// @param Det the determinate of Jacobian with N in 9 Gauss Points
    /// @param weight the weightness in 9 Gauss Points
    void get_Fint_s(double F[9], const double N[9 * 9], const double T_fld[9], const double Tdx_fld[9 * 2], 
                 const double Ndx[9 * 9 * 2], const double Det[9], const double weight[9]);

    /// @brief calculate the element stiffness matrix K {solid element}
    /// @param K the element stiffness matrix K
    /// @param N shape function
    /// @param T_fld T field value in 9 Gauss Point
    /// @param Tdx_fld the derivative of T field value in 9 Gauss Point
    /// @param Ndx the derivative of shape function
    /// @param Det the determinate of Jacobian with N in 9 Gauss Points
    /// @param weight the weightness in 9 Gauss Points
    void get_K_s(double K[9][9], const double N[9 * 9], const double T_fld[9], const double Tdx_fld[9 * 2], 
                 const double Ndx[9 * 9 * 2], const double Det[9], const double weight[9]);

    /// @brief calculate the element stiffness matrix K and Load vector F {fluid element}
    /// @param K the element stiffness matrix K
    /// @param F the element load vector
    void KF_f(double K[31][31], double F[31]);

    /// @brief calculate the element stiffness matrix K and Load vector F {solid element}
    /// @param K the element stiffness matrix K
    /// @param F the element load vector
    void KF_s(double K[9][9], double F[9]);

    /// @brief an interface to access the value of Nodetag
    /// @param index the i-th element of Nodetag
    /// @return value
    int at_Nodetag(int index);

    ~P9SF(){};

};

#endif