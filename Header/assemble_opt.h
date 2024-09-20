#include <iostream>
#include <string.h>
#include <vector>
#include <mkl.h>
#include <mkl_spblas.h>
#include <GFE_API.h>
#include <GFE_Struct/GFE_Outp.h>
#include "assemble.h"
#include "solver.h"

#ifndef ASSEMBLE_OPT_H
#define ASSEMBLE_OPT_H

class asb_opt_manager: public asb_manager{
public:
    SparseMatrix* KdE_mat;
    double* udE;
    
// public:
    asb_opt_manager(std::string name):asb_manager(name){
        KdE_mat = nullptr;
        udE = nullptr;
    };

    /// @brief initial K_{ij} \frac{\partial K_{ij}}{\partial E_k}
    void init_KE(int* invcolume = nullptr, int* invcolume_E = nullptr, bool is_write = false);

    /// @brief get the stifness matrix of the element, \frac{\partial K_{ij}}{\partial E_k} and assemble
    /// @param ele the c3d4 type element
    /// @param invcol the (i,j)th element's Stiffness Matrix cell position in value array or colume array
    void asb_KE(c3d4 &ele, int* invcol = nullptr, int* invcol_E = nullptr);

    /// @brief initial work matrix (K_{ij}, F_i, \frac{\partial K_{ij}}{\partial E_k}, ...)
    void initialize(int mode = 1);

    /// @brief solve the displacement get the matrix \frac{\partial u_i}{\partial E_j} and 
    /// \frac{\partial u_i}{\partial E_j \partial E_k}
    /// @param get_udE if true, get the matrix \frac{\partial u_i}{\partial E_j}
    /// @param get_udE2 if true, get the matrix \frac{\partial u_i}{\partial E_j \partial E_k}
    /// @param clear_kde if true, clear the KdE_mat when done this procedure
    void solve(bool &alloc_err, bool get_udE = false);

    /// @brief reduced matrix
    /// @param F if true; reduced Fout
    void sub_mat_vec(bool F = true);

    /// @brief reduced matrix or vector by measure nodes
    void sub_measure();

    /// @brief modify the K_{ij} to add the boundary conditions and push back position to sub_list
    /// @param mode method to add boundary conditions, default 1
    /// 1 - by multiple large numbers; 2 - direct substitution;
    void addboundry(int mode);

    /// @brief optimize the Young's Modules to minimize the point function
    /// @param method method (1-Gauss-Newton; 2-Lagrange Multiplier)
    void opt_val(int method = 1);

    void wirte(std::ofstream &file_stream);
    void wirte(std::ofstream &file_stream, const std::vector<int> &node_list);

    /// @brief initial the symbolic Stiffness and the \frac{\partial K_{ij}}{\partial E_k}
    void init_KE_symbolic();

    /// @brief write answer to *.db file
    void write_db();

    /// @brief write geometry information to *.db file
    /// @param db .db file io
    /// @return if done, return true.
    bool write_geo2db(std::shared_ptr<GFE::DB> db);

    /// @brief write displacement information to *.db file
    /// @param db .db file io
    /// @param frame frame number
    /// @return if done, return true.
    bool write_disp2db(std::shared_ptr<GFE::DB> db, int frame);

    ~asb_opt_manager(){
        delete[] udE; udE = nullptr;
        for(int i = 0; i < Mater_lib.size(); i++){
            KdE_mat[i].del();
        }
        delete[] KdE_mat; KdE_mat = nullptr;
    }
};

#endif