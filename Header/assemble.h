#ifndef ASSEMBLE_H
#define ASSEMBLE_H
#include <vector>
#include <iostream>
#include <string.h>
#include <mkl.h>
#include "solver.h"
#include "read.h"
#include "c3d4_ele.h"

class asb_manager
{
public:
    // problem information
    std::string inp_file_name;
    // work matrix
    SparseMatrix K_mat;
    double* Fout;
    // work information
    std::vector<c3d4> elements;
    std::vector<double*> xyz_coord;
    std::vector<set> Nset;
    std::vector<set> ELset;
    std::vector<std::pair<std::string,std::string>> mat_map;
    std::vector<Material> Mater_lib;
    std::vector<Boundary<std::string,double,int>> bound_set_map;
    std::vector<Boundary<set*,double,int>> bound_set;
    std::vector<Boundary<int,double,int>> bound_node;

    std::vector<Boundary<std::string,double,int>> load_set_map;
    std::vector<Boundary<set*,double,int>> load_set;
    std::vector<Boundary<int,double,int>> load_node;

    std::vector<Boundary<std::string,double,double>> dload_set_map;
    std::vector<Boundary<set*,double,double>> dload_set;
    std::vector<Boundary<int,double,double>> dload_ele;
    // solved answer
    double* uvw_ans;

    std::vector<std::pair<int,double>> sub_list;
    std::vector<int> fnode_list;
    std::vector<int> probe_list;
    std::vector<double> u_real_probe;

// public:

    asb_manager(std::string file_name){
        inp_file_name = file_name;
        uvw_ans = NULL;
        Fout = NULL;
    }

    /// @brief initial elements informations
    void init_ele();

    /// @brief initial the boudary conditions
    void init_bnd();

    /// @brief initial K_{ij}
    void init_K();

    /// @brief initial work matrix (K_{ij}, F_i, ...)
    void initialize(int mode = 1);

    /// @brief read *.inp file
    /// @param file_stream *.inp file stream
    void read_manager(std::ifstream &file_stream);

    /// @brief get the stifness matrix of the element and assemble
    /// @param ele the c3d4 type element
    /// @param invcol the (i,j)th element's Stiffness Matrix cell position in value array or colume array
    void asb_K(c3d4 &ele, int* invcol);

    void getFout();

    /// @brief find the Nset's or ELset's position whose name is ${name}
    /// @param name Nset's or ELset's name
    /// @param mode Fine Nset -> "N" or ELset -> "EL", default "N"
    /// @return set's position (set*)
    set* whereset(const std::string &name, const char* mode);

    /// @brief find the Material's position whose name is ${name}
    /// @param name Material's name
    /// @return Material's position (Material*)
    Material* whereset(const std::string &name);

    /// @brief modify the K_{ij} to add the boundary conditions
    /// @param mode method to add boundary conditions, default 1
    /// 1 - by multiple large numbers; 2 - direct substitution;
    void addboundry(int mode = 1);

    void solve();

    void sub_mat_vec(bool F = true);

    /// @brief initial the symbolic Stiffness
    void init_K_symbolic();

    ~asb_manager(){
        for(int i = 0; i < xyz_coord.size(); i++){
            delete[] xyz_coord.at(i);
        }
        K_mat.del();
        delete[] Fout; Fout = NULL;
        delete[] uvw_ans; uvw_ans = NULL;
    };
};

#endif