#include <vector>
#include <iostream>
#include <string.h>
#include <mkl.h>
#include <fstream>
#include <stdio.h>
#include "solver.h"
#include "P9SF_ele.h"
#include "mesh.h"
#include "read.h"

#ifndef ASSEMBLE_H
#define ASSEMBLE_H

class asb_manager
{
public:
    // problem information
    // std::string inp_file_name;
    double C2 = 8e7;
    double C1 = 2E5;
    // work matrix
    SparseMatrix K_mat;
    double* Fint;
    double* Fout;
    // work information
    double P_in;
    std::vector<P9SF> elements;
    double* xy_coord;
    std::vector<Material> Mater_lib;
    int typenode_num[3];

    std::vector<set> Nset;
    std::vector<set> ELset;

    std::vector<std::pair<int,int>> inlet_ele;
    std::vector<std::pair<int,int>> outlet_ele;
    std::vector<int> wall_node;
    std::vector<int> face_node;
    // solved answer
    double* uvtp_ans;

    double fsinterval;
    std::vector<double> ifw;
	std::vector<double> ofw;

// public:

    asb_manager(/* std::string file_name */){
        // inp_file_name = file_name;
        P_in = 0.01;
        xy_coord = nullptr;
        Fint = nullptr;
        uvtp_ans = nullptr;
        Fout = nullptr;
    }

    /// @brief read *.inp file
    /// @param file_stream *.inp file stream
    void read_manager(std::ifstream &file_stream);

    /// @brief initial elements informations
    void init_ele();

    void init_mesh();

    /// @brief initial the boudary conditions
    // void init_bnd();

    /// @brief initial K_{ij}
    void init_KF();

    /// @brief initial work matrix (K_{ij}, F_i, ...)
    void initialize(int mode = 1);

    /// @brief get the stifness matrix of the element and assemble
    /// @param ele the c3d4 type element
    /// @param invcol the (i,j)th element's Stiffness Matrix cell position in value array or colume array
    void asb_KF_f(P9SF &ele, int* invcol = nullptr);
    void asb_KF_s(P9SF &ele, int* invcol = nullptr);

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

    /// @brief get the row/col of this point's U/V/T/P value in the global stiffness matrix
    /// @param nodetag point tag
    /// @param uvtp 0-U value; 1-V value; 2-T value; 3-P value
    /// @return row/col
    int wherend(int nodetag, int uvtp);

    /// @brief modify the K_{ij} to add the boundary conditions
    void addboundry();

    /// @brief solve the answer with single step.
    /// @return if true, the answer astringency to the real answer
    bool solve();

    /// @brief initial the symbolic Stiffness
    void init_K_symbolic();

    void get_face(int m, int x_numele, int y_numele);

    /// @brief write answer to *.out file
    void write(std::string file_name = "demo");

    ~asb_manager(){
        delete[] xy_coord; xy_coord = nullptr;
        K_mat.del();
        delete[] Fout; Fout = nullptr;
        delete[] uvtp_ans; uvtp_ans = nullptr;
        delete[] Fint; Fint = nullptr;
    };
};

#endif