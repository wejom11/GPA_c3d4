#include <string.h>
#include <fstream>
#include <vector>
#include <math.h>
#include "c3d4_ele.h"
#ifndef READ_H
#define READ_H

struct set{
    std::string name;
    std::vector<int> sets;
};

template <typename T1, typename T2, typename T3>
struct Boundary
{
    T1 nodes;
    std::vector<T3> dims;
    T2 value;
    bool is_zero = false;
    bool is_alldim = false;
};


/// @brief read all nodes' coordinates
/// @param file_stream *.inp file name
/// @return last line string where program exit
std::string read_coord(std::ifstream &file_stream, std::vector<double*> &xyz_c);

/// @brief read elements' node tags 
/// @param file_stream *.inp file name
/// @param eles c3d4 elements sets
/// @return last line string where program exit
std::string read_Element(std::ifstream &file_stream, std::vector<c3d4> &eles);

/// @brief read node sets
/// @param file_stream *.inp file name
/// @param ELset node sets
/// @param name node set name
/// @return last line string where program exit
std::string read_Nset(std::ifstream &file_stream, std::vector<set> &Nset, std::string name);

/// @brief read element sets
/// @param file_stream *.inp file name
/// @param ELset element sets
/// @param name element set name
/// @return last line string where program exit
std::string read_Elset(std::ifstream &file_stream, std::vector<set> &ELset, std::string name);

/// @brief read materials
/// @param file_stream *.inp file name
/// @param mat_lib material library
/// @param name material name
/// @return last line string where program exit
std::string read_mater(std::ifstream &file_stream, std::vector<Material> &mat_lib, std::string name);

/// @brief read boundary conditions
/// @param file_stream *.inp file name
/// @param bsetmap set boundary mapping vetor
/// @param bnode node boundary conditions
/// @return last line string where program exit
std::string read_boundary(std::ifstream &file_stream, std::vector<Boundary<std::string,double,int>> &bsetmap, std::vector<Boundary<int,double,int>> &bnode);

/// @brief read load conditions
/// @param file_stream *.inp file name
/// @param lsetmap set load mapping vetor
/// @param lnode node load conditions
/// @return last line string where program exit
std::string read_cload(std::ifstream &file_stream, std::vector<Boundary<std::string,double,int>> &lsetmap, std::vector<Boundary<int,double,int>> &lnode);

/// @brief read load conditions
/// @param file_stream *.inp file name
/// @param dlsetmap set load mapping vetor
/// @param dlele node load conditions
/// @return last line string where program exit
std::string read_dload(std::ifstream &file_stream, std::vector<Boundary<std::string,double,double>> &dlsetmap, std::vector<Boundary<int,double,double>> &dlele);

/// @brief read node measure value
/// @param file_stream *.inp file name
/// @param pb_list measured node tag
/// @param u_pb_val measured node displacement
/// @return last line string where program exit
std::string read_measure(std::ifstream &file_stream, std::vector<int> &pb_list, std::vector<double> &u_pb_val);

/// @brief read a line of the specified file, skip the comment.
/// @param file_stream specified file stream
/// @param str line string
/// @return is reach end of file
bool getlmsg(std::ifstream &file_stream, std::string &str);

/// @brief read the first word of a string
/// @param str string without first word (output)
/// @param word first word (output)
void first_wd(std::string &str, std::string &word);

#endif