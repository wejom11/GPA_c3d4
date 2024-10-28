#ifndef SOLVER_H
#define SOLVER_H
#include <stdio.h>
#include <vector>
#include <math.h>

class SparseMatrix{
public:
    int type;
    double* val;
    int* col;
    int* row_st;

    SparseMatrix(int mtype = 2){
        type = mtype;
        val = nullptr;
        col = nullptr;
        row_st = nullptr;
    }

    void del(){
        delete[] val; val = nullptr;
        delete[] col; col = nullptr;
        delete[] row_st; row_st = nullptr;
    }

    /// @brief show this SparseMatrix
    /// @param row total rows of this matrix
    void show(int row);

    ~SparseMatrix(){
        if(!(val == nullptr && col == nullptr && row_st == nullptr)){
            this -> del();
        }
    }
};

class pardiso_cfg{
public:
    int maxfct;
    int mnum;
    int perm;
    int iparm[64];
    int msglvl;
    int error;

    pardiso_cfg(int &mtype){
        initial(mtype);
    };

    /// @brief intial the parameter of Intel Pardiso solver
    void initial(int &type);

    ~pardiso_cfg(){};
};

/// @brief find (i,j) Matrix element's position in Sparse Matrix storage format 'val' array
/// @param row the i-th row dimension
/// @param col the j-th col dimension
/// @attention the row/col number and SparseMatrix is one-based indexing;
/// @return the position where this element stored in SparseMatrix.val Array.
int match(const int row, const int col,  SparseMatrix &SPM);

/// @brief fine val's position in the array, if doesn't exist, return -1.
/// @param val val
/// @param array array
/// @return position
/// @attention the array should be monotonically increasing array(Mathematically,
///            this array must be a strictly defined set)
int check(const int val, const std::vector<int> &array);

/// @brief get the normal vector of line AB which point to the right side of AB
/// @param ptA point A
/// @param ptB point B
/// @return normal vector
std::vector<double> get_normal(double* ptA, double* ptB);

#endif