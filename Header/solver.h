#ifndef SOLVER_H
#define SOLVER_H
#include <stdio.h>

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

#endif