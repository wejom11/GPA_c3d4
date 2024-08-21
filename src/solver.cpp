#include <stdio.h>
#include "solver.h"

void SparseMatrix::show(int row){
    int i, j, k, tiscol, lstcol;

    for(i = 0; i < row; i++){
        printf("{");
        lstcol = 1;
        for(j = this->row_st[i]-1; j < this->row_st[i+1]-1; j++){
            tiscol = this->col[j];
            for(k = lstcol; k < tiscol; k++){
                printf("%.2f, ", 0.);
            }
            if(tiscol == row){
                printf("%.8f},", this->val[j]);
            }
            else{
                printf("%.8f, ", this->val[j]);
            }
            lstcol = tiscol + 1;
        }

        for(k = lstcol - 1; k < row; k++){
            if(k == row - 1){
                printf("%.2f}, ",0.);
            }
            else{
                printf("%.2f, ",0.);
            }
        }
    }
};

void pardiso_cfg::initial(int &type){
    int i = 0;
    maxfct = 1;
    mnum = 1;
    msglvl = 0;
    for(i = 0; i < 64; i++){
        iparm[i] = 0;
    }
    if(type == 2){
        iparm[0] = 1;
        iparm[1] = 2;
        iparm[9] = 13;
    }
    else if(type == 11 || type == 1){
        iparm[0] = 1;
        iparm[1] = 2;
        iparm[9] = 13;
    }
    else{
        printf("matrix type %i haven't supported yet!", type);
    }
}

int match(const int row, const int col,  SparseMatrix &SPM){
    int start = SPM.row_st[row-1] - 1;
    int end = SPM.row_st[row] - 1;
    int mid = (start + end) / 2;
    int ptr = mid;
    bool is_done = false;

    while(!is_done){
        if(SPM.col[ptr] < col){
            start = mid;
            mid = (start + end) / 2;
            if(mid == ptr){
                is_done = true;
            }else{
                ptr = mid;
            }
        }
        else if(SPM.col[ptr] > col){
            end = mid;
            mid = (start + end) / 2;
            if(mid == ptr){
                is_done = true;
            }
            else{
                ptr = mid;
            }
        }
        else{
            return ptr;
        }
    }
    
    return -1;
};