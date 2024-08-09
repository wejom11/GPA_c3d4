#include <stdio.h>
#include "solver.h"

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