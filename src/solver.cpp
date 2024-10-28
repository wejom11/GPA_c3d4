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

int match(const int row, const int col, SparseMatrix &SPM){
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

int check(const int val, const std::vector<int> &array){
    bool done = false;
    int front_ptr = 0;
    int back_ptr = array.size() - 1;
    int ptr = (front_ptr + back_ptr) / 2;

    if(array.back() == val){
        return back_ptr;
    }
    
    while(!done){
        if(array.at(ptr) == val){
            return ptr;
        }
        else if(array.at(ptr) > val){
            back_ptr = ptr;
            ptr = (front_ptr + back_ptr) / 2;
            if(ptr == back_ptr){
                break;
            }
        }
        else{
            front_ptr = ptr;
            ptr = (front_ptr + back_ptr) / 2;
            if(ptr == front_ptr){
                break;
            }
        }
    }
    return -1;
    
};

std::vector<double> get_normal(double* ptA, double* ptB){
    double x = 0., y = 0., len = 0.0;

    x = ptB[0] - ptA[0];
    y = ptB[1] - ptA[1];
    len = sqrt(pow(x,2) + pow(y,2));
    x = x / len;
    y = y / len;

    std::vector<double> norm = {y, -x, len};
    return norm;

};