#include <iostream>
#include <string.h>
#include <vector>
#include <mkl.h>
#include <mkl_spblas.h>
#include "assemble.h"
#include "solver.h"
#include "assemble_opt.h"

void asb_opt_manager::init_KE(int* invcolume, int* invcolume_E, bool is_write){
    int dof = 3 * xy_coord.size();
    int opt_var = Mater_lib.size();
    int now_p = 1;
    int to_ele = elements.size();
    int ROW_I, ROW_E;
    int* invcol = nullptr;
    int* invcol_E = nullptr;

    if(K_mat.val == nullptr && KdE_mat == nullptr){
        init_KE_symbolic();
    }
    else{
        for(int i = 0; i < K_mat.row_st[dof]-1; i++){
            K_mat.val[i] = 0.;
        }
        for(int j = 0; j < opt_var; j++){
            for(int i = 0; i < KdE_mat[j].row_st[dof]-1; i++){
                KdE_mat[j].val[i] = 0.;
            }
        }
    }

    if(invcolume != nullptr && invcolume_E != nullptr){
        invcol = invcolume;
        invcol_E = invcolume_E;
    }
    else{
        try{
            invcol = new int[dof * dof];
            invcol_E = new int[opt_var * dof * dof];
        }
        catch(const std::bad_alloc e){
            std::cerr << e.what() << '\n';
            delete[] invcol; invcol = nullptr;
            delete[] invcol_E; invcol_E = nullptr;
        }

        if(invcol != nullptr && invcol_E != nullptr){
            for(int i = 0; i < dof; i++){
                ROW_I = i * dof;
                for(int j = K_mat.row_st[i] - 1; j < K_mat.row_st[i+1] - 1; j++){
                    invcol[ROW_I + K_mat.col[j] - 1] = j;
                }
            }
            for(int k = 0; k < opt_var; k++){
                ROW_E = k * dof * dof;
                for(int i = 0; i < dof; i++){
                    ROW_I = i * dof;
                    for(int j = KdE_mat[k].row_st[i] - 1; j < KdE_mat[k].row_st[i+1] - 1; j++){
                        invcol_E[ROW_E + ROW_I + KdE_mat[k].col[j] - 1] = j;
                    }
                }
            }
        }
    }
    
    for(std::vector<c3d4>::iterator itele = elements.begin(); itele != elements.end(); itele++){
        if(itele->Nodetag.size() == 0){
            continue;
        }
        asb_KE(*itele, invcol, invcol_E);
        printf("\r assembling K: %.2f%%", now_p*100.0/to_ele);
        now_p++;
    }
    printf("\n");

    if(is_write){
        invcolume = invcol;
        invcolume_E = invcol_E;
    }
    else{
        delete[] invcol; invcol = nullptr;
        delete[] invcol_E; invcol_E = nullptr;
    }
};

void asb_opt_manager::asb_KE(c3d4 &ele, int* invcol, int* invcol_E){
    double K[12][12];
    double KdE[12][12];
    int ROW_I, ROW_E;
    int dof = 3 * xy_coord.size();
    ele.getK_e(K,KdE);

    if((K_mat.type == 2 || K_mat.type == -2) && (invcol != nullptr && invcol_E != nullptr)){
        int I, J;
        ROW_E = ele.mat->where * dof * dof;
        for (int i = 0; i < 4; i++){
            for (int j = i; j < 4; j++){
                I = ele.Nodetag[i];
                J = ele.Nodetag[j];
                I = 3*I - 3;
                J = 3*J - 3;
                if(I > J){
                    for(int m = 0; m < 3; m++){
                        ROW_I = (J+m) * dof;;
                        for(int n = 0; n < 3; n++){
                            K_mat.val[invcol[ROW_I+I+n]] += K[3*j + m][3*i + n];
                            KdE_mat[ele.mat->where].val[invcol_E[ROW_E + ROW_I + I+n]] += KdE[3*j + m][3*i + n];
                        }
                    }
                }
                else if(I < J){
                    for(int m = 0; m < 3; m++){
                        ROW_I = (I+m) * dof;
                        for(int n = 0; n < 3; n++){
                            K_mat.val[invcol[ROW_I+J+n]] += K[3*i + m][3*j + n];
                            KdE_mat[ele.mat->where].val[invcol_E[ROW_E + ROW_I + J+n]] += KdE[3*i + m][3*j + n];
                        }
                    }
                }
                else{
                    for(int m = 0; m < 3; m++){
                        ROW_I = (I+m) * dof;
                        for(int n = m; n < 3; n++){
                            K_mat.val[invcol[ROW_I+I+n]] += K[3*i + m][3*i + n];
                            KdE_mat[ele.mat->where].val[invcol_E[ROW_E + ROW_I + I+n]] += KdE[3*i + m][3*i + n];
                        }
                    }
                }
            }
        }
    }
    else if((K_mat.type == 2 || K_mat.type == -2) && (invcol == nullptr || invcol_E == nullptr)){
        int I, J;
        for (int i = 0; i < 4; i++){
            for (int j = i; j < 4; j++){
                I = ele.Nodetag[i];
                J = ele.Nodetag[j];
                I = 3*I - 3;
                J = 3*J - 3;
                if(I > J){
                    for(int m = 0; m < 3; m++){
                        ROW_I = J + m + 1;
                        for(int n = 0; n < 3; n++){
                            K_mat.val[match(ROW_I, I+n+1, K_mat)] += K[3*j + m][3*i + n];
                            KdE_mat[ele.mat->where].val[match(ROW_I, I+n+1, KdE_mat[ele.mat->where])] += KdE[3*j + m][3*i + n];
                        }
                    }
                }
                else if(I < J){
                    for(int m = 0; m < 3; m++){
                        ROW_I = I + m + 1;
                        for(int n = 0; n < 3; n++){
                            K_mat.val[match(ROW_I, J+n+1, K_mat)] += K[3*i + m][3*j + n];
                            KdE_mat[ele.mat->where].val[match(ROW_I, J+n+1, KdE_mat[ele.mat->where])] += KdE[3*i + m][3*j + n];
                        }
                    }
                }
                else{
                    for(int m = 0; m < 3; m++){
                        ROW_I = I + m + 1;
                        for(int n = m; n < 3; n++){
                            K_mat.val[match(ROW_I, I+n+1, K_mat)] += K[3*i + m][3*i + n];
                            KdE_mat[ele.mat->where].val[match(ROW_I, I+n+1, KdE_mat[ele.mat->where])] += KdE[3*i + m][3*i + n];
                        }
                    }
                }
            }
        }
    }
};

void asb_opt_manager::initialize(int mode){
    std::ifstream inp_file(inp_file_name);
    read_manager(inp_file);
    init_ele();
    init_bnd();
    init_KE();
    printf("done assemble K \n");
    // for(int i = 0; i < KdE_mat.at(0).outerSize(); i++){
    //     printf("col: %d ", i);
    //     for(Eigen::SparseMatrix<double>::InnerIterator it(KdE_mat.at(0),i); it; ++it){
    //         printf("%.4f ",it.value());
    //     }
    //     printf("\n");
    // }
    getFout();
    printf("done get Fout \n");
    // for(Eigen::VectorXd::iterator it = Fout.begin(); it != Fout.end(); it++){
    //     printf("%.4f ", *it);
    // }
    // printf("\n");
    addboundry(mode);
    printf("done add boundry \n");
};

void asb_opt_manager::solve(bool &alloc_err, bool get_udE){
    const int dof = (fnode_list.size() != 0) ? fnode_list.size() : 3*xy_coord.size();
    int optvar_num = Mater_lib.size();
    int* row_fst = nullptr;
    int* row_lst = nullptr;
    sparse_matrix_t A;
    try{
        uvw_ans = new double[dof];
        udE = new double[optvar_num * dof];
    }
    catch(std::bad_alloc& e){
        std::cerr << e.what() << '\n';
        delete[] uvw_ans; uvw_ans = nullptr;
        delete[] udE; udE = nullptr;
    }
    if(uvw_ans == nullptr || udE == nullptr){alloc_err = true;return;};

    void* ptsolver[64];
    for(int i = 0; i < 64; i++){
        ptsolver[i] = 0; 
    }
    int nrhs = 1;
    pardiso_cfg para(K_mat.type);
    double ddum;

    printf("start solving ... \n");
    int phase = 11;
    pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(K_mat.type), &phase, &dof, K_mat.val, 
            K_mat.row_st, K_mat.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), &ddum, &ddum, &(para.error));
    if(para.error != 0 ){
        printf ("ERROR during symbolic factorization: %i \n", para.error);
    }
    
    phase = 22;
    pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(K_mat.type), &phase, &dof, K_mat.val, 
            K_mat.row_st, K_mat.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), &ddum, &ddum, &(para.error));
    if(para.error != 0 ){
        printf ("ERROR during numerical factorization: %i \n", para.error);
    }

    phase = 33;
    pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(K_mat.type), &phase, &dof, K_mat.val, 
            K_mat.row_st, K_mat.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), Fout, uvw_ans, &(para.error));
    if(para.error != 0 ){
        printf ("ERROR during solution: %i \n", para.error);
    }
    printf("done solve u \n");

    printf("solving udE ... \n");
    double* KudE = nullptr;
    double* KudE_e = nullptr;
    if(!get_udE){
        phase = -11;
        pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(K_mat.type), &phase, &dof, K_mat.val, 
                K_mat.row_st, K_mat.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), Fout, &ddum, &(para.error));
        return;
    };
    try{
        KudE = new double[dof * optvar_num];
        KudE_e = new double[dof];
    }
    catch(std::bad_alloc &e){
        std::cerr<< e.what() << std::endl;
        delete[] KudE; KudE = nullptr;
        delete[] KudE_e; KudE_e = nullptr;
    }
    if(KudE == nullptr || KudE_e == nullptr){
        alloc_err = true;
        delete[] uvw_ans; uvw_ans = nullptr;
        delete[] udE; udE = nullptr;
        phase = -11;
        pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(K_mat.type), &phase, &dof, K_mat.val, 
                K_mat.row_st, K_mat.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), Fout, &ddum, &(para.error));
        return;
    }

    matrix_descr KdE_des = {SPARSE_MATRIX_TYPE_SYMMETRIC, SPARSE_FILL_MODE_UPPER, SPARSE_DIAG_NON_UNIT};
    row_fst = new int[dof];
    row_lst = new int[dof];
    for(int i = 0; i < optvar_num; i++){
        for(int j = 0; j < dof; j++){
            row_fst[j] = KdE_mat[i].row_st[j];
            row_lst[j] = KdE_mat[i].row_st[j+1];
        }
        mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ONE, dof, dof, row_fst, row_lst,
                                KdE_mat[i].col, KdE_mat[i].val);
        //product
        mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, KdE_des, uvw_ans, 0, KudE_e);
        for(int k = 0; k < dof; k++){
            KudE[i * dof + k] = - KudE_e[k];
        }
        mkl_sparse_destroy(A);
    }
    delete[] row_fst; row_fst = nullptr;
    delete[] row_lst; row_lst = nullptr;
    delete[] KudE_e; KudE_e = nullptr;

    pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(K_mat.type), &phase, &dof, K_mat.val, 
            K_mat.row_st, K_mat.col, &(para.perm), &optvar_num, para.iparm, &(para.msglvl), KudE, udE, &(para.error));
    if(para.error != 0 ){
        printf ("ERROR during solution: %i \n", para.error);
    }
    printf("done solve udE\n");

    phase = -11;
    pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(K_mat.type), &phase, &dof, K_mat.val, 
            K_mat.row_st, K_mat.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), Fout, &ddum, &(para.error));
    delete[] KudE; KudE = nullptr;
    // for(int i = 0; i < dof; i++){
    //     for(int j = 0; j < optvar_num; j++){
    //         printf("%.10f ", udE(i,j));
    //     }
    //     printf("\n");
    // }

};

void asb_opt_manager::sub_mat_vec(bool F){
    // int subsize = (fnode_list.size() != 0) ? fnode_list.size() : 3*xy_coord.size();
    // int opt_p = 0;
    // int total_opt = Mater_lib.size();
    // double val = 0.0;
    // Eigen::SparseMatrix<double> Wk_mat;
    // Eigen::VectorXd F_wk;
    // if(F){F_wk.resize(subsize);}
    // Wk_mat.resize(subsize,subsize);

    // int i = 0;
    // for(std::vector<int>::iterator itfi = fnode_list.begin(); itfi!= fnode_list.end(); itfi++){
    //     if(F){
    //         val = 0.0;
    //         for(std::vector<std::pair<int,double>>::reverse_iterator itsub = sub_list.rbegin(); itsub < sub_list.rend(); ++itsub){
    //             val += K_mat.coeff(*itfi, itsub->first) * itsub->second;
    //         }
    //         F_wk.coeffRef(i) = Fout.coeffRef(*itfi) - val;
    //     }
        
    //     int j = 0;
    //     for(std::vector<int>::iterator itfj = fnode_list.begin(); itfj != fnode_list.end(); itfj++){
    //         val = K_mat.coeff(*itfi, *itfj);
    //         if(abs(val) > 1E-12){
    //             Wk_mat.coeffRef(i,j) = val;
    //         }
    //         j++;
    //     }
    //     i++;
    //     printf("\r K processing: %.2f%%",i*100.0/subsize);
    // }
    // K_mat.resize(subsize,subsize);
    // K_mat = Wk_mat;
    // if(F){
    //     Fout.resize(subsize);
    //     Fout = F_wk;
    // }
    // printf("\n");

    // for(std::vector<Eigen::SparseMatrix<double>>::iterator itkde = KdE_mat.begin(); itkde != KdE_mat.end(); itkde++){
    //     Wk_mat.setZero();
    //     int i = 0;
    //     for(std::vector<int>::iterator itfi = fnode_list.begin(); itfi!= fnode_list.end(); itfi++){
    //         int j = 0;
    //         for(std::vector<int>::iterator itfj = fnode_list.begin(); itfj!= fnode_list.end(); itfj++){
    //             val = itkde->coeffRef(*itfi, *itfj);
    //             if(abs(val) > 1E-12){
    //                 Wk_mat.coeffRef(i,j) = val;
    //             }
    //             j++;
    //         }
    //         i++;
    //         printf("\r K_mat processing: %.2f%%",i*100.0/subsize/total_opt + opt_p*100.0/total_opt);
    //     }
    //     itkde->resize(subsize, subsize);
    //     *itkde = Wk_mat;
    //     opt_p++;
    // }

    // for(Eigen::VectorXd::iterator uit = uans_sub.begin(); uit != uans_sub.end(); uit++){
    //     printf("%.4f ",*uit);
    // }
    // printf("\n");
};

void asb_opt_manager::sub_measure(){
    double* wk_vec = nullptr;
    double* wk_mat = nullptr;

    int moni_size = probe_list.size();
    int optvar_num = Mater_lib.size();
    int dof = 3 * xy_coord.size();
    int ROW_I = 0, ROW_J = 0;
    wk_vec = new double[moni_size];
    wk_mat = new double[moni_size * optvar_num];

    if(fnode_list.size() != 0){
        // rfnode_list.resize(3 * xy_coord.size(), -1);
        // int i = 0;
        // for(std::vector<int>::iterator itp = fnode_list.begin(); itp != fnode_list.end(); itp++){
        //     rfnode_list.at(*itp) = i;
        //     i++;
        // }

        // for(int i = 0; i < optvar_num; i++){
        //     int j = 0;
        //     for(std::vector<int>::iterator itp = probe_list.begin(); itp != probe_list.end(); itp++){
        //         wk_mat.coeffRef(j,i) = udE.coeff(rfnode_list.at(*itp),i);
        //         wk_vec.coeffRef(j) = uvw_ans.coeff(rfnode_list.at(*itp));
        //         j++;
        //     }
        // }
    }
    else{
        int j = 0;
        for(std::vector<int>::iterator itp = probe_list.begin(); itp != probe_list.end(); itp++){
            wk_vec[j] = uvw_ans[*itp];
            j++;
        }
        for(int i = 0; i < optvar_num; i++){
            ROW_I = i * moni_size;
            ROW_J = i * dof;
            j = 0;
            for(std::vector<int>::iterator itp = probe_list.begin(); itp != probe_list.end(); itp++){
                wk_mat[ROW_I + j] = udE[ROW_J + *itp];
                j++;
                printf("\rsub measure : %.2f %%", (j * 100.0 / moni_size / optvar_num + i * 100.0 / optvar_num));
            }
        }
    }

    delete[] uvw_ans; uvw_ans = wk_vec;
    wk_vec = nullptr;
    delete[] udE; udE = wk_mat;
    wk_mat = nullptr;
    printf("\n");

};

void asb_opt_manager::addboundry(int mode){
    if(mode == 1){
        int posi = 0;
        bool isz = false;

        for(std::vector<Boundary<set*,double,int>>::iterator itset = bound_set.begin(); itset != bound_set.end(); itset++){
            isz = itset->is_zero;
            if(itset -> is_alldim){
                for(std::vector<int>::iterator it = itset->nodes->sets.begin(); it != itset->nodes->sets.end(); it++){
                    posi = 3 * (*it - 1);
                    for(int i = 0; i < 3; i++){
                        K_mat.val[K_mat.row_st[posi]-1] = 1E12;
                        Fout[posi] = isz ? 0.0 : ((1E12) * itset->value);
                        posi++;
                    }
                }
            }
            else{
                for(std::vector<int>::iterator it = itset->nodes->sets.begin(); it != itset->nodes->sets.end(); it++){
                    posi = 3 * (*it - 1);
                    for(std::vector<int>::iterator itdim = itset->dims.begin(); itdim != itset->dims.end(); itdim++){
                        K_mat.val[K_mat.row_st[posi+*itdim-1]-1] = 1E12;
                        Fout[posi+*itdim-1] = isz ? 0.0 : ((1E12) * itset->value);
                    }
                }
            }
        }

        for(std::vector<Boundary<int,double,int>>::iterator itnode = bound_node.begin(); itnode != bound_node.end(); itnode++){
            isz = itnode->is_zero;
            if(itnode->is_alldim){
                posi = 3*(itnode->nodes-1);
                for(int i = 0; i < 3; i++){
                    K_mat.val[K_mat.row_st[posi]-1] = 1E12;
                    Fout[posi] = isz ? 0.0 : ((1E12) * itnode->value);
                    posi++;
                }
            }
            else{
                posi = 3*(itnode->nodes-1);
                for(std::vector<int>::iterator itdim = itnode->dims.begin(); itdim != itnode->dims.end(); itdim++){
                    K_mat.val[K_mat.row_st[posi+*itdim-1]-1] = 1E12;
                    Fout[posi+*itdim-1] = isz ? 0.0 : ((1E12) * itnode->value);
                }
            }
        }
    }
    else if(mode == 2){
        printf("ADD_BOUNDARY ERROR: direct substitution method (mode = 2) not available presently!\n\
                    try multiple large numbers method (mode = 1)");
        exit(1);

        int posi = 0;
        bool isz = false;

        fnode_list.clear();
        sub_list.clear();

        for(std::vector<Boundary<set*,double,int>>::iterator itset = bound_set.begin(); itset != bound_set.end(); itset++){
            isz = itset->is_zero;
            if(itset -> is_alldim){
                for(std::vector<int>::iterator it = itset->nodes->sets.begin(); it != itset->nodes->sets.end(); it++){
                    posi = 3 * (*it - 1);
                    for(int i = 0; i < 3; i++){
                        sub_list.push_back({posi, itset->value});
                        posi++;
                    }
                }
            }
            else{
                for(std::vector<int>::iterator it = itset->nodes->sets.begin(); it != itset->nodes->sets.end(); it++){
                    posi = 3 * (*it - 1);
                    for(std::vector<int>::iterator itdim = itset->dims.begin(); itdim != itset->dims.end(); itdim++){
                        sub_list.push_back({posi+*itdim-1, itset->value});
                    }
                }
            }
        }

        for(std::vector<Boundary<int,double,int>>::iterator itnode = bound_node.begin(); itnode != bound_node.end(); itnode++){
            isz = itnode->is_zero;
            if(itnode->is_alldim){
                posi = 3*(itnode->nodes-1);
                for(int i = 0; i < 3; i++){
                    sub_list.push_back({posi, itnode->value});
                    posi++;
                }
            }
            else{
                posi = 3*(itnode->nodes-1);
                for(std::vector<int>::iterator itdim = itnode->dims.begin(); itdim != itnode->dims.end(); itdim++){
                    sub_list.push_back({posi+*itdim-1, itnode->value});
                }
            }
        }

        for(int i = 0; i < 3*xy_coord.size(); i++){
            fnode_list.push_back(i);
        }
        for(std::vector<std::pair<int,double>>::reverse_iterator itsub = sub_list.rbegin(); itsub < sub_list.rend(); ++itsub){
            fnode_list.erase(fnode_list.begin() + itsub->first);
        }
        printf("start sub Matrix \n");
        this->sub_mat_vec();
        printf("done sub Matrix \n");

        // for(int i = 0; i < K_mat.outerSize(); i++){
        //     printf("col: %d ", i);
        //     for(Eigen::SparseMatrix<double>::InnerIterator it(K_mat,i); it; ++it){
        //         printf("%.4f ",it.value());
        //     }
        //     printf("\n");
        // }
    }
    else{
        printf("mode %d is unsupported",mode);
    }
}

void asb_opt_manager::opt_val(int method){
    int max_iteration = 50;
    double start = clock();

    if(method == 1){ // Gauss-Newton
        this->initialize(1);
        int dof = probe_list.size();
        int optvar_num = Mater_lib.size();
        int* ipiv = nullptr;
        double* FE0it = nullptr;
        double* fdE0 = nullptr;
        double* JTJ = nullptr;
        //allocate memory
        try{
            ipiv = new int[optvar_num] {0};
            FE0it = new double[optvar_num] {0};
            fdE0 = new double[dof] {0};
            JTJ = new double[optvar_num * optvar_num] {0};
        }
        catch(const std::bad_alloc& e){
            std::cerr << e.what() << '\n';
            delete[] ipiv; ipiv = nullptr;
            delete[] FE0it; FE0it = nullptr;
            delete[] fdE0; fdE0 = nullptr;
            delete[] JTJ; JTJ = nullptr;
        }
        if(ipiv == nullptr || FE0it == nullptr || fdE0 == nullptr || JTJ == nullptr){return;}

        bool is_done = false;
        bool alloc_err = false;
        double normf = 0.0;
        double val = 0.0;
        int iter_time = 0;
        while(!is_done && iter_time < max_iteration){
            normf = 0.0;
            this->solve(alloc_err, true);
            this->sub_measure();

            for(int i = 0; i < dof; i++){
                fdE0[i] = uvw_ans[i] - u_real_probe.at(i);
                normf += fdE0[i] * fdE0[i];
            }
            if(normf < 1e-12){
                is_done = true;
                for(std::vector<Material>::iterator itE = Mater_lib.begin(); itE != Mater_lib.end(); itE++){
                    printf("%.2f ",itE->E);
                }
                printf("\n");
            }
            else{
                // calculated {J^T J}
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, optvar_num, optvar_num, dof, 1, udE, dof,
                            udE, dof, 0, JTJ, optvar_num);
                // calculated {J^T fdE0}
                cblas_dgemv(CblasRowMajor, CblasNoTrans, optvar_num, dof, 1, udE, dof, fdE0, 1, 0, FE0it, 1);
                // calculated E_{n+1}
                LAPACKE_dgetrf(LAPACK_ROW_MAJOR, optvar_num, optvar_num, JTJ, optvar_num, ipiv);
                LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', optvar_num, 1, JTJ, optvar_num, ipiv, FE0it, 1);

                printf("%.8f \n",normf);
                for(std::vector<Material>::iterator itE = Mater_lib.begin(); itE != Mater_lib.end(); itE++){
                    val = itE->E - FE0it[itE->where];
                    itE->E = (val > 100.0) ? val : 100;
                    printf("%.2f ",itE->E);
                }
                printf("\n");

                this->init_KE();
                // this->sub_mat_vec(false);
                this->addboundry(1);
            }
            iter_time ++;
            delete[] uvw_ans; uvw_ans = nullptr;
            delete[] udE; udE = nullptr;
        }

        delete[] fdE0; fdE0 = nullptr;
        delete[] JTJ; JTJ = nullptr;
        delete[] ipiv; ipiv = nullptr;
        delete[] FE0it; FE0it = nullptr;
    }
    else if(method == 2){ // Lagrange Multiplier Method
        // this->initialize(1);
        // int fsize = fnode_list.size();
        // int dof = (fsize != 0) ? fsize : 3*xy_coord.size();
        // int optvar_num = Mater_lib.size();

        // bool is_done = false;
        // int iter_t = 0;
        // int posi = -1;
        // double normf = 0.0;
        // double val = 0.0;
        // std::vector<int> rprobe_list;
    
        // rprobe_list.resize(dof, -1);
        // {
        //     int i = 0;
        //     for(std::vector<int>::iterator itpb = probe_list.begin(); itpb != probe_list.end(); itpb++){
        //         rprobe_list.at(*itpb) = i;
        //         i++;
        //     }
        // }

        // double* f_eqn = nullptr, *allvar = nullptr;
        // SparseMatrix K_eqn; 
        // try{
        //     fsize = 2 * dof + optvar_num;
        //     f_eqn = new double[fsize]{0};
        //     allvar = new double[fsize]{0};
        // }
        // catch(std::bad_alloc &e){
        //     std::cerr << e.what() << std::endl;
        //     delete[] f_eqn; f_eqn = nullptr;
        //     delete[] allvar; allvar = nullptr;
        // }
        // if(f_eqn == nullptr || allvar == nullptr){exit(1);}

        // for(int i = 0; i < dof; i++){
        //     posi = (fsize != 0) ? rprobe_list.at(fnode_list.at(i)) : rprobe_list.at(i);
        //     allvar.coeffRef(i) = (posi != -1) ? u_real_probe.at(posi) : 0.00;
        // }
        // for(int i = dof; i < 2 * dof; i++){
        //     allvar.coeffRef(i) = 0.0;
        // }
        // for(int i = 2*dof, j = 0; i < 2*dof + optvar_num; i++, j++){
        //     allvar.coeffRef(i) = Mater_lib.at(j).E;
        // }
        // // for(std::vector<int>::iterator itsub = sub_list.begin(); itsub != sub_list.end(); itsub++){
        // //     allvar.coeffRef(dof + *itsub) = 0.;
        // //     allvar.coeffRef(*itsub) = u_real_probe(*itsub);
        // // }
        
        // while(!is_done && iter_t < max_iteration){
        //     for(int i = 0; i < dof; i++){
        //         val = 0.0;
        //         for(int j = dof, rj = 0; j < 2*dof; j++, rj++){
        //             val += allvar.coeff(j) * K_mat.coeff(rj,i);
        //         }
        //         posi = (fsize != 0) ? rprobe_list.at(fnode_list.at(i)) : rprobe_list.at(i);
        //         f_eqn.coeffRef(i) = (posi != -1) ? \
        //             (2 * (allvar.coeff(i) - u_real_probe.at(posi)) + val) : val;
        //         printf("\rget f: %.2f%%", i * 100.0 / (2*dof + optvar_num));
        //     }
        //     for(int i = dof, ri = 0; i < 2 * dof; i++, ri++){
        //         val = 0.0;
        //         for(int j = 0; j < dof; j++){
        //             val += allvar.coeff(j) * K_mat.coeff(ri,j);
        //         }
        //         f_eqn.coeffRef(i) = val - Fout.coeff(ri);
        //         printf("\rget f: %.2f%%", i * 100.0 / (2*dof + optvar_num));
        //     }
        //     for(int i = 2*dof, ri = 0; i < 2*dof + optvar_num; i++, ri++){
        //         val = 0.0;
        //         for(int j = dof, rj = 0; j < 2 * dof; j++, rj++){
        //             for(int k = 0; k < dof; k++){
        //                 val += KdE_mat.at(ri).coeff(rj, k) * allvar.coeff(j) * allvar.coeff(k);
        //             }
        //         }
        //         f_eqn.coeffRef(i) = val;
        //         printf("\rget f: %.2f%%", i * 100.0 / (2*dof + optvar_num));
        //     }
        //     printf("\n");
            
        //     normf = 0.0;
        //     for(Eigen::VectorXd::iterator itf = f_eqn.begin(); itf != f_eqn.end(); itf++){
        //         normf += (*itf) * (*itf);
        //     }

        //     if(normf < 1E-8){
        //         is_done = true;
        //         for(int i = 2*dof; i < 2*dof + optvar_num; i++){
        //             printf("%.2f ", allvar.coeff(i));
        //         }
        //         printf("\n");
        //     }
        //     else{
        //         printf("%.6f \n",normf);
        //         K_eqn.setZero();
        //         for(int i = 0; i < dof; i++){ // block11 R:1-dof C:1-dof
        //             posi = (fsize != 0) ? rprobe_list.at(fnode_list.at(i)) : rprobe_list.at(i);
        //             if(posi != -1){
        //                 K_eqn.coeffRef(i,i) = 2.0;
        //             }

        //             for(int j = dof, rj = 0; j < 2 * dof; j++, rj++){ // block12,21 R:1-dof C:(dof+1)-2dof
        //                 K_eqn.coeffRef(i,j) = K_mat.coeff(rj, i);
        //                 K_eqn.coeffRef(j,i) = K_mat.coeff(rj, i);
        //             }

        //             for(int j = 2*dof, rj = 0; j < 2 * dof + optvar_num; j++, rj++){
        //                 val = 0.0;
        //                 for(int k = 0, rk = dof; k < dof; k++, rk++){
        //                     val += allvar.coeff(rk) * KdE_mat.at(rj).coeff(k,i);
        //                 }
        //                 K_eqn.coeffRef(i,j) = val;
        //                 K_eqn.coeffRef(j,i) = val;
        //             }
        //             printf("\rget K: %.2f%%", i * 100.0 / (2*dof));
        //         }
        //         for(int i = dof, ri = 0; i < 2 * dof; i++, ri++){ 
        //             for(int j = 2 * dof, rj = 0; j < 2 * dof + optvar_num; j++, rj++){
        //                 val = 0.0;
        //                 for(int k = 0; k < dof; k++){
        //                     val += KdE_mat.at(rj).coeff(ri,k) * allvar.coeff(k);
        //                 }
        //                 K_eqn.coeffRef(i,j) = val;
        //                 K_eqn.coeffRef(j,i) = val;
        //             }             
        //             printf("\rget K: %.2f%%", i * 100.0 / (2*dof));       
        //         }
        //         printf("\n");

        //         printf("itertion solve, times: %d",iter_t+1);
        //         Eigen::SparseLU<Eigen::SparseMatrix<double>> K_solver;
        //         K_solver.analyzePattern(K_eqn);
        //         K_solver.factorize(K_eqn);
        //         allvar -= K_solver.solve(f_eqn);

        //         for(int i = 2*dof, ri = 0; i < 2*dof + optvar_num; i++, ri++){
        //             if(allvar.coeff(i) < 0){
        //                 allvar.coeffRef(i) = 100;
        //             }
        //             printf("%.2f ", allvar.coeff(i));
        //             Mater_lib.at(ri).E = allvar.coeff(i);
        //         }
        //         printf("\n");
        //         this->init_KE();
        //         // this->sub_mat_vec(false);
        //         this->addboundry(1);
        //     }
        //     iter_t ++;
        // }
    
        printf("OPTIMIZE ERROR: Lagrange method (= 2) doesn't available presently\n\
                try Gauss-Newton method (= 1) again!");
        exit(2);
    }
    else{
        printf("method %d do not supported!", method);
    }
    double end = clock();
    printf("time: %.4f \n", end - start);
}

void asb_opt_manager::wirte(std::ofstream &file_stream){
    int fsize = fnode_list.size();
    int i = 0;

    printf("start writing to file ... \n");
    file_stream.setf(file_stream.scientific);
    file_stream.precision(10);
    file_stream << "*measure\n";
    if(fsize != 0){
        for(i = 0; i < fsize; i++){
            // snprintf(buffer, sizeof(buffer),"%d, %d, %.10f\n",fnode_list.at(i)/3 + 1,fnode_list.at(i) % 3 + 1,*itu);
            file_stream << fnode_list.at(i)/3 + 1 <<", " << fnode_list.at(i) % 3 + 1 << ", " << uvw_ans[i] << "\n";
        }
    }
    else{
        for(i = 0; i < 3*xy_coord.size(); i++){
            // snprintf(buffer, sizeof(buffer),"%d, %d, %.10f\n",i/3 + 1,i % 3 + 1,*itu);
            file_stream << i/3 + 1 <<", " << i % 3 + 1 <<", " << uvw_ans[i] << "\n";
        }
    }
    printf("done wirte \n");
}

void asb_opt_manager::wirte(std::ofstream &file_stream, const std::vector<int> &node_list){
    int fsize = fnode_list.size();
    int i = 0;
    std::vector<int> rfnode_list;

    printf("start writing to file ... \n");
    file_stream.setf(file_stream.scientific);
    file_stream.precision(10);
    file_stream << "*measure\n";
    if(fsize != 0){
        rfnode_list.resize(3 * xy_coord.size(), -1);
        i = 0;
        for(std::vector<int>::iterator itp = fnode_list.begin(); itp != fnode_list.end(); itp++){
            rfnode_list.at(*itp) = i;
            i++;
        }

        for(std::vector<int>::const_iterator itu = node_list.begin(); itu != node_list.end(); itu++){
            // snprintf(buffer, sizeof(buffer),"%d, %d, %.10f\n",fnode_list.at(i)/3 + 1,fnode_list.at(i) % 3 + 1,*itu);
            for(int j = 0; j < 3; j++){
                file_stream << *itu <<", " << j+1 << ", " << uvw_ans[rfnode_list.at(3*(*itu-1)+j)] << "\n";
            }
        }
    }
    else{
        for(std::vector<int>::const_iterator itu = node_list.begin(); itu != node_list.end(); itu++){
            // snprintf(buffer, sizeof(buffer),"%d, %d, %.10f\n",fnode_list.at(i)/3 + 1,fnode_list.at(i) % 3 + 1,*itu);
            for(int j = 0; j < 3; j++){
                file_stream << *itu <<", " << j+1 << ", " << uvw_ans[3*(*itu-1)+j] << "\n";
            }
        }
    }
    printf("done wirte \n");
}

void asb_opt_manager::init_KE_symbolic(){
    const int dof = xy_coord.size();
    const int opt_var = Mater_lib.size();
    int nznv = 0;
    int dignzn = 0;
    int nznv_E = 0;
    int dignzn_E = 0;
    int ROW_I, I3;
    bool* syb_K_mat = new bool[dof * dof];
    bool* syb_KdE_mat = new bool[dof * dof];
    
    if(K_mat.type == 2 || K_mat.type == -2){
        for(int i = 0; i < dof; i++){
            ROW_I = i * dof;
            for(int j = i; j < dof; j++){
                syb_K_mat[ROW_I + j] = false;
                syb_KdE_mat[ROW_I + j] = false;
            }
        }

        for(std::vector<c3d4>::iterator itele = elements.begin(); itele != elements.end(); itele++){
            if(itele->Nodetag.size() == 0){
                continue;
            }
            for(int i = 0; i < 4; i++){
                ROW_I = (itele->Nodetag[i] - 1) * dof;
                for(int j = i; j < 4; j++){
                    if(itele->Nodetag[i] > itele->Nodetag[j]){
                        syb_K_mat[(itele->Nodetag[j]-1)*dof + itele->Nodetag[i]-1] = true;
                    }
                    else{
                        syb_K_mat[ROW_I + itele->Nodetag[j]-1] = true;
                    }
                }
            }
        }

        for(int i = 0; i < dof; i++){
            ROW_I = i * dof;
            for(int j = i; j < dof; j++){
                if(syb_K_mat[ROW_I + j]){
                    nznv++;
                    if(i == j){
                        dignzn++;
                    }
                }
            }
        }
        int icol = 0;

        const int nzn = 9*nznv - 3*dignzn;
        K_mat.val = new double[nzn];
        K_mat.col = new int[nzn];
        K_mat.row_st = new int[3*dof+1];
        K_mat.row_st[0] = 1;

        bool is_dig;
        for(int i = 0; i < dof; i++){
            I3 = 3 * i;
            ROW_I = i * dof;
            is_dig = false;
            
            for(int j = i; j < dof; j++){
                if(syb_K_mat[ROW_I + j]){

                    K_mat.col[icol] = 3*j + 1;
                    K_mat.col[icol+1] = 3*j + 2;
                    K_mat.col[icol+2] = 3*j + 3;
                    K_mat.val[icol] = 0.;
                    K_mat.val[icol+1] = 0.;
                    K_mat.val[icol+2] = 0.;

                    if(i == j){
                        is_dig = true;
                    }
                    icol += 3;
                }
            }
            K_mat.row_st[I3 + 1] = icol + 1;

            if(is_dig){
                for(int j = 1; j < 3; j++){
                    for(int k = K_mat.row_st[I3] + j - 1; k < K_mat.row_st[I3 + 1] - 1; k++, icol++){
                        K_mat.col[icol] = K_mat.col[k];
                        K_mat.val[icol] = 0.;
                    }
                    K_mat.row_st[I3 + j + 1] = icol + 1;
                }
            }
            else{
                for(int j = 1; j < 3; j++){
                    for(int k = K_mat.row_st[I3] - 1; k < K_mat.row_st[I3 + 1] - 1; k++, icol++){
                        K_mat.col[icol] = K_mat.col[k];
                        K_mat.val[icol] = 0.;
                    }
                    K_mat.row_st[I3 + j + 1] = icol + 1;
                }
            }
        }

        // initial the KdE symbolic
        KdE_mat = new SparseMatrix[opt_var];
        for(int var = 0; var < opt_var; var++){
            for (int i = 0; i < dof; i++) {
                ROW_I = i * dof;
                for (int j = i; j < dof; j++) {
                    syb_KdE_mat[ROW_I + j] = false;
                }
            }

            for(std::vector<c3d4>::iterator itele = elements.begin(); itele != elements.end(); itele++){
                if(itele->Nodetag.size() == 0 || itele->mat->where != var){
                    continue;
                }
                for(int i = 0; i < 4; i++){
                    ROW_I = (itele->Nodetag[i] - 1) * dof;
                    for(int j = i; j < 4; j++){
                        if(itele->Nodetag[i] > itele->Nodetag[j]){
                            syb_KdE_mat[(itele->Nodetag[j]-1)*dof + itele->Nodetag[i]-1] = true;
                        }
                        else{
                            syb_KdE_mat[ROW_I + itele->Nodetag[j]-1] = true;
                        }
                    }
                }
            }

            nznv_E = 0;
            dignzn_E = 0;
            for(int i = 0; i < dof; i++){
                ROW_I = i * dof;
                for(int j = i; j < dof; j++){
                    if(syb_KdE_mat[ROW_I + j]){
                        nznv_E++;
                        if(i == j){
                            dignzn_E++;
                        }
                    }
                }
            }
            int icol_E = 0;

            const int nzn_E = 9*nznv_E - 3*dignzn_E;
            try{
                KdE_mat[var].val = new double[nzn_E];
                KdE_mat[var].col = new int[nzn_E];
                KdE_mat[var].row_st = new int[3 * dof + 1];
            }
            catch (const std::bad_alloc& e){
                std::cerr << e.what() << std::endl;
                exit(1);
            }
            KdE_mat[var].row_st[0] = 1;

            bool is_dig_E;
            for(int i = 0; i < dof; i++){
                I3 = 3 * i;
                ROW_I = i * dof;
                is_dig_E = false;
                
                for(int j = i; j < dof; j++){
                    if(syb_KdE_mat[ROW_I + j]){

                        KdE_mat[var].col[icol_E] = 3*j + 1;
                        KdE_mat[var].col[icol_E+1] = 3*j + 2;
                        KdE_mat[var].col[icol_E+2] = 3*j + 3;
                        KdE_mat[var].val[icol_E] = 0.;
                        KdE_mat[var].val[icol_E+1] = 0.;
                        KdE_mat[var].val[icol_E+2] = 0.;

                        if(i == j){
                            is_dig_E = true;
                        }
                        icol_E += 3;
                    }
                }
                KdE_mat[var].row_st[I3 + 1] = icol_E + 1;

                if(is_dig_E){
                    for(int j = 1; j < 3; j++){
                        for(int k = KdE_mat[var].row_st[I3] + j - 1; k < KdE_mat[var].row_st[I3 + 1] - 1; k++, icol_E++){
                            KdE_mat[var].col[icol_E] = KdE_mat[var].col[k];
                            KdE_mat[var].val[icol_E] = 0.;
                        }
                        KdE_mat[var].row_st[I3 + j + 1] = icol_E + 1;
                    }
                }
                else{
                    for(int j = 1; j < 3; j++){
                        for(int k = KdE_mat[var].row_st[I3] - 1; k < KdE_mat[var].row_st[I3 + 1] - 1; k++, icol_E++){
                            KdE_mat[var].col[icol_E] = KdE_mat[var].col[k];
                            KdE_mat[var].val[icol_E] = 0.;
                        }
                        KdE_mat[var].row_st[I3 + j + 1] = icol_E + 1;
                    }
                }
            }
        }
    }
    else{
        printf("matrix type %i haven't supported yet",K_mat.type);
    }
    delete[] syb_K_mat;
    delete[] syb_KdE_mat;
}