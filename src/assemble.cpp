#include "assemble.h"

void asb_manager::init_ele(){
    double H = 0.004, B = 1;
    int y_numele = 1, x_numele = 800;
    Material mat1;
    Mater_lib.push_back(mat1);

    double dH = H / y_numele, dB = B / x_numele;
    int y_numnode = 2*y_numele + 1, x_numnode = 2*x_numele + 1,
        i , j, ROW, ROW1, ROW2, ROW3, ROW4;
    double x_i, y_i;

    xy_coord = new double[2 * y_numnode * x_numnode]{0};
    typenode_num[0] = (y_numele + 1) * (x_numele + 1);
    typenode_num[1] = y_numnode * x_numnode - typenode_num[0];
    typenode_num[2] = 0;
    int dof = 4*typenode_num[0] + 3*typenode_num[1] + typenode_num[2];
    uvtp_ans = new double[dof]{0.};

    for(i = 0; i < typenode_num[0]; i++){
        uvtp_ans[4*i + 2] = 25;
    }
    for(i = typenode_num[0]; i < typenode_num[0] + typenode_num[1]; i++){
        uvtp_ans[3*i + typenode_num[0] + 2] = 25;
    }

    inlet_node.resize(y_numnode);
    outlet_node.resize(y_numnode);
    for(i = 0; i < x_numele + 1; i++){
        ROW = 2 * i * (y_numele + 1);
        x_i = dB * i;
        for(j = 0; j < y_numele + 1; j++){
            xy_coord[ROW + 2*j] = x_i;
            xy_coord[ROW + 2*j + 1] = dH * j;
        }
        if(i == 0){
            for(j = 0; j < y_numele + 1; j++){
                inlet_node.at(2*j) = j + 1;
            }
        }
        else if(i == x_numele){
            int itmp = x_numele*(y_numele+1);
            for(j = 0; j < y_numele + 1; j++){
                outlet_node.at(2*j) = itmp + j + 1;
            }
        }
    }

    double* xy_nstart = xy_coord + 2*typenode_num[0];
    for(i = 0; i < x_numele + 1; i++){
        ROW = 2 * i * y_numele;
        x_i = dB * i;
        for(j = 0; j < y_numele; j++){
            xy_nstart[ROW + 2*j] = x_i;
            xy_nstart[ROW + 2*j + 1] = dH * j + dH / 2;
        }

        if(i == 0){
            for(j = 0; j < y_numele; j++){
                inlet_node.at(2*j+1) = typenode_num[0] + j + 1;
            }
        }
        else if(i == x_numele){
            int itmp = typenode_num[0] + x_numele * y_numele;
            for(j = 0; j < y_numele; j++){
                outlet_node.at(2*j+1) = itmp + j + 1;
            }
        }
    }

    xy_nstart += 2 * (x_numele + 1) * y_numele;
    for(i = 0; i < x_numele; i++){
        ROW = 2 * i * (2*y_numele + 1);
        x_i = dB * i + dB / 2;

        for(j = 0; j < 2*y_numele + 1; j++){
            xy_nstart[ROW + 2*j] = x_i;
            xy_nstart[ROW + 2*j + 1] = dH * j / 2;
        }
    }
    xy_nstart = nullptr;

    elements.resize(x_numele * y_numele);
    P9SF m_ele(&Mater_lib.at(0));
    int nt[9]{0}; 
    for(i = 0; i < x_numele; i++){
        ROW = i * y_numele;
        ROW1 = i * (y_numele + 1);
        ROW2 = (i+1) * (y_numele+1);
        ROW3 = (i+1) * y_numele;
        ROW4 = i * (2*y_numele + 1) + (x_numele + 1) * (2*y_numele + 1);

        for (j = 0; j < y_numele; j++) {
            nt[0] = ROW1 + j + 1;
            nt[1] = ROW2 + j + 1;
            nt[2] = ROW2 + j + 2;
            nt[3] = ROW1 + j + 2;

            nt[4] = ROW4 + 2 * j + 1;
            nt[5] = ROW3 + typenode_num[0] + j + 1;
            nt[6] = ROW4 + 2 * j + 3;
            nt[7] = ROW + typenode_num[0] + j + 1;

            nt[8] = ROW4 + 2*j + 2;

            m_ele.set_nval(nt, xy_coord, uvtp_ans, typenode_num);
            elements.at(ROW + j) = m_ele;
        }
    }

    wall_node.resize(2*x_numnode);
    for(i = 0; i < x_numele+1; i++){
        wall_node[2*i] = i * (y_numele+1) + 1;
        wall_node[2*i + 1] = (i+1) * (y_numele+1);
    }
    for(i = 0; i < x_numele; i++){
        wall_node[2*(x_numele+1)+2*i] = (x_numele+1)*y_numnode + i * y_numnode + 1;
        wall_node[2*(x_numele+1)+2*i+1] = (x_numele+1)*y_numnode + (i+1) * y_numnode;
    }

};

// void asb_manager::init_bnd(){
//     for(std::vector<Boundary<std::string,double,int>>::iterator itbm = bound_set_map.begin(); itbm < bound_set_map.end(); itbm++){
//         Boundary<set*,double,int> bdele;
//         bdele.dims = itbm->dims; bdele.is_alldim = itbm->is_alldim;
//         bdele.is_zero = itbm->is_zero; bdele.value = itbm->value;
//         bdele.nodes = whereset(itbm->nodes,"N");
//         bound_set.push_back(bdele);
//     }
//     bound_set_map.clear();

//     for(std::vector<Boundary<std::string,double,int>>::iterator itlm = load_set_map.begin(); itlm < load_set_map.end(); itlm++){
//         Boundary<set*,double,int> ldele;
//         ldele.dims = itlm->dims; ldele.is_alldim = itlm->is_alldim;
//         ldele.is_zero = itlm->is_zero; ldele.value = itlm->value;
//         ldele.nodes = whereset(itlm->nodes,"N");
//         load_set.push_back(ldele);
//     }
//     load_set_map.clear();

//     for(std::vector<Boundary<std::string,double,double>>::iterator itdlm = dload_set_map.begin(); itdlm != dload_set_map.end(); itdlm++){
//         Boundary<set*,double,double> dldele;
//         dldele.dims = itdlm->dims; dldele.value = itdlm->value;
//         if(itdlm->nodes.size() == 0){
//             for(std::vector<set>::iterator itele = ELset.begin(); itele != ELset.end(); itele++){
//                 dldele.nodes = whereset(itele->name,"EL");
//                 dload_set.push_back(dldele);
//             }
//         }
//         else{
//             dldele.nodes = whereset(itdlm->nodes, "EL");
//             dload_set.push_back(dldele);
//         }
//     }
//     dload_set_map.clear();
// }

void asb_manager::init_KF(){
    const int dof = 4*typenode_num[0] + 3*typenode_num[1] + typenode_num[2];
    int ROW_I;

    if(K_mat.val == nullptr && Fint == nullptr){
        K_mat.type = 11;
        init_K_symbolic();
        Fint = new double[dof]{0.};
    }
    else{
        for(int i = 0; i < K_mat.row_st[dof]-1; i++){
            K_mat.val[i] = 0.;
        }
        for(int i = 0; i < dof; i++){
            Fint[i] = 0.;
        }
    }

    for(std::vector<P9SF>::iterator itele = elements.begin(); itele != elements.end(); itele++){
        if(!itele->is_valid()){
            continue;
        }
        if(itele->is_f()){
            asb_KF_f(*itele);
        }
        else{
            asb_KF_s(*itele);
        }
    }
};

void asb_manager::initialize(int mode){
    // std::ifstream inp_file(inp_file_name);
    // read_manager(inp_file);
    init_ele();
    // init_bnd();
    init_KF();  

    getFout();
    // for(Eigen::VectorXd::iterator it = Fout.begin(); it != Fout.end(); it++){
    //     printf("%.4f ", *it);
    // }
    // printf("\n");
    addboundry();

    //int dof = 4*typenode_num[0] + 3*typenode_num[1] + typenode_num[2];
    //K_mat.show(dof);
     //for (int i = 0; i < 31; i++) {
     //    printf("{");
     //    for (int j = 0; j < 30; j++) {
     //        printf("%.8f, ", K_mat.val[i*31+j]);
     //    }
     //    printf("%.8f},", K_mat.val[i*31+30]);
     //}
};

void asb_manager::asb_KF_f(P9SF &ele, int* invcol){
    const int dof = 4*typenode_num[0] + 3*typenode_num[1] + typenode_num[2];
    int dof1 = typenode_num[0];
    int ROW_I, I4, J4, J3 ,I3;
    double K[31][31];
    double F[31];
    ele.KF_f(K, F);

    if(K_mat.type == 1 || K_mat.type == 11){
        int I, J;
        for (int i = 0; i < 4; i++){
            I4 = 4 * i;
            I = ele.at_Nodetag(i);
            I = 4*I - 4;

            for (int j = 0; j < 4; j++){
                J = ele.at_Nodetag(j);
                J = 4*J - 4;
                J4 = 4 * j;
                
                for(int m = 0; m < 4; m++){
                    ROW_I = I + m + 1;
                    for(int n = 0; n < 4; n++){
                        K_mat.val[match(ROW_I, J+n+1, K_mat)] += K[I4 + m][J4 + n];
                    }
                }
            }

            for(int j = 4; j < 9; j++){
                J = ele.at_Nodetag(j);
                J = 3*J - 3 + dof1;
                J3 = j * 3 + 4;

                for(int m = 0; m < 4; m++){
                    ROW_I = I + m + 1;
                    for(int n = 0; n < 3; n++){
                        K_mat.val[match(ROW_I, J+n+1, K_mat)] += K[I4 + m][J3 + n];
                    }
                }
            }

            for(int m = 0; m < 4; m++){
                Fint[I + m] += F[I4 + m];
            }
        }

        for(int i = 4; i < 9; i++){
            I = ele.at_Nodetag(i);
            I = 3*I + dof1 - 3;
            I3 = 3 * i + 4;

            for(int j = 0; j < 4; j++){
                J = ele.at_Nodetag(j);
                J = 4*J - 4;
                J4 = 4 * j;

                for(int m = 0; m < 3; m++){
                    ROW_I = I + m + 1;
                    for(int n = 0; n < 4; n++){
                        K_mat.val[match(ROW_I, J+n+1, K_mat)] += K[I3 + m][J4 + n];
                    }
                }
            }

            for(int j = 4; j < 9; j++){
                J = ele.at_Nodetag(j);
                J = 3*J - 3 + dof1;
                J3 = j * 3 + 4;

                for(int m = 0; m < 3; m++){
                    ROW_I = I + m + 1;
                    for(int n = 0; n < 3; n++){
                        K_mat.val[match(ROW_I, J+n+1, K_mat)] += K[I3 + m][J3 + n];
                    }
                }
            }

            for(int m = 0; m < 3; m++){
                Fint[I + m] += F[I3 + m];
            }
        }
    }
    
};

void asb_manager::asb_KF_s(P9SF &ele, int* invcol){
    double K[9][9];
    double F[9];
    ele.KF_s(K, F);

    if(K_mat.type == 1){
        int I, J;
        for (int i = 0; i < 9; i++){
            I = ele.at_Nodetag(i);
            I = wherend(I, 2);

            for (int j = 0; j < 9; j++){
                J = ele.at_Nodetag(j);
                J = wherend(J, 2);
                
                K_mat.val[match(I+1, J+1, K_mat)] += K[i][j];
            }

            Fint[I] += F[i];
        }
    }
    
};

void asb_manager::getFout(){
    int dof = 4*typenode_num[0] + 3*typenode_num[1] + typenode_num[2];

    int where = wherend(inlet_node.at(inlet_node.size()-2), 0);

    Fout = new double[dof]{0};
    Fout[where] = 2e4;

};

// set* asb_manager::whereset(const std::string &name, const char* mode){
//     if(!_strcmpi(mode,"N")){
//         for(std::vector<set>::iterator itn = Nset.begin(); itn < Nset.end(); itn ++){
//             if(!strcmp(itn->name.data(),name.data())){
//                 return &*itn;
//             }
//         }
//     }
//     else if(!_strcmpi(mode,"EL")){
//         for(std::vector<set>::iterator itel = ELset.begin(); itel < ELset.end(); itel ++){
//             if(!strcmp(itel->name.data(),name.data())){
//                 return &*itel;
//             }
//         }
//     }
//     else{
//         printf("option: %s is not supported",mode);
//     }
//     return nullptr;
// };

// Material* asb_manager::whereset(const std::string &name){
//     for(std::vector<Material>::iterator itma = Mater_lib.begin(); itma < Mater_lib.end(); itma ++){
//         if(!strcmp(itma->name.data(),name.data())){
//             return &*itma;
//         }
//     }
//     return nullptr;
// }

int asb_manager::wherend(int nodetag, int uvtp){
    if(nodetag <= typenode_num[0]){
        return 4 * nodetag - 4 + uvtp;
    }
    else if(nodetag > typenode_num[0] && nodetag <= typenode_num[0] + typenode_num[1]){
        return 3 * nodetag + typenode_num[0] - 3 + uvtp;
    }
    else{
        if(uvtp != 2){
            std::cerr << "ERROR: this node is third type node!" << std::endl;
        }
        return typenode_num[0] * 3 + typenode_num[1] * 2 + nodetag - 1;
    }
};

void asb_manager::solve(){
    this->initialize();
    const int dof = 4*typenode_num[0] + 3*typenode_num[1] + typenode_num[2];
    const int max_iteration = 50;
    int mtype = K_mat.type;
    int nrhs = 1, itr_time = 0;
    int phase;
    double* uvtp_itr = new double[dof]{0};
    void* ptsolver[64];
    pardiso_cfg para = pardiso_cfg(mtype);
    double ddum, norm;

    for(int i = 0; i < 64; i++){
        ptsolver[i] = 0; 
    }

    do{
        norm = 0.;
        for(int i = 0; i < dof; i++){
            Fint[i] = Fout[i] - Fint[i];
            norm += Fint[i] * Fint[i];
        }

        if(norm < 1e-10){
            printf("done solving !\n");
            for (int i = 0; i < dof; i++) {
                printf("%.4f ", uvtp_ans[i]);
            }
            printf("\n");
            break;
        }
        else{
            printf("itertion time: %i, Residual: %.8f\n", itr_time+1, norm);

            phase = 11;
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
                    K_mat.row_st, K_mat.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), Fint, uvtp_itr, &(para.error));
            if(para.error != 0 ){
                printf ("ERROR during solution: %i \n", para.error);
            }
            printf("iteration solve completed ...\n");

            for(int i = 0; i < dof; i++){
                uvtp_ans[i] += uvtp_itr[i];
            }

            phase = -11;
            pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(K_mat.type), &phase, &dof, K_mat.val, 
                    K_mat.row_st, K_mat.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), Fint, &ddum, &(para.error));

            this->init_KF();
            this->addboundry();
        }

        itr_time++;
    }while(itr_time < max_iteration);

    delete[] uvtp_itr; uvtp_itr = nullptr;

};

void asb_manager::addboundry(){
    int i = 0, j = 0, k = 0,
        ROW_I, ROW_Ip1, COL,
        Tag[3];

    double Gauss_point[3] = {-sqrt(0.6), 0, sqrt(0.6)};
    double Gauss_point2[2] = {-1/sqrt(3), 1/sqrt(3)};
    double weight[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
    double y[3], p[2], u[3];
    double N[3][3], Np[2][2];
    double Ndl[3][3], Npdl[2][2];
    double Jaco, p_temp, u_temp;

    for(i = 0; i < 3; i++){
        double val = Gauss_point[i];
        N[i][0] = val * (val-1) / 2;
        N[i][1] = 1 - val * val;
        N[i][2] = val * (val+1) / 2;

        Ndl[i][0] = val - 0.5;
        Ndl[i][1] = - 2 * val;
        Ndl[i][2] = val + 0.5;
    }

    for(i = 0; i < 2; i++){
        double val = Gauss_point2[i];
        Np[i][0] = (1-val) / 2;
        Np[i][1] = (1+val) / 2;

        Npdl[i][0] = -0.5;
        Npdl[i][1] = 0.5;
    }

    //for(i = 1; i < inlet_node.size()-2; i++){
    //    ROW_I = wherend(inlet_node.at(i), 0);
    //    ROW_Ip1 = wherend(inlet_node.at(i+1), 0);
    //    for(j = K_mat.row_st[ROW_I] - 1; j < K_mat.row_st[ROW_I + 1] - 1; j++){
    //        COL = match(ROW_Ip1+1, K_mat.col[j], K_mat);
    //        if(COL != -1){
    //            K_mat.val[j] -= K_mat.val[COL];
    //        }
    //    }

    //    Fint[ROW_I] -= Fint[ROW_Ip1];
    //}

    //ROW_I = wherend(inlet_node.at(inlet_node.size()-2), 0);
    //for(i = K_mat.row_st[ROW_I] - 1; i < K_mat.row_st[ROW_I + 1] - 1 ; i++){
    //    K_mat.val[i] = 0;
    //}
    //Fint[ROW_I] = 0;

    //for(j = 0; j < 2; j++){
    //    for(i = 0; i < (inlet_node.size()-1) / 2; i++){
    //        for(k = 0; k < 2; k++){
    //            Tag[k] = inlet_node.at(2*i+2*k);
    //            y[k] = xy_coord[(Tag[k]-1)*2 + 1];
    //            p[k] = uvtp_ans[wherend(Tag[k], 3)];
    //        }
    //        Jaco = abs(Npdl[j][0]*y[0] + Npdl[j][1]*y[1]);
    //        p_temp = Np[j][0]*p[0] + Np[j][1]*p[1];

    //        K_mat.val[match(ROW_I+1, wherend(Tag[0], 3)+1, K_mat)] += Jaco * Np[j][0] / (y[2] - y[0]);
    //        K_mat.val[match(ROW_I+1, wherend(Tag[1], 3)+1, K_mat)] += Jaco * Np[j][1] / (y[2] - y[0]);

    //        Fint[ROW_I] += p_temp / (y[2] - y[0]) * Jaco;
    //    }

    //    for(i = 0; i < (outlet_node.size()-1) / 2; i++){
    //        for(k = 0; k < 2; k++){
    //            Tag[k] = outlet_node.at(2*i+2*k);
    //            y[k] = xy_coord[(Tag[k]-1)*2 + 1];
    //            p[k] = uvtp_ans[wherend(Tag[k], 3)];
    //        }
    //        Jaco = abs(Npdl[j][0]*y[0] + Npdl[j][1]*y[1]);
    //        p_temp = Np[j][0]*p[0] + Np[j][1]*p[1];

    //        K_mat.val[match(ROW_I+1, wherend(Tag[0], 3)+1, K_mat)] -= Jaco * Np[j][0] / (y[2] - y[0]);
    //        K_mat.val[match(ROW_I+1, wherend(Tag[1], 3)+1, K_mat)] -= Jaco * Np[j][1] / (y[2] - y[0]);

    //        Fint[ROW_I] -= p_temp * Jaco / (y[2] - y[0]);
    //    }
    //}

    //for(j = 0; j < 3; j++){
    //    for(i = 0; i < (inlet_node.size()-1) / 2; i++){
    //        for(k = 0; k < 3; k++){
    //            Tag[k] = inlet_node.at(2*i+k);
    //            y[k] = xy_coord[(Tag[k]-1)*2 + 1];
    //            u[k] = uvtp_ans[wherend(Tag[k], 0)];
    //        }

    //        Jaco = abs(Ndl[j][0]*y[0] + Ndl[j][1]*y[1] + Ndl[j][2]*y[2]);
    //        if (i == 0 && i != (inlet_node.size() - 3) / 2) {
    //            u_temp = N[j][1] * u[1] + N[j][2] * u[2];
    //        }
    //        else if (i == (inlet_node.size() - 3) / 2 && i != 0) {
    //            u_temp = N[j][0] * u[0] + N[j][1] * u[1];
    //        }
    //        else if (i == 0 && i == (inlet_node.size() - 3) / 2) {
    //            u_temp = N[j][1]*u[1];
    //        }
    //        else {
    //            u_temp = N[j][0] * u[0] + N[j][1] * u[1] + N[j][2] * u[2];
    //        }

    //        if(i != 0){
    //            K_mat.val[match(ROW_I + 1, wherend(Tag[0], 0) + 1, K_mat)] += 8e7 * weight[j] * Jaco * N[j][0];
    //        }
    //        K_mat.val[match(ROW_I+1, wherend(Tag[1], 0)+1, K_mat)] += 8e7 * weight[j] * Jaco * N[j][1];
    //        if(i != (inlet_node.size()-3) / 2){
    //            K_mat.val[match(ROW_I+1, wherend(Tag[2], 0)+1, K_mat)] += 8e7 * weight[j] * Jaco * N[j][2];
    //        }
    //        

    //        Fint[ROW_I] += 8e7 * u_temp * Jaco * weight[j];
    //    }
    //}

    for (i = 1; i < inlet_node.size() - 2; i++) {
        ROW_I = wherend(inlet_node.at(i), 0);
        ROW_Ip1 = wherend(inlet_node.at(i + 1), 0);
        for (j = K_mat.row_st[ROW_I] - 1; j < K_mat.row_st[ROW_I + 1] - 1; j++) {
            COL = match(ROW_Ip1 + 1, K_mat.col[j], K_mat);
            if (COL != -1) {
                K_mat.val[j] -= K_mat.val[COL];
            }
        }

        Fint[ROW_I] -= Fint[ROW_Ip1];
    }

    ROW_I = wherend(inlet_node.at(inlet_node.size() - 2), 0);
    for (i = K_mat.row_st[ROW_I] - 1; i < K_mat.row_st[ROW_I + 1] - 1; i++) {
        K_mat.val[i] = 0;
    }
    Fint[ROW_I] = 0;

    p_temp = 0.;
    double pnode = (inlet_node.size() + 1) / 2;
    for (i = 0; i < inlet_node.size(); i += 2) {
        Tag[0] = inlet_node.at(i);
        p_temp += uvtp_ans[wherend(Tag[0], 3)];

        K_mat.val[match(ROW_I + 1, wherend(Tag[0], 3) + 1, K_mat)] += 1.0 / pnode;
    }
    Fint[ROW_I] += p_temp / pnode;

    p_temp = 0.;
    pnode = (outlet_node.size() + 1) / 2;
    for (i = 0; i < outlet_node.size(); i+=2) {
        Tag[0] = outlet_node.at(i);
        p_temp += uvtp_ans[wherend(Tag[0], 3)];

        K_mat.val[match(ROW_I + 1, wherend(Tag[0], 3) + 1, K_mat)] += -1.0 / pnode;
    }
    Fint[ROW_I] += - p_temp / pnode;

    u_temp = 0.;
    for (i = 1; i < inlet_node.size() - 1; i++) {
        Tag[0] = inlet_node.at(i - 1);
        Tag[1] = inlet_node.at(i);
        Tag[2] = inlet_node.at(i + 1);
        y[0] = xy_coord[(Tag[0] - 1) * 2 + 1];
        y[1] = xy_coord[(Tag[1] - 1) * 2 + 1];
        y[2] = xy_coord[(Tag[2] - 1) * 2 + 1];

        double dh = abs(y[2]-y[0]) / 2.0;

        u_temp += uvtp_ans[wherend(Tag[1], 0)] * dh;

        K_mat.val[match(ROW_I + 1, wherend(Tag[1], 0) + 1, K_mat)] += 8e7 * dh;
    }
    Fint[ROW_I] += u_temp * 8e7;

    for(std::vector<int>::iterator itwx = wall_node.begin(); itwx != wall_node.end(); itwx++){
        int where_x = wherend(*itwx, 0);
        int where_y = wherend(*itwx, 1);
        int where_t = wherend(*itwx, 2);

        K_mat.val[match(where_x+1, where_x+1, K_mat)] = 1E10;
        Fint[where_x] = 0;
        K_mat.val[match(where_y + 1, where_y + 1, K_mat)] = 1E10;
        Fint[where_y] = 0;
        K_mat.val[match(where_t + 1, where_t + 1, K_mat)] = 1E10;
        Fint[where_t] = 0;
    }

}

void asb_manager::init_K_symbolic(){
    const int dof = typenode_num[0] + typenode_num[1] + typenode_num[2];
    int nzn44 = 0, nzn43 = 0, nzn41 = 0,
                   nzn33 = 0, nzn31 = 0,
                              nzn11 = 0;
    int ROW_I, ROW_Is1, I3, I4, I;
    bool* syb_K_mat = new bool[dof * dof];
    
    if(K_mat.type == 1 || K_mat.type == 11){
        for(int i = 0; i < dof; i++){
            ROW_I = i * dof;
            for(int j = 0; j < dof; j++){
                syb_K_mat[ROW_I + j] = false;
            }
        }

        for(std::vector<P9SF>::iterator itele = elements.begin(); itele != elements.end(); itele++){
            if(!itele->is_valid()){
                continue;
            }
            for(int i = 0; i < 9; i++){
                ROW_I = (itele->at_Nodetag(i) - 1) * dof;
                for(int j = 0; j < 9; j++){
                    syb_K_mat[ROW_I + itele->at_Nodetag(j)-1] = true;
                }
            }
        }

        //for (int i = 0; i < dof; i++) {
        //    ROW_I = i * dof;
        //    printf("\nROW %i: ", i+1);
        //    for (int j = 0; j < dof; j++) {
        //        printf("%i ",syb_K_mat[ROW_I + j]);
        //    }
        //}

        for(int i = 2; i < inlet_node.size()-1; i++){
            ROW_I = (inlet_node[i] - 1) * dof;
            ROW_Is1 = (inlet_node[i - 1] - 1) * dof;
            for(int j = 0; j < dof; j++){
                if(syb_K_mat[ROW_I + j]){
                    syb_K_mat[ROW_Is1 + j] = true;
                }
            }
        }

        // ROW_I = (inlet_node.at(inlet_node.size()-2) - 1) * dof;
        ROW_I = (inlet_node.at(inlet_node.size()-2) - 1) * dof;
        for(int i = 0; i < inlet_node.size(); i++){
            syb_K_mat[ROW_I + inlet_node.at(i) - 1] = true;
        }
        for(int i = 0; i < outlet_node.size(); i += 2){
            syb_K_mat[ROW_I + outlet_node.at(i) - 1] = true;
        }

        for(int i = 0; i < typenode_num[0]; i++){
            ROW_I = i * dof;
            for(int j = 0; j < typenode_num[0]; j++){
                if(syb_K_mat[ROW_I + j]){
                     nzn44++;
                }
            }

            for(int j = typenode_num[0]; j < typenode_num[1] + typenode_num[0]; j++){
                if(syb_K_mat[ROW_I + j]){
                    nzn43 ++;
                }
            }

            for(int j = typenode_num[1] + typenode_num[0]; j < dof; j++){
                if(syb_K_mat[ROW_I + j]){
                    nzn41 ++;
                }
            }
        }

        for(int i = typenode_num[0]; i < typenode_num[1] + typenode_num[0]; i++){
            ROW_I = dof * i;
            
            for (int j = 0; j < typenode_num[0]; j++) {
                if (syb_K_mat[ROW_I + j]) {
                    nzn43++;
                }
            }

            for (int j = typenode_num[0]; j < typenode_num[1] + typenode_num[0]; j++) {
                if (syb_K_mat[ROW_I + j]) {
                    nzn33++;
                }
            }

            for (int j = typenode_num[1] + typenode_num[0]; j < dof; j++) {
                if (syb_K_mat[ROW_I + j]) {
                    nzn31++;
                }
            }
        }

        for(int i = typenode_num[0] + typenode_num[1]; i < dof; i++){
            ROW_I = i * dof;

            for (int j = 0; j < typenode_num[0]; j++) {
                if (syb_K_mat[ROW_I + j]) {
                    nzn41++;
                }
            }

            for (int j = typenode_num[0]; j < typenode_num[1] + typenode_num[0]; j++) {
                if (syb_K_mat[ROW_I + j]) {
                    nzn31++;
                }
            }

            for (int j = typenode_num[1] + typenode_num[0]; j < dof; j++) {
                if (syb_K_mat[ROW_I + j]) {
                    nzn11++;
                }
            }
        }

        const int nzn = 16 * nzn44 + 12 * nzn43 + 4 * nzn41 + 9 * nzn33 + 3 * nzn31 + nzn11;
        K_mat.val = new double[nzn]{0};
        K_mat.col = new int[nzn];
        K_mat.row_st = new int[4 * typenode_num[0] + 3 * typenode_num[1] + typenode_num[2] + 1] {0};

        int icol = 0;
        K_mat.row_st[0] = 1;
        for(int i = 0; i < typenode_num[0]; i++){
            I4 = 4 * i;
            ROW_I = i * dof;
            
            for(int j = 0; j < typenode_num[0]; j++){
                if(syb_K_mat[ROW_I + j]){
                    K_mat.col[icol] = 4*j + 1;
                    K_mat.col[icol+1] = 4*j + 2;
                    K_mat.col[icol+2] = 4*j + 3;
                    K_mat.col[icol+3] = 4*j + 4;
                    icol += 4;
                }
            }

            for(int j = typenode_num[0]; j < typenode_num[0] + typenode_num[1]; j++){
                if(syb_K_mat[ROW_I + j]){
                    K_mat.col[icol] = 3*j + typenode_num[0] + 1;
                    K_mat.col[icol+1] = 3*j + typenode_num[0] + 2;
                    K_mat.col[icol+2] = 3*j + typenode_num[0] + 3;
                    icol += 3;
                }
            }

            for(int j = typenode_num[0] + typenode_num[1]; j < dof; j++){
                if(syb_K_mat[ROW_I + j]){
                    K_mat.col[icol] = j + 3*typenode_num[0] + 2*typenode_num[1] + 1;
                    icol++;
                }
            }
            K_mat.row_st[I4 + 1] = icol + 1;

            for(int j = 1; j < 4; j++){
                for(int k = K_mat.row_st[I4] - 1; k < K_mat.row_st[I4 + 1] - 1; k++, icol++){
                    K_mat.col[icol] = K_mat.col[k];
                }
                K_mat.row_st[I4 + j + 1] = icol + 1;
            }
        }

        for(int i = typenode_num[0]; i < typenode_num[0] + typenode_num[1]; i++){
            I3 = 3 * i + typenode_num[0];
            ROW_I = i * dof;
            
            for(int j = 0; j < typenode_num[0]; j++){
                if(syb_K_mat[ROW_I + j]){
                    K_mat.col[icol] = 4*j + 1;
                    K_mat.col[icol+1] = 4*j + 2;
                    K_mat.col[icol+2] = 4*j + 3;
                    K_mat.col[icol+3] = 4*j + 4;
                    icol += 4;
                }
            }

            for(int j = typenode_num[0]; j < typenode_num[0] + typenode_num[1]; j++){
                if(syb_K_mat[ROW_I + j]){
                    K_mat.col[icol] = 3*j + typenode_num[0] + 1;
                    K_mat.col[icol+1] = 3*j + typenode_num[0] + 2;
                    K_mat.col[icol+2] = 3*j + typenode_num[0] + 3;
                    icol += 3;
                }
            }

            for(int j = typenode_num[0] + typenode_num[1]; j < dof; j++){
                if(syb_K_mat[ROW_I + j]){
                    K_mat.col[icol] = j + 3*typenode_num[0] + 2*typenode_num[1] + 1;
                    icol++;
                }
            }
            K_mat.row_st[I3 + 1] = icol + 1;

            for(int j = 1; j < 3; j++){
                for(int k = K_mat.row_st[I3] - 1; k < K_mat.row_st[I3 + 1] - 1; k++, icol++){
                    K_mat.col[icol] = K_mat.col[k];
                }
                K_mat.row_st[I3 + j + 1] = icol + 1;
            }
        }

        for(int i = typenode_num[0] + typenode_num[1]; i < dof; i++){
            I = i + 3*typenode_num[0] + 2*typenode_num[1];
            ROW_I = i * dof;
            
            for(int j = 0; j < typenode_num[0]; j++){
                if(syb_K_mat[ROW_I + j]){
                    K_mat.col[icol] = 4*j + 1;
                    K_mat.col[icol+1] = 4*j + 2;
                    K_mat.col[icol+2] = 4*j + 3;
                    K_mat.col[icol+3] = 4*j + 4;
                    icol += 4;
                }
            }

            for(int j = typenode_num[0]; j < typenode_num[0] + typenode_num[1]; j++){
                if(syb_K_mat[ROW_I + j]){
                    K_mat.col[icol] = 3*j + typenode_num[0] + 1;
                    K_mat.col[icol+1] = 3*j + typenode_num[0] + 2;
                    K_mat.col[icol+2] = 3*j + typenode_num[0] + 3;
                    icol += 3;
                }
            }

            for(int j = typenode_num[0] + typenode_num[1]; j < dof; j++){
                if(syb_K_mat[ROW_I + j]){
                    K_mat.col[icol] = j + 3*typenode_num[0] + 2*typenode_num[1] + 1;
                    icol++;
                }
            }
            K_mat.row_st[I + 1] = icol + 1;
        }
    }
    else{
        printf("matrix type %i haven't supported yet",K_mat.type);
    }

    //printf("\n");
    //for (int i = 0; i < K_mat.row_st[4 * typenode_num[0] + 3 * typenode_num[1]]-1; i++) {
    //    printf("%i ", K_mat.col[i]);
    //}
    //printf("\n");

    delete[] syb_K_mat; syb_K_mat = nullptr;

}