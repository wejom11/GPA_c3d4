#include "assemble.h"

void asb_manager::init_ele(){
    for(std::vector<c3d4>::iterator it = elements.begin(); it < elements.end(); it++){
        for (std::vector<int>::iterator its = it->Nodetag.begin(); its < it->Nodetag.end(); its++){
            it->xyz0.push_back(xyz_coord.at(*its - 1));
        }
        if(it->Nodetag.size() != 0){
            it->vol = it->getvolume();
        }
    }
    std::vector<std::pair<set*, Material*>> map_p;
    for(int i = 0; i < mat_map.size(); i++){
        map_p.push_back({whereset(mat_map.at(i).first,"EL"), whereset(mat_map.at(i).second)});
    }
    for(int i = 0; i < map_p.size(); i++){
        for(int j = 0; j < map_p.at(i).first->sets.size(); j++){
            elements[map_p.at(i).first->sets.at(j)-1].mat = map_p.at(i).second;
        }
    }
    mat_map.clear();
};

void asb_manager::init_bnd(){
    for(std::vector<Boundary<std::string,double,int>>::iterator itbm = bound_set_map.begin(); itbm < bound_set_map.end(); itbm++){
        Boundary<set*,double,int> bdele;
        bdele.dims = itbm->dims; bdele.is_alldim = itbm->is_alldim;
        bdele.is_zero = itbm->is_zero; bdele.value = itbm->value;
        bdele.nodes = whereset(itbm->nodes,"N");
        bound_set.push_back(bdele);
    }
    bound_set_map.clear();

    for(std::vector<Boundary<std::string,double,int>>::iterator itlm = load_set_map.begin(); itlm < load_set_map.end(); itlm++){
        Boundary<set*,double,int> ldele;
        ldele.dims = itlm->dims; ldele.is_alldim = itlm->is_alldim;
        ldele.is_zero = itlm->is_zero; ldele.value = itlm->value;
        ldele.nodes = whereset(itlm->nodes,"N");
        load_set.push_back(ldele);
    }
    load_set_map.clear();

    for(std::vector<Boundary<std::string,double,double>>::iterator itdlm = dload_set_map.begin(); itdlm != dload_set_map.end(); itdlm++){
        Boundary<set*,double,double> dldele;
        dldele.dims = itdlm->dims; dldele.value = itdlm->value;
        if(itdlm->nodes.size() == 0){
            for(std::vector<set>::iterator itele = ELset.begin(); itele != ELset.end(); itele++){
                dldele.nodes = whereset(itele->name,"EL");
                dload_set.push_back(dldele);
            }
        }
        else{
            dldele.nodes = whereset(itdlm->nodes, "EL");
            dload_set.push_back(dldele);
        }
    }
    dload_set_map.clear();
}

void asb_manager::init_K(){
    const int dof = 3 * xyz_coord.size();
    int* invcolume = new int[dof * dof];
    int ROW_I;

    if(K_mat.val == NULL){
        init_K_symbolic();
    }
    else{
        for(int i = 0; i < K_mat.row_st[dof+1]-2; i++){
            K_mat.val[i] = 0.;
        }
    }

    for(int i = 0; i < dof; i++){
        ROW_I = i * dof;
        for(int j = K_mat.row_st[i]-1; j < K_mat.row_st[i+1]-1; j++){
            invcolume[ROW_I + K_mat.col[j]-1] = j;
        }
    }

    for(std::vector<c3d4>::iterator itele = elements.begin(); itele != elements.end(); itele++){
        if(itele->Nodetag.size() == 0){
            continue;
        }
        asb_K(*itele, invcolume);
    }
    delete[] invcolume;
};

void asb_manager::initialize(int mode){
    std::ifstream inp_file(inp_file_name);
    read_manager(inp_file);
    init_ele();
    init_bnd();
    init_K();
    // for(int i = 0; i < KdE_mat.at(0).outerSize(); i++){
    //     printf("col: %d ", i);
    //     for(Eigen::SparseMatrix<double>::InnerIterator it(KdE_mat.at(0),i); it; ++it){
    //         printf("%.4f ",it.value());
    //     }
    //     printf("\n");
    // }
    getFout();
    // for(Eigen::VectorXd::iterator it = Fout.begin(); it != Fout.end(); it++){
    //     printf("%.4f ", *it);
    // }
    // printf("\n");
    addboundry(mode);
};

void asb_manager::read_manager(std::ifstream &file_stream){
    std::string wkstr;
    std::string wkwd;
    bool is_done = false;
    bool is_end = false;
    do{
        is_end = getlmsg(file_stream,wkstr);
        if(is_end){
            printf("reach the end of file \n");
        }
    } while (wkstr.front() != '*' && !is_end);

    while(!is_done && !is_end){
        if(wkstr.front() == '*'){
            wkstr.erase(wkstr.begin());
        }
        else if(wkstr.front() == ','){
            is_end = getlmsg(file_stream, wkstr);
            continue;
        }
        else{
            printf("expect a '*' before keyword \n");
        }
        first_wd(wkstr,wkwd);
        if(!_strcmpi(wkwd.data(),"Heading")){
            printf("start reading ... \n");
            is_end = getlmsg(file_stream, wkstr);
        }
        else if(!_strcmpi(wkwd.data(),"Node")){
            printf("reading node ... \n");
            wkstr = read_coord(file_stream, xyz_coord);
        }
        else if(!_strcmpi(wkwd.data(),"ENDSTEP")){
            printf("end read\n");
            is_done = true;
        }
        else if(!_strcmpi(wkwd.data(),"NSET")){
            printf("reading Node Sets ... \n");
            first_wd(wkstr,wkwd);
            int posi = wkwd.find('=');
            wkstr = read_Nset(file_stream, Nset, wkwd.substr(posi+1,wkwd.size()-posi));
        }
        else if(!_strcmpi(wkwd.data(),"ELSET")){
            printf("reading Element Sets ... \n");
            first_wd(wkstr,wkwd);
            int posi = wkwd.find('=');
            wkstr = read_Elset(file_stream, ELset, wkwd.substr(posi+1,wkwd.size()-posi));
        }
        else if(!_strcmpi(wkwd.data(),"Element")){
            printf("reading Elements ... \n");
            wkstr = read_Element(file_stream, elements);
        }
        else if(!_strcmpi(wkwd.data(),"SolidSection")){
            printf("reading Solid Section ... \n");
            first_wd(wkstr,wkwd);
            int posi = wkwd.find("=") + 1;
            std::pair<std::string,std::string> map;
            map.first = wkwd.substr(posi, wkwd.size()-posi);

            first_wd(wkstr,wkwd);
            posi = wkwd.find("=") + 1;
            map.second = wkwd.substr(posi, wkwd.size()-posi);
            is_end = getlmsg(file_stream,wkstr);

            mat_map.push_back(map);
        }
        else if(!_strcmpi(wkwd.data(),"Material")){
            printf("reading Material ... \n");
            first_wd(wkstr,wkwd);
            int posi = wkwd.find("=") + 1;
            wkstr = read_mater(file_stream, Mater_lib, wkwd.substr(posi, wkwd.size()-posi));
        }
        else if(!_strcmpi(wkwd.data(),"BOUNDARY")){
            printf("reading Boundary ... \n");
            wkstr = read_boundary(file_stream, bound_set_map, bound_node);
        }
        else if(!_strcmpi(wkwd.data(),"Step") || !_strcmpi(wkwd.data(),"OutPut")){
            is_end = getlmsg(file_stream, wkstr);
        }
        else if(!_strcmpi(wkwd.data(),"Static")){
            is_end = getlmsg(file_stream, wkstr);
            is_end = getlmsg(file_stream, wkstr);
        }
        else if(!_strcmpi(wkwd.data(), "Measure")){
            printf("reading measure ... \n");
            wkstr = read_measure(file_stream, probe_list, u_real_probe);
        }
        else if(!_strcmpi(wkwd.data(),"CLOAD")){
            printf("reading centre load ... \n");
            wkstr = read_cload(file_stream, load_set_map, load_node);
        }
        else if(!_strcmpi(wkwd.data(),"DLOAD")){
            printf("reading distribution load ... \n");
            wkstr = read_dload(file_stream, dload_set_map, dload_ele);
        }
        else{
            printf("syntax error: no such keyword <%s> \n", wkwd.data());
            is_end = getlmsg(file_stream,wkstr);
            is_done = true;
        }
        if(is_end){printf("reach end of file \n");};
    }
    file_stream.close();
};

void asb_manager::asb_K(c3d4 &ele, int* invcol){
    const int dof = 3 * xyz_coord.size();
    int ROW_I;
    double K[12][12];
    ele.getK_e(K);

    if(K_mat.type == 2 || K_mat.type == -2){
        int I, J;
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
                            K_mat.val[/**(invcol+ROW_I+I+n)*/invcol[ROW_I+I+n]] += K[3*j + m][3*i + n];
                        }
                    }
                }
                else if(I < J){
                    for(int m = 0; m < 3; m++){
                        ROW_I = (I+m) * dof;
                        for(int n = 0; n < 3; n++){
                            K_mat.val[/**(invcol+ROW_I+J+n)*/invcol[ROW_I+J+n]] += K[3*i + m][3*j + n];
                        }
                    }
                }
                else{
                    for(int m = 0; m < 3; m++){
                        ROW_I = (I+m) * dof;
                        for(int n = m; n < 3; n++){
                            K_mat.val[/**(invcol+ROW_I+I+n)*/invcol[ROW_I+I+n]] += K[3*i + m][3*i + n];
                        }
                    }
                }
            }
        }
    }
    
};

void asb_manager::getFout(){    
    const int dof = 3 * xyz_coord.size();
    double grav = 0.0;
    double rho = 0.0;
    double vol = 0.0;
    double G = 0.0;
    double norm = 0.0;
    double n_vec[3];
    Fout = new double[dof];

    for(int i = 0; i < dof; i++){
        Fout[i] = 0.;
    }

    int posi = 0;
    bool isz = false;
    for(std::vector<Boundary<set*,double,int>>::iterator itset = load_set.begin(); itset != load_set.end(); itset++){
        if(itset -> is_alldim){
            for(std::vector<int>::iterator it = itset->nodes->sets.begin(); it != itset->nodes->sets.end(); it++){
                posi = 3 * (*it - 1);
                for(int i = 0; i < 3; i++){
                    Fout[posi] += itset->value;
                    posi++;
                }
            }
        }
        else{
            for(std::vector<int>::iterator it = itset->nodes->sets.begin(); it != itset->nodes->sets.end(); it++){
                posi = 3 * (*it - 1);
                for(std::vector<int>::iterator itdim = itset->dims.begin(); itdim != itset->dims.end(); itdim++){
                    Fout[posi+*itdim-1] += itset->value;
                }
            }
        }
    }

    for(std::vector<Boundary<int,double,int>>::iterator itnode = load_node.begin(); itnode != load_node.end(); itnode++){
        if(itnode->is_alldim){
            posi = 3*(itnode->nodes-1);
            for(int i = 0; i < 3; i++){
                Fout[posi] = itnode->value;
                posi++;
            }
        }
        else{
            posi = 3*(itnode->nodes-1);
            for(std::vector<int>::iterator itdim = itnode->dims.begin(); itdim != itnode->dims.end(); itdim++){
                Fout[posi+*itdim-1] = itnode->value;
            }
        }
    }

    for(std::vector<Boundary<set*,double,double>>::iterator itdset = dload_set.begin(); itdset != dload_set.end(); itdset++){
        grav = itdset->value;
        norm = 0.0;
        for(int i = 0; i < 3; i++){
            norm += itdset->dims.at(i) * itdset->dims.at(i);
        }
        norm = sqrt(norm);
        if(norm < 1E-12){throw;};
        for(int i = 0; i < 3; i++){
            n_vec[i] = itdset->dims.at(i) / norm;
        }

        for(std::vector<int>::iterator itele = itdset->nodes->sets.begin(); itele != itdset->nodes->sets.end(); itele++){
            vol = elements.at(*itele-1).vol;
            rho = elements.at(*itele-1).mat->rho;
            G = vol * rho * grav / 4;
            for(std::vector<int>::iterator itelenode = elements.at(*itele-1).Nodetag.begin(); itelenode != elements.at(*itele-1).Nodetag.end(); itelenode++){
                posi = 3 * (*itelenode - 1);
                for(int i = 0; i < 3; i++){
                    Fout[posi] += G * n_vec[i];
                    posi++;
                }
            }
        }   
    }

    for(std::vector<Boundary<int,double,double>>::iterator itdele = dload_ele.begin(); itdele != dload_ele.end(); itdele++){
        grav = itdele->value;
        for(int i = 0; i < 3; i++){
            norm += itdele->dims.at(i) * itdele->dims.at(i);
        }
        norm = sqrt(norm);
        if(norm < 1E-12){throw;};
        for(int i = 0; i < 3; i++){
            n_vec[i] = itdele->dims.at(i) / norm;
        }

        vol = elements.at(itdele->nodes-1).vol;
        rho = elements.at(itdele->nodes-1).mat->rho;
        G = vol * rho * grav / 4.0;
        for(std::vector<int>::iterator itelenode = elements.at(itdele->nodes-1).Nodetag.begin(); itelenode != elements.at(itdele->nodes-1).Nodetag.end(); itelenode++){
            posi = 3 * (*itelenode - 1);
            for(int i = 0; i < 3; i++){
                Fout[posi] += G * n_vec[i];
                posi++;
            }
        }
    }
};

set* asb_manager::whereset(const std::string &name, const char* mode){
    if(!_strcmpi(mode,"N")){
        for(std::vector<set>::iterator itn = Nset.begin(); itn < Nset.end(); itn ++){
            if(!strcmp(itn->name.data(),name.data())){
                return &*itn;
            }
        }
    }
    else if(!_strcmpi(mode,"EL")){
        for(std::vector<set>::iterator itel = ELset.begin(); itel < ELset.end(); itel ++){
            if(!strcmp(itel->name.data(),name.data())){
                return &*itel;
            }
        }
    }
    else{
        printf("option: %s is not supported",mode);
    }
    return nullptr;
};

Material* asb_manager::whereset(const std::string &name){
    for(std::vector<Material>::iterator itma = Mater_lib.begin(); itma < Mater_lib.end(); itma ++){
        if(!strcmp(itma->name.data(),name.data())){
            return &*itma;
        }
    }
    return nullptr;
}

void asb_manager::solve(){
    const int dof = (fnode_list.size() != 0) ? fnode_list.size() : 3*xyz_coord.size();
    int mtype = K_mat.type;
    int nrhs = 1;
    uvw_ans = new double[dof];
    void* ptsolver[64];
    pardiso_cfg para = pardiso_cfg(mtype);
    double ddum;

    for(int i = 0; i < 64; i++){
        ptsolver[i] = 0; 
    }

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
    printf("solve completed ...\n");

    phase = -11;
    pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(K_mat.type), &phase, &dof, K_mat.val, 
            K_mat.row_st, K_mat.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), Fout, &ddum, &(para.error));
};

void asb_manager::addboundry(int mode){
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

        for(int i = 0; i < 3*xyz_coord.size(); i++){
            fnode_list.push_back(i);
        }
        for(std::vector<std::pair<int,double>>::reverse_iterator itsub = sub_list.rbegin(); itsub < sub_list.rend(); ++itsub){
            fnode_list.erase(fnode_list.begin() + itsub->first);
        }

        this->sub_mat_vec();

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

void asb_manager::sub_mat_vec(bool F){
    // int subsize = 0;
    // double val = 0.0;
    // Eigen::SparseMatrix<double> K_wkmat;
    // Eigen::VectorXd F_wk;
    // subsize = fnode_list.size();
    // K_wkmat.resize(subsize,subsize);
    // if(F){F_wk.resize(subsize);}

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
    //             K_wkmat.coeffRef(i,j) = val;
    //         }
    //         j++;
    //     }
    //     i++;
    // }
    // K_mat.resize(subsize,subsize);
    // K_mat = K_wkmat;
    // if(F){
    //     Fout.resize(subsize);
    //     Fout = F_wk;
    // }
}

void asb_manager::init_K_symbolic(){
    const int dof = xyz_coord.size();
    int nznv = 0;
    int dignzn = 0;
    int ROW_I, I3;
    bool* syb_K_mat = new bool[dof * dof];
    
    if(K_mat.type == 2 || K_mat.type == -2){
        for(int i = 0; i < dof; i++){
            ROW_I = i * dof;
            for(int j = i; j < dof; j++){
                syb_K_mat[ROW_I + j] = false;
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
        const int nzn = 9*nznv - 3*dignzn;
        K_mat.val = new double[nzn];
        K_mat.col = new int[nzn];
        K_mat.row_st = new int[3*dof+1];

        double* value = K_mat.val;
        int* colume = K_mat.col;
        int* rowindex = K_mat.row_st;
        int icol = 0;
        rowindex[0] = 1;
        bool is_dig;
        for(int i = 0; i < dof; i++){
            I3 = 3 * i;
            ROW_I = i * dof;
            is_dig = false;
            
            for(int j = i; j < dof; j++){
                if(syb_K_mat[ROW_I + j]){
                    colume[icol] = 3*j + 1;
                    colume[icol+1] = 3*j + 2;
                    colume[icol+2] = 3*j + 3;

                    value[icol] = 0.;
                    value[icol+1] = 0.;
                    value[icol+2] = 0.;

                    if(i == j){
                        is_dig = true;
                    }
                    icol += 3;
                }
            }
            rowindex[I3 + 1] = icol + 1;

            if(is_dig){
                for(int j = 1; j < 3; j++){
                    for(int k = rowindex[I3] + j - 1; k < rowindex[I3 + 1] - 1; k++, icol++){
                        colume[icol] = colume[k];
                        value[icol] = 0.;
                    }
                    rowindex[I3 + j + 1] = icol + 1;
                }
            }
            else{
                for(int j = 1; j < 3; j++){
                    for(int k = rowindex[I3] - 1; k < rowindex[I3 + 1] - 1; k++, icol++){
                        colume[icol] = colume[k];
                        value[icol] = 0.;
                    }
                    rowindex[I3 + j + 1] = icol + 1;
                }
            }
        }
    }
    else{
        printf("matrix type %i haven't supported yet",K_mat.type);
    }
    delete[] syb_K_mat;

}