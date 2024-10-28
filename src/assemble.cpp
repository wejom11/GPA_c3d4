#include <stdio.h>
#include "assemble.h"

void asb_manager::read_manager(std::ifstream &file_stream){
    std::string wkstr;
    std::string wkwd;
    bool is_done = false;
    bool is_end = false;

    std::vector<double*> xyz_coord;
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
            is_end = getlmsg(file_stream, wkstr);
            continue;
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
        // else if(!_strcmpi(wkwd.data(),"NSET")){
        //     printf("reading Node Sets ... \n");
        //     first_wd(wkstr,wkwd);
        //     int posi = wkwd.find('=');
        //     wkstr = read_Nset(file_stream, Nset, wkwd.substr(posi+1,wkwd.size()-posi));
        // }
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
        // else if(!_strcmpi(wkwd.data(),"SolidSection")){
        //     printf("reading Solid Section ... \n");
        //     first_wd(wkstr,wkwd);
        //     int posi = wkwd.find("=") + 1;
        //     std::pair<std::string,std::string> map;
        //     map.first = wkwd.substr(posi, wkwd.size()-posi);

        //     first_wd(wkstr,wkwd);
        //     posi = wkwd.find("=") + 1;
        //     map.second = wkwd.substr(posi, wkwd.size()-posi);
        //     is_end = getlmsg(file_stream,wkstr);

        //     mat_map.push_back(map);
        // }
        // else if(!_strcmpi(wkwd.data(),"Material")){
        //     printf("reading Material ... \n");
        //     first_wd(wkstr,wkwd);
        //     int posi = wkwd.find("=") + 1;
        //     wkstr = read_mater(file_stream, Mater_lib, wkwd.substr(posi, wkwd.size()-posi));
        // }
        // else if(!_strcmpi(wkwd.data(),"BOUNDARY")){
        //     printf("reading Boundary ... \n");
        //     wkstr = read_boundary(file_stream, bound_set_map, bound_node);
        // }
        else if(!_strcmpi(wkwd.data(),"Step") || !_strcmpi(wkwd.data(),"OutPut")){
            is_end = getlmsg(file_stream, wkstr);
        }
        else if(!_strcmpi(wkwd.data(),"Static")){
            is_end = getlmsg(file_stream, wkstr);
            is_end = getlmsg(file_stream, wkstr);
        }
        // else if(!_strcmpi(wkwd.data(), "Measure")){
        //     printf("reading measure ... \n");
        //     wkstr = read_measure(file_stream, probe_list, u_real_probe);
        // }
        // else if(!_strcmpi(wkwd.data(),"CLOAD")){
        //     printf("reading centre load ... \n");
        //     wkstr = read_cload(file_stream, load_set_map, load_node);
        // }
        // else if(!_strcmpi(wkwd.data(),"DLOAD")){
        //     printf("reading distribution load ... \n");
        //     wkstr = read_dload(file_stream, dload_set_map, dload_ele);
        // }
        else{
            printf("syntax error: no such keyword <%s> \n", wkwd.data());
            is_end = getlmsg(file_stream,wkstr);
            // is_done = true;
        }
        if(is_end){printf("reach end of file \n");};
    }
    file_stream.close();

    Mater_lib.resize(1);
    Mater_lib.front() = Material();

    int nnum = xyz_coord.size();
    int elenum = elements.size();
    int idum;
    double mid_x, mid_y;
    xy_coord = new double[2*(nnum + elenum)]{0.};
    for(int i = 0; i < nnum; i++){
        xy_coord[2*i] = xyz_coord.at(i)[0];
        xy_coord[2*i + 1] = xyz_coord.at(i)[1];

        delete[] xyz_coord[i]; xyz_coord[i] = nullptr;
    }
    xyz_coord.clear(); xyz_coord.shrink_to_fit();

    for(int i = nnum; i < nnum + elenum; i++){
        mid_x = 0.;
        mid_y = 0.;
        for(int j = 0; j < 4; j++){
            idum = 2 * (elements.at(i-nnum).Nodetag[j] - 1);
            mid_x += xy_coord[idum];
            mid_y += xy_coord[idum + 1];
        }
        xy_coord[2*i] = mid_x;
        xy_coord[2*i + 1] = mid_y;
        elements.at(i-nnum).Nodetag[8] = i + 1;
        elements.at(i-nnum).mat = &Mater_lib.front();
    }

    int* map_node = new int[nnum + elenum]{-1};
    short int* node_type = new short int[nnum + elenum]{0};
    int i, j;
    set* f_eleset = whereset("fluid", "EL");
    for(std::vector<int>::iterator itfele = f_eleset->sets.begin(); itfele != f_eleset->sets.end(); itfele++){
        elements.at(*itfele - 1).set_fs(true);
        for(i = 0; i < 4; i++){
            node_type[elements.at(*itfele-1).Nodetag[i] - 1] = 3;
        }
        for(i = 4; i < 9; i++){
            if(node_type[elements.at(*itfele-1).Nodetag[i] - 1] < 2){
                node_type[elements.at(*itfele-1).Nodetag[i] - 1] = 2;
            }
        }
    }

    set* s_eleset = whereset("solid", "EL");
    for(std::vector<int>::iterator itsele = s_eleset->sets.begin(); itsele != s_eleset->sets.end(); itsele++){
        elements.at(*itsele - 1).set_fs(false);
        for(i = 0; i < 9; i++){
            if(node_type[elements.at(*itsele-1).Nodetag[i] - 1] < 2){
                node_type[elements.at(*itsele-1).Nodetag[i] - 1] = 1;
            }
            else{
                wall_node.push_back(elements.at(*itsele-1).Nodetag[i]);
            }
        }
    }
    wall_node.shrink_to_fit();
    f_eleset->sets.clear(); f_eleset->sets.shrink_to_fit();
    s_eleset->sets.clear(); s_eleset->sets.shrink_to_fit();

    set* hs1 = whereset("Source1","EL");
    set* hs2 = whereset("Source2","EL");
    for(std::vector<int>::iterator iths = hs1->sets.begin(); iths != hs1->sets.end(); iths++){
        elements.at(*iths - 1).set_hs(true);
    }
    for(std::vector<int>::iterator iths = hs2->sets.begin(); iths != hs2->sets.end(); iths++){
        elements.at(*iths - 1).set_hs(true);
    }
    hs1->sets.clear(); hs1->sets.shrink_to_fit();
    hs2->sets.clear(); hs2->sets.shrink_to_fit();

    int ticker = 0;
    for(i = 0; i < nnum + elenum; i++){
        if(node_type[i] == 3){
            map_node[i] = ticker;
            ticker++;
        }
    }
    typenode_num[0] = ticker;
    for(i = 0; i < nnum + elenum; i++){
        if(node_type[i] == 2){
            map_node[i] = ticker;
            ticker++;
        }
    }
    typenode_num[1] = ticker - typenode_num[0];
    for(i = 0; i < nnum + elenum; i++){
        if(node_type[i] == 1){
            map_node[i] = ticker;
            ticker++;
        }
    }
    typenode_num[2] = ticker - typenode_num[0] - typenode_num[1];
    delete[] node_type; node_type = nullptr;

    uvtp_ans = new double[4*typenode_num[0] + 3*typenode_num[1] + typenode_num[2]]{0.};
    for(i = 1; i < nnum + elenum + 1; i++){
        uvtp_ans[wherend(i, 2)] = 25.0;
    }

    double* xy_copy = new double[2 * (nnum + elenum)]{0};
    for(i = 0; i < nnum + elenum; i++){
        xy_copy[2*map_node[i]] = xy_coord[2*i];
        xy_copy[2*map_node[i] + 1] = xy_coord[2*i + 1];
    }
    delete[] xy_coord; xy_coord = nullptr;
    xy_coord = xy_copy;
    xy_copy = nullptr;

    int n_tag[9]{0};
    for(i = 0; i < elenum; i++){
        for(j = 0; j < 9; j++){
            n_tag[j] = map_node[elements.at(i).Nodetag[j] - 1] + 1;
        }
        elements.at(i).set_nval(n_tag, xy_coord, uvtp_ans, typenode_num);
    }

    int mnt = 0;
    for(std::vector<int>::iterator itwn = wall_node.begin(); itwn != wall_node.end(); itwn++){
        mnt = map_node[*itwn - 1] + 1;
        *itwn = mnt;
    }

    inlet_ele = {{3,4},{5,4},{3894,2},{1,4},{3897,2}};
    outlet_ele = {{4,2},{6,2},{3893,4},{2,2},{3896,4}};

    delete[] map_node; map_node = nullptr;    
};

void asb_manager::init_ele(){
    double H = 0.04, B = 10;
    int y_numele = 4, x_numele = 1000;
    Material mat1;
    Mater_lib.push_back(mat1);

    double dH = H / y_numele, dB = B / x_numele;
    int y_numnode = 2*y_numele + 1, x_numnode = 2*x_numele + 1,
        i , j, ROW, ROW1, ROW2, ROW3, ROW4;
    double x_i;

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

    inlet_ele.resize(y_numele);
    outlet_ele.resize(y_numele);
    face_node.resize(y_numnode);
    for(i = 0; i < x_numele + 1; i++){
        ROW = 2 * i * (y_numele + 1);
        x_i = dB * i;
        for(j = 0; j < y_numele + 1; j++){
            xy_coord[ROW + 2*j] = x_i;
            xy_coord[ROW + 2*j + 1] = dH * j;
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

    elements.resize(x_numele * y_numele, P9SF(&Mater_lib.front()));
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

            elements.at(ROW + j).set_nval(nt, xy_coord, uvtp_ans, typenode_num);
        }
    }

    for(i = 0; i < y_numele; i++){
        inlet_ele.at(i) = {i + 1, 4};
        outlet_ele.at(i) = {y_numele * (x_numele - 1) + i + 1, 2};
    }

    for(i = 0; i < y_numele + 1; i++){
        uvtp_ans[4*i + 3] = P_in;
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

    //  int m = (inlet_node.size() + 1) / 2-1;
    //  for(int i = 1; i < inlet_node.size()-1; i++){
    //    uvtp_ans[wherend(inlet_node.at(i), 0)] = 0.2 * (-pow(i - m, 2) / pow(m, 2) + 1.0);
    //  }

    //  for(int i = 1; i < outlet_node.size()-1; i++){
    //    uvtp_ans[wherend(outlet_node.at(i), 0)] = 0.2 * (-pow(i - m, 2) / pow(m, 2) + 1.0);
    //  }

    face_node = {2, 28, 5, 33, 8, 38, 11, 43, 14};
};

void asb_manager::init_mesh(){
    int i = 0, j = 0;
    int* map_node = nullptr;

    fsinterval = 1.0;
    ifw.resize(14);
    ifw = {4.0,2.0,3.0,3.0,4.0,4.0,4.0,3.0,4.0,4.0,4.0,4.0,4.0,4.0};
    ofw.resize(16);
    ofw = {4.0,4.0,4.0,2.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0};

    Module md(fsinterval, ifw, ofw);
    md.Init_geo();
    md.Get_num_element_xy();
    md.Get_nodes();
    md.Get_elements();

    Material mat_f;
    Material mat_s;
    Mater_lib.resize(2);
    Mater_lib.at(0) = mat_f;
    Mater_lib.at(1) = mat_s;

    int node_dof = md.nodes.size();
    typenode_num[0] = md.number_nodes_P.size();
    typenode_num[2] = md.number_nodes_T.size();
    typenode_num[1] = node_dof -typenode_num[0] - typenode_num[2];

    try{
        map_node = new int[node_dof]{-1};
        xy_coord = new double[2*node_dof]{0.};
    }
    catch(const std::bad_alloc &e){
        std::cerr << e.what() << '\n';
        map_node = nullptr;
        exit(1);
    }
    
    if(map_node != nullptr){
        j = 0;
        for(std::vector<int>::iterator itn = md.number_nodes_P.begin(); itn != md.number_nodes_P.end(); itn++){
            map_node[*itn] = j;
            j++;
        }

        j = typenode_num[0];
        for(i = 0; i < node_dof; i++){
            if((check(i, md.number_nodes_P) == -1)&&(check(i, md.number_nodes_T) == -1)){
                map_node[i] = j;
                j++;
            }
        }

        for(std::vector<int>::iterator itn = md.number_nodes_T.begin(); itn != md.number_nodes_T.end(); itn++){
            map_node[*itn] = j;
            j++;
        }
    }

    for(i = 0; i < node_dof; i++){
        xy_coord[2*map_node[i]] = md.nodes.at(i).x;
        xy_coord[2*map_node[i] + 1] = md.nodes.at(i).y;
    }

    this->elements.resize(md.elements.size());
    int node_tag[9]{0};
    for(std::vector<Element>::iterator itele = md.elements.begin(); itele != md.elements.end(); i++){
        P9SF ele(&(itele->type == 0 ? mat_f:mat_s), !itele->type);
        for(i = 0; i < 9; i++){
            node_tag[i] = map_node[itele->nodes[i]] + 1;
        }

        ele.set_nval(node_tag, xy_coord, uvtp_ans, typenode_num);
        i++;
    }

    delete[] map_node; map_node = nullptr;
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

        if(itele->is_hs){
            asb_KF_s(*itele);
        }
    }
};

void asb_manager::initialize(int mode){
    // std::ifstream inp_file(inp_file_name);
    // read_manager(inp_file);
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
    int id[3]{0}, face[3]{0};
    if(Fout){
        for(int i = 0; i < dof; i++){
            Fout[i] = 0;
        }
    }
    else{
        Fout = new double[dof]{0.0};
    }

    for(std::vector<std::pair<int,int>>::iterator innd = inlet_ele.begin(); innd != inlet_ele.end(); innd++){
        face[0] = innd->second;
        face[1] = innd->second + 4;
        face[2] = innd->second % 4 + 1;
        for(int i = 0; i < 3; i++){
            id[i] = wherend(elements.at(innd->first - 1).at_Nodetag(face[i] - 1), 0);
        }
        elements.at(innd->first - 1).P.at(face[0]-1)[0] = P_in;
        elements.at(innd->first - 1).P.at(face[2]-1)[0] = P_in;
        elements.at(innd->first - 1).get_Fout(Fout, innd->second, id);
    }

    for(std::vector<std::pair<int,int>>::iterator outnd = outlet_ele.begin(); outnd != outlet_ele.end(); outnd++){
        face[0] = outnd->second;
        face[1] = outnd->second + 4;
        face[2] = outnd->second % 4 + 1;
        for(int i = 0; i < 3; i++){
            id[i] = wherend(elements.at(outnd->first - 1).at_Nodetag(face[i] - 1), 0);
        }
        elements.at(outnd->first - 1).P.at(face[0] - 1)[0] = 0;
        elements.at(outnd->first - 1).P.at(face[2] - 1)[0] = 0;
        elements.at(outnd->first - 1).get_Fout(Fout, outnd->second, id);
    }

    //for(int i = 1; i < inlet_node.size()-1; i++){
    //    Fout[wherend(inlet_node.at(i), 0)] = 0.01;
    //}
    //for (int i = 1; i < outlet_node.size() - 1; i++) {
    //    Fout[wherend(outlet_node.at(i), 0)] = 0.01;
    //}

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

int asb_manager::wherend(int nodetag, int uvtp){
    if(nodetag <= typenode_num[0]){
        return 4 * nodetag - 4 + uvtp;
    }
    else if(nodetag > typenode_num[0] && nodetag <= typenode_num[0] + typenode_num[1]){
        if(uvtp > 2){
            return -1;
        }
        return 3 * nodetag + typenode_num[0] - 3 + uvtp;
    }
    else{
        if(uvtp != 2){
            // std::cerr << "ERROR: this node is third type node!" << std::endl;
            return -1;
        }
        return typenode_num[0] * 3 + typenode_num[1] * 2 + nodetag - 1;
    }
};

bool asb_manager::solve(){
    this->init_KF();
    const int dof = 4*typenode_num[0] + 3*typenode_num[1] + typenode_num[2];
    const int max_iteration = 15;
    int mtype = K_mat.type;
    int nrhs = 1, itr_time = 0;
    int phase;
    double* uvtp_itr = new double[dof]{0};
    void* ptsolver[64];
    pardiso_cfg para = pardiso_cfg(mtype);
    double ddum, norm, norm0;
    bool solved = false;

    for(int i = 0; i < 64; i++){
        ptsolver[i] = 0; 
    }

    // std::ofstream file("F:/work_data/competition_fuild_opt/GPA_c3d4/answer.txt");
    // for (int i = 0; i < dof; i++) {
    //     file << K_mat.val[match(i + 1, i + 1, K_mat)] << "    ";
    // }
    // file.close();

    phase = 11;
    pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(K_mat.type), &phase, &dof, K_mat.val, 
            K_mat.row_st, K_mat.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), &ddum, &ddum, &(para.error));
    if(para.error != 0 ){
        printf ("ERROR during symbolic factorization: %i \n", para.error);
    }

    do{
        norm = 0.;
        //for (int j = 1; j < inlet_node.size() - 1; j++) {
        //    printf("%.4f %.4f %.4f\n", Fint[wherend(inlet_node[j], 0)], Fint[wherend(inlet_node[j], 1)], Fint[wherend(inlet_node[j], 2)]);
        //}
        this->getFout();

        for(int i = 0; i < dof; i++){
            Fint[i] = Fout[i] - Fint[i];
        }
        this->addboundry();
        for(int i = 0; i < dof; i++){
            norm += Fint[i] * Fint[i];
        }

        if (itr_time == 0){ norm0 = norm / 1000 < 1E-6 ? norm / 1000 : 1E-6; }
        if(norm < norm0){
            printf("itertion time: %i, Residual: %.8f\n", itr_time+1, norm);
            printf("done solving !\n");
            solved = true;
            break;
        }
        else if(norm > 1E15){
            printf("itertion time: %i, Residual: %.8f\n", itr_time + 1, norm);
            printf("Residual too large!\n");
            solved = false;
            break;
        }
        else{
            printf("itertion time: %i, Residual: %.8f\n", itr_time+1, norm);

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
        }

        this->init_KF();
        itr_time++;
    }while(itr_time < max_iteration);

    phase = -11;
    pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(K_mat.type), &phase, &dof, K_mat.val, 
            K_mat.row_st, K_mat.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), Fint, &ddum, &(para.error));

    delete[] uvtp_itr; uvtp_itr = nullptr;

    return solved;
};

void asb_manager::addboundry(){
    int i = 0, j = 0, k = 0;

    for(std::vector<int>::iterator itwx = wall_node.begin(); itwx != wall_node.end(); itwx++){
        int where_x = wherend(*itwx, 0);
        int where_y = wherend(*itwx, 1);
        if (where_x == -1 || where_y == -1) {
            printf("there\n");
        }

        K_mat.val[match(where_x+1, where_x+1, K_mat)] = 1E10;
        Fint[where_x] = 0;
        K_mat.val[match(where_y + 1, where_y + 1, K_mat)] = 1E10;
        Fint[where_y] = 0;
    }

    // for(i = 1; i < inlet_node.size()-1; i++){
        // int where_x = wherend(inlet_node.at(i), 1);
        //int where_y = wherend(inlet_node.at(i), 1);

        // K_mat.val[match(where_x+1, where_x+1, K_mat)] = 1E10;
        // Fint[where_x] = 0;
        //K_mat.val[match(where_y + 1, where_y + 1, K_mat)] = 1E10;
        //Fint[where_y] = 0;
    // }
    // for(i = 1; i < outlet_node.size()-1; i++){
        // int where_x = wherend(outlet_node.at(i), 1);
        //int where_y = wherend(outlet_node.at(i), 1);

        // K_mat.val[match(where_x+1, where_x+1, K_mat)] = 1E10;
        // Fint[where_x] = 0;
        // K_mat.val[match(where_y + 1, where_y + 1, K_mat)] = 1E10;
        // Fint[where_y] = 0;
    // }

}

void asb_manager::init_K_symbolic(){
    const int dof = typenode_num[0] + typenode_num[1] + typenode_num[2];
    int nzn44 = 0, nzn43 = 0, nzn41 = 0,
                   nzn33 = 0, nzn31 = 0,
                              nzn11 = 0;
    int ROW_I, I3, I4, I;
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

        // for(int i = 2; i < inlet_node.size()-1; i++){
        //     ROW_I = (inlet_node[i] - 1) * dof;
        //     ROW_Is1 = (inlet_node[i - 1] - 1) * dof;
        //     for(int j = 0; j < dof; j++){
        //         if(syb_K_mat[ROW_I + j]){
        //             syb_K_mat[ROW_Is1 + j] = true;
        //         }
        //     }
        // }

        // ROW_I = (inlet_node.at(inlet_node.size()-2) - 1) * dof;
        // ROW_I = (inlet_node.at(inlet_node.size()-2) - 1) * dof;
        // for(int i = 0; i < inlet_node.size(); i++){
        //     syb_K_mat[ROW_I + inlet_node.at(i) - 1] = true;
        // }
        // for(int i = 0; i < outlet_node.size(); i += 2){
        //     syb_K_mat[ROW_I + outlet_node.at(i) - 1] = true;
        // }

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

void asb_manager::get_face(int m, int x_numele, int y_numele) {
    int i;

    for (i = 0; i < y_numele + 1; i++) {
        face_node.at(2 * i) = m * (y_numele + 1) + i + 1;
    }
    for (i = 0; i < y_numele; i++) {
        face_node.at(2 * i + 1) = typenode_num[0] + m * y_numele + i + 1;
    }
}

void asb_manager::write(std::string file_name){
    std::ofstream out_file_io; 
    std::string fn = "E:/personal data/mesh_generation/FEMT_practice/Tutorial Program/PP_process/" + file_name + ".out";
    out_file_io.open(fn, std::ios_base::out);
    
    out_file_io.setf(out_file_io.scientific);
    out_file_io.precision(15);
    out_file_io << "Geometry" << std::endl;
    out_file_io << std::endl;

    // node
    int nnum = typenode_num[0] + typenode_num[1] + typenode_num[2];
    int elenum = elements.size();
    int i = 0, j = 0, p = 0;
    double val = 0.;
    out_file_io << " COORDINATES " << nnum << std::endl;
    out_file_io << "! NODE  X   Y   Z" << std::endl;
    for(i = 0; i < nnum; i++){
        out_file_io << "    " << i+1 << "    " << xy_coord[2*i] << "    " << xy_coord[2*i+1] << "    " << 0.0 << std::endl;
    }
    out_file_io << std::endl;

    //elements
    out_file_io << " ELEMENTS " << elenum << std::endl;
    out_file_io << "    ELEMENT_NODES" << std::endl;
    for(i = 0; i < elenum; i++){
        out_file_io << "    " << i+1;
        for(j = 0; j < 9; j++){
            out_file_io << "    " << elements.at(i).at_Nodetag(j);
        }
        out_file_io << std::endl;
    }
    out_file_io << std::endl;
    out_file_io << "    ELEMENT_MATERIAL" << std::endl;
    out_file_io << std::endl;

    //answer
    out_file_io << " *** NODE_VALUE ***" << std::endl;
    out_file_io << "NODE    U    V    T    P" << std::endl;
    for(i = 0; i < nnum; i++){
        out_file_io << "    " << i+1;
        for(j = 0; j < 4; j++){
            p = wherend(i+1, j);
            val = p > -1 ? uvtp_ans[p] : 0.0;
            out_file_io << "    " << val;
        }
        out_file_io << std::endl;
    }
    out_file_io << std::endl;
    out_file_io << "*END" << std::endl;
}