#include <vector>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "c3d4_ele.h"
#include "assemble.h"
#include "assemble_opt.h"
#include "read.h"
#include "solver.h"

void read_nodelist(std::ifstream &file_stream, std::vector<int> &nlist){
    bool is_done = false;
    bool is_end = false;
    int start,end,step;
    std::string work_str;
    std::string word;
    std::vector<std::string> word_list;
    while(!is_end){
        is_end = getlmsg(file_stream, work_str);
        is_done = false;
        word_list.clear();
        do{
            if(work_str.size() != 0){
                first_wd(work_str,word);
            }
            else{
                is_done = true;
                break;
            }

            if(word.size() != 0){
                word_list.push_back(word);
            }
        }while(!is_done);

        if(word_list.size() == 3){
            start = std::stoi(word_list.at(0));
            end = std::stoi(word_list.at(1));
            step = std::stoi(word_list.at(2));
            for(int i = start; i <= end; i += step){
                nlist.push_back(i);
            }
        }
        else if(word_list.size() == 2){
            start = std::stoi(word_list.at(0));
            end = std::stoi(word_list.at(1));
            for(int i = start; i <= end; i ++){
                nlist.push_back(i);
            }
        }
        else if(word_list.size() == 1){
            nlist.push_back(std::stoi(word_list.front()));
        }
    }
};

int main(){
    // Material mat0 = {3000,0.3};
    // std::vector<double *> xyz0 = {new double[3]{-1.0, 0.0, 0.0},new double[3]{1.0, 0.0, 0.0},new double[3]{0.0, 1.0, 0.0},new double[3]{0.0, 0.5, 0.5}};
    // std::vector<int> nodet = {1,2,3,4};
    // c3d4 ele(xyz0,nodet, & mat0);
    // asb_manager asbma;
    // asbma.getK_asb(ele);
    // delete[] xyz0[0];
    

    // std::cout << wkmsg << std::endl;
    // inp_f.close();
    // first_wd(wkmsg,word);
    // std::cout << word << std::endl;
    // std::cout << wkmsg << std::endl;
    // printf("%.d",1);
    // std::vector<double *> xyz0;
    // {
    // double *x =new double[3]{0.75,1.0,2.0};
    // xyz0.push_back(x);}
    // delete xyz0[0];
    // xyz0.erase(xyz0.begin());
    // printf("%d",xyz0.size());   

    int mode;
    bool is_wri;
    std::string name;
    printf("file name:\n");
    std::cin >> name;
    asb_opt_manager asb("E:/personal data/Summer internship/work_dir/nonliner_FEM_basic/assemble/"+name);
    // asb.initialize();
    // asb.solve();
    printf("is_write, mode:\n");
    std::cin >> is_wri >> mode;

    if(!is_write){
        asb.opt_val(mode);
    }
    else{
        bool err;
        std::vector<int> nl;
        std::ifstream vec_file("E:/personal data/Summer internship/work_dir/nonliner_FEM_basic/assemble/c++/node_list.txt");
        read_nodelist(vec_file, nl);

        asb.initialize(mode);
        std::ofstream out_inp_file;
        out_inp_file.open("../vol_ans.inp",std::ios::app);
        asb.solve(err);
        if(err){
            printf("Excessive memory");
        }
        if(nl.size() == 0){
            asb.wirte(out_inp_file);
        }
        else{
            asb.wirte(out_inp_file, nl);
        }
    }
    // for (int i = 0; i < 4; i++){
    //     for (int j = 0; j < 3; j++){
    //         printf("%.4f ",asb.elements.at(2).xyz0[i][j]);
    //     }
    //     printf("\n");
    // }
    // for(int i = 0; i < 3; i++){
    //     printf("%d ",asb.bound_set.at(0).nodes->sets.at(i));
    // }
    // printf("\n");

    // std::string name;
    // printf("file name: \n");
    // std::cin >> name;
    // asb_manager asb(name.size() == 0 ? "E:\\personal data\\Summer internship\\work_dir\\nonliner_FEM_basic\\assemble\\vol_mesh.inp"
    //                 : name);
    // asb.initialize(1);
    // asb.solve();
    // for(int i = 3*54; i < 3*55; i++){
    //     printf("%.6f ",asb.Fout[i]);
    // }

    // double depsilon[3][3][12];
    // for(int i = 0; i < 3; i++){
    //     for(int j = 0; j < 3; j++){
    //         double temp[12] = {1.0,1.0,1.0,1.0,1.1,1.,2.0,1.,1.,1.,1.,1.};
    //         for(int k = 0; k < 12; k++){
    //             depsilon[i][j][k] = temp[k];
    //         }
    //         // *depsilon[i][j] = *temp;
    //     }
    // }
    // double*** pdep = (double ***) depsilon;
    // printf("%.2f",*((double *)pdep+6+10));
}