#include <string.h>
#include <fstream>
#include <vector>
#include <math.h>
#include "c3d4_ele.h"
#include "read.h"

std::string read_coord(std::ifstream &file_stream, std::vector<double*> &xyz_c){
    std::string work_str;
    std::string word;
    int index;
    double x, y, z;
    bool is_done = false;
    do{
        getlmsg(file_stream, work_str);
        if (work_str.front() == '*'){
            is_done = true;
        }
        else{
            first_wd(work_str,word);
            index = std::stoi(word) - 1;
            first_wd(work_str,word);
            x = std::stod(word);
            first_wd(work_str,word);
            y = std::stod(word);
            first_wd(work_str,word);
            z = std::stod(word);
            if(index == xyz_c.size()){
                xyz_c.push_back(new double[3]{x,y,z});
            }
            else if(index > xyz_c.size()){
                xyz_c.resize(index);
                xyz_c.push_back(new double[3]{x,y,z});
            }
            else{
                if(xyz_c.at(index) == NULL){
                    xyz_c[index] = new double[3]{x,y,z};
                }
                else{
                    throw "Duplicate node";
                }
            }
        }
    } while (! is_done);
    return work_str;
}

std::string read_Element(std::ifstream &file_stream, std::vector<c3d4> &eles){
    std::string work_str;
    std::string word;
    int it_tag[5];
    bool is_done = false;
    do{
        getlmsg(file_stream,work_str);
        if(work_str.front() == '*'){
            is_done =true;
        }
        else{
            for (int i = 0; i < 5; i++){
                first_wd(work_str,word);
                it_tag[i] = std::stoi(word);
            }
            it_tag[0]--;
            if(it_tag[0] == eles.size()){
                eles.push_back(c3d4({},{it_tag[1],it_tag[2],it_tag[3],it_tag[4]}));
            }
            else if(it_tag[0] > eles.size()){
                eles.resize(it_tag[0]);
                eles.push_back(c3d4({},{it_tag[1],it_tag[2],it_tag[3],it_tag[4]}));
            }
            else{
                if(eles.at(it_tag[0]).Nodetag.size() == 0){
                    eles[it_tag[0]] = c3d4({},{it_tag[1],it_tag[2],it_tag[3],it_tag[4]});
                }
                else{
                    throw "Duplicate node";
                }
            }
        }
    } while (!is_done);
    return work_str;
};

std::string read_Nset(std::ifstream &file_stream, std::vector<set> &Nset, std::string name){
    std::string work_str;
    std::string word;
    bool is_done = false;
    Nset.push_back({name,{}});
    do{
        getlmsg(file_stream,work_str);
        if(work_str.front() == '*'){
            is_done =true;
        }
        else{
            while(work_str.size() != 0){
                first_wd(work_str,word);
                if(word.size() == 0){break;};
                Nset.back().sets.push_back(std::stoi(word));
            }
        }
    } while (!is_done);
    return work_str;
};

std::string read_Elset(std::ifstream &file_stream, std::vector<set> &ELset, std::string name){
    std::string work_str;
    std::string word;
    bool is_done = false;
    ELset.push_back({name,{}});
    do{
        getlmsg(file_stream,work_str);
        if(work_str.front() == '*'){
            is_done =true;
        }
        else{
            while(work_str.size() != 0){
                first_wd(work_str,word);
                if(word.size() == 0){break;};
                ELset.back().sets.push_back(std::stoi(word));
            }
        }
    } while (!is_done);
    return work_str;
};

std::string read_mater(std::ifstream &file_stream, std::vector<Material> &mat_lib, std::string name){
    std::string work_str;
    std::string work_str_m;
    std::string word;
    Material mater;
    mater.name = name;
    bool is_done = false;
    getlmsg(file_stream, work_str);
    do{
        if(work_str.front() == '*'){
            work_str_m = work_str.substr(1,work_str.size()-1);
        }
        first_wd(work_str_m,word);
        if(!_strcmpi(word.data(), "elastic")){
            getlmsg(file_stream, work_str);
            first_wd(work_str,word);
            mater.E = std::stod(word);
            
            first_wd(work_str,word);
            mater.mu = std::stod(word);
            getlmsg(file_stream,work_str);
        }
        else if(!_strcmpi(word.data(),"Density")){
            getlmsg(file_stream,work_str);
            first_wd(work_str, word);
            mater.rho = std::stod(word);
            getlmsg(file_stream,work_str);
        }
        else{
            is_done = true;
        }
    }while(!is_done);
    mater.where = mat_lib.size();
    mat_lib.push_back(mater);
    return work_str;
};

std::string read_boundary(std::ifstream &file_stream, std::vector<Boundary<std::string,double,int>> &bsetmap, std::vector<Boundary<int,double,int>> &bnode){
    std::string work_str;
    std::string word;
    std::vector<std::string> word_list;
    bool is_done = false;
    int tag, start, end, step;
    double val;
    do{
        getlmsg(file_stream,work_str);
        if(work_str.front() == '*'){
            is_done = true;
        }
        else{
            first_wd(work_str, word);
            word_list.clear();
            do{
                word_list.push_back(word);
                if(work_str.size() > 0){
                    first_wd(work_str,word);
                }
                else{
                    break;
                }
            }while(true);

            if(word_list.front().front() < 58 && word_list.front().front() > 47){
                if(word_list.size() == 3){
                    tag = std::stoi(word_list[0]);
                    start = std::stoi(word_list[1]);
                    val = std::stod(word_list[2]);
                    bnode.push_back(Boundary<int,double,int>{tag,{start},val});
                }
                else if(word_list.size() == 4){
                    tag = std::stoi(word_list[0]);
                    start = std::stoi(word_list[1]);
                    end = std::stoi(word_list[2]);
                    val = std::stod(word_list[3]);
                    if(start == 1 && end == 3){
                        if(abs(val) < 1e-20){
                            Boundary<int,double,int> bd;
                            bd.is_alldim = true; bd.is_zero = true; bd.nodes = tag;
                            bnode.push_back(bd);
                        }
                        else{
                            bnode.push_back(Boundary<int,double,int> {tag,{},val,false,true});
                        }
                    }
                    else{
                        bnode.push_back(Boundary<int,double,int>{tag,{start,end},val});
                    }
                }
                else if(word_list.size() == 5){
                    tag = std::stoi(word_list[0]);
                    start = std::stoi(word_list[1]);
                    end = std::stoi(word_list[2]);
                    step = std::stoi(word_list[3]);
                    val = std::stod(word_list[4]);
                    std::vector<int> ddim;
                    for(int id = start; id < end; id += step){
                        ddim.push_back(id);
                    }
                    if(ddim.size() == 3){
                        bnode.push_back(Boundary<int,double,int>{tag,{},val,false,true});
                    }
                    else{
                        bnode.push_back(Boundary<int,double,int>{tag,ddim,val});
                    }
                }
                else{
                    printf("syntax error: parameter too much or less");
                }
            }
            else{
                if(word_list.size() == 2){
                    if(word_list.at(1).front() < 58 && word_list.at(1).front() > 47){
                        val = std::stod(word_list[1]);
                        bsetmap.push_back(Boundary<std::string, double,int>{word_list[0],{},val,false,true});
                    }
                    else{
                        if(!_strcmpi(word_list.at(1).data(),"ENCASTRE")){
                            bsetmap.push_back(Boundary<std::string, double,int>{word_list[0],{},0,true,true});
                        }
                        else{
                            printf("syntax error: no such option");
                        }
                    }
                }
                else if(word_list.size() == 3){
                    start = std::stoi(word_list[1]);
                    val = std::stod(word_list[2]);
                    bsetmap.push_back(Boundary<std::string, double,int>{word_list[0],{start},val});
                }
                else if(word_list.size() == 4){
                    start = std::stoi(word_list.at(1));
                    end = std::stoi(word_list.at(2));
                    val = std::stod(word_list.at(3));
                    bsetmap.push_back(Boundary<std::string, double,int>{word_list.front(), {start,end},val});
                }
            }
        }
    }while(!is_done);
    return work_str;
};

std::string read_cload(std::ifstream &file_stream, std::vector<Boundary<std::string,double,int>> &lsetmap, std::vector<Boundary<int,double,int>> &lnode){
    std::string work_str;
    std::string word;
    std::vector<std::string> word_list;
    bool is_done = false;
    int tag, start, end, step;
    double val;
    do{
        getlmsg(file_stream,work_str);
        if(work_str.front() == '*'){
            is_done = true;
        }
        else{
            first_wd(work_str, word);
            word_list.clear();
            do{
                word_list.push_back(word);
                if(work_str.size() > 0){
                    first_wd(work_str,word);
                }
                else{
                    break;
                }
            }while(true);

            if(word_list.front().front() < 58 && word_list.front().front() > 47){
                if(word_list.size() == 3){
                    tag = std::stoi(word_list[0]);
                    start = std::stoi(word_list[1]);
                    val = std::stod(word_list[2]);
                    lnode.push_back(Boundary<int,double,int>{tag,{start},val});
                }
                else if(word_list.size() == 4){
                    tag = std::stoi(word_list[0]);
                    start = std::stoi(word_list[1]);
                    end = std::stoi(word_list[2]);
                    val = std::stod(word_list[3]);
                    if(start == 1 && end == 3){
                        lnode.push_back(Boundary<int,double,int> {tag,{},val,false,true});
                    }
                    else{
                        lnode.push_back(Boundary<int,double,int>{tag,{start,end},val});
                    }
                }
                else if(word_list.size() == 5){
                    tag = std::stoi(word_list[0]);
                    start = std::stoi(word_list[1]);
                    end = std::stoi(word_list[2]);
                    step = std::stoi(word_list[3]);
                    val = std::stod(word_list[4]);
                    std::vector<int> ddim;
                    for(int id = start; id < end; id += step){
                        ddim.push_back(id);
                    }
                    if(ddim.size() == 3){
                        lnode.push_back(Boundary<int,double,int>{tag,{},val,false,true});
                    }
                    else{
                        lnode.push_back(Boundary<int,double,int>{tag,ddim,val});
                    }
                }
                else{
                    printf("syntax error: parameter too much or less");
                }
            }
            else{
                if(word_list.size() == 2){
                    if(word_list.at(1).front() < 58 && word_list.at(1).front() > 47){
                        val = std::stod(word_list[1]);
                        lsetmap.push_back(Boundary<std::string, double,int>{word_list[0],{},val,false,true});
                    }
                    else{
                        printf("syntax error: no such option");
                        
                    }
                }
                else if(word_list.size() == 3){
                    start = std::stoi(word_list[1]);
                    val = std::stod(word_list[2]);
                    lsetmap.push_back(Boundary<std::string, double,int>{word_list[0],{start},val});
                }
            }
        }
    }while(!is_done);
    return work_str;
};

std::string read_dload(std::ifstream &file_stream, std::vector<Boundary<std::string,double,double>> &dlsetmap, std::vector<Boundary<int,double,double>> &dlele){
    std::string work_str;
    std::string word;
    std::vector<std::string> word_list;
    bool is_done = false;
    int tag;
    double vec[3];
    double val;
    do{
        getlmsg(file_stream,work_str);
        if(work_str.front() == '*'){
            is_done = true;
        }
        else{
            first_wd(work_str, word);
            word_list.clear();
            do{
                word_list.push_back(word);
                if(work_str.size() > 0){
                    first_wd(work_str,word);
                }
                else{
                    break;
                }
            }while(true);

            if(word_list.front().size() == 0){
                val = std::stod(word_list.at(2));
                for(int i = 0; i < 3; i++){
                    vec[i] = std::stod(word_list.at(i+3));
                }
                dlsetmap.push_back(Boundary<std::string,double,double>{word_list.at(0),{vec[0],vec[1],vec[2]},val});
            }
            else{
                if(word_list.front().front() < 58 && word_list.front().front() > 47){
                    if(word_list.size() == 6){
                        tag = std::stoi(word_list[0]);
                        if(!_strcmpi(word_list.at(1).data(), "GRAV")){
                            val = std::stod(word_list.at(2));
                            for(int i = 0; i < 3; i++){
                                vec[i] = std::stod(word_list.at(i+3));
                            }
                            dlele.push_back(Boundary<int,double,double>{tag,{vec[0],vec[1],vec[2]},val});
                        }
                    }
                    else{
                        printf("syntax error: parameter too much or less");
                    }
                }
                else{
                    if(word_list.size() == 6){
                        if(!_strcmpi(word_list.at(1).data(), "GRAV")){
                            val = std::stod(word_list.at(2));
                            for(int i = 0; i < 3; i++){
                                vec[i] = std::stod(word_list.at(i+3));
                            }
                            dlsetmap.push_back(Boundary<std::string,double,double>{word_list.at(0),{vec[0],vec[1],vec[2]},val});
                        }
                    }
                    else{
                        printf("syntax error: parameter too much or less");
                    }
                }
            }
        }
    }while(!is_done);
    return work_str;
};

std::string read_measure(std::ifstream &file_stream, std::vector<int> &pb_list, std::vector<double> &u_pb_val){
    std::string work_str;
    std::string word;
    std::vector<std::string> word_list;

    int tag, dim;
    double val;
    bool is_done = false;
    do{
        getlmsg(file_stream, work_str);
        if(work_str.front() == '*'){
            is_done = true;
        }
        else{
            first_wd(work_str,word);
            word_list.clear();
            do{
                word_list.push_back(word);
                if(work_str.size() > 0){
                    first_wd(work_str,word);
                }
                else{
                    break;
                }
            }while(true);
            tag = std::stoi(word_list.at(0));
            dim = std::stoi(word_list.at(1));
            val = std::stod(word_list.at(2));

            pb_list.push_back(3*(tag-1)+dim-1);
            u_pb_val.push_back(val);
        }

    }while(!is_done);
    return work_str;
};

bool getlmsg(std::ifstream &file_stream, std::string &str){
    std::getline(file_stream, str);
    while (true){
        if (!(str.length() == 0 || (str.at(0) == '*' && str.at(1) == '*'))){
            break;
        }
        else{
            if(file_stream.eof()){
                return true;
            }
            std::getline(file_stream, str);
        }
    }
    return false;
}

void first_wd(std::string &str, std::string &word){
    word.clear();
    std::vector<std::string::iterator> del_list;
    for(std::string::iterator i = str.begin(); i != str.end(); i++){
        if(*i == ','){
            del_list.push_back(i);
            break;
        }
        else if(*i == ' '){
            del_list.push_back(i);
        }
        else{
            word.push_back(*i);
            del_list.push_back(i);
        }
    }
    for (std::vector<std::string::iterator>::reverse_iterator i = del_list.rbegin(); i != del_list.rend(); i++){
        str.erase(*i);
    }
};