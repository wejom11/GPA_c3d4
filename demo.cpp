#include <stdio.h>
#include <mkl.h>
// #include <Eigen/Dense>
#include "P9SF_ele.h"
#include "assemble.h"
#include "mesh.h"

int main(){
	// Eigen::Matrix2d a;
	asb_manager asb;
    /*asb.init_ele();
    asb.getFout();
    asb.solve(asb.Fout);*/
    //asb.multi_solve_rho(100);
	std::ifstream inp_file("./Job-1.inp");
	asb.read_manager(inp_file);
	//asb.init_ele();
	asb.solve();
	asb.write();

   	// std::vector<double>ifw = { 4,2,3,3,4,4,4,3,4,4,4,4,4,4};
	//std::vector<double> ofw{ 4,4,4,2,4,4,4,4,4,4,4,4,4,4,4,4};
	//int nif = 3;
	//int nof = 4;
	//double fsinterval = 1;
	//Module module(fsinterval,ifw, ofw);
	//module.Init_geo();
	//module.Get_num_element_xy();
	//module.Get_nodes();
	//module.Get_elements();
	//module.Write_data();
	//module.Write_heatsourcedata();
}