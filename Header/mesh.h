#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <new>

#ifndef MESH_H
#define MESH_H

const int numoutfork = 4;   //number of outfork
const int numinfork=3;      //number of intfork

bool Isinregion(const std::vector<double>contour,double* point);

bool Isinando(const std::vector<double>contour1, const std::vector<double> contour2,double* point);

bool Iscollinear(double* p1, double* p2, double* p3);

bool Isonboundary(const std::vector<double>contour,double* point);

class XYnode {
public:
	double xy;
	double Dxy_para[31]{0.};
	XYnode(double XY, double* DXY_para){	
		xy = XY;
		for(int i=0;i<31;i++)Dxy_para[i] = DXY_para[i];
	}

	bool operator==(const XYnode& other) const {
		return xy == other.xy;
	}
};

bool Compare(const XYnode&a, const XYnode& b);

class Geo{
public:
	int numpara;
	int dim;
	double fsinterval;//The interval between fluid and solid
	std::vector<double> iffw{};
	std::vector<double> offw{};
	std::vector<double>bdp{};
	std::vector<double>ifibdp{};
	std::vector<double>ifobdp{};
	std::vector<double>ofibdp{};
	std::vector<double>ofobdp{};
	std::vector<double>ifcrack{};
	std::vector<double>ofcrack{};
	std::vector<double>Difi_para{};
	std::vector<double>Difo_para{};
	std::vector<double>Dofi_para{};
	std::vector<double>Dofo_para{};
	std::vector<double>Difcrack_para{};
	std::vector<double>Dofcrack_para{};
	std::vector<double>sobdp{};
	std::vector<double>Dsobdp_para{};
	std::vector<std::vector<std::vector<XYnode>>>xynodes{};//1-zones,2-dim,3-points(.xy,.Dxy_para)
	std::vector<double>pumpsize{ 4,4 };//pump width,pump height
	std::vector<double>heatsourcesize{ 20,5,5,10 };//hs1 width,hs1 height,hs2 width,hs2 height
	double pumpxy[2*2]{27,14,27,8};
	double heatsourcexy[2 * 2]{70,7.5,139.5,57};

	Geo() {};

	Geo(double fs, std::vector<double> ifw, std::vector<double> ofw){
		numpara = 31;
		dim =2;
		fsinterval = fs;
		iffw = ifw;
		offw = ofw;
		bdp.resize(static_cast<std::vector<double, std::allocator<double>>::size_type>(12) * 2);
		ifibdp.resize(static_cast<std::vector<double, std::allocator<double>>::size_type>(14) * 2);
		ifobdp.resize(static_cast<std::vector<double, std::allocator<double>>::size_type>(14) * 2);
		ofibdp.resize(static_cast<std::vector<double, std::allocator<double>>::size_type>(16) * 2);
		ofobdp.resize(static_cast<std::vector<double, std::allocator<double>>::size_type>(16) * 2);
		ifcrack.resize(static_cast<std::vector<double, std::allocator<double>>::size_type>((numinfork - 1) * 6) * 2);
		ofcrack.resize(static_cast<std::vector<double, std::allocator<double>>::size_type>((numoutfork - 1) * 6) * 2);
		Difcrack_para.resize(static_cast<std::vector<double, std::allocator<double>>::size_type>((numinfork - 1) * 6 * 2) * 31);
		Dofcrack_para.resize(static_cast<std::vector<double, std::allocator<double>>::size_type>((numoutfork - 1) * 6 * 2) * 31);
		sobdp.resize(static_cast<std::vector<double, std::allocator<double>>::size_type>((72 + (numinfork + numoutfork - 2) * 6)) * 2);
		Dsobdp_para.resize(static_cast<std::vector<double, std::allocator<double>>::size_type>((72 + (numinfork + numoutfork - 2) * 6) * 2) * 31);
		Difi_para.resize(static_cast<std::vector<double, std::allocator<double>>::size_type>(14 * 2) * 31);
		Difo_para.resize(static_cast<std::vector<double, std::allocator<double>>::size_type>(14 * 2) * 31);
		Dofi_para.resize(static_cast<std::vector<double, std::allocator<double>>::size_type>(16 * 2) * 31);
		Dofo_para.resize(static_cast<std::vector<double, std::allocator<double>>::size_type>(16 * 2) * 31);
		xynodes.resize(5);
		for (auto& dim_vector : xynodes) {
			dim_vector.resize(dim);
		}
	}
	
    void Init_bdp();

	void Init_ifibdp();
    
    void Get_Difi_para();
	
	void Init_ifobdp();
	
	void Get_Difo_para();

	void Init_ofibdp();
	
	void Get_Dofi_para();

	void Init_ofobdp();

	void Get_Dofo_para();

	void Init_ifcrack();

	void Get_Difcrack_para();
	
	void Init_ofcrack();
	
	void Get_Dofcrack_para();
	
	void Get_sobdp();

	void Init_xynodes();

	void Denser();
};

class Element{
public:
	int number;
	int nodes[9];
	int type = 0;         //0-fluid,1-solid;

	Element(int num,int* ns){
		number = num;
		for (int i = 0; i < 9; i++)nodes[i] = ns[i];
	}
};

class Node{
public:
	int number;
	double x;
	double y;
	double DXY_para[2 * 31];

	Node(int num,double xc,double yc,double* DXY_parac){
		number = num;
		x = xc;
		y = yc;
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 31; j++)
			{
				DXY_para[31*i+j] = DXY_parac[31*i+j];
			}
		}
	}
};

class Module{
public:
	Geo geo;
	std::vector<Element>elements{};
	std::vector<Node>nodes{};
	std::vector<int>number_nodes_T{};
	std::vector<int>number_nodes_P{};
	std::vector<int>number_element_hs{};
	int numnode = 0, numelement = 0;
	std::vector<int> num_node_zone{};
	std::vector<int> num_node_zone_xy{};
	std::vector<int> num_element_xy{};

	Module(double fs, std::vector<double> ifw, std::vector<double> ofw){
		geo = Geo(fs, ifw,ofw);
		num_element_xy.resize(static_cast<std::vector<int, std::allocator<int>>::size_type>(5) * 2);
		num_node_zone.resize(5, 0);
		num_node_zone_xy.resize(static_cast<std::vector<int, std::allocator<int>>::size_type>(5) * 2, 0);
	}

	void Init_geo();

	void Get_num_element_xy();

	void Get_nodes();

	int Get_number_node(int index, int jndex, int zonenum);

	bool Issolid(double* centerxy);

	void Get_elements();

	void Write_data();

	void Write_heatsourcedata();

	void Get_number_node_T();

};

#endif