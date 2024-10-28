#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <new>
#include "mesh.h"

bool Isinregion(const std::vector<double>contour,double* point){
	std::vector<double>cont = contour;
	cont.push_back(cont[0]); cont.push_back(cont[1]);
	bool result = false;
	double p1[2]{}, p2[2]{};
	double ycross = 0;
	double xcross = 0;
	int n = contour.size() / 2;
	for (int i = 0; i < n; i++)
	{
		p1[0] = cont[2 * i + 0]; p1[1] = cont[2 * i + 1];
		p2[0] = cont[2 * i + 2]; p2[1] = cont[2 * i + 3];
		double maxx = std::max(p1[0], p2[0]);
		double minx = std::min(p1[0], p2[0]);
		double maxy = std::max(p1[1], p2[1]);
		double miny = std::min(p1[1], p2[1]);
		if (point[0] == minx)continue;
		if (p2[0] == p1[0])continue;
		ycross = (p2[1] - p1[1]) * (point[0] - p1[0]) / (p2[0] - p1[0]) + p1[1];
		if (p2[1] == p1[1])xcross = point[0];
		else xcross = (p2[0] - p1[0]) * (point[1] - p1[1]) / (p2[1] - p1[1]) + p1[0];
		if(ycross>=miny&&ycross<= maxy&&xcross>=minx&&xcross<=maxx&&ycross>point[1])result = !result;
	}
	return result;
}

bool Isinando(const std::vector<double>contour1, const std::vector<double> contour2,double* point){
	//Determine whether the point is in the dual connected domain
	bool r1=Isinregion(contour1,point);
	bool r2 = Isinregion(contour2,point);
	return (r1 != r2);
}

bool Iscollinear(double* p1, double* p2, double* p3){
	double p1p2[2]{p1[0] - p2[0],p1[1] - p2[1]};
	double p1p3[2]{p1[0] - p3[0],p1[1] - p3[1]};
	if (p1p2[0] * p1p3[1] - p1p2[1] * p1p3[0] == 0)return true;
	else return false;
}

bool Isonboundary(const std::vector<double>contour,double* point){
	std::vector<std::vector<double>>cont(contour.size() / 2);
	int n = contour.size() / 2;
	for (int i = 0; i < n; i++)
	{
		cont[i][0] = contour[2 * i + 0]; cont[i][1] = contour[2 * i + 1];
	}
	cont.push_back(cont[0]);
	for (int i = 0; i < n; i++)
	{
		double p1[2]{ cont[i][0],cont[i][1]};
		double p2[2]{ cont[i+1][0],cont[i+1][1]};
		if (Iscollinear(p1 , p2, point))return true;
		else{
            return false;
        }
	}
    return false;
}

bool Compare(const XYnode&a, const XYnode& b){
	return a.xy < b.xy;
}

void Geo::Init_bdp(){
    bdp[2 * 0 + 0]=0; bdp[2 * 0 + 1] = 0;
    bdp[2 * 1 + 0] = 150; bdp[2 * 1 + 1]=0;
    bdp[2 * 2 + 0] = 150; bdp[2 * 2 + 1]=120;
    bdp[2 * 3 + 0] =130; bdp[2 * 3 + 1] = 120;
    bdp[2 * 4 + 0] =130; bdp[2 * 4 + 1] = 20;
    bdp[2 * 5 + 0] =110; bdp[2 * 5 + 1] = 20;
    bdp[2 * 6 + 0] = 110; bdp[2 * 6 + 1]=  70;
    bdp[2 * 7 + 0] = 90; bdp[2 * 7 + 1] = 70;
    bdp[2 * 8 + 0] = 90; bdp[2 * 8 + 1] = 20;
    bdp[2 * 9 + 0] = 70; bdp[2 * 9 + 1] = 20;
    bdp[2 * 10 + 0] = 70; bdp[2 * 10 + 1] = 100;
    bdp[2 * 11 + 0] = 0; bdp[2 * 11 + 1] = 100;
}

void Geo::Init_ifibdp(){
    //0
    ifibdp[2 * 0 + 0] = 2 * fsinterval + iffw[13] + offw[15]; ifibdp[2 * 0 + 1] = pumpxy[2 * 0 + 1] + pumpsize[1] / 2;
    //1
    ifibdp[2 * 1 + 0] = bdp[2 * 9 + 0] - (numinfork + numoutfork + 1) * fsinterval - numinfork * iffw[11] - numoutfork * offw[13]; ifibdp[2 * 1 + 1] = ifibdp[2*0+1];
    //2
    ifibdp[2 * 2 + 0] = ifibdp[2*1+0]; ifibdp[2 * 2 + 1] = iffw[2] + iffw[4] + 2 * fsinterval;
    //3
    ifibdp[2 * 3 + 0] = bdp[2 * 1 + 0] - 2 * fsinterval - iffw[3] - offw[5]; ifibdp[2 * 3 + 1] = ifibdp[2*2+1];
    //4
    ifibdp[2 * 4 + 0] = ifibdp[2 * 3 + 0]; ifibdp[2 * 4 + 1] = bdp[2 * 2 + 1] - 2 * fsinterval - iffw[4] - offw[6];
    //5
    ifibdp[2 * 5 + 0] = bdp[2 * 3 + 0] + 2 * fsinterval + iffw[5] + offw[7]; ifibdp[2 * 5 + 1] = ifibdp[2*4+1];
    //6
    ifibdp[2 * 6 + 0] = ifibdp[2 * 5 + 0]; ifibdp[2 * 6 + 1] = bdp[2 * 4 + 1] - 2 * fsinterval - iffw[6] - offw[8];
    //7
    ifibdp[2 * 7 + 0] = bdp[2 * 5 + 0] - 2 * fsinterval - iffw[7] - offw[9]; ifibdp[2 * 7 + 1] = ifibdp[2*6+1];
    //8
    ifibdp[2 * 8 + 0] = ifibdp[2 * 7 + 0]; ifibdp[2 * 8 + 1] = bdp[2 * 6 + 1] - 2 * fsinterval - iffw[8] - offw[10];
    //9
    ifibdp[2 * 9 + 0] = bdp[2 * 7 + 0] + 2 * fsinterval + iffw[9] + offw[11]; ifibdp[2 * 9 + 1] = ifibdp[2*8+1];
    //10
    ifibdp[2 * 10 + 0] = ifibdp[2 * 9 + 0]; ifibdp[2 * 10 + 1] = bdp[2 * 8 + 1] - 2 * fsinterval - iffw[10] - offw[12];
    //11
    ifibdp[2 * 11 + 0] = ifibdp[2 * 1 + 0] + fsinterval; ifibdp[2 * 11 + 1] = ifibdp[2*10+1];
    //12
    ifibdp[2 * 12 + 0] = ifibdp[2 * 11 + 0]; ifibdp[2 * 12 + 1] = bdp[2 * 10 + 1] - numinfork * iffw[12] - numoutfork * offw[14] - (numinfork + numoutfork) * fsinterval;
    //13
    ifibdp[2 * 13 + 0] = ifibdp[2 * 0 + 0]; ifibdp[2 * 13 + 1] = ifibdp[2*12+1];
}

void Geo::Get_Difi_para(){
    //0
    Difi_para[2*31*0+31*0+13] = 1; Difi_para[2 * 31 * 0 + 31 * 0 + 29] = 1; Difi_para[2 * 31 * 0 + 31 * 0 + 30] = 2;
    //1
    Difi_para[2*31*1+31*0+11] = -numinfork; Difi_para[2 * 31 * 1 + 31 * 0 + 27] = -numoutfork; Difi_para[2 * 31 * 1 + 31 * 0 + 30] = -(numinfork + numoutfork + 1);
    for (int i = 0; i < 31; i++)Difi_para[2 * 31 * 1 + 31 * 1 + i] = Difi_para[2 * 31 * 0 + 31 * 1 + i];
    //2
    for(int i=0;i<31;i++)Difi_para[2*31*2+31*0+i] = Difi_para[2*31*1+31*0+i];
    Difi_para[2 * 31 * 2+31*1+2] = 1; Difi_para[2 * 31 * 2+31*1+4] = 1; Difi_para[2 * 31 * 2+31*1+30] = 2;
    //3
    Difi_para[2*31*3+31*0+3] = -1; Difi_para[2*31*3+31*0+19] = -1; Difi_para[2 * 31 * 3+31*0+30] = -2;
    for (int i = 0; i < 31; i++)Difi_para[2 * 31 * 3 + 31 * 1 + i] = Difi_para[2 * 31 * 2 + 31 * 1 + i]; 
    //4
    for (int i = 0; i < 31; i++)Difi_para[2 * 31 * 4 + 31 * 0 + i] = Difi_para[2 * 31 * 3 + 31 * 0 + i];
    Difi_para[2 * 31 * 4+31*1+4] = -1; Difi_para[2 * 31 * 4+31*1+20] = -1; Difi_para[2 * 31 * 4+31*1+30] = -2;
    //5
    Difi_para[2 * 31 * 5+31*0+5] = 1; Difi_para[2 * 31 * 5+31*0+21] = 1; Difi_para[2 * 31 * 5+31*0+30] = 2;
    for (int i = 0; i < 31; i++)Difi_para[2 * 31 * 5 + 31 * 1 + i] = Difi_para[2 * 31 * 4 + 31 * 1 + i];
    //6
    for (int i = 0; i < 31; i++)Difi_para[2 * 31 * 6 + 31 * 0 + i] = Difi_para[2 * 31 * 5 + 31 * 0 + i];
    Difi_para[2 * 31 * 6+31*1+6] = -1; Difi_para[2 * 31 * 6+31*1+22] = -1; Difi_para[2 * 31 * 6+31*1+30] = -2;
    //7
    Difi_para[2 * 31 * 7+31*0+7] = -1; Difi_para[2 * 31 * 7+31*0+23] = -1; Difi_para[2 * 31 * 7+31*0+30] = -2;
    for (int i = 0; i < 31; i++)Difi_para[2 * 31 * 7 + 31 * 1 + i] = Difi_para[2 * 31 * 6 + 31 * 1 + i];
    //8
    for (int i = 0; i < 31; i++)Difi_para[2 * 31 * 8 + 31 * 0 + i] = Difi_para[2 * 31 * 7 + 31 * 0 + i];
    Difi_para[2 * 31 * 8+31*1+8] = -1; Difi_para[2 * 31 * 8+31*1+24] = -1; Difi_para[2 * 31 * 8+31*1+30] = -2;
    //9
    Difi_para[2 * 31 * 9+31*0+9] = 1; Difi_para[2 * 31 * 9+31*0+25] = 1; Difi_para[2 * 31 * 9+31*0+30] = 2; 
    for (int i = 0; i < 31; i++)Difi_para[2 * 31 * 9 + 31 * 1 + i] = Difi_para[2 * 31 * 8 + 31 * 1 + i];
    //10
    for (int i = 0; i < 31; i++)Difi_para[2 * 31 * 10 + 31 * 0 + i] = Difi_para[2 * 31 * 9 + 31 * 0 + i];
    Difi_para[2 * 31 * 10+31*1+10] = -1; Difi_para[2 * 31 * 10+31*1+26] = -1; Difi_para[2 * 31 * 10+31*1+30] = -2;
    //11
    for (int i = 0; i < 31; i++)Difi_para[2 * 31 * 11 + 31 * 0 + i] = Difi_para[2 * 31 * 1 + 31 * 0 + i];
    Difi_para[2 * 31 * 11+31*0+30] += 1;
    for (int i = 0; i < 31; i++)Difi_para[2 * 31 * 11 + 31 * 1 + i] = Difi_para[2 * 31 * 10 + 31 * 1 + i];
    //12
    for (int i = 0; i < 31; i++)Difi_para[2 * 31 * 12 + 31 * 0 + i] = Difi_para[2 * 31 * 11 + 31 * 0 + i];
    Difi_para[2 * 31 * 12+31*1+12] = -numinfork; Difi_para[2 * 31 * 12+31*1+28] = -numoutfork; Difi_para[2 * 31 * 12+31*1+30] = -(numinfork + numoutfork);
    //13
    for (int i = 0; i < 31; i++)Difi_para[2 * 31 * 13 + 31 * 0 + i] = Difi_para[2 * 31 * 0 + 31 * 0 + i];
    for (int i = 0; i < 31; i++)Difi_para[2 * 31 * 13 + 31 * 1 + i] = Difi_para[2 * 31 * 12 + 31 * 1 + i];
}

void Geo::Init_ifobdp(){    
    //0
    ifobdp[2 * 0 + 0] = ifibdp[2*0+0] - iffw[13]; ifobdp[2 * 0 + 1] = ifibdp[2*0+1] - iffw[0];
    //1
    ifobdp[2 * 1 + 0] = ifibdp[2*1+0] - iffw[1]; ifobdp[2 * 1 + 1] = ifobdp[2*0+1];
    //2
    ifobdp[2 * 2 + 0] = ifobdp[2*1+0]; ifobdp[2 * 2 + 1] = ifibdp[2*2+1] - iffw[2];
    //3
    ifobdp[2 * 3 + 0] = ifibdp[2*3+0] + iffw[3]; ifobdp[2 * 3 + 1] = ifobdp[2*2+1];
    //4
    ifobdp[2 * 4 + 0] = ifobdp[2*3+0]; ifobdp[2 * 4 + 1] = ifibdp[2*4+1] + iffw[4];
    //5
    ifobdp[2 * 5 + 0] = ifibdp[2*5+0] - iffw[5]; ifobdp[2 * 5 + 1] = ifobdp[2*4+1];
    //6
    ifobdp[2 * 6 + 0] = ifobdp[2*5+0]; ifobdp[2 * 6 + 1] = ifibdp[2*6+1] + iffw[6];
    //7
    ifobdp[2 * 7 + 0] = ifibdp[2*7+0] + iffw[7]; ifobdp[2 * 7 + 1] = ifobdp[2*6+1];
    //8
    ifobdp[2 * 8 + 0] = ifobdp[2*7+0]; ifobdp[2 * 8 + 1] = ifibdp[2*8+1] + iffw[8];
    //9
    ifobdp[2 * 9 + 0] = ifibdp[2*9+0] - iffw[9]; ifobdp[2 * 9 + 1] = ifobdp[2*8+1];
    //10
    ifobdp[2 * 10 + 0] = ifobdp[2*9+0]; ifobdp[2 * 10 + 1] = ifibdp[2*10+1] + iffw[10];
    //11
    ifobdp[2 * 11 + 0] = ifibdp[2*11+0] + numinfork * iffw[11] + (numinfork - 1) * fsinterval; ifobdp[2 * 11 + 1] = ifobdp[2*10+1];
    //12
    ifobdp[2 * 12 + 0] = ifobdp[2*11+0]; ifobdp[2 * 12 + 1] = ifibdp[2*12+1] + numinfork * iffw[12] + (numinfork - 1) * fsinterval;
    //13
    ifobdp[2 * 13 + 0] = ifibdp[2*13+0] - iffw[13]; ifobdp[2 * 13 + 1] = ifobdp[2*12+1];
}

void Geo::Get_Difo_para(){
    //0
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 0 + 31 * 0 + i] = Difi_para[2 * 31 * 0 + 31 * 0 + i];
    Difo_para[2*31*0+31*0+13] += -1;
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 0 + 31 * 1 + i] = Difi_para[2 * 31 * 0 + 31 * 1 + i];
    Difo_para[2 * 31 * 0+31*1+0] += -1;
    //1
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 1 + 31 * 0 + i] = Difi_para[2 * 31 * 1 + 31 * 0 + i];
    Difo_para[2 * 31 * 1+31 * 0+1] += -1; 
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 1 + 31 * 1 + i] = Difo_para[2 * 31 * 0 + 31 * 1 + i];
    //2
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 2 + 31 * 0 + i] = Difo_para[2 * 31 * 1 + 31 * 0 + i];
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 2 + 31 * 1 + i] = Difi_para[2 * 31 * 2 + 31 * 1 + i];
    Difo_para[2 * 31 * 2+31 * 1+2] += -1;
    //3
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 3 + 31 * 0 + i] = Difi_para[2 * 31 * 3 + 31 * 0 + i];
    Difo_para[2 * 31 * 3+31 * 0+3] += 1; 
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 3 + 31 * 1 + i] = Difo_para[2 * 31 * 2 + 31 * 1 + i];
    //4
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 4 + 31 * 0 + i] = Difo_para[2 * 31 * 3 + 31 * 0 + i];
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 4 + 31 * 1 + i] = Difi_para[2 * 31 * 4 + 31 * 1 + i];
    Difo_para[2 * 31 * 4+31 * 1+4] += 1;
    //5
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 5 + 31 * 0 + i] = Difi_para[2 * 31 * 5 + 31 * 0 + i];
    Difo_para[2 * 31 * 5+31 * 0+5] += -1;
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 5 + 31 * 1 + i] = Difo_para[2 * 31 * 4 + 31 * 1 + i];
    //6
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 6 + 31 * 0 + i] = Difo_para[2 * 31 * 5 + 31 * 0 + i];
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 6 + 31 * 1 + i] = Difi_para[2 * 31 * 6 + 31 * 1 + i];
    Difo_para[2 * 31 * 6+31 * 1+6] += 1;
    //7
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 7 + 31 * 0 + i] = Difi_para[2 * 31 * 7 + 31 * 0 + i];
    Difo_para[2 * 31 * 7+31 * 0+7] += 1;
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 7 + 31 * 1 + i] = Difo_para[2 * 31 * 6 + 31 * 1 + i];
    //8
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 8 + 31 * 0 + i] = Difo_para[2 * 31 * 7 + 31 * 0 + i];
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 8 + 31 * 1 + i] = Difi_para[2 * 31 * 8 + 31 * 1 + i];
    Difo_para[2 * 31 * 8+31 * 1+8] += 1;
    //9
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 9 + 31 * 0 + i] = Difi_para[2 * 31 * 9 + 31 * 0 + i];
    Difo_para[2 * 31 * 9+31 * 0+9] += -1;
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 9 + 31 * 1 + i] = Difo_para[2 * 31 * 8 + 31 * 1 + i];
    //10
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 10 + 31 * 0 + i] = Difo_para[2 * 31 * 9 + 31 * 0 + i];
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 10 + 31 * 1 + i] = Difi_para[2 * 31 * 10 + 31 * 1 + i];
    Difo_para[2 * 31 * 10+31 * 1+10] += 1;
    //11
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 11 + 31 * 0 + i] = Difi_para[2 * 31 * 11 + 31 * 0 + i];
    Difo_para[2 * 31 * 11+31 * 0+11] += numinfork;
    Difo_para[2 * 31 * 11+31 * 0+30] += numinfork - 1;
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 11 + 31 * 1 + i] = Difo_para[2 * 31 * 10 + 31 * 1 + i];
    //12
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 12 + 31 * 0 + i] = Difo_para[2 * 31 * 11 + 31 * 0 + i];
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 12 + 31 * 1 + i] = Difi_para[2 * 31 * 12 + 31 * 1 + i];
    Difo_para[2 * 31 * 12+31 * 1+12] += numinfork;
    Difo_para[2 * 31 * 12+31 * 1+30] += numinfork - 1;
    //13
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 13 + 31 * 0 + i] = Difi_para[2 * 31 * 13 + 31 * 0 + i];
    Difo_para[2 * 31 * 13+31*0+13] += -1; 
    for (int i = 0; i < 31; i++)Difo_para[2 * 31 * 13 + 31 * 1 + i] = Difo_para[2 * 31 * 12 + 31 * 1 + i];
}

void Geo::Init_ofibdp(){
    //0
    ofibdp[2 * 0 + 0] = fsinterval + offw[15]; ofibdp[2 * 0 + 1] = fsinterval + offw[0];
    //1
    ofibdp[2 * 1 + 0] = 10; ofibdp[2 * 1 + 1] = ofibdp[2*0+1];
    //2
    ofibdp[2 * 2 + 0] = ofibdp[2*1+0]; ofibdp[2 * 2 + 1] = 2 * fsinterval + offw[0] + offw[2];
    //3
    ofibdp[2 * 3 + 0] = ifobdp[2*1+0] - fsinterval; ofibdp[2 * 3 + 1] = ofibdp[2*2+1];
    //4
    ofibdp[2 * 4 + 0] = ofibdp[2*3+0]; ofibdp[2 * 4 + 1] = fsinterval + offw[4];
    //5
    ofibdp[2 * 5 + 0] = bdp[2 * 1 + 0] - fsinterval - offw[5]; ofibdp[2 * 5 + 1] = ofibdp[2*4+1];
    //6
    ofibdp[2 * 6 + 0] = ofibdp[2*5+0]; ofibdp[2 * 6 + 1] = bdp[2 * 2 + 1] - fsinterval - offw[6];
    //7
    ofibdp[2 * 7 + 0] = bdp[2 * 3 + 0] + fsinterval + offw[7]; ofibdp[2 * 7 + 1] = ofibdp[2*6+1];
    //8
    ofibdp[2 * 8 + 0] = ofibdp[2*7+0]; ofibdp[2 * 8 + 1] = bdp[2 * 4 + 1] - fsinterval - offw[8];
    //9
    ofibdp[2 * 9 + 0] = bdp[2 * 5 + 0] - fsinterval - offw[9]; ofibdp[2 * 9 + 1] = ofibdp[2*8+1];
    //10
    ofibdp[2 * 10 + 0] = ofibdp[2*9+0]; ofibdp[2 * 10 + 1] = bdp[2 * 6 + 1] - fsinterval - offw[10];
    //11
    ofibdp[2 * 11 + 0] = bdp[2 * 7 + 0] + fsinterval + offw[11]; ofibdp[2 * 11 + 1] = ofibdp[2*10+1];
    //12
    ofibdp[2 * 12 + 0] = ofibdp[2*11+0]; ofibdp[2 * 12 + 1] = bdp[2 * 8 + 1] - fsinterval - offw[12];
    //13
    ofibdp[2 * 13 + 0] = bdp[2 * 9 + 0] - numoutfork * (fsinterval + offw[13]); ofibdp[2 * 13 + 1] = ofibdp[2*12+1];
    //14
    ofibdp[2 * 14 + 0] = ofibdp[2*13+0]; ofibdp[2 * 14 + 1] = bdp[2 * 10 + 1] - numoutfork * (fsinterval + offw[14]);
    //15
    ofibdp[2 * 15 + 0] = bdp[2 * 11 + 0] + fsinterval + offw[15]; ofibdp[2 * 15 + 1] = ofibdp[2*14+1];
}

void Geo::Get_Dofi_para(){
    //0
    Dofi_para[2*31*0+31*0+29] = 1; Dofi_para[2 * 31 * 0+31 * 0+30] = 1; Dofi_para[2 * 31 * 0+31 * 1+14] = 1; Dofi_para[2 * 31 * 0+31 * 1+30] = 1;
    //1
    for (int i = 0; i < 31; i++)Dofi_para[2 * 31 * 1 + 31 * 1 + i] = Dofi_para[2 * 31 * 0 + 31 * 1 + i];
    //2
    for (int i = 0; i < 31; i++)Dofi_para[2 * 31 * 2 + 31 * 0 + i] = Dofi_para[2 * 31 * 1 + 31 * 0 + i];
    Dofi_para[2 * 31 * 2+31 * 1+14] = 1;Dofi_para[2 * 31 * 2+31 * 1+16] = 1;Dofi_para[2 * 31 * 2+31 * 1+30] = 2;
    //3
    for (int i = 0; i < 31; i++)Dofi_para[2 * 31 * 3 + 31 * 0 + i] = Difo_para[2 * 31 * 1 + 31 * 0 + i];
    Dofi_para[2 * 31 * 3+31 * 0+30] += -1;
    for (int i = 0; i < 31; i++)Dofi_para[2 * 31 * 3 + 31 * 1 + i] = Dofi_para[2 * 31 * 2 + 31 * 1 + i];
    //4
    for (int i = 0; i < 31; i++)Dofi_para[2 * 31 * 4 + 31 * 0 + i] = Dofi_para[2 * 31 * 3 + 31 * 0 + i];
    Dofi_para[2 * 31 * 4+31 * 1+18] = 1; 
    Dofi_para[2 * 31 * 4+31 * 1+30] = 1;
    //5
    Dofi_para[2 * 31 * 5+31 * 0+19] = -1;
    Dofi_para[2 * 31 * 5+31 * 0+30] = -1;
    for (int i = 0; i < 31; i++)Dofi_para[2 * 31 * 5 + 31 * 1 + i] = Dofi_para[2 * 31 * 4 + 31 * 1 + i];
    //6
    for (int i = 0; i < 31; i++)Dofi_para[2 * 31 * 6 + 31 * 0 + i] = Dofi_para[2 * 31 * 5 + 31 * 0 + i];
    Dofi_para[2 * 31 * 6+31 * 1+20] = -1; 
    Dofi_para[2 * 31 * 6+31 * 1+30] = -1;
    //7
    Dofi_para[2 * 31 * 7+31 * 0+21] = 1;
    Dofi_para[2 * 31 * 7+31 * 0+30] = 1;
    for (int i = 0; i < 31; i++)Dofi_para[2 * 31 * 7 + 31 * 1 + i] = Dofi_para[2 * 31 * 6 + 31 * 1 + i];
    //8
    for (int i = 0; i < 31; i++)Dofi_para[2 * 31 * 8 + 31 * 0 + i] = Dofi_para[2 * 31 * 7 + 31 * 0 + i];
    Dofi_para[2 * 31 * 8+31 * 1+22] = -1;
    Dofi_para[2 * 31 * 8+31 * 1+30] = -1;
    //9
    Dofi_para[2 * 31 * 9+31 * 0+23] = -1;
    Dofi_para[2 * 31 * 9+31 * 0+30] = -1;
    for (int i = 0; i < 31; i++)Dofi_para[2 * 31 * 9 + 31 * 1 + i] = Dofi_para[2 * 31 * 8 + 31 * 1 + i];
    //10
    for (int i = 0; i < 31; i++)Dofi_para[2 * 31 * 10 + 31 * 0 + i] = Dofi_para[2 * 31 * 9 + 31 * 0 + i];
    Dofi_para[2 * 31 * 10+31 * 1+24] = -1;
    Dofi_para[2 * 31 * 10+31 * 1+30] = -1;
    //11
    Dofi_para[2 * 31 * 11+31 * 0+25] = 1;
    Dofi_para[2 * 31 * 11+31 * 0+30] = 1; 
    for (int i = 0; i < 31; i++)Dofi_para[2 * 31 * 11 + 31 * 1 + i] = Dofi_para[2 * 31 * 10 + 31 * 1 + i];
    //12
    for (int i = 0; i < 31; i++)Dofi_para[2 * 31 * 12 + 31 * 0 + i] = Dofi_para[2 * 31 * 11 + 31 * 0 + i];
    Dofi_para[2 * 31 * 12+31 * 1+26] = -1;
    Dofi_para[2 * 31 * 12+31 * 1+30] = -1;
    //13
    Dofi_para[2 * 31 * 13+31 * 0+27] = -numoutfork;
    Dofi_para[2 * 31 * 13+31 * 0+30] = -numoutfork; 
    for (int i = 0; i < 31; i++)Dofi_para[2 * 31 * 13 + 31 * 1 + i] = Dofi_para[2 * 31 * 12 + 31 * 1 + i];
    //14
    for (int i = 0; i < 31; i++)Dofi_para[2 * 31 * 14 + 31 * 0 + i] = Dofi_para[2 * 31 * 13 + 31 * 0 + i];
    Dofi_para[2 * 31 * 14+31 * 1+28] = -numoutfork; 
    Dofi_para[2*31*14+31 * 1+30] = -numoutfork;
    //15
    Dofi_para[2 * 31 * 15+31 * 0+29] = 1;
    Dofi_para[2 * 31 * 15+31 * 0+30] = 1; 
    for (int i = 0; i < 31; i++)Dofi_para[2 * 31 * 15 + 31 * 1 + i] = Dofi_para[2 * 31 * 14 + 31 * 1 + i];
}

void Geo::Init_ofobdp(){
    //0
    ofobdp[2 * 0 + 0] = ofibdp[2*0+0] - offw[15]; ofobdp[2 * 0 + 1] = ofibdp[2*0+1] - offw[0];
    //1
    ofobdp[2 * 1 + 0] = ofibdp[2*1+0] + offw[1]; ofobdp[2 * 1 + 1] = ofobdp[2*0+1];
    //2
    ofobdp[2 * 2 + 0] = ofobdp[2*1+0]; ofobdp[2 * 2 + 1] = ofibdp[2*2+1] - offw[2];
    //3
    ofobdp[2 * 3 + 0] = ofibdp[2*3+0] - offw[3]; ofobdp[2 * 3 + 1] = ofobdp[2*2+1];
    //4
    ofobdp[2 * 4 + 0] = ofobdp[2*3+0]; ofobdp[2 * 4 + 1] = ofibdp[2*4+1] - offw[4];
    //5
    ofobdp[2 * 5 + 0] = ofibdp[2*5+0] + offw[5]; ofobdp[2 * 5 + 1] = ofobdp[2*4+1];
    //6
    ofobdp[2 * 6 + 0] = ofobdp[2*5+0]; ofobdp[2 * 6 + 1] = ofibdp[2*6+1] + offw[6];
    //7
    ofobdp[2 * 7 + 0] = ofibdp[2*7+0] - offw[7]; ofobdp[2 * 7 + 1] = ofobdp[2*6+1];
    //8
    ofobdp[2 * 8 + 0] = ofobdp[2*7+0]; ofobdp[2 * 8 + 1] = ofibdp[2*8+1] + offw[8];
    //9
    ofobdp[2 * 9 + 0] = ofibdp[2*9+0] + offw[9]; ofobdp[2 * 9 + 1] = ofobdp[2*8+1];
    //10
    ofobdp[2 * 10 + 0] = ofobdp[2*9+0]; ofobdp[2 * 10 + 1] = ofibdp[2*10+1] + offw[10];
    //11
    ofobdp[2 * 11 + 0] = ofibdp[2*11+0] - offw[11]; ofobdp[2 * 11 + 1] = ofobdp[2*10+1];
    //12
    ofobdp[2 * 12 + 0] = ofobdp[2*11+0]; ofobdp[2 * 12 + 1] = ofibdp[2*12+1] + offw[12];
    //13
    ofobdp[2 * 13 + 0] = ofibdp[2*13+0] + numoutfork * offw[13] + (numoutfork - 1) * fsinterval; ofobdp[2 * 13 + 1] = ofobdp[2*12+1];
    //14
    ofobdp[2 * 14 + 0] = ofobdp[2*13+0]; ofobdp[2 * 14 + 1] = ofibdp[2*14+1] + numoutfork * offw[14] + (numoutfork - 1) * fsinterval;
    //15
    ofobdp[2 * 15 + 0] = ofibdp[2*15+0] - offw[15]; ofobdp[2 * 15 + 1] = ofobdp[2*14+1];
}

void Geo::Get_Dofo_para(){
    //0
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 0 + 31 * 0 + i] = Dofi_para[2 * 31 * 0 + 31 * 0 + i];
    Dofo_para[2*31*0+31*0+29] += -1; 
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 0 + 31 * 1 + i] = Dofi_para[2 * 31 * 0 + 31 * 1 + i];
    Dofo_para[2 * 31 * 0+31 * 1+14] += -1;
    //1
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 1 + 31 * 0 + i] = Dofi_para[2 * 31 * 1 + 31 * 0 + i];
    Dofo_para[2 * 31 * 1+31 * 0+15] += 1;
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 1 + 31 * 1 + i] = Dofo_para[2 * 31 * 0 + 31 * 1 + i];
    //2
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 2 + 31 * 0 + i] = Dofo_para[2 * 31 * 1 + 31 * 0 + i];
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 2 + 31 * 1 + i] = Dofi_para[2 * 31 * 2 + 31 * 1 + i];
    Dofo_para[2 * 31 * 2+31 * 1+16] += -1;
    //3
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 3 + 31 * 0 + i] = Dofi_para[2 * 31 * 3 + 31 * 0 + i];
    Dofo_para[2 * 31 * 3+31 * 0+17] += -1;
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 3 + 31 * 1 + i] = Dofo_para[2 * 31 * 2 + 31 * 1 + i]; 
    //4
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 4 + 31 * 0 + i] = Dofo_para[2 * 31 * 3 + 31 * 0 + i];
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 4 + 31 * 1 + i] = Dofi_para[2 * 31 * 4 + 31 * 1 + i];
    Dofo_para[2 * 31 * 4+31 * 1+18] += -1;
    //5
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 5 + 31 * 0 + i] = Dofi_para[2 * 31 * 5 + 31 * 0 + i];
    Dofo_para[2 * 31 * 5+31 * 0+19] += 1;
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 5 + 31 * 1 + i] = Dofo_para[2 * 31 * 4 + 31 * 1 + i];
    //6
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 6 + 31 * 0 + i] = Dofo_para[2 * 31 * 5 + 31 * 0 + i];
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 6 + 31 * 1 + i] = Dofi_para[2 * 31 * 6 + 31 * 1 + i];
    Dofo_para[2 * 31 * 6+31 * 1+20] += 1;
    //7
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 7 + 31 * 0 + i] = Dofi_para[2 * 31 * 7 + 31 * 0 + i];
    Dofo_para[2 * 31 * 7+31 * 0+21] += -1;
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 7 + 31 * 1 + i] = Dofo_para[2 * 31 * 6 + 31 * 1 + i];
    //8
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 8 + 31 * 0 + i] = Dofo_para[2 * 31 * 7 + 31 * 0 + i];
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 8 + 31 * 1 + i] = Dofi_para[2 * 31 * 8 + 31 * 1 + i];
    Dofo_para[2 * 31 * 8+31 * 1+22] += 1;
    //9
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 9 + 31 * 0 + i] = Dofi_para[2 * 31 * 9 + 31 * 0 + i];
    Dofo_para[2 * 31 * 9+31 * 0+23] += 1;
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 9 + 31 * 1 + i] = Dofo_para[2 * 31 * 8 + 31 * 1 + i];
    //10
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 10 + 31 * 0 + i] = Dofo_para[2 * 31 * 9 + 31 * 0 + i];
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 10 + 31 * 1 + i] = Dofi_para[2 * 31 * 10 + 31 * 1 + i];
    Dofo_para[2 * 31 * 10+31 * 1+24] += 1;
    //11
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 11 + 31 * 0 + i] = Dofi_para[2 * 31 * 11 + 31 * 0 + i];
    Dofo_para[2 * 31 * 11+31 * 0+25] += -1;
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 11 + 31 * 1 + i] = Dofo_para[2 * 31 * 10 + 31 * 1 + i];
    //12
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 12 + 31 * 0 + i] = Dofo_para[2 * 31 * 11 + 31 * 0 + i];
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 12 + 31 * 1 + i] = Dofi_para[2 * 31 * 12 + 31 * 1 + i];
    Dofo_para[2 * 31 * 12+31 * 1+26] += 1;
    //13
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 13 + 31 * 0 + i] = Dofi_para[2 * 31 * 13 + 31 * 0 + i];
    Dofo_para[2 * 31 * 13+31 * 0+27] += numoutfork;
    Dofo_para[2 * 31 * 13+31 * 0+30] += numoutfork - 1;
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 13 + 31 * 1 + i] = Dofo_para[2 * 31 * 12 + 31 * 1 + i];
    //14
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 14 + 31 * 0 + i] = Dofo_para[2 * 31 * 13 + 31 * 0 + i];
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 14 + 31 * 1 + i] = Dofi_para[2 * 31 * 14 + 31 * 1 + i];
    Dofo_para[2 * 31 * 14+31 * 1+28] += numoutfork;
    Dofo_para[2 * 31 * 14+31 * 1+30] += numoutfork - 1;
    //15
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 15 + 31 * 0 + i] = Dofi_para[2 * 31 * 15 + 31 * 0 + i];
    Dofo_para[2 * 31 * 15+31 * 0+29] += -1;
    for (int i = 0; i < 31; i++)Dofo_para[2 * 31 * 15 + 31 * 1 + i] = Dofo_para[2 * 31 * 14 + 31 * 1 + i];
}

void Geo::Init_ifcrack(){
    for (int i=0;i<numinfork-1;i++)
    {
        //0
        ifcrack[6 * 2 * i + 2 * 0 + 0] = 2 * fsinterval + iffw[13] + offw[15]; ifcrack[6 * 2 * i + 2 * 0 + 1] = bdp[2 * 11 + 1] - (numinfork + numoutfork - i) * fsinterval - (numinfork - i - 1) * iffw[12] - numoutfork * offw[14];
        //1
        ifcrack[6 * 2 * i + 2 * 1 + 0] = bdp[2 * 10 + 0] - (numinfork + numoutfork - i) * fsinterval - (numinfork - i - 1) * iffw[11] - numoutfork * offw[13]; ifcrack[6 * 2 * i + 2 * 1 + 1]= ifcrack[6*2*i+2*0+1];
        //2
        ifcrack[6 * 2 * i + 2 * 2 + 0] = ifcrack[6*2*i+2*1+0]; ifcrack[6 * 2 * i + 2 * 2 + 1]= ifibdp[2 * 11 + 1] + iffw[10];
        //3
        ifcrack[6 * 2 * i + 2 * 3 + 0]= ifcrack[6 * 2 * i+2*2+0] + fsinterval; ifcrack[6 * 2 * i + 2 * 3 + 1]= ifcrack[6 * 2 * i+2*2+1];
        //4
        ifcrack[6 * 2 * i + 2 * 4 + 0]= ifcrack[6 * 2 * i+2*1+0] + fsinterval; ifcrack[6 * 2 * i + 2 * 4 + 1]= ifcrack[6 * 2 * i+2*1+1] + fsinterval;
        //5
        ifcrack[6 * 2 * i + 2 * 5 + 0]= ifcrack[6 * 2 * i+2*0+0]; ifcrack[6 * 2 * i + 2 * 5 + 1]= ifcrack[6 * 2 * i+2*0+1] + fsinterval;
    }
}

void Geo::Get_Difcrack_para(){
    for (int i = 0; i < numinfork - 1; i++)
    {
        //0
        Difcrack_para[6 * 2 * 31 * i+2 * 31 * 0+31 * 0+13] = 1;Difcrack_para[6 * 2 * 31 * i+2 * 31 * 0+31 * 0+29] = 1;Difcrack_para[6 * 2 * 31 * i+2 * 31 * 0+31 * 0+30] = 2;
        Difcrack_para[6 * 2 * 31 * i+2 * 31 * 0+31 * 1+12] = -(numinfork - i - 1);Difcrack_para[6 * 2 * 31 * i+2 * 31 * 0+31 * 1+28] = -numoutfork;Difcrack_para[6 * 2 * 31 * i+2 * 31 * 0+31 * 1+30] = -(numinfork + numoutfork - i);
        //1
        Difcrack_para[6 * 2 * 31 * i + 2 * 31 * 1 + 31 * 0 + 11] = -(numinfork - i - 1);
        Difcrack_para[6 * 2 * 31 * i+2 * 31 * 1+31 * 0+27] = -numoutfork;
        Difcrack_para[6 * 2 * 31 * i+2 * 31 * 1+31 * 0+30] = -(numinfork + numoutfork - i);
        for (int j = 0; j < 31; j++)Difcrack_para[6 * 2 * 31 * i + 2 * 31 * 1 + 31 * 1+j] = Difcrack_para[6 * 2 * 31 * i + 2 * 31 * 0+31 * 1+j];
        //2
        for (int j = 0; j < 31; j++)Difcrack_para[6 * 2 * 31 * i + 2 * 31 * 2 + 31 * 0+j] = Difcrack_para[6 * 2 * 31 * i + 2 * 31 * 1 + 31 * 0+j];
        for (int j = 0; j < 31; j++)Difcrack_para[6 * 2 * 31 * i + 2 * 31 * 2 + 31 * 1+j] = Difi_para[2*31*11+31*1+j];
        Difcrack_para[6 * 2 * 31 * i+2 * 31 * 2+31 * 1+10] += 1;
        //3
        for (int j = 0; j < 31; j++)Difcrack_para[6 * 2 * 31 * i + 2 * 31 * 3 + 31 * 0+j] = Difcrack_para[6 * 2 * 31 * i+2 * 31 * 2+31 * 0+j];
        Difcrack_para[6 * 2 * 31 * i+2 * 31 * 3+31 * 0+30] += 1;
        for (int j = 0; j < 31; j++)Difcrack_para[6 * 2 * 31 * i + 2 * 31 * 3 + 31 * 1+j] = Difcrack_para[6 * 2 * 31 * i+2 * 31 * 2+31 * 1+j];
        //4
        for (int j = 0; j < 31; j++)Difcrack_para[6 * 2 * 31 * i + 2 * 31 * 4 + 31 * 0+j] = Difcrack_para[6 * 2 * 31 * i+2 * 31 * 1+31 * 0+j];
        Difcrack_para[6 * 2 * 31 * i+2 * 31 * 4+31 * 0+30] += 1;
        for (int j = 0; j < 31; j++)Difcrack_para[6 * 2 * 31 * i + 2 * 31 * 4 + 31 * 1+j] = Difcrack_para[6 * 2 * 31 * i+2 * 31 * 1+31 * 1+j];
        Difcrack_para[6 * 2 * 31 * i+2 * 31 * 4+31 * 1+30] += 1;
        //5
        for (int j = 0; j < 31; j++)Difcrack_para[6 * 2 * 31 * i + 2 * 31 * 5 + 31 * 0+j] = Difcrack_para[6 * 2 * 31 * i+2 * 31 * 0+31 * 0+j];
        for(int j=0;j<31;j++)Difcrack_para[6 * 2 * 31 * i + 2 * 31 * 5 + 31 * 1+j] = Difcrack_para[6 * 2 * 31 * i+2 * 31 * 0+31 * 1+j];
        Difcrack_para[6 * 2 * 31 * i+2 * 31 * 5+31 * 1+30] += 1;
    }
}

void Geo::Init_ofcrack(){
    for (int i = 0; i < numoutfork - 1; i++)
    {
        //0
        ofcrack[6 * 2 * i + 2 * 0 + 0] = fsinterval + offw[13]; ofcrack[6 * 2 * i + 2 * 0 + 1] = bdp[2 * 11 + 1] - (numoutfork - i) * fsinterval - (numoutfork - 1 - i) * offw[12];
        //1
        ofcrack[6 * 2 * i + 2 * 1 + 0] = bdp[2 * 10 + 0] - (numoutfork - i) * fsinterval - (numoutfork - 1 - i) * offw[11]; ofcrack[6 * 2 * i + 2 * 1 + 1] = ofcrack[6*2*i+2*0+1];
        //2
        ofcrack[6 * 2 * i + 2 * 2 + 0] = ofcrack[6 * 2 * i+2*1+0]; ofcrack[6 * 2 * i + 2 * 2 + 1] = ofibdp[2*13+1] + offw[12];
        //3
        ofcrack[6 * 2 * i + 2 * 3 + 0] = ofcrack[6 * 2 * i+2*2+0] + fsinterval; ofcrack[6 * 2 * i + 2 * 3 + 1] = ofcrack[6 * 2 * i+2*2+1];
        //4
        ofcrack[6 * 2 * i + 2 * 4 + 0] = ofcrack[6 * 2 * i+2*1+0] + fsinterval; ofcrack[6 * 2 * i + 2 * 4 + 1] = ofcrack[6 * 2 * i+2*1+1] + fsinterval;
        //5
        ofcrack[6 * 2 * i + 2 * 5 + 0] = ofcrack[6 * 2 * i+2*0+0]; ofcrack[6 * 2 * i + 2 * 5 + 1] = ofcrack[6 * 2 * i+2*0+1] + fsinterval;
    }
}

void Geo::Get_Dofcrack_para(){
    for (int i = 0; i < numoutfork - 1; i++)
    {
        //0
        Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 0+31 * 0+27] = 1;Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 0+31 * 0+30] = 1;
        Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 0+31 * 1+26]=-(numoutfork-1-i);Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 0+31 * 1+30] = -(numoutfork - i);
        //1
        Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 1+31 * 0+25] = -(numoutfork - 1 - i);Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 1+31 * 0+30] = -(numoutfork - i);
        for (int j = 0; j < 31; j++)Dofcrack_para[6 * 2 * 31 * i + 2 * 31 * 1 + 31 * 1+j] = Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 0+31 * 1+j];
        //2
        for (int j = 0; j < 31; j++)Dofcrack_para[6 * 2 * 31 * i + 2 * 31 * 2 + 31 * 0+j] = Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 1+31 * 0+j];
        Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 2+31 * 1]= Dofi_para[2*31*13+31*1];
        Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 2+31 * 1+26] += 1;
        //3
        for (int j = 0; j < 31; j++)Dofcrack_para[6 * 2 * 31 * i + 2 * 31 * 3 + 31 * 0+j] = Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 2+31 * 0+j];
        Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 3+31 * 0+30] += 1;
        for (int j = 0; j < 31; j++)Dofcrack_para[6 * 2 * 31 * i + 2 * 31 * 3 + 31 * 1+j] = Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 2+31 * 1+j];
        //4
        for (int j = 0; j < 31; j++)Dofcrack_para[6 * 2 * 31 * i + 2 * 31 * 4 + 31 * 0 + j] = Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 1+31 * 0+j];
        Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 4+31 * 0+30] += 1;
        for (int j = 0; j < 31; j++)Dofcrack_para[6 * 2 * 31 * i + 2 * 31 * 4 + 31 * 1+j] = Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 1+31 * 1+j];
        Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 4+31 * 1+30] += 1;
        //5
        for (int j = 0; j < 31; j++)Dofcrack_para[6 * 2 * 31 * i + 2 * 31 * 5 + 31 * 0+j] = Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 0+31 * 0+j];
        for (int j = 0; j < 31; j++)Dofcrack_para[6 * 2 * 31 * i + 2 * 31 * 5 + 31 * 1+j] = Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 0+31 * 1+j];
        Dofcrack_para[6 * 2 * 31 * i+2 * 31 * 5+31 * 1+30] += 1;
    }
}

void Geo::Get_sobdp(){
    for (int i = 0; i < 12; i++)
    {
        sobdp[2 * i + 0] = bdp[2 * i + 0]; sobdp[2 * i + 1] = bdp[2 * i + 1];
    }
    for (int i = 0; i < 14; i++)
    { 
        for (int j = 0; j < 2; j++)
        {
            sobdp[24 + 2 * i + j] = ifibdp[2 * i + j];
            for (int k = 0; k < 31; k++)
            {
                Dsobdp_para[2 * 31 * (12 + i) + 31 * j + k] = Difi_para[2 * 31 * i + 31 * j + k];
            }
        }
    }
    for (int i = 0; i < 14; i++)
    { 
        for (int j = 0; j < 2; j++)
        {
            sobdp[52 + 2 * i + j] = ifobdp[2 * i + j];
            for (int k = 0; k < 31; k++)
            {
                Dsobdp_para[2 * 31 * (26 + i) + 31 * j + k] = Difo_para[2 * 31 * i + 31 * j + k];
            }
        }
    }
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            sobdp[80 + 2 * i + j] = ofibdp[2 * i + j];
            for (int k = 0; k < 31; k++)
            {
                Dsobdp_para[2 * 31 * (40 + i) + 31 * j + k] = Dofi_para[2 * 31 * i + 31 * j + k];
            }
        }
    }
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            sobdp[112 + 2 * i + j] = ofobdp[2 * i + j];
            for (int k = 0; k < 31; k++)
            {
                Dsobdp_para[2 * 31 * (56 + i) + 31 * j + k] = Dofo_para[2 * 31 * i + 31 * j + k];
            }
        }
    }
    for (int i = 0; i < numinfork - 1; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                sobdp[144 + 12 * i + 2 * j + k] = ifcrack[6 * 2 * i + 2 * j + k];
                for (int m = 0; m < 31; m++)
                {
                    Dsobdp_para[2 * 31 * (72 + 6 * i + j) + 31 * k + m] = Difcrack_para[6 * 2 * 31 * i + 2 * 31 * j + 31 * k + m];
                }
            }
        }
    }
    for (int i = 0; i < numoutfork - 1; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                sobdp[144 + 12 * (numinfork - 1) + 12 * i + 2 * j + k] = ofcrack[6 * 2 * i + 2 * j + k];
                for (int m = 0; m < 31; m++)
                {
                    Dsobdp_para[2 * 31 * (72 + 6 * (numinfork - 1 + i) + j) + 31 * k + m] = Dofcrack_para[6 * 2 * 31 * i + 2 * 31 * j + 31 * k + m];
                }
            }
        }
    }
}

void Geo::Init_xynodes(){
    std::vector<double>ip((72 + (numinfork + numoutfork - 2) * 6 + 3) * 2,0);
    std::vector<double>Dip_para((72 + (numinfork + numoutfork - 2) * 6 + 3) * 2 * 31,0);
    int n = 72 + (numinfork + numoutfork - 2) * 6;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            ip[2 * i + j] = sobdp[2 * i + j];
            for (int k = 0; k < 31; k++)
            {
                Dip_para[2 * 31 * i + 31 * j + k] = Dsobdp_para[2 * 31 * i + 31 * j + k];
            }
        }
    }
    ip[2 * (n + 0) + 0] = heatsourcexy[2 * 0 + 0] + heatsourcesize[0] / 2; ip[2 * (n + 0) + 1] = heatsourcexy[2 * 0 + 1] + heatsourcesize[1] / 2;
    ip[2 * (n + 1) + 0] = heatsourcexy[2 * 1 + 0] + heatsourcesize[2] / 2; ip[2 * (n + 1) + 1] = heatsourcexy[2 * 1 + 1] - heatsourcesize[3] / 2;
    ip[2 * (n + 2) + 0] = heatsourcexy[2 * 1 + 0] - heatsourcesize[2] / 2; ip[2 * (n + 2) + 1] = heatsourcexy[2 * 1 + 1] + heatsourcesize[3] / 2;
    n += 3;
    double temp[31]{};
    for (int i = 0; i < n; i++)
    {
        if (ip[2*i+0] <= bdp[2*9+0])
        {
            for (int j = 0; j < 31; j++)temp[j] = Dip_para[2 * 31 * i + 31 * 0 + j];
            xynodes[0][0].push_back(XYnode(ip[2*i+0],temp));
            for (int j = 0; j < 31; j++)temp[j] = Dip_para[2 * 31 * i + 31 * 1 + j];
            xynodes[0][1].push_back(XYnode(ip[2*i+1],temp));
            if (ip[2 * i + 1] <= bdp[2 * 9 + 1])
            {
                for (int j = 1; j < 5; j++)
                {
                    for (int k = 0; k < 31; k++)temp[k] = Dip_para[2 * 31 * i + 31 * 1 + k];
                    xynodes[j][1].push_back(XYnode(ip[2 * i + 1], temp));
                }
            }
        }
        if (ip[2*i+0] >= bdp[2*9+0] && ip[2*i+0] <= bdp[2*8+0]&& ip[2*i+1] <= bdp[2*9+1] && ip[2*i+1] <= bdp[2*8+1])
        {
            for (int j = 0; j < 31; j++)temp[j] = Dip_para[2 * 31 * i + 31 * 0 + j];
            xynodes[1][0].push_back(XYnode(ip[2*i+0], temp));
            for (int j = 0; j < 31; j++)temp[j] = Dip_para[2 * 31 * i + 31 * 1 + j];
            xynodes[1][1].push_back(XYnode(ip[2*i+1], temp));
        }
        if (ip[2*i+0] >= bdp[2*8+0] && ip[2*i+0] <= bdp[2*5+0])
        {
            for (int j = 0; j < 31; j++)temp[j] = Dip_para[2 * 31 * i + 31 * 0 + j];
            xynodes[2][0].push_back(XYnode(ip[2*i+0], temp));
            for (int j = 0; j < 31; j++)temp[j] = Dip_para[2 * 31 * i + 31 * 1 + j];
            xynodes[2][1].push_back(XYnode(ip[2*i+1], temp));
        }
        if (ip[2*i+0] >= bdp[2*5+0] && ip[2*i+0] <= bdp[2*4+0]&& ip[2*i+1] <= bdp[2*5+1] && ip[2*i+1] <= bdp[2*4+1])
        {
            for (int j = 0; j < 31; j++)temp[j] = Dip_para[2 * 31 * i + 31 * 0 + j];
            xynodes[3][0].push_back(XYnode(ip[2*i+0],temp));
            for (int j = 0; j < 31; j++)temp[j] = Dip_para[2 * 31 * i + 31 * 1 + j];
            xynodes[3][1].push_back(XYnode(ip[2*i+1], temp));
        }
        if(ip[2 * i + 0] >= bdp[2 * 4 + 0])
        {
            for (int j = 0; j < 31; j++)temp[j] = Dip_para[2 * 31 * i + 31 * 0 + j];
            xynodes[4][0].push_back(XYnode(ip[2*i+0], temp));
            for (int j = 0; j < 31; j++)temp[j] = Dip_para[2 * 31 * i + 31 * 1 + j];
            xynodes[4][1].push_back(XYnode(ip[2*i+1], temp));
        }
    }
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            std::sort(xynodes[i][j].begin(), xynodes[i][j].end(),Compare);
            auto last = std::unique(xynodes[i][j].begin(), xynodes[i][j].end());
            xynodes[i][j].erase(last, xynodes[i][j].end());
        }
    }
}

void Geo::Denser(){
    double xy;
    double Dxy_para[31]{};
    for (int nz = 0; nz < 5; nz++)
    {
        for (int dim = 0; dim < 2; dim++)
        {
            int n = xynodes[nz][dim].size();
            for (int i = 0; i < n - 1; i++)
            {
                double delta = xynodes[nz][dim][i + 1].xy - xynodes[nz][dim][i].xy;
                if ( delta > fsinterval)
                {
                    int nd = int(delta/fsinterval)-1;
                    for (int j = 1; j <= nd; j++)
                    {
                        xy = (xynodes[nz][dim][i].xy * j + xynodes[nz][dim][i + 1].xy * (nd + 1 - j)) / (nd + 1);
                        for (int k = 0; k < 31; k++)
                        {
                            Dxy_para[k] = (xynodes[nz][dim][i].Dxy_para[k] * j + xynodes[nz][dim][i + 1].Dxy_para[k] * (nd + 1 - j)) / (nd + 1);
                        }
                        xynodes[nz][dim].push_back(XYnode(xy,Dxy_para));
                    }
                }
            }
            std::sort(xynodes[nz][dim].begin(),xynodes[nz][dim].end(),Compare);
        }
    }
}

void Module::Init_geo(){
    geo.Init_bdp();//Initialize the overall boundary
    geo.Init_ifibdp();//Initialize the inner boundary of the internal bifurcation structure
    geo.Init_ifobdp();//Initialize the outer boundary of the internal bifurcation structure
    geo.Init_ofibdp();//Initialize the inner boundary of the outer fork structure
    geo.Init_ofobdp();//Initialize the outer boundary of the outer fork structure
    geo.Init_ifcrack();//Initialize the gap of the internal branching structure
    geo.Init_ofcrack();//Initialize the gap of the external fork structure
    geo.Get_Difi_para();
    geo.Get_Difo_para();
    geo.Get_Dofi_para();
    geo.Get_Dofo_para();
    geo.Get_Difcrack_para();
    geo.Get_Dofcrack_para();
    geo.Get_sobdp();
    geo.Init_xynodes();
    //geo.Denser();
}

void Module::Get_num_element_xy(){
    for (int i = 0; i < 5; i++)
    {
        num_element_xy[2 * i + 0] = int(geo.xynodes[i][0].size()) - 1;
        num_element_xy[2 * i + 1] = int(geo.xynodes[i][1].size()) - 1;
    }
}

void Module::Get_nodes(){
    std::vector<double>xnode{};
    std::vector<double>ynode{};
    int xnum, ynum;
    double nodex, nodey;
    double Dxy_para[2*31]{};
    for (int nz = 0; nz < 5; nz++)
    {
        xnum = int(2 * geo.xynodes[nz][0].size()) - 1;
        ynum = int(2 * geo.xynodes[nz][1].size()) - 1;
        if (nz == 0 || nz == 2 || nz == 4)
        {
            num_node_zone_xy[2 * nz + 0] = xnum;
            num_node_zone_xy[2 * nz + 1] = ynum;
            for (int i = 0; i < xnum; i++)
            {
                for (int j = 0; j < ynum; j++)
                {
                    if (i % 2 == 0 && j % 2 == 0)number_nodes_P.push_back(numnode);
                    if (i % 2 != 0)
                    {
                        nodex = (geo.xynodes[nz][0][(i - 1) / 2].xy + geo.xynodes[nz][0][(i + 1) / 2].xy) / 2;
                        for(int k=0;k<31;k++)Dxy_para[31*0+k] = (geo.xynodes[nz][0][(i - 1) / 2].Dxy_para[k] + geo.xynodes[nz][0][(i + 1) / 2].Dxy_para[k]) / 2;
                    }
                    else 
                    {
                        nodex = geo.xynodes[nz][0][i / 2].xy;
                        for (int k = 0; k < 31; k++)Dxy_para[31*0+k] = geo.xynodes[nz][0][i / 2].Dxy_para[k];
                    } 
                    if (j % 2 != 0)
                    {
                        nodey = (geo.xynodes[nz][1][(j - 1) / 2].xy + geo.xynodes[nz][1][(j + 1) / 2].xy) / 2;
                        for (int k = 0; k < 31; k++)Dxy_para[31*1+k] = (geo.xynodes[nz][1][(j - 1) / 2].Dxy_para[k] + geo.xynodes[nz][1][(j + 1) / 2].Dxy_para[k]) / 2;
                    }
                    else
                    {
                        nodey = geo.xynodes[nz][1][j / 2].xy;
                        for (int k = 0; k < 31; k++)Dxy_para[31*1+k] = geo.xynodes[nz][1][j / 2].Dxy_para[k];
                    }
                    nodes.push_back(Node(numnode, nodex, nodey,Dxy_para));
                    numnode++;
                }
            }
            num_node_zone[nz]=numnode;
        }
        else
        {
            num_node_zone_xy[2 * nz + 0] = xnum - 2;
            num_node_zone_xy[2 * nz + 1] = ynum;
            for (int i = 1; i < xnum - 1; i++)
            {
                for (int j = 0; j < ynum; j++)
                {
                    if (i % 2 != 0 && j % 2 == 0)number_nodes_P.push_back(numnode);
                    if (i % 2 != 0)
                    {
                        nodex = (geo.xynodes[nz][0][(i - 1) / 2].xy + geo.xynodes[nz][0][(i + 1) / 2].xy) / 2;
                        for (int k = 0; k < 31; k++)Dxy_para[31*0+k] = (geo.xynodes[nz][0][(i - 1) / 2].Dxy_para[k] + geo.xynodes[nz][0][(i + 1) / 2].Dxy_para[k]) / 2;
                    }
                    else
                    {
                        nodex = geo.xynodes[nz][0][i / 2].xy;
                        for (int k = 0; k < 31; k++)Dxy_para[31*0+k] = geo.xynodes[nz][0][i / 2].Dxy_para[k];
                    }
                    if (j % 2 != 0)
                    {
                        nodey = (geo.xynodes[nz][1][(j - 1) / 2].xy + geo.xynodes[nz][1][(j + 1) / 2].xy) / 2;
                        for (int k = 0; k < 31; k++)Dxy_para[31*1+k] = (geo.xynodes[nz][1][(j - 1) / 2].Dxy_para[k] + geo.xynodes[nz][1][(j + 1) / 2].Dxy_para[k]) / 2;
                    }
                    else 
                    {
                        nodey = geo.xynodes[nz][1][j / 2].xy;
                        for (int k = 0; k < 31; k++)Dxy_para[31*1+k] = geo.xynodes[nz][1][j / 2].Dxy_para[k];
                    } 
                    nodes.push_back(Node(numnode, nodex, nodey,Dxy_para));
                    numnode++;
                }
            }
            num_node_zone[nz]=numnode;
        }
    }
}

int Module::Get_number_node(int index, int jndex, int zonenum){
    int num = 0;
    if(zonenum==0)
    {
        num = index * num_node_zone_xy[2 * zonenum + 1] + jndex;
    }
    else
    {
        if(index==0)num = num_node_zone[zonenum-1]- num_node_zone_xy[2 * (zonenum-1) + 1] + index * num_node_zone_xy[2 * zonenum + 1] + jndex;
        else num = num_node_zone[zonenum - 1] + index * num_node_zone_xy[2 * zonenum + 1] + jndex;
    }
    return num;
}

bool Module::Issolid(double* centerxy){
    std::vector<double> ofcracki(12,0);
    std::vector<double> ifcracki(12,0);
    if (Isinando(geo.bdp,geo.ofobdp,centerxy))return true;
    for (int i = 0; i < numoutfork - 1; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            ofcracki[2 * j + 0] = geo.ofcrack[6 * 2 * i +2 * j + 0];
            ofcracki[2 * j + 1] = geo.ofcrack[6 * 2 * i + 2 * j + 1];
        }
        if (Isinregion(ofcracki,centerxy))return true;
    }
    for (int i = 0; i < numinfork - 1; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            ifcracki[2*j+0] = geo.ifcrack[6 * 2 * i + 2*j+0];
            ifcracki[2 * j + 1] = geo.ifcrack[6 * 2 * i + 2 * j + 1];
        }
        if (Isinregion(ifcracki,centerxy))return true;
    }
    if (Isinando(geo.ofibdp,geo.ifobdp,centerxy))return true;
    if (Isinregion(geo.ifibdp,centerxy))return true;
    return false;
}

void Module::Get_elements(){
    for (int numzone = 0; numzone < 5; numzone++)
    {
        int xnum = num_element_xy[2 * numzone+0];
        int ynum = num_element_xy[2 * numzone+1];
        for (int i = 0; i < xnum; i++)
        {
            for (int j = 0; j < ynum; j++)
            {
                int i1 = Get_number_node(2 * i, 2 * j, numzone);
                int i2 = Get_number_node(2 * i + 2, 2 * j, numzone);
                int i3 = Get_number_node(2 * i + 2, 2 * j + 2, numzone);
                int i4 = Get_number_node(2 * i, 2 * j + 2, numzone);
                int temp[9]{i1,i2,i3,i4,(i1 + i2) / 2 ,(i2 + i3) / 2 , (i3 + i4) / 2 ,(i1 + i4) / 2,(i1 + i2 + i3 + i4) / 4 };
                double nodex = nodes[(i1 + i2 + i3 + i4) / 4].x;
                double nodey = nodes[(i1 + i2 + i3 + i4) / 4].y;
                if ((nodex > 60 && nodex < 80 && nodey > 5 && nodey < 10) || (nodex > 137 && nodex < 142 && nodey > 52 && nodey < 62))
                {
                    number_element_hs.push_back(numelement);
                }
                elements.push_back(Element(numelement, temp));
                double cxy[2]{nodes[(i1 + i2 + i3 + i4) / 4].x,nodes[(i1 + i2 + i3 + i4) / 4].y};
                if (Issolid(cxy))
                {
                    elements[numelement].type = 1;
                }
                numelement++;
            }
        }
    }
}

void Module::Write_data(){
    std::ofstream outFile("./data.txt");
    if (outFile.is_open())
    {
        outFile << "xy and number of nodes" << std::endl;
        for (int i = 0; i < numnode; i++)
        {
            outFile << nodes[i].number << "  " << nodes[i].x << "  " << nodes[i].y << std::endl;
        }
        outFile << "nodes and number of elements" << std::endl;
        for (int i = 0; i < numelement; i++)
        {
            outFile << elements[i].number;
            for (int j = 0; j < 9; j++)
            {
                outFile << "  " << elements[i].nodes[j];
            }
            outFile << std::endl;
        }
        outFile << "Dxy_paras of nodes" << std::endl;
        for (int i = 0; i < numnode; i++)
        {
            for(int j=0;j<31;j++)
            {
                outFile << nodes[i].DXY_para[j]<< "  ";
            }
        }
    }
    outFile.close();
}

void Module::Write_heatsourcedata(){
    std::ofstream outFile("./hsdata.txt");
    int index;
    int n = int(number_element_hs.size());
    if (outFile.is_open())
    {
        outFile << "number,number of nodes,xy of nodes of heatsource elements" << std::endl;
        for (int i = 0; i < n; i++)
        {
            index = number_element_hs[i];
            outFile << index << std::endl;
            for (int j = 0; j < 9; j++)
            {
                outFile << elements[index].nodes[j] << " , "<< nodes[elements[index].nodes[j]].x << " , " << nodes[elements[index].nodes[j]].y << " ; ";
            }
            outFile << std::endl;
        }
    }
    outFile.close();
}

void Module::Get_number_node_T(){
    double nodexy[2]{};
    for (int i = 0; i < numelement; i++)
    {
        if (elements[i].type == 1)
        {
            for (int j = 0; j < 9; j++)
            {
                nodexy[0] = nodes[elements[i].nodes[j]].x;
                nodexy[1] = nodes[elements[i].nodes[j]].y;
                if (!Isonboundary(geo.sobdp,nodexy))number_nodes_T.push_back(elements[i].nodes[j]);
            }
        }
    }
}
