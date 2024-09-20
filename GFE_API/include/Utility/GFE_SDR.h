#pragma once
#include <vector>
#include <string>
#include "../api_a.h"
#include "../GFE_Struct/GFE_Outp.h"

class GFE_API LayerDisplacementAngleCalc
{
public:
    LayerDisplacementAngleCalc(std::shared_ptr<GFE::DB> db,std::vector<std::string> const &nodesetName);
    ~LayerDisplacementAngleCalc();
    bool prepare(int Z = 2);
    bool startCalculation(int frameno,std::vector<double>& dispx,std::vector<double>& dispy, int Z = 2);
    const std::vector<std::pair<int, int>>& MaxNodesX() { return max_nodes[0]; }
    const std::vector<std::pair<int, int>>& MaxNodesY() { return max_nodes[1]; }
    const auto& InvalidLayers() { return invalid_layers; }

private:

    enum Dir {
        X, Y, Z
    };

    //集合
    std::vector<std::string>        m_nodesetName;
    std::shared_ptr<GFE::DB>        m_db;
    struct storeyPair
    {
        std::vector<int> floordownnode;
        std::vector<int> floorupnode;
    };
    //每一层的点对应关系
    std::vector<storeyPair> m_storeynodes;
    //每一层多组对应关系
    std::vector<std::vector<std::pair<int, int>>> m_storeynodes2;
    std::vector<float> nodecoords;
    std::vector<std::pair<int, int>> max_nodes[2];
    std::vector<int> invalid_layers;
};
