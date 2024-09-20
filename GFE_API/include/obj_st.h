#pragma once
#include "api_a.h"

namespace GFE
{
//struct Quantity_ {
//    string name;
//    int type1;
//    int type2;
//    string description;
//};

//struct Compontent_ {
//    string name;
//    int quantity;
//    string description;
//};

struct Variable_ {
    string quantity;
    string component;
    char type;
};

struct Element_ {
    int id;
    string name;
    int node_num;
    int face_num;
    int max_face_node_num;
    vector<short> face_order;
};

struct ElementSubType_ {
    int id;
    string name;
    int parent;
};

// 若新增单元, 高阶单元请放在低阶之后, 否则getElementFaceNode会出错
enum CT{
    CT_VERTEX = -2,         // 节点
    CT_UNKNOWN = -1,
    // Linear Element
    CT_POINT,
    CT_LINE,
    CT_TRIANGLE,
    CT_QUAD,
    CT_TETRA,
    CT_HEXA,
    CT_WEDGE,
    CT_PYRAMID,
    // Quadratic Element
    CT_TETRA_2,
    CT_NUM
};

// 后处理中的分类
enum VT{
    VT_UNKNOWN = -1,
    VT_NODE,
    VT_ELEMENT,
    VT_SCALAR,
    VT_NUM
};

// 处理中的分类
enum VT2 {
    VT2_UNKNOWN = -1,
    VT2_NODE,
    VT2_ELEMENT,
    VT2_ENERGY,
    VT2_CONTACT,
    VT2_INTEGRATED,
    VT2_OTHER,
    VT2_NUM
};

// 若新增单元, 高阶单元请放在低阶之后, 否则getElementFaceNode会出错
namespace [[deprecated("use enum CT_XXXX instead")]] CellType {
enum{
    VERTEX = -2,         // 节点
    UNKNOWN = -1,
    // Linear Element
    POINT,
    LINE,
    TRIANGLE,
    QUAD,
    TETRA,
    HEXA,
    WEDGE,
    PYRAMID,
    // Quadratic Element
    TETRA_2,
    TYPE_NUM,
};
};

// 后处理中的分类
namespace [[deprecated("use enum VT_XXXX instead")]] VariableType {
enum {
    UNKNOWN = -1,
    NODE,
    ELEMENT,
    SCALAR
};
}

// 前处理中的分类
namespace [[deprecated("use enum VT2_XXXX instead")]] VariableType2 {
enum {
    UNKNOWN = -1,
    NODE,
    ELEMENT,
    ENERGY,
    CONTACT,
    INTEGRATED = 4
};
}

}
