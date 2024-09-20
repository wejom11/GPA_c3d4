#pragma once
#include "obj_st.h"
#include <memory>
#include <map>
#include "GFE_Struct/GFE_Function.h"
#include "GFE_Struct/GFE_BC.h"
#include "GFE_Struct/GFE_Outp.h"
#include "GFE_Struct/GFE_XYData.h"
#include "GFE_Struct/GFE_Material.h"
#include "GFE_Struct/GFE_Connector.h"

namespace GFE
{
struct Node {
    // TODO: 考虑inp中存在多个part时, label要与part name组合才构成节点的唯一索引
    id_t id;
    id_t label;
    data_t x;
    data_t y;
    data_t z;
};

struct Element {
    // TODO: 考虑inp中存在多个part时, label要与part name组合才构成单元的唯一索引
    id_t id;
    id_t label;
    int type;
    int material;
    vector<id_t> nodes;
};

struct NodeSet {
    string name;
    vector<id_t> nodes;
};

struct ElementSet {
    string name;
    vector<id_t> elements;
};

struct Surface {
    string name;
    vector<id_t> elements;
    vector<char> faces;
    vector<id_t> nodes;
};

struct SurfacePair {
    enum {
        Tie
    };
    int type;
    string name;
    string surf_1;
    string surf_2;
    int n_param = 0;
    vector<double> params;
    // Tie: Position_Tolerance
};

struct FieldOutput {
    int frame;
    string variable;
    vector<data_t> data;
};

struct FieldOutputRequest {
    string variable;
    string region;
    int type;
};

struct FieldOutputFrame {
    int frame;
    double time;
    bool has_output;
};

struct HistoryOutput {
    id_t id;
    id_t point_id;
    int position;
    string region;
    string variable;
    vector<data_t> data;
    vector<double> frequency;
};

struct HistoryOutputRequest {
    int id;
    string region;
    int type;
    string variables;
    vector<double> frequency;
};

struct Property {
    Property() : type(Null) {}
    explicit Property(int _t) : type(_t) {}
    virtual ~Property(){}
    enum {
        Null = -2,
        Unknown,
        Solid,
        Bush,
        Beam,
        Shell,
        Membrane,
        BeamGeneral,
        Connector
    };
    static const auto& TypeStrMap() {
        static std::map<int, string> aMap = {
            {Null, "Null"},
            {Unknown, "Unknown"},
            {Solid, "Solid"},
            {Bush, "Bush"},
            {Beam, "Beam"},
            {Shell, "Shell"},
            {Membrane, "Membrane"},
            {BeamGeneral, "BeamGeneral"},
            {Connector, "Connector"}
        };
        return aMap;
    }
    decltype(auto) TypeStr() const { return TypeStrMap().at(type); }
    virtual void SetMaterial(const std::string& mat) {};
    virtual string Material() { return ""; }
    id_t id = -1;
    string name;
    string elset_name;
    vector<double> ecc = vector<double>(4, 0);
    int type = Solid;
    int subtable_row = -1;
};

struct PropertySolid : public Property {
    PropertySolid() : Property(Property::Solid) {}
    void SetMaterial(const std::string& mat) override { mat_name = mat; }
    string Material() override { return mat_name; }
    int mat_id;
    string mat_name;
    bool has_thickness = false;
    double thickness = 1;
};

struct PropertyBush : public Property {
    PropertyBush() : Property(Property::Bush) {}
    // TODO: bush type
    // TODO: orientation
    vector<double> params;
    // 接地弹簧: kx, ky, kz, krx, kry, krz, cx, cy, ...

};

struct PropertyBeamGeneral : public Property {
    PropertyBeamGeneral() : Property(Property::BeamGeneral) {}
    double density, poisson = 0;
//    string section = "GENERAL";
    vector<double> param1;      // A, I11, I12, I22, J, gamma0, gammaW
    vector<double> axis;
    vector<double> param2;      // E, G, ...extension
};

struct PropertyBeam : public Property {
    PropertyBeam() : Property(Property::Beam) {}
    void SetMaterial(const std::string& mat) override { mat_name = mat; }
    string Material() override { return mat_name; }
    enum {
        RECT,
        BOX,
        I,
        CIRCULAR,
        L,
        PIPE,
        THICKPIPE,
        ARB
    };
    int shape = RECT;
    int mat_id;
    string mat_name;
    int fiber_num = 1;
    vector<double> shape_params;
    vector<double> params;
    // 梁单元：ah*fiber_num, x*fiber_num, y*fiber_num, A, Asy, Asz, Ixx, Iy, Iz, Gxy, Gyz, shape
    vector<double> direction;
    vector<double> shear;
};

struct GFE_API RebarLayer {
	RebarLayer() = default;
	explicit RebarLayer(const char*);
    std::vector<char> ToBinary() const;
    string TypeStr() const { return "Rebar"; }

    int rebar_geometry = 0; // 0: constant 1: lift-equation
    string orientation_name = "";
    int rebar_num = 0;
    vector<string> layer_name;
    vector<string> mat_name;
    vector<double> params;
    // rebar_num * (Area, Spacing, Position, Angle, Direction, Ext_Ratio, Radius)
};

struct PropertyShell : public Property {
    PropertyShell() : Property(Property::Shell) {}
    void SetMaterial(const std::string& mat) override { mat_name = mat; }
    string Material() override { return mat_name; }
    int mat_id;
    string mat_name;
    double thickness = 0;
    int integral_point = 1;
    int layer_num = 1;
    vector<double> params;
    bool has_rebar = false;
    RebarLayer rebar;
    // 每层厚度
};

struct PropertyMembrane : public Property {
    PropertyMembrane() : Property(Property::Membrane) {}
    void SetMaterial(const std::string& mat) override { mat_name = mat; }
    string Material() override { return mat_name; }
    int mat_id;
    string mat_name;
    double thickness = 0;
    bool has_rebar = false;
    RebarLayer rebar;
};

struct GFE_API PropertyConnector : public Property {
    PropertyConnector() : Property(Property::Connector) {}
    enum ConnectorType {
        CARTESIAN = 0,
        JOIN,
        LINK,
        ROTATION
    };
    string behavior;
    string orientation;
    int translational_type = -1;
    int rotational_type = -1;

    static std::string TypeStr(int type);
};

struct RigidBody {
    int id;
    int type = 0;		// 0 element rigid, 1 surface rigid
    string name;
    string set_name;
    double thickness;
    double density;
    int ref_node;
    std::string ref_set;
};

struct MPC {
    int id;
    int n_param;
    vector<int> params;
};

struct Embed {
    int id;
    string name;
    string host_name;
    double roundoff_tolerance;
    double exterior_tolerance;
    vector<string> embedded_names;
};

struct IncidentWaveProperty {
    enum {
        AirBlast,
        SurfaceBlast
    };
    int id;
    int def = AirBlast;
    string name;
    vector<double> data;
};

struct IncidentWave {
    int id;
    string name;
    bool is_node_set;
    string set_name;
    int node_id;
    string surf_name;
    double time_detonation;     //CONWEP
    double mag_scale_factor;    //CONWEP
    string prop_name;
};

struct ArtBoundary {
    int id;
    string name;
    vector<double> struct_centre = vector<double>(3, 0);
    string surf_name;
};

struct SoilLayer {
    int id;
    string name;
    int n_layer = 0;
    vector<double> layer_thickness;
    vector<string> layer_mat;
    string bedrockMat;
};

struct VibLoad {
    int id;
    string name;
    std::unordered_map<string, string> param;
};

struct SpringDashpot {
    enum SDType{
        spring1 =  0b000001,
        spring2 =  0b000010,
        springA =  0b000100,
        dashpot1 = 0b001000,
        dashpot2 = 0b010000,
        dashpotA = 0b100000
    };
    static const auto& TypeStrMap() {
        static std::map<int, string> aMap = {
            {spring1, "Spring1"},
            {spring2, "Spring2"},
            {springA, "SpringA"},
            {dashpot1, "Dashpot1"},
            {dashpot2, "Dashpot2"},
            {dashpotA, "DashpotA"},
        };
        return aMap;
    }
    decltype(auto) TypeStr() const { return TypeStrMap().at(type); }
    int id;
    int type = 0;
    string name;
    // internal elset name
    string nset;					// spring, dash
    vector<string> nset2;		// spring2, springA
    vector<int> nodes = vector<int>(2, 0);           // 用于后处理，弹簧端点id
    double dof = 0;
    double dof2 = 0;
    double stiffness = 0;
    double coefficient = 0;
    std::string orientation;
};

}
