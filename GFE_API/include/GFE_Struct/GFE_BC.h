#pragma once
#pragma warning(disable: 4251)
#include "../api_a.h"
#include <map>

namespace GFE {
struct GFE_API GenericBC {
    enum BCType {
        CONSTRAIN,
        DISPLACEMENT,
        INITIAL_VELOCITY,
        VELOCITY,
        ACCELERATION,
        FORCE,
        SURFACE,
        GRAVITY,
        LINELOAD,
        TRAINLOAD,
        MASS,
        ACOUSTIC,

        //! Optistruct
        RLOAD,

        //! 没用到
        PRESSURE,
        TEMPERATURE,

        //! 水压力
        WATER_PRESSURE,
        PORE_PRESSURE,
        C_PORE_FLOW,		// concentrated pore fluid
        S_PORE_FLOW,		// surface pore fluid
        BODY_FORCE,
        AC_DISP,			// acoustic dof 8
        SIZE
    };
    GenericBC() = default;
    GenericBC(bool _is_ref, int _type, string _name) : is_ref(_is_ref), type(_type), name(_name) {}
    virtual ~GenericBC() {}
    static const std::map<int, string>& TypeStrMap();
    string TypeStr() const;
    static const std::map<int, string>& TypeStrMap_Short();
    string TypeStr_Short() const;

    bool is_ref = false;        // 标识该边界条件是否仅用于被其他边界条件引用
    int type = -1;
    string name;
};
using BCType = GenericBC::BCType;

struct GFE_API RegularBC : public GenericBC {
    enum DistType {
        DIST_UNIFORM = 0,
        DIST_X = 1,     // 沿x轴分布
        DIST_Y = 2,
        DIST_Z = 3,
        DIST_USER = 4,

        DIST_END_ABS = 0xffff,
        DIST_RELATIVE = 0x10000		// 相对的空间分布

    };
    virtual ~RegularBC() {}

    bool is_node_set = true;        //列车荷载时，表示Region是否为网格节点集
    vector<int> node_id;            //列车荷载时，为起点终点号
    string set_name;
    int valid_dof;  // 六位0/1代表六个自由度 mass时为Units
    vector<double> values;      //列车荷载时，前三位为ScaleFactor,Velocity,StartTime，3-5为方向，6为选择的算法，后面为列车数据
    string func_name;
    vector<double> values_im;
    string func_name_im;
    string step_name;
    vector<int> track_id;
    vector<double> track_coord;

    //! 边界条件在空间上不均匀分布的相关设置
    string dist_name;               // 荷载的空间分布函数, 引用AMP
    [[deprecated("use field instead")]] int dist_type = 0;              // 空间分布的函数类型, 见DistType
    [[deprecated("use field instead")]] std::vector<double> dist_dir = {0, 0, 0};		// 用户空间分布方向 DistType = DIST_USER

    //! 列车荷载考虑惯性力
    bool has_if = false;
    int if_mode = 0;   // 0: 轨道实测法, 1: 轨道不平顺法
    string if_acce_path;
    vector<double> if_acce;
    vector<double> if_param;    // 轮对质量，转向架质量，车厢质量，阻尼1，阻尼2，刚度1，刚度2
    double if_itv;      // 仅用于轨道实测法
    double if_tot;      // 仅用于轨道不平顺法
    int if_grade = 0;  // 仅用于轨道不平顺法
    vector<double> if_force;    // 仅用于惯性力直接输入
};
using Boundary = RegularBC;

struct GFE_API WaterPressure : public GenericBC
{
    WaterPressure() : GenericBC() { type = WATER_PRESSURE; }
    virtual ~WaterPressure() {}

    enum eModel
    {
        Westgarrd,
        Housner,
        eModel_Size
    };

    int model;              // 计算模型
    string surf_name;
    double density;
    double lvl_height;
    double sink_wide;       // 水槽宽度，适用于Housner模型
};

namespace OptiStruct {
struct Rload : public GenericBC {
    Rload() : GenericBC() { type = RLOAD; }
    virtual ~Rload() {}

    int sub_type = -1;			// 0 RLOAD1, 1 RLOAD2, 2 ACSRCE->SLOAD
    double scale = 1;
    double phase = 0;
    double t = 0;
    string func_name;
    string func_name_im;
    string ref_name;		// reference bc name
};
}
namespace Opti = OptiStruct;
using BC_Rload = Opti::Rload;

namespace BC {
GFE_API void Add(shared_ptr<DB>, const GenericBC*);
GFE_API shared_ptr<GenericBC> Find(shared_ptr<DB>, const std::string& name, bool nocase = false);
GFE_API vector<string> AllName(shared_ptr<DB>, int type = -1);
GFE_API vector<shared_ptr<GenericBC>> All(shared_ptr<DB>, int type = -1);
GFE_API vector<shared_ptr<GenericBC>> AllNotRef(shared_ptr<DB>);


template <class T,
         typename = typename std::enable_if<std::is_base_of<GenericBC, T>::value>::type>
shared_ptr<T> Find(shared_ptr<DB> db, const std::string& name, bool nocase = false) {
    auto p = Find(db, name, nocase);
    if(!p) return nullptr;
    return std::dynamic_pointer_cast<T>(p);
}

template <class T,
         typename = typename std::enable_if<std::is_base_of<GenericBC, T>::value>::type,
         typename = typename std::enable_if<!std::is_same<T, GenericBC>::value>::type>
vector<shared_ptr<T>> All(shared_ptr<DB> db) {
    vector<shared_ptr<T>> ret;
    if(std::is_same<T, Opti::Rload>::value) {
        auto tmp = All(db, BCType::RLOAD);
        ret.reserve(tmp.size());
        for(auto p : tmp) ret.push_back(std::static_pointer_cast<T>(p));
    }
    else {
       auto tmp = All(db);
        ret.reserve(tmp.size());
       for(auto p : tmp) {
            auto d = std::dynamic_pointer_cast<T>(p);
           if(d) ret.push_back(d);
       }
    }
    return ret;
}
}

};

