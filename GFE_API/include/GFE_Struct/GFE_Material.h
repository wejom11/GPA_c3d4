#pragma once
#include "../api_a.h"
#include "GFE_Function.h"
#include <map>
#include <vector>
#include "math.h"

namespace GFE {
struct Mat_General {
    explicit Mat_General(int _t) : type(_t) {}
    virtual ~Mat_General(){}
    enum EntryType {
        Density,
        Elastic,
        Plastic,
        HyperFoam,
        HyperElastic,
        Damping,
        ViscoElastic,
        ConcreteDamaged,
        MohrCoulomb,
        User,
        TestData,
        Creep,
        Permeability,
        PorousBulkModuli,
        Sorption,
        Expansion,
        BedCoefficient,			// 基床系数
        RateDependent,
        General_Type_Size
    };
    enum {
        Long_term,
        Instantaneous
    };
    inline int Type() { return type; }
    static const auto& TypeStrMap() {
        static std::map<int, string> aMap = {
            {Density, "Density"},
            {Elastic, "Elastic"},
            {Plastic, "Plastic"},
            {HyperFoam, "HyperFoam"},
            {HyperElastic, "HyperElastic"},
            {Damping, "Damping"},
            {ViscoElastic, "ViscoElastic"},
            {ConcreteDamaged, "ConcreteDamaged"},
            {MohrCoulomb, "MohrCoulomb"},
            {User, "User"},
            {TestData, "TestData"},
            {Creep, "Creep"},
            {Permeability, "Permeability"},
            {PorousBulkModuli, "PorousBulkModuli"},
            {Sorption, "Sorption"},
            {Expansion, "Expansion"},
            {Creep, "Creep"},
            {BedCoefficient, "BedCoeff"},
            {RateDependent, "RateDep"}
        };
        return aMap;
    }
    auto TypeStr() { return TypeStrMap().at(type); }
    virtual Mat_General* Clone() = 0;
private:
    int type;
};

template <class T>
struct MatGenTemp : public Mat_General {
    explicit MatGenTemp(EntryType _t) : Mat_General(_t) {}
    virtual Mat_General* Clone() override {
        return new T(*static_cast<T*>(this));
    }
};

struct Mat_Density : public MatGenTemp<Mat_Density> {
    Mat_Density() : MatGenTemp<Mat_Density>(Mat_General::Density) {}

    bool temp_dp = false;       // temperature dependence
    int n_param = 1;            // params存了一张表, n_param是列数
    vector<double> params;
};

struct Mat_Elastic : public MatGenTemp<Mat_Elastic> {
    Mat_Elastic() : MatGenTemp<Mat_Elastic>(Mat_General::Elastic) {}
    enum {
        Isotropic
    };

    bool temp_dp = false;
    int n_param = 2;
    int type = Isotropic;
    int moduli_time_scale = Long_term;
    bool compression = false;
    bool tension = false;
    vector<double> params;     // 一般0号位是Young's Modulus, 1是Poisson's Ratio, 其余是temp_dp, n_param扩展
};

struct Mat_Damping : public MatGenTemp<Mat_Damping> {
    Mat_Damping() : MatGenTemp<Mat_Damping>(Mat_General::Damping) {}

    int n_param = 2;
    vector<double> params;
};

struct Mat_Plastic : public MatGenTemp<Mat_Plastic> {
    Mat_Plastic() : MatGenTemp<Mat_Plastic>(Mat_General::Plastic) {}
    enum {
        Isotropic,
        Johnson_cook
    };
    int harden_type = Isotropic;
    bool rate_dp = false;
    bool temp_dp = false;
    vector<double> params;     // 前两位为shearfailure Isotropic时为stress strain (rate) (temp), Johnson_cook时为6个参数
};

struct Mat_HyperFoam : public MatGenTemp<Mat_HyperFoam> {
    Mat_HyperFoam() : MatGenTemp<Mat_HyperFoam>(Mat_General::HyperFoam) {}
    bool test_data = false;
    int N = 1;
    bool temp_dp = false;
    int moduli_time_scale = Long_term;
    // 直接定义参数
    vector<double> params;    // mu1 alpha1 ... muN alphaN nu1 ... nuN (temp...)
    // 使用Test Data
    vector<double> uniaxial, biaxial, simple_shear, planar, volumetric;
};

struct Mat_HyperElastic : public MatGenTemp<Mat_HyperElastic> {
    Mat_HyperElastic() : MatGenTemp<Mat_HyperElastic>(Mat_General::HyperElastic) {}
    enum TestDataType {
        MOONEY_RIVLIN,
        OGDEN,
        YEOH,
        POLY,
        MARLOW,
        REPOLY
    };
    int he_type = MOONEY_RIVLIN;
    bool test_data = false;
    int moduli_time_scale = Long_term;
    bool has_poisson = false;
    double poisson = 0;
    int N = 1;      // Marlow时无意义
    bool temp_dp = false;
    // 直接定义参数
    vector<double> params;     // para1 ...paraN (temp)
    // 使用Test Data
    vector<double> uniaxial, biaxial, planar, volumetric;
};

struct Mat_ViscoElastic : public MatGenTemp<Mat_ViscoElastic> {
    Mat_ViscoElastic() : MatGenTemp<Mat_ViscoElastic>(Mat_General::ViscoElastic) {}
    int type = -1;  // 0 Frequency; 1 Time; -1 Null
    int n_param = 0;
    // type=0: 5个为一组(omegaRg, omegaIg, omegaRk, omegaIk, freq)
    // type=1: 3个为一组(g, k, tau)
    vector<double> params;
};

struct Mat_ConcreteDamaged : public MatGenTemp<Mat_ConcreteDamaged> {
    Mat_ConcreteDamaged() : MatGenTemp<Mat_ConcreteDamaged>(Mat_General::ConcreteDamaged) {}
    int n_plasticity = 1;
    vector<double> plasticity;
    int n_comp_harden = 0;
    vector<double> comp_harden;
    int n_comp_damage = 0;
    vector<double> comp_damage;
    int n_tens_stiff = 0;
    vector<double> tens_stiff;
    int n_tens_damage = 0;
    vector<double> tens_damage;
};

struct Mat_MohrCoulomb : public MatGenTemp<Mat_MohrCoulomb> {
    Mat_MohrCoulomb() : MatGenTemp<Mat_MohrCoulomb>(Mat_General::MohrCoulomb) {}
    int n_plasticity = 2;
    vector<double> plasticity;
    int n_cohesion = 0;
    vector<double> cohesion;
};

struct Mat_User : public MatGenTemp<Mat_User> {
    Mat_User() : MatGenTemp<Mat_User>(Mat_General::User) {}
    enum {
        Normal,
        OneDN,
        Davi
    };
    int user_type = Normal;
    int n_constants = 0;
    vector<double> constants;
};

struct Mat_TestData : public MatGenTemp<Mat_TestData> {
    Mat_TestData() : MatGenTemp<Mat_TestData>(Mat_General::TestData) {}
    int n_test_data = 0; // 行数
    vector<double> test_data;
};

struct Mat_Creep : public MatGenTemp<Mat_Creep> {
    Mat_Creep() : MatGenTemp<Mat_Creep>(Mat_General::Creep) {}

    enum LawEnum {
        Strain = 1,
        Time,
        HyperB,
        User,
        Anand,
        Darveaux,
        DoublePower
    };
    enum TimeEnum {
        TIME_Creep = 1,
        TIME_Total
    };

    int law = 0;
    int time = 0;
    int nRow = 1;       // 跟温度相关时data里可能包含了多行数据, 目前应该用不到
    vector<double> data;
};

class Mat_Permeability : public MatGenTemp<Mat_Permeability> {
public:
	Mat_Permeability() : MatGenTemp<Mat_Permeability>(Mat_General::Permeability) {}

	enum Type {
		ISOTROPIC = 0,
		ORTHOTRPIC,
		ANISOTROPIC,
		SATURATION,
		VELOCITY
	};
	struct Internal {
		int type = 0;
		double specific = 0;
		std::vector<double> data;
	};
	std::vector<Internal> entities;
};

class Mat_PorousBulkModuli : public MatGenTemp<Mat_PorousBulkModuli> {
public:
	Mat_PorousBulkModuli() : MatGenTemp<Mat_PorousBulkModuli>(Mat_General::PorousBulkModuli) {}

	double solid_grains = 0;
	double permeating_fluid = 0;
};

class Mat_Sorption : public MatGenTemp<Mat_Sorption> {
public:
	Mat_Sorption() : MatGenTemp<Mat_Sorption>(Mat_General::Sorption) {}

	enum Law { LOG = 0, TABULAR };
	enum Type { ABSORPTION = 0, EXSORPTION, SCANNING };

	struct Internal {
		int type = 0;
		int law = 1;
		std::vector<double> data;
	};
	std::vector<Internal> entities;
};

class Mat_Expansion : public MatGenTemp<Mat_Expansion> {
public:
	Mat_Expansion() : MatGenTemp<Mat_Expansion>(Mat_General::Expansion) {}

	enum Type {
		ISOTROPIC = 0,
		ORTHOTRPIC,
		ANISOTROPIC,
	};
	int sub_type = 0;
	std::vector<double> value;
};

struct Mat_BedCoeff : public MatGenTemp<Mat_BedCoeff> {
	Mat_BedCoeff() : MatGenTemp<Mat_BedCoeff>(Mat_General::BedCoefficient) {}

	double kh = 0;
	double kv = 0;
};

struct Mat_RateDep : public MatGenTemp<Mat_RateDep> {
	Mat_RateDep() : MatGenTemp<Mat_RateDep>(Mat_General::RateDependent) {}

	enum Type {
		POWER_LAW = 0,
		JOHNSON_COOK,
		YIELD_RATIO
	};

	int sub_type = 0;
	std::vector<double> value;
};

struct GFE_API Material {
    Material() = default;
    Material(const Material&);
    Material(Material&&);

    Material& operator= (const Material& rhs) { DeepCopy(rhs, *this); return *this; }
    Material& operator= (Material&&);

    static void DeepCopy(const Material& src, Material& tg);
    static void ShallowCopy(const Material& src, Material& tg);

    id_t id;
    string name;
    vector<shared_ptr<Mat_General>> entries;
    void AsElastic(double rou, double e, double nu);

#define GFE_MAT_GETTER(T) \
    shared_ptr<Mat_##T> Get##T() const {            \
            for(const auto& i : entries)                \
            if(i->Type() == Mat_General::T)    \
            return std::dynamic_pointer_cast<Mat_##T>(i);    \
            return nullptr; \
    }

    GFE_MAT_GETTER(Density);
    GFE_MAT_GETTER(Elastic);
    GFE_MAT_GETTER(Damping);
    GFE_MAT_GETTER(Plastic);
    GFE_MAT_GETTER(HyperFoam);
    GFE_MAT_GETTER(HyperElastic);
    GFE_MAT_GETTER(ViscoElastic);
    GFE_MAT_GETTER(ConcreteDamaged);
    GFE_MAT_GETTER(MohrCoulomb);
    GFE_MAT_GETTER(User);
    GFE_MAT_GETTER(TestData);
};
}

//! Material2
//! 辅助数据结构，解决Material取数据麻烦的问题
namespace GFE {
class GFE_API Material2 {
public:
    Material2() { Reset(); }
    Material2(const Material& aMat) { Reset(aMat); }
    Material2& operator =(const Material& aMat) { Reset(aMat); return *this; }

    void Reset();
    void Reset(const Material&);
//    long long Read(const std::string& path);
//    long long Write(const std::string& path);

    const std::string& Name() { return _mat.name; }
    double GetDensity() { return _density; }
    const auto& GetElastic() { return _elastic; }
    const auto& GetDamping() { return _damping; }

    struct Para_SSA   // 地下结构常用参数, SSA: Sub-Structure Analysis
    {
        double Rou = NAN;   // 密度
        double E = NAN;     // 弹性模量
        double V = NAN;     // 泊松比
        double G = NAN;     // 剪切模量(二阶拉梅常数)
        double K = NAN;     // 一阶拉梅常数
        double CS = NAN;    // 剪切波速
        double CP = NAN;    // 压缩波速
        double DA = 0;      // α阻尼
        GFE::Function2 TestGR, TestDR;       // 非线性：应变-剪切模量、应变-阻尼比曲线
    };
    const Para_SSA& GetPara_SSA() { return _ssa; }

private:
    using D2 = std::pair<double, double>;
    Material _mat;
    double _density;
    D2 _elastic;
    D2 _damping;
    Para_SSA _ssa;
};
}

