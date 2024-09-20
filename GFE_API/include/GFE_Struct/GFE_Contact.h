#pragma once

#include <optional>
#include <array>
#include "api_a.h"

namespace GFE {

using std::vector;
using std::array;
using std::optional;

struct GFE_API CohesiveBehavior
{
    enum {
        ELI_ORIGIN,
        ELI_CURRENT,
        ELI_SPECIFIED
    };

    int eligibility = ELI_ORIGIN;
    vector<array<double, 8>> datas;

    void Clear();
};

struct GFE_API DamageEvolution
{
    enum {
        TP_DISPLACEMENT,
        TP_ENERGY,
        TP_HYSTERESIS
    };
    enum {
        SO_LINEAR,
        SO_EXPONENTIAL,
        SO_TABULAR
    };
    enum {
        MMB_TABULAR,
        MMB_POWER,
        MMB_BK
    };

    int type = TP_ENERGY;
    int softening = SO_LINEAR;
    int mixed_mode_behavior = MMB_POWER;
    double power = 2.0;
    vector<array<double, 8>> datas;

    void Clear();
};

struct GFE_API DamageInitiation
{
    enum {
        CRI_QUADS
    };

    int criterion = CRI_QUADS;
    vector<array<double, 8>> datas;

    void Clear();
};

struct GFE_API SurfaceInteraction
{
    string name;
    optional<CohesiveBehavior> cohesive_behavior;
    optional<DamageEvolution> damage_evolution;
    optional<DamageInitiation> damage_initiation;
};

struct GFE_API ContPropAssign
{
    std::string step;       // 这个变量暂时固定是Initial, 没用
    std::vector<std::string> datas;

    void Clear();
};

struct ContPropAssign_Partition
{
    std::string step;
    std::string surface1;
    std::string surface2;
    std::string property;
    int partition_id = -1;
};

GFE_API void AddSurfInt(shared_ptr<DB>, const SurfaceInteraction&);
GFE_API vector<SurfaceInteraction> GetAllSurfInt(shared_ptr<DB>);
GFE_API std::optional<SurfaceInteraction> GetSurfInt(shared_ptr<DB>, const string& name);

GFE_API void AddContPropAssign(shared_ptr<DB>, const ContPropAssign&);
GFE_API vector<ContPropAssign> GetAllContPropAssign(shared_ptr<DB>);
GFE_API vector<ContPropAssign_Partition> GetAllContPropAssign_Partition(shared_ptr<DB>);
}

