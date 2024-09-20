#pragma once
#include "../api_a.h"
#include <optional>
#ifdef __linux
#include <math.h>
#endif

namespace GFE
{
namespace Step
{
// 空数据结构，暂时没接具体数据
struct Visco {
    double init_inc = 0;
    double period = 0;
    double min_inc = 0;
    double max_inc = 0;
    double allsdtol = NAN;
    double cetol = NAN;
    double stabilize = NAN;
};
GFE_API void SetVisco(const std::shared_ptr<DB>&, const Visco&);
GFE_API std::optional<Visco> GetVisco(const std::shared_ptr<DB>&);

struct MassScaling {
    int frequency = 100;
    double dt = 5e-5;
};
GFE_API void SetMassScaling(const std::shared_ptr<DB>&, const MassScaling&);
GFE_API std::optional<MassScaling> GetMassScaling(const std::shared_ptr<DB>&);

struct Static {
    double init_inc = 0;
    double period = 0;
    double min_inc = 0;
    double max_inc = 0;
};
GFE_API void SetStatic(const std::shared_ptr<DB>&, const Static&);
GFE_API std::optional<Static> GetStatic(const std::shared_ptr<DB>&);

struct Geostatic {
    double init_inc = 0;
    double period = 0;
    double min_inc = 0;
    double max_inc = 0;
};
GFE_API void SetGeostatic(const std::shared_ptr<DB>&, const Geostatic&);
GFE_API void AddGeostaticIdx(const std::shared_ptr<DB>&, int);
GFE_API std::optional<Geostatic> GetGeostatic(const std::shared_ptr<DB>&);
GFE_API std::vector<int> GetGeostaticIdx(const std::shared_ptr<DB>&);

struct Frequency {
    int eigen_num = NAN;
    double min_freq = NAN;
    double max_freq = NAN;
};
GFE_API void SetFrequency(const std::shared_ptr<DB>&, const Frequency&);
GFE_API std::optional<Frequency> GetFrequency(const std::shared_ptr<DB>&);


struct Dynamic {
    bool direct = false;
    bool explicit_ = true;
    double init_inc = 0;
    double period = 0;
    double min_inc = 0;
    double max_inc = 0;
};
GFE_API void SetDynaminc(const std::shared_ptr<DB>&, const Dynamic&);
GFE_API std::optional<Dynamic> GetDynamic(const std::shared_ptr<DB>&);

struct Soils {

	double cetol = 0;
	double utol = 0;
	int end = 0;
	bool is_consolidation = false;

	double init_inc;
	double period;
	double min_inc;
	double max_inc;
};
GFE_API void SetSoils(const std::shared_ptr<DB>&, const Soils&);
GFE_API std::optional<Soils> GetSoils(const std::shared_ptr<DB>&);

enum eAnalysisType {
    T_Unknown,
    T_Freq,
    T_Static,
    T_DynaExp,
    T_ResponseSpectrum,
    T_ModalDyna,
    T_SteadyStateDynamic,
    T_DynaImp,
    T_Visco,
    T_GeoStatic,
    T_Soil
};

GFE_API int GetAnalysisType(std::shared_ptr<DB>);

}
}



