#pragma once
#include "../api_a.h"
#include <array>
#include <optional>

namespace GFE
{
// 这个数据实际存在Parameter表中，这里的接口将数据转存到该数据结构中
struct ModalDamping
{
    enum eType {
        Direct,
        Composite,
        Rayleigh,
        Structural
    };
    enum eDef {
        Mode,
        Freq
    };
    enum eField {
        All,
        Mechanical,
        Acoustic
    };

    int type = Direct;
    int definition = Mode;
    int field = All;
    vector<std::array<double, 4>> data;
};

GFE_API void SetModalDamping(shared_ptr<DB>, const ModalDamping&);
GFE_API std::optional<ModalDamping> GetModalDamping(shared_ptr<DB>);

struct GlobalDamping
{
    enum eField {
        All,
        Mechanical,
        Acoustic
    };

    double alpha = 0;
    double beta = 0;
    int field = All;
    double structual = 0;
};

GFE_API void SetGlobalDamping(shared_ptr<DB>, const GlobalDamping&);
GFE_API std::optional<GlobalDamping> GetGlobalDamping(shared_ptr<DB>);

struct ModelChange
{
    int addOrRemove = 0;    // 0 无效; 1 add; 2 remove
    vector<string> data;
};

struct Spectrum
{
    string name;
    int type = 3;       // 1 位移; 2 速度; 3 加速度
    vector<double> data;    // {mag1, frq1, damp1, mag2, frq2, ...}
};

GFE_API void SetSpectrum(shared_ptr<DB>, const Spectrum&);
GFE_API std::vector<Spectrum> GetAllSpectrum(shared_ptr<DB>);

struct ResponseSpectrum
{
    enum SumType {
        Unknown,
        ABS,
        CQC,
        DSC,
        GRP,
        NRL,
        SRSS,
        TENP
    };

    static string SumTypeStr(int i) {
        static vector<string> strMap = {"Unknown", "ABS", "CQC", "DSC", "GRP", "NRL", "SRSS", "TENP"};
        if(i >= 0 && i < strMap.size()) return strMap[i];
        return "Unknown";
    }

    string name;
    int sum = ABS;
    string spectrum;
    std::vector<double> data;       // {xcos, ycos, zcos, factor, 后续扩展, ...}
};

GFE_API void SetResponseSpectrum(shared_ptr<DB>, const ResponseSpectrum&);
GFE_API std::optional<ResponseSpectrum> GetResponseSpectrum(shared_ptr<DB>);

struct ModalDynamic {
	bool cont = false;		// continue
	double time_increment = 0;
	double time_period = 0;
};

GFE_API void SetModalDynamic(const std::shared_ptr<DB>&, const ModalDynamic&);
GFE_API std::optional<ModalDynamic> GetModalDynamic(const std::shared_ptr<DB>&);

struct GFE_API SteadyStateDynamic {
    enum IntervalEnum {
        EIGENFREQUENCY,
        RANGE,
        SPREAD
    };
    enum ScaleEnum {
        LOG,
        LINEAR
    };

    bool direct = false;
    int interval = EIGENFREQUENCY;
    int scale = LOG;
    vector<std::array<double, 6>> data;      // 对应Inp里的每一行，一行最多6个数

    vector<double> GetSinglePoints(const std::vector<float>& modFrq);
};

GFE_API void SetSteadyStateDynamic(const std::shared_ptr<DB>&, const SteadyStateDynamic&);
GFE_API std::optional<SteadyStateDynamic> GetSteadyStateDynamic(const std::shared_ptr<DB>&);

}

