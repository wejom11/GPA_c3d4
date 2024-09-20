#pragma once
#include "Vibration/EERASolver.h"

//class ERASolver;
class GFE_API ARSolver {
public:
    struct ARInput {
        void Clear();

        double StrcTopD = NAN;
        double StrcBotD = NAN;
        bool UpdateERA = false;
        std::shared_ptr<ERASolver> eraSolver;
    } Input_;
    struct Output {
        void Clear();

        double Time = NAN;      // 结构顶底点位移差最大时刻
        GFE::Function2 StrcTopU, StrcBotU;
        GFE::Function2 A,V,U,S,A_S;

        bool IsEmpty = true;
        std::string ErrStr = "";
    } Output_;

    ARSolver();
    bool Calc();

    //! 临时

    //! 辅助数据结构, 因为GFE::Material的结构复杂，取数据不方便
    struct MatInfo
    {
        double Density = NAN;
        std::string Name;
    };
    std::vector<std::shared_ptr<MatInfo>> GetMatInfo(const std::vector<GFE::Material>& input);
    void Densify(std::vector<std::shared_ptr<MatInfo>>& inMat, const std::vector<double>& depths, const std::vector<double>& acD);
};

