#pragma once
#include "ERASolver.h"
#include <GFE_Struct/GFE_Material.h>
#ifdef __linux
#include <cmath>
#endif

//! 场地反应分析--积分方法
//! 不支持非线性材料
class EERASolver;
class GFE_API IERASolver : public ERASolver
{
public:
    struct Param : public ERASolver::Param {
        int WaveType = 0;               // 0 S波 1 P波
    } Param_;
    struct Input : public ERASolver::Input {
        void Clear();

        std::shared_ptr<EERASolver> EERA;       // 上游的EERA求解器，如果设置了，Param中部分参数将沿用上游设置
        GFE::Function2 A;                       // 基岩处加速度, 用于封闭算法
        GFE::Function2 V;                       // 基岩露头速度, 用于开放算法
    } Input_;
    struct Output : public ERASolver::Output{
        virtual void Clear() override;

        GFE::Function2 MaxA, MaxV, MaxU, MaxE, MaxS;
        GFE::Function2 MaxA_R, MaxV_R, MaxU_R;      // 相对底部
        std::string LastErr;

        friend class IERASolver;

    } Output_;

    IERASolver();
    virtual bool Calc() override;
    virtual ERASolver::Param* GetBaseParam() override { return &Param_; }
    virtual ERASolver::Input* GetBaseInput() override { return &Input_; }
    virtual ERASolver::Output* GetBaseOutput() override { return &Output_; }

private:
    void PrePare();
    bool Calc_CDM();            // 中心差分法，刚性底
    bool Calc_CDM_Free();       // 中心差分，自由底
    bool Calc_Newmark();        // Newmark法，刚性底，垂直入射。未实现

    using MatPtr = std::shared_ptr<GFE::Material2::Para_SSA>;
};

