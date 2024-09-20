#pragma once
#include <Math/GFE_Interp.h>
#include "ERASolver.h"

class GFE_API EERASolver : public ERASolver
{
public:
    EERASolver();
    void Init();
    virtual bool Calc() override;
    void ReverseOutput();
    //! @brief 传给求解器需要做的步骤, 包括:
    //! 1.将结果改为从低到高排序
    //! 2.删去观测点的时程数据
    void ToGFEX();

    virtual ERASolver::Param* GetBaseParam() override { return &Param; }
    virtual ERASolver::Input* GetBaseInput() override { return &Input; }
    virtual ERASolver::Output* GetBaseOutput() override { return &Result; }

    struct EERAParam : public ERASolver::Param {
        bool WaveType = 0;          // 0 S波 1 P波
        int N = 4096;               // 地震波插值点个数, 因为是基于均匀采样做FFT
        int MaxIter = 10;           // 最大迭代次数
        double Tol = 1e-7;          // 达到收敛的最小容差
        double Rr = 0.5;            // 有效应变系数
        double DefKsi = NAN;
        double DampConvOrder = 0;       // 瑞利阻尼转阻尼比所用的模态阶数
        double DampScale = 1;

        std::map<std::string, std::string> Others;
        // Detrend_LeastSquare: 使用最小二乘法对积分结果消除趋势项

    } Param;

    struct EERAInput : public ERASolver::Input {
        bool IsOutcrop = true;                  // Outcrop=true时, A_Layer无效
        int A_Layer = 1;                        // 输入加速度所在层, 从高到低, 1表示地表. 暂时只支持地表跟基岩处
        GFE::Material BrMat;                   // 基岩材料
        GFE::Function A;                        // 输入加速度, (1, N)
    } Input;

    struct GFE_API EERAResult : public ERASolver::Output
    {
        virtual void Clear() override;

        std::string ErrStr;

        GFE::Function A_oc, V_oc, U_oc;
        std::vector<GFE::Interp<double>> GRatio, KXI, MaxE, MaxS, MaxU, MaxV, MaxA;
        std::vector<GFE::Interp<double>> MaxU_R, MaxV_R, MaxA_R;         // 相对底部
        std::vector<double> ErrG, ErrKXI;
        // 土体模态
        std::vector<std::vector<double>> SoilModal;
        std::vector<double> SoilModalFrq;
        // 材料等效线性化
        std::vector<double> ERatio, DampA;
        std::vector<GFE::Material> EqMats;

        friend class EERASolver;
    } Result;

protected:
    struct MatInfo {
        enum MatType {
            UnkownNoElastic = -2,
            Other,
            Elastic = 0,
            Davidenkov
        };
        std::string Name;
        int Type = -2;
        std::string TypeName;
        double G;
        double aR;
        double Den;
        int idOfInMats = -1;
        GFE::Interp<double> mat_test_g, mat_test_damp;
    };

    static std::vector<std::shared_ptr<MatInfo>> GetMatInfo(const std::vector<GFE::Material>&, bool isPWave);

    static std::tuple<std::vector<std::vector<double>>, std::vector<double>>
    Calc_SoilModal(const std::vector<double> depths, const std::vector<std::shared_ptr<MatInfo>>& matInfos);

    // 计算土体等效线性化参数
    void Calc_MatEL(const std::vector<std::shared_ptr<EERASolver::MatInfo>>& matInfos);
};

