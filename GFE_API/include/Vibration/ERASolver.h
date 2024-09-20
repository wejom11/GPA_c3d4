#pragma once
#include <vector>
#include <map>
#include <GFE_API.h>

template <class T>
using sp3 = std::array<std::shared_ptr<T>, 3>;

class GFE_API ERASolver {
public:
    struct Param {
        double SubLayerHeight = 1;      // 细分层层高
        double TimeInterval = 0.02;     // 输出频率
    };
    struct Input {
        virtual void Clear();

        std::vector<GFE::Material> Mats;        // 各层材料
        std::vector<double> Depths;             // 各层厚度, 从高到低
        std::map<std::string, double> ObsPts;       // 观测点
    };
    struct Output {
        virtual void Clear();
        const std::vector<double>& AcD() const { return _acD; }
        const std::vector<double>& MidAcD() const { return _midAcD; }
        const std::vector<double>& SubT() const { return _subT; }
        const std::map<std::string, int>& ObsIdx() const { return _obsIdx; }

        std::vector<GFE::Function2> A,V,U,S;
        bool IsEmpty = true;

    protected:
        //! @variable m_acD: 堆叠深度, 前面会多一个0, 长度比midAcD和subT多1
        std::vector<double> _acD;
        //! @variable m_midAcD: 层中间点深度
        std::vector<double> _midAcD;
        //! @variable m_subT: 各层厚度
        std::vector<double> _subT;
        //! @variable m_obsIdx: 观测点下标
        std::map<std::string, int> _obsIdx;

        friend class ERASolver;
    };

    ERASolver();

    virtual bool Calc() { return false; };
    virtual Param* GetBaseParam() { return nullptr; }
    virtual Input* GetBaseInput() { return nullptr; }
    virtual Output* GetBaseOutput() { return nullptr; }

protected:
    /**
     * @brief 土层细化
     * @arg inProp 附加于层的一种属性, 比如, 材料
     */
    template <class T>
    void Densify(std::vector<T>& inProp);
};

#include <Math/GFE_Accum.h>
template <class T>
void ERASolver::Densify(std::vector<T>& inProp)
{
    using std::vector;
    using std::shared_ptr;
    using std::string;

    //! @variable acOriD: 原始层深度
    auto input_ = GetBaseInput();
    auto acOriD = input_->Depths;
    GFE::FloatAcc(acOriD.data(), acOriD.size());

    //! @variable: acSubD: 细分层深度
    vector<double> acSubD;
    auto tt = acOriD.back();
    auto& ml = GetBaseParam()->SubLayerHeight;
    {
        int i = 1;
        double cc = i*ml;
        while(cc < tt) {
            acSubD.push_back(cc);
            i++;
            cc = i*ml;
        }
        if(cc-tt < ml/2)
            acSubD.push_back(tt);
        else {
            acSubD.pop_back();
            acSubD.push_back(tt);
        }
    }

    //! 插入观测点, 记录观测点位置
    //! @variable acD: 深度
    //! @variable obsIdx: 记录观测点下标
    auto output_ = GetBaseOutput();
    auto& obsIdx = output_->_obsIdx;
    obsIdx.clear();
    auto& acD = output_->_acD;
    acD.clear();
    {
        auto& obsPts = input_->ObsPts;
        using Tp2 = std::pair<double, string>;
        vector<Tp2> temp;
        temp.reserve(acSubD.size()+obsPts.size());
        for(auto d : acSubD)
            temp.emplace_back(d, "");
        string zeroObsPt;       // 因为后面会在acD前面补一个0，所以深度为0的观测点需要特殊处理，避免出现层高为0
        for (auto d : obsPts)
            if (d.second == 0)
                zeroObsPt = d.first;
            else
                temp.emplace_back(d.second, d.first);

        std::sort(temp.begin(), temp.end(), [](const Tp2& a, const Tp2& b) {
            if(std::abs(std::get<0>(a)-std::get<0>(b)) < 1e-12)
                return std::get<1>(a).size() > std::get<1>(b).size();
            else
                return std::get<0>(a) < std::get<0>(b);
        });
        temp.erase(std::unique(temp.begin(), temp.end(), [](const Tp2& a, const Tp2& b) {
                       return std::abs(std::get<0>(a)-std::get<0>(b)) < 1e-12;
                   }), temp.end());

        for (std::size_t i = 0; i < temp.size(); i++)
            if (!temp[i].second.empty())
                obsIdx[temp[i].second] = i + 1; // +1是因为后面会在acD前面补一个0
        if (!zeroObsPt.empty())
            obsIdx[zeroObsPt] = 0;

        acD.reserve(temp.size());
        for(const auto& tp2 : temp)
            acD.push_back(tp2.first);
    }

    //! @variable outProp: 各层属性
    auto nAcD = acD.size();
    vector<T> outProp(nAcD);
    {
        std::size_t j = 0;
        for(std::size_t i = 0; i < acOriD.size(); i++)
            for(; j < nAcD; j++)
                if(i != acOriD.size()-1 && acD[j] > acOriD[i]) break;
                else outProp[j] = inProp[i];
    }

    //! @variable midAcD: 层中间点深度
    auto& midAcD = output_->_midAcD;
    midAcD.clear();
    midAcD.reserve(nAcD);
    midAcD.push_back(acD.at(0)/2);
    for(std::size_t i = 0; i < nAcD-1; i++)
        midAcD.push_back((acD.at(i+1)+acD.at(i))/2);

    output_->_subT.reserve(nAcD);
    output_->_subT.push_back(acD.at(0));
    for(std::size_t i = 0; i < nAcD-1; i++)
        output_->_subT.push_back(acD.at(i+1)-acD.at(i));

    acD.insert(acD.begin(), 0);

    inProp = std::move(outProp);
}

