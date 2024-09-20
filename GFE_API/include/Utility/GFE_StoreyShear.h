#pragma once
#include <memory>
#include <vector>
#include <api_a.h>
#include <optional>
#include <obj_st.h>

namespace GFE {
class DB;
struct ElementSet;
}

class GFE_API GFE_StoreyShear
{
public:
    enum ErrType {
        Success,
        NoDB,
        NoVar       // 缺少必要的物理量
    };

    GFE_StoreyShear() {}

	inline void SetElementType(const std::vector<int> types) {
		_elementType = types;
	}
    //! @param dir x:0 y:1 z:2
    //! @param dim 3d:3 2d:2
    bool Perform(std::shared_ptr<GFE::DB> db, const std::vector<std::string>& setNames, int dir, int dim);
    const std::vector<double>& GetOneFrameResultA(int iFrame) { return _resultOfEachFrameA.at(iFrame); }
    const std::vector<double>& GetOneFrameResultB(int iFrame) { return _resultOfEachFrameB.at(iFrame); }
    const std::vector<double>& GetEnvelopeA() { return _resultOfEnvelopeA; }
    const std::vector<double>& GetEnvelopeB() { return _resultOfEnvelopeB; }
    std::vector<double> GetHeight();
    std::vector<int> GetStorey() { return _storey; }
    int GetLastError() { return _lastError; }

    // 记录进度的两个变量
    int Progress1 = 0;                 // 1:预处理-计算区域 2/4：计算每帧结果 3/5：计算包络
    std::pair<int,int> Progress2;      // 子过程进度, n/m

private:
    using RegionTuple = std::tuple<float, std::vector<int>, std::vector<short>>;    // 高度，单元，单元节点flag
    std::optional<RegionTuple> GetRegionOfOneStorey(std::unique_ptr<GFE::ElementSet> pSet);

    bool CalcShearForceOfOneFrame(int iFrame, int dir, std::vector<std::vector<double>>& toVec);
    void CalcEnvelop(const std::vector<std::vector<double>>& eachFrame, std::vector<double>& toVec);

    std::shared_ptr<GFE::DB> _db;
    int _dir = 2;
    int _dim = 3;
    std::vector<float> _nodeCoordZ;
    std::vector<RegionTuple> _regionOfEachStorey;
    std::vector<std::vector<double>> _resultOfEachFrameA, _resultOfEachFrameB;
    std::vector<double> _resultOfEnvelopeA, _resultOfEnvelopeB;
    int _lastError = 0;
    std::vector<int> _storey;
    std::vector<int> _elementType = {GFE::CT_LINE, GFE::CT_TRIANGLE, GFE::CT_QUAD};
};

