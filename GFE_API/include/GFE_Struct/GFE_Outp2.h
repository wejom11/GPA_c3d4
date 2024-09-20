#pragma once

#include <api_a.h>
#ifdef __linux
#include <cmath>
#include <algorithm>
#endif
#include <GFE_Struct/GFE_Mesh.h>
#include <GFE_API_2.h>

//! @header
//! 考虑多积分点

namespace GFE
{

//! @brief 场输出
namespace FO
{

//! @brief 得到节点上的结果数据，如果原始数据是在积分点上，则会先经过外插
//! @return
//! {结果数组，{最小值，最大值}}。注意：因为取平均算法特殊，一个节点上不一定只有一个值，所以结果数组的最值范围并不准确，要另外返回
GFE_API GridData GetData_Node(shared_ptr<DB>, const string& variable, int frame, bool inc_sm = false);

//! @brief 得到节点上的结果包络数据
GFE_API GridData GetData_Node_Envelope(shared_ptr<DB>, const string& variable, const vector<int>& frames,
                                       const string& envelope_type, bool inc_sm = false);

//! @brief 得到包络结果（不进行外推）
GFE_API vector<data_t> GetData_Envelope(shared_ptr<DB>, const string& variable, const vector<int>& frames,
                                        const string& envelope_type);

//! @brief 得到单元上的结果数据。不变量在API内进行计算
GFE_API vector<data_t> GetData2(shared_ptr<DB>, const string& variable, int frame);

GFE_API vector<GFE::FO::Data> GetData2(shared_ptr<DB>, const string& variable, bool fill_empty_frame,
                                       const vector<int>& frames, const vector<int>& subset);

//! @brief 设置后处理平均阈值
GFE_API void SetAverageThreshold(double);
//! @brief 获取后处理平均阈值。默认75%
GFE_API double GetAverageThreshold();

//! @brief 更新多积分点数据
GFE_API void SetIntegData(shared_ptr<DB>, int frame, const string& var, const vector<data_t>& data);
//! @brief 更新多积分点数据，带region的版本，odb转db用到
GFE_API void SetIntegData(shared_ptr<DB>, int frame, const string& var, const string& region, const vector<data_t>& data);
//! @brief 更新多积分点数据，带region的版本，odb转db用到
GFE_API void SetIntegData2(shared_ptr<DB>, int frame, const string& var, const string& region, const vector<data_t>& data);


}

namespace HO {
//! @brief 临时函数
//! 获取历史输出，如果variable是不变量，则通过分量进行计算
//! 关于不变量，目前仅硬编码单独处理了e,le,s的principal
GFE_API vector<HistoryOutput> GetData_Test(shared_ptr<DB>, const string& region, const string& variable);

GFE_API void AddFilter(shared_ptr<DB>, const OutpFilter&);
GFE_API std::optional<OutpFilter> GetFilter(shared_ptr<DB>, const string& name);
}
}
