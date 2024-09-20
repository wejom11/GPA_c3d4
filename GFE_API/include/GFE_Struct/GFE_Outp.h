#pragma once
#include "../api_a.h"
#include <map>
#include <set>
#include <optional>

namespace GFE
{

enum PrsReq
{
    PR_Modal,
    PR_SSD,         // Steady State Dynamic Analysis
    PR_Static
};

struct Parameter {
    string name;
    string content;
    int step_id = -1;
};

struct OutpReq
{
    string name;
    int type;           // 0 field; 1 history; 2 envelope 目前只有field. history后续重构
    string region;
    string vars;        // 物理量名称. 以‘;’分隔
    vector<double> frq;
    int inc = -1;          // 考虑隐式分析的时候新增的变量. 当指定输出频率的方式为指定每x个增量步输出一次时用到。
    int method;         // 0 AbsMax; 1 AbsMin; 2 Max; 3 Min. 仅在type=2时生效
    std::string filter;
};

struct OutpFilter
{
    string name;
    double cutoff;
    int N;
};


GFE_API void ClearOutpReq(shared_ptr<DB>);

//! @brief 添加输出请求
//! OutpReq.name为空时会设置默认名称, 传参的name会修改为最终写入db的名称
//! 因为设置默认名称需要查询整个表, 如果需要大量插入, 请事先设置好名称
GFE_API void SetOutpReq(shared_ptr<DB>, OutpReq*);

GFE_API vector<OutpReq> GetAllReqs(shared_ptr<DB>);
GFE_API vector<int> GetIncrement(shared_ptr<DB>, int type = -1);
GFE_API vector<double> GetFrequency(shared_ptr<DB>, int type = -1);


//! @brief 获取预设的输出请求
//! 返回值的frq为空，需要手动赋值
GFE_API vector<OutpReq> GetPresReq(PrsReq);

GFE_API shared_ptr<vector<id_t>> GetRegionIDs(shared_ptr<DB>, const std::string& region_union, int vt2);
//! 无序，不去重版本
GFE_API shared_ptr<vector<id_t>> GetRegionIDs2(shared_ptr<DB>, const std::string& region_union, int vt2);


namespace FO {

struct Frame
{
    int frame;
    double frq;
    bool has_output;            // 该帧已否输出
};
struct Region
{
    string set_names;           // 支持多个集合, 以‘;’分隔
};
struct Data
{
    int frame;
    string var;
    vector<data_t> data;
    vector<data_t> data_im;
    int region_id;
    vector<data_t> data_integ;  // 多积分点数据
};

//! @brief 将输出请求转成场输出的一些前置信息，包括：
//! 1. 对所有场输出请求的时间轴进行合并、排序、去重，得到一个全局的时间轴、帧与时刻间的一一对应关系，存到FO_Frame表中
//! 2. 计算{物理量，帧} -> 输出区域的映射，存到FO_Data中，数据列为空。
//! 3. 第2步中计算的输出区域存到FO_Region中，FO_Data只保留一个region_id，避免大量的字符串的重复记录
GFE_API void ConvReq(shared_ptr<DB>);

GFE_API shared_ptr<OutpReq> GetReq(shared_ptr<DB>, const string& name);
GFE_API vector<OutpReq> GetReq(shared_ptr<DB>);

GFE_API void SetFrame(shared_ptr<DB>, const vector<double>& frq);
GFE_API void AddFrame(shared_ptr<DB>, double time, bool sta);
GFE_API void UpdateFrame(shared_ptr<DB>, int frame, bool sta = true);
GFE_API vector<double> GetFrame(shared_ptr<DB>, bool has_output);
GFE_API vector<FO::Frame> GetFrame(shared_ptr<DB>);

GFE_API std::map<string, string> GetVar2Region(shared_ptr<DB>, int frame);

//! @brief GetVar2Region_Implict
//! 接隐式分析时新增的函数
//! @param inc 表示第几个增量步
//! @param time 该增量步对应的时刻
//! @return map:var->region, region可以是多个集合的并集, 集合间以;分隔
GFE_API std::map<string, std::string> GetVar2Region(shared_ptr<DB>, int inc, double time);

GFE_API shared_ptr<vector<id_t>> GetRegionIDs(shared_ptr<DB>, int frame, const string& var);
GFE_API shared_ptr<vector<id_t>> GetRegionIDs(shared_ptr<DB>, int region_id, int vt2);

GFE_API void SetData(shared_ptr<DB>, int frame, const string& var, const vector<data_t>& data, bool isImg = false);
GFE_API void SetData(shared_ptr<DB>, int frame, const string& var, const string& region, const vector<data_t>& data, const vector<data_t>& data_im = {});
GFE_API void SetData(shared_ptr<DB>, int frame, std::string&& var, std::string&& region, vector<data_t>&& data, vector<data_t>&& data_im = {});
//! region无序版本
GFE_API void SetData2(const shared_ptr<DB>&, int frame, const string& var, const string& region, const vector<data_t>& data, const vector<data_t>& data_im = {});

//! @brief getFOData
//! @param fill_empty_frame, 对于无数据的帧补上全NAN的数据; 为假时不返回该帧
GFE_API vector<FO::Data> GetData(shared_ptr<DB>, const string& var, bool fill_empty_frame, const vector<int>& frames = {-1}, vector<id_t> subset = {-1});
inline vector<data_t> GetData(shared_ptr<DB> db, const string& var, bool fill_empty_frame, int frame) {
    auto ret = GetData(db, var, fill_empty_frame, vector<int>{frame});
    return ret.empty() ? vector<data_t>() : ret[0].data;
}
inline vector<data_t> GetImgData(shared_ptr<DB> db, const string& var, bool fill_empty_frame, int frame) {
    auto ret = GetData(db, var, fill_empty_frame, vector<int>{frame});
    return ret.empty() ? vector<data_t>() : ret[0].data_im;
}


//! @brief 只获取有输出的数据，不填NAN
GFE_API vector<FO::Data> GetPureData(shared_ptr<DB>, const string& var, const vector<int>& frames = {-1});
GFE_API vector<int> GetPureDataNum(shared_ptr<DB>, const string& var, const vector<int>& frames = {-1}, bool isImg = false);

GFE_API std::set<string> GetVars(shared_ptr<DB>, bool has_output);
}

namespace EO {
    struct Region
    {
        string set_names;           // 支持多个集合, 以‘;’分隔
    };
    struct Data
    {
        int method;
        string var;
        vector<data_t> data;
        int region_id;
        vector<int> frame;
    };

    GFE_API void ConvReq(shared_ptr<DB>);
    GFE_API shared_ptr<OutpReq> GetReq(shared_ptr<DB>, const string& name);
    GFE_API vector<OutpReq> GetReq(shared_ptr<DB>);
    GFE_API std::set<string> GetVars(shared_ptr<DB>);
    GFE_API std::vector<std::tuple<int, string>> GetMethodVars(shared_ptr<DB>);
    GFE_API void SetData(shared_ptr<DB>, int method, const string& var, const vector<data_t>& data, const vector<int>& frame);
    GFE_API void SetData(shared_ptr<DB>, int method, const string& var, const string& region, const vector<data_t>& data, const vector<int>& frame);
    GFE_API vector<EO::Data> GetAllData(shared_ptr<DB>);
    GFE_API shared_ptr<EO::Data> GetData(shared_ptr<DB>, int method, const string& var);
    GFE_API vector<std::tuple<int, string>> GetConvReq(shared_ptr<DB>);
    GFE_API shared_ptr<vector<id_t>> GetRegionIDs(shared_ptr<DB>, int method, const string& var);
    GFE_API shared_ptr<vector<id_t>> GetRegionIDs(shared_ptr<DB>, int region_id, int vt2);
    GFE_API std::string GetRegionName(shared_ptr<DB>, int region_id);
}
}
