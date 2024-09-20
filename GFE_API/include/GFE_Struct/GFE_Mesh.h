#pragma once
#include <GFE_API_2.h>
#include <api_a.h>
#include <array>
#ifdef __linux
#include <GFE_API.h>
#include <algorithm>
#include <cmath>
#endif

namespace GFE
{
using std::array;
using std::vector;

//! @brief 每种单元子类型的名称，下标i的元素对应类型i
GFE_API vector<string> GetElementSubtypeName();

//! @brief 每种单元类型的节点数，下标i的元素对应类型i
GFE_API vector<int> GetElementTypeNodeNum();

//! @brief 下标i的元素对应id为i的单元的类型
GFE_API vector<int> GetElementType(shared_ptr<DB>);

//! @brief 下标i的元素对应id为i的单元的子类型
GFE_API vector<int> GetElementSubtype(shared_ptr<DB>);

//! @brief 下标i的元素对应id为i的单元的节点
//! 这个函数相对耗时，次数较多的循环里，建议不要重复调用
[[deprecated("use ModelInfo::elem2node instead")]] GFE_API vector<array<int, 10>> GetElementNode(shared_ptr<DB>);

/**
 * @brief GetElemIntgPtNum
 * @return 各单元的积分点个数，长度等于单元数
 */
GFE_API vector<int> GetElemIntgPtNum(shared_ptr<DB>);

struct ModelInfo;

//! @brief 可以存每个单元的各节点值，也可以存每个节点的各单元值
class GFE_API GridData
{
public:
    GridData()
    {
        Clear();
    }

    //! @brief 以一个数组初始化GridData, 每个Gird只有一个数据，Neighbor是它自身
    //! 拷贝构造版本
    GridData(const std::vector<data_t>& data, int n_grid);

    //! @brief 以一个数组初始化GridData, 每个Gird只有一个数据，Neighbor是它自身
    //! 移动构造版本
    GridData(std::vector<data_t>&& data, int n_grid);

    //! @brief 插入1个单元/节点的数据
    virtual void Append(int grid_id, const data_t* data, const int* neighbor_id, size_t length)
    {
        _data.insert(_data.end(), data, data + length);
        _neighbor_id.insert(_neighbor_id.end(), neighbor_id, neighbor_id + length);
        _row_offset.push_back(_data.size());
        _grid_id.push_back(grid_id);
    }

    //! @brief 传入下标，获取1个单元/节点的数据
    std::pair<const data_t*, int> GetData(size_t pos) const
    {
        auto begin = _row_offset.at(pos);
        auto end = _row_offset.at(pos + 1);
        return {_data.data() + begin, end - begin};
    }

    //! @brief 传入下标，获取1个单元/节点的数据
    std::pair<data_t*, int> GetData(size_t pos)
    {
        auto begin = _row_offset.at(pos);
        auto end = _row_offset.at(pos + 1);
        return {_data.data() + begin, end - begin};
    }

    //! @brief 传入下标，获取1个单元/节点的相邻节点/单元
    std::pair<const int*, int> GetNeighbor(size_t pos) const
    {
        auto begin = _row_offset.at(pos);
        auto end = _row_offset.at(pos + 1);
        return {_neighbor_id.data() + begin, end - begin};
    }

    //! @brief 传入下标，获取1个单元/节点的相邻节点/单元
    std::pair<int*, int> GetNeighbor(size_t pos)
    {
        auto begin = _row_offset.at(pos);
        auto end = _row_offset.at(pos + 1);
        return {_neighbor_id.data() + begin, end - begin};
    }

    //! @brief 预先分配空间，配合SetData可以更高效地赋值
    void Allocate(const vector<std::size_t>& row_offset)
    {
        if (row_offset.empty() || row_offset[0] != 0)
            throw std::runtime_error("GFE::GridData: Invalid row offset");
        _row_offset = row_offset;
        _data.resize(row_offset.back(), NAN);
        _neighbor_id.resize(row_offset.back(), -1);
        _grid_id.resize(row_offset.size() - 1, -1);
    }

    //! @brief 配合Allocate可以更高效地赋值
    void SetData(size_t pos, int grid_id, const data_t* data, const int* neighbor_id, size_t length)
    {
        size_t offset = _row_offset.at(pos);
        std::copy_n(data, length, _data.begin() + offset);
        std::copy_n(neighbor_id, length, _neighbor_id.begin() + offset);
        _grid_id.at(pos) = grid_id;
    }

    int GetGridID(size_t pos) const
    {
        return _grid_id.at(pos);
    }

    bool IsEmpty() const
    {
        return _data.empty();
    }
    virtual void Clear()
    {
        _data.clear();
        _neighbor_id.clear();
        _row_offset = {0};
        _grid_id.clear();
    }
    int GetGridSize() const
    {
        return _row_offset.size() - 1;
    }

    //! @brief 将网格信息进行翻转，比如原来是记录每个单元各节点值，转换为记录每个节点的各单元值
    GridData Reverse() const;

    //! @brief 与compare做比较，取包络。注意：会修改当前对象
    void CalcEnvelope(const GridData& compare, const char* oper_type);

    const vector<data_t>& RawData() const
    {
        return _data;
    }

    //! @return {最小值pos，最大值pos，最小值，最大值}
    std::tuple<int, int, data_t, data_t> GetRange() const;

    //! @param use_neighbor_id = true，则selection里的是neighbor id，否则是grid id
    //! @return {最小值pos，最大值pos，最小值，最大值}
    std::tuple<int, int, data_t, data_t> GetRange(const std::vector<int>& selection, bool use_neighbor_id) const
    {
        auto r = GetRange2(selection, use_neighbor_id);
        return {r.first[0], r.first[1], r.second[0], r.second[1]};
    }

    //! @param use_neighbor_id = true，则selection里的是neighbor id，否则是grid id
    //! @return {最小值pos，最大值pos，绝对值最小值pos，绝对值最大值pos} {最小值，最大值，绝对值最小值，绝对值最大值}
    std::pair<array<int, 4>, array<data_t, 4>> GetRange2(const std::vector<int>& selection, bool use_neighbor_id) const;

    //! @param fill 给定总长度，grid data中不存在数据的下标填NAN
    //! @return 返回一个长度等于fill的数组，每个位置上取那个grid的第一位非NAN数据
    std::vector<data_t> Flatten(int fill) const
    {
        vector<data_t> result;
        result.resize(fill, NAN);
        for (size_t i = 0; i < _grid_id.size(); ++i)
        {

            auto begin = _row_offset.at(i);
            auto end = _row_offset.at(i + 1);
            auto values = _data.data() + begin;
            auto length = end - begin;

            if (length == 0)
                throw std::exception();

            result[_grid_id.at(i)] = values[0];
            for (size_t j = 0; j < length; ++j)
            {
                if (!std::isnan(values[j]))
                    result[_grid_id.at(i)] = values[j];
            }
        }
        return result;
    }

    //! @param fill 给定总长度，grid data中不存在数据的下标填NAN
    //! @param selection 填neighbor id
    //! @return 返回一个长度等于fill的数组，每个位置上优先取selection中的neighbor的位置
    std::vector<data_t> Flatten(const std::vector<int>& selection, int fill) const;

    //! @todo 目前Average里集成了Reverse这一步，后面应该分离出来
    GridData Average(double threshold, const ModelInfo& model_info) const;

    //! @brief
    //! 返回连通域，如果this是单元-节点，reverse是节点-单元，返回的就是单元的最大连通域。反过来就是节点的最大连通域
    static vector<vector<int>> GetConnectedDoamin(const GridData& grid, const GridData& reverse);

    //! @attention 返回某个单元的所有相邻单元
    static std::vector<int> GetNeighborElemIds(int elem_id, const GFE::GridData& elem2node,
                                               const GFE::GridData& node2elem, const std::vector<int>& elem_grid_id2pos,
                                               const std::vector<int>& node_grid_id2pos);

    //! @attention 返回GridData的id到pos的映射
    std::vector<int> GetID2Pos();

    //! @attention 膜壳单元划分平均区域
    //! @param region_sm_elem_ids 截面属性区域中的膜壳单元ID数组
    static std::vector<std::vector<int>> ShellMembRigionDivision(
        std::shared_ptr<GFE::DB> db, GFE::ModelInfo& model_info, const std::vector<GFE::id_t>& region_sm_elem_ids,
        const std::vector<GFE::Element>& gfe_elems, const std::vector<GFE::data_t>& gfe_nodes_coord,
        const std::vector<int>& elem_grid_id2pos, const std::vector<int>& node_grid_id2pos);

    //! @attention 包含膜壳单元取平均
    //! @param inc_sm 是否包含膜壳单元特征边
    GridData Average(double threshold, const ModelInfo& model_info, bool inc_sm) const;

protected:
    vector<data_t> _data;
    vector<int> _neighbor_id;
    vector<std::size_t> _row_offset;
    vector<int> _grid_id;
};

struct ModelInfo
{
    int n_node;
    int n_elem;
    vector<int> intg_ee;               // 各单元积分点个数
    vector<int> type;                  // 各单元类型
    vector<int> subtype;               // 各单元子类型
    vector<std::array<int, 10>> enode; // 单元节点数组
    vector<int> element_type_nnode;    // 各单元类型对应的节点个数

    GridData elem2node, node2elem;
    vector<vector<int>> average_region;    // 平均区域，用户后处理，不考虑膜壳单元特征边
    vector<vector<int>> sm_average_region; // 膜壳单元平均区域
    vector<vector<int>> s_average_region;  // 实体单元平均区域
};

GFE_API const ModelInfo& GetModelInfo(const shared_ptr<DB>&);

//! @brief 有时要多次调用GFE::open创建DB对象，避免ModelInfo重复计算，可以调用此函数设置DB的ModelInfo成员
GFE_API void SetModelInfo(const shared_ptr<DB>&, const ModelInfo&);
} // namespace GFE
