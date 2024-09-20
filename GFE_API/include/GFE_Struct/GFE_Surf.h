#pragma once

#include <vector>
#include <array>
#include <memory>
#include <map>
#include <api_a.h>

namespace GFE {
class DB;

struct GFE_API SurfPartition
{
    string origin;
    int partition_id;
    string splited_surf;
};

namespace Surf {

//!
//! \brief GetTopoNodes
//! \param surf_name
//! \param label   true返回节点label, false返回节点id
//! \param middle_point   高阶单元, 如二阶四面体, true返回值包含中间节点, false不包含
//! \return 每个小单元面的节点, 节点顺序按右手法则定义面法向. 线surf叉乘(0,0,1)为法向
//!
GFE_API vector<std::array<int, 6>> GetTopoNodes(shared_ptr<DB> db, const std::string& surf_name, bool label = false, bool middle_point = false);

GFE_API void AddPartition(shared_ptr<DB>, const SurfPartition&);

//! \brief GetPart 获取分区后的surface
GFE_API string GetPart(shared_ptr<DB>, const string& origin, int partition_id);
//! \brief GetPart 获取分区后的surface
GFE_API std::map<int, string> GetPart(shared_ptr<DB>, const string& origin);

}
}

