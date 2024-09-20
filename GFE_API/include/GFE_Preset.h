#pragma once
#include "api_a.h"
#include "obj_st.h"
#include <set>

namespace GFE {

GFE_API std::array<int, 3> ApiVersion();
GFE_API std::string ApiVerStr();

//! @brief 物理量相关预设
namespace Var {

GFE_API vector<GFE::Variable_> GetAll();
GFE_API std::set<string> GetAllStr(bool no_component);

//! @brief 物理量符号描述
GFE_API string GetDCRP(const string&);

GFE_API int GetType(const string&);
GFE_API int GetType2(const string&);

//! @brief 输入主量返回所有分量; 输入分量返回自身
GFE_API vector<string> Extend(const string&);

}

}

