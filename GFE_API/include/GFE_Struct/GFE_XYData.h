#pragma once
#include "../api_a.h"
#include <optional>

namespace GFE
{
struct GFE_API XYData {
    string name, xlabel, ylabel;
    vector<double> x, y;
    vector<double> x_img, y_img;

    void ToTxt(std::ostream&, bool noHeader) const;
//    static std::tuple<XYData, int> FromTxt(std::istream&);    // 未完成
};

GFE_API std::optional<XYData> GetXYData(shared_ptr<DB>, const string& name);
GFE_API void SetXYData(shared_ptr<DB>, const XYData&);
GFE_API vector<string> GetAllXYDataName(shared_ptr<DB>);
GFE_API void RemoveXYData(shared_ptr<DB>, const string& name);
}
