#pragma once
#include "obj_st.h"
#include "GFE_Preset.h"
#include <map>
/**
 * @brief Function of pre-defined settings.
 * e.g. variables' name, elements' node order
 */

namespace GFE {

GFE_API vector<Variable_> getAllVariables();
/**
 * @brief getQuantityDCRP
 * @return Description of quantity inputted
 */
GFE_API string getQuantityDCRP(const string&);
GFE_API int getQuantityType(const string&);
GFE_API int getVariableType(const string&);
GFE_API int getQuantityType2(const string&);
GFE_API int getVariableType2(const string&);
GFE_API vector<string> getComponentByQuantity(const string& quantity);   // deprecated, use getVariable
GFE_API vector<string> getVariable(const string& q_or_c);


GFE_API string getElementTypeName(int id);
GFE_API int getElementTypeId(const string& name);
GFE_API string getElementSubTypeName(int id);
GFE_API const vector<string> getElementSubTypeName();
GFE_API int getElementSubTypeId(const string& name);
GFE_API const vector<int>& getElementSubTypeParent();
GFE_API const vector<int>& getElementNodeNum();

GFE_API vector<string> getPreVar(GFE::VT2, int field_or_hist);
GFE_API void clearStaticData();

namespace Surf {
    GFE_API const std::unordered_map<int, std::string>& LabelHash();
    GFE_API const std::unordered_map<std::string, int>& LabelHashInv();
    GFE_API const std::map<std::pair<int, int>, vector<int>>& FaceMap(bool nl = false);      // {element subtype id, face label id} -> node order
}


}
