#ifndef DEPRECATED_H
#define DEPRECATED_H

#include <GFE_API.h>
/**
 * @brief getElementFaceNode
 * @param type
 * @param face_id
 * @param nl = flase: 不包含中间节点
 * @return
 */
namespace GFE {
GFE_API vector<int> getElementFaceNode(int type, int face_id, bool nl = false);
GFE_API const vector<short>& getElementFaceNum();
GFE_API const std::map<int, vector<vector<short>>>& getElmFaceNodeMap(bool nl = false);
}


#endif // DEPRECATED_H
