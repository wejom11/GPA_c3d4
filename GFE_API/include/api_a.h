#pragma once

/**
 * @brief Generic
 */

#ifdef _MSC_VER
#if defined(GFE_EXPORT)
#define GFE_API __declspec(dllexport)
#else
#define GFE_API __declspec(dllimport)
#endif
#elif defined __GNUC__
#define GFE_API __attribute__((__visibility__("default")))
#endif

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

namespace GFE {
using std::vector;
using std::string;
using id_t = int;
using data_t = float;
class DB;
using std::shared_ptr;
}
