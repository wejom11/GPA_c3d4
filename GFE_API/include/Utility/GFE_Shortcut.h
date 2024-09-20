#pragma once
#include "../api_a.h"
#include "../Vibration/EERASolver.h"
#include "../Vibration/IERASolver.h"

#include <memory>
#include <string>

namespace GFE {
class DB;
}

namespace Shortcut
{
/**
 * @return ErrCode
 * 0: success;
 * 1: no db;
 * 2: dim dir err;
 * 3: output path err;
 * -1: unknown err;
 */
GFE_API int ExportSDR(std::shared_ptr<GFE::DB>, int hDir, int dim);

/**
 * @return ErrCode
 * 0: success;
 * 1: no db;
 * 2: no vibload;
 * 3: output path err;
 * -1: unknown err;
 */
GFE_API int ExportEERAResult(std::shared_ptr<GFE::DB> db, const std::string& path = "");
}
