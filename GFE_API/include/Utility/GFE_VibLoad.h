#pragma once
#include "../api_a.h"
#include "../Vibration/ERASolver.h"

namespace GFE {

/**
 * @return ErrCode
 * 0: success;
 * 1: no db;
 * 2: no vibload;
 * 3: output path err;
 * -1: unknown err;
 */
GFE_API int GetERASolversFromDB(std::shared_ptr<GFE::DB> db, sp3<ERASolver>& solvers);

}
