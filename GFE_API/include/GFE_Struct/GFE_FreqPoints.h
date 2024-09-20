#ifndef GFE_FREQPOINTS_H
#define GFE_FREQPOINTS_H

#include "../api_a.h"

namespace GFE {

// 扫频点定义，仅用于频响分析
#ifdef __linux
struct GFE_API FreqPoints {         // linux上用下面的写法会报错
#else
struct GFE_API [[deprecated]] FreqPoints {
#endif

	enum TYPE {
		EIGENFREQUENCY,
		RANGE,
        SPREAD,
	};
	enum SCALE {
		LOG,
		LINEAR,
	};
	bool direct = false;
	int type = 0;
	int scale_type = 0;
	std::vector<double> lower;
	std::vector<double> upper;
	std::vector<int> nPoint;
	std::vector<double> bias;
	std::vector<double> scale;
	std::vector<double> spread;
	std::vector<double> singlePoint;
};

}




#endif // GFE_FREQPOINTS_H
