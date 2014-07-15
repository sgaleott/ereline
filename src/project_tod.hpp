#ifndef PROJECT_TOD_HPP
#define PROJECT_TOD_HPP

#include "datatypes.hpp"

#include <cstdlib>
#include <vector>

struct Projected_tod_t {
    int pointing_id;
    int nside;
    std::vector<size_t> pixel_idx;
    std::vector<double> value;
    std::vector<size_t> hits;
};

void project_tod(int quality_flag,
		 const Range_t<uint64_t> & obt_range,
		 const PointingData & pnt_data,
		 const DifferencedData & diff_data,
		 Projected_tod_t & proj_tod);

#endif

