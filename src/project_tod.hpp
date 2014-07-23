#ifndef PROJECT_TOD_HPP
#define PROJECT_TOD_HPP

#include "datatypes.hpp"

#include <cstdlib>
#include <vector>
#include <stdexcept>

struct Projected_tod_t {
    int pointing_id;
    int nside;
    std::vector<size_t> pixel_idx;
    std::vector<double> sum;
    std::vector<size_t> hits;

    // This operator can only be used with two Projected_tod_t objects
    // that have followed the same scanning strategy (i.e., all the
    // indexes in the two pixel_idx must be equal -- this is not
    // enforced here).
    Projected_tod_t & operator+=(Projected_tod_t & proj_tod) {
        if(pointing_id != proj_tod.pointing_id
           || nside != proj_tod.nside
           || pixel_idx.size() != proj_tod.pixel_idx.size())
            throw std::domain_error("The two projected TODs are not compatible "
                                    "for a += operation");

        for(size_t idx = 0; idx < pixel_idx.size(); ++idx) {
            sum[idx] += proj_tod.sum.at(idx);
            hits[idx] += proj_tod.hits.at(idx);
        }
    }
};

void project_tod(int quality_flag,
                 const Range_t<uint64_t> & obt_range,
                 const PointingData & pnt_data,
                 const DifferencedData & diff_data,
                 Range_t<double> & extrema,
                 Projected_tod_t & proj_tod);

#endif

