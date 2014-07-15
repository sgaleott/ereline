#include "project_tod.hpp"
#include "logging.hpp"
#include "misc.hpp"

#include <algorithm>
#include <stdexcept>

extern "C" {
#include "chealpix.h"
}

////////////////////////////////////////////////////////////////////////////////

/* The variable "proj_tod" must have its "nside" and "pointing_id"
 * fields already initialized before calling this function (the latter
 * field is only used for diagnostic messages). */
void project_tod(int quality_flag,
		 const Range_t<uint64_t> & obt_range,
		 const PointingData & pnt_data,
		 const DifferencedData & diff_data,
		 Projected_tod_t & proj_tod)
{
    Logger * log = Logger::get_instance();
    
    int start_idx, end_idx;
    find_range_of_indexes(pnt_data.obt_time,
			  obt_range.start, obt_range.end,
			  start_idx, end_idx);
    if(start_idx < 0 || end_idx < 0 || end_idx < start_idx) {
	auto msg = boost::format("no data found for pointing ID %1%") 
	    % proj_tod.pointing_id;
	log->error(msg);
	throw std::runtime_error(msg.str());
    }

    const int num_of_pixels = 12 * proj_tod.nside * proj_tod.nside;
    std::vector<size_t> hit_map(num_of_pixels, 0);
    std::vector<double> data_map(num_of_pixels, 0.0);
    size_t num_of_seen_pixels = 0;
    for(int idx = start_idx; idx <= end_idx; ++idx) {
	if((diff_data.flags[idx] & quality_flag) != 0)
	    continue;

	long pixel_num = 0;
	ang2pix_nest(proj_tod.nside, 
		     pnt_data.theta[idx], 
		     pnt_data.phi[idx], 
		     &pixel_num);

	auto & datum = data_map.at(pixel_num);
	datum += diff_data.sky_load.at(idx);

	auto & hits = hit_map.at(pixel_num);
	++hits;
	if(hits == 1)
	    ++num_of_seen_pixels;
    }

    proj_tod.pixel_idx.resize(num_of_seen_pixels);
    proj_tod.value.resize(num_of_seen_pixels);
    proj_tod.hits.resize(num_of_seen_pixels);
    size_t proj_idx = 0;
    for(int idx = 0; idx < data_map.size(); ++idx) {
	if(hit_map[idx] == 0)
	    continue;

	proj_tod.pixel_idx[proj_idx] = idx;
	proj_tod.value[proj_idx] = data_map[idx];
	proj_tod.hits[proj_idx] = hit_map[idx];
	++proj_idx;
    }
}
