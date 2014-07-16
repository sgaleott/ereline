#include "project_tod.hpp"
#include "logging.hpp"
#include "misc.hpp"

#include <algorithm>
#include <stdexcept>

extern "C" {
#include "chealpix.h"
}

////////////////////////////////////////////////////////////////////////////////

void
compress_map(const std::vector<double> & sum_map,
	     const std::vector<size_t> & hit_map,
	     size_t num_of_seen_pixels,
	     Projected_tod_t & proj_tod)
{
    proj_tod.pixel_idx.resize(num_of_seen_pixels);
    proj_tod.sum.resize(num_of_seen_pixels);
    proj_tod.hits.resize(num_of_seen_pixels);

    if(num_of_seen_pixels == 0)
	return;

    size_t proj_idx = 0;
    for(size_t idx = 0; idx < sum_map.size(); ++idx) {
	if(hit_map[idx] == 0)
	    continue;

	proj_tod.pixel_idx[proj_idx] = idx;
	proj_tod.sum[proj_idx] = sum_map[idx];
	proj_tod.hits[proj_idx] = hit_map[idx];
	++proj_idx;
    }
}

////////////////////////////////////////////////////////////////////////////////

/* The variable "proj_tod" must have its "nside" and "pointing_id"
 * fields already initialized before calling this function (the latter
 * field is only used for diagnostic messages). */
void project_tod(int quality_flag,
		 const Range_t<uint64_t> & obt_range,
		 const PointingData & pnt_data,
		 const DifferencedData & diff_data,
		 Range_t<double> & extrema,
		 Projected_tod_t & proj_tod)
{
    int start_idx, end_idx;
    find_range_of_indexes(pnt_data.obt_time,
			  obt_range.start, obt_range.end,
			  start_idx, end_idx);
    if(start_idx < 0 || end_idx < 0 || end_idx < start_idx) {
	auto msg = boost::format("no data found for pointing ID %1%") 
	    % proj_tod.pointing_id;
	Logger::get_instance()->error(msg);
	throw std::runtime_error(msg.str());
    }

    extrema.start = extrema.end =-1.6e+30;

    // Project the TOD on a (nested) Healpix map
    const int num_of_pixels = 12 * proj_tod.nside * proj_tod.nside;
    std::vector<size_t> hit_map(num_of_pixels, 0);
    std::vector<double> sum_map(num_of_pixels, 0.0);
    size_t num_of_seen_pixels = 0;
    for(int idx = start_idx; idx <= end_idx; ++idx) {
	if((diff_data.flags[idx] & quality_flag) != 0)
	    continue;

	long pixel_num = 0;
	ang2pix_nest(proj_tod.nside, 
		     pnt_data.theta[idx], 
		     pnt_data.phi[idx], 
		     &pixel_num);

	auto & pixel_sum = sum_map.at(pixel_num);
	auto datum = diff_data.sky_load.at(idx);
	pixel_sum += datum;

	// Keep track of the minimum/maximum value seen in "datum"
	if(extrema.start < -1e+30 || datum < extrema.start)
	    extrema.start = datum;
	if(extrema.end < -1e+30 || datum > extrema.end)
	    extrema.end = datum;

	auto & hits = hit_map.at(pixel_num);
	++hits;
	if(hits == 1)
	    ++num_of_seen_pixels;
    }

    // Compress the map to save space
    compress_map(sum_map, hit_map, num_of_seen_pixels, proj_tod);
}
