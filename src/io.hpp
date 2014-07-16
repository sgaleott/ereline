#ifndef IO_HPP
#define IO_HPP

#include "ahf_info.hpp"
#include "healpix_map.hpp"
#include "fits_object.hpp"

#include <string>
#include <vector>
#include <boost/format.hpp>

#include <sqlite3.h>

class PlanckVelocity;
class SQLite3Connection;
class gainTable;
class LfiRadiometer;

template<typename T>
void load_map(const std::string & file_name, 
	      int column, 
	      Healpix::Map_t<T> & map)
{
    FitsObject file;
    size_t num_of_pixels;
    file.openTable(file_name);

    file.getKey("NAXIS2", num_of_pixels);
    map.pixels.resize(num_of_pixels);

    std::string ordering_str;
    file.getKey("ORDERING", ordering_str);

    if(ordering_str == "NESTED")
	map.ordering = Healpix::Ordering_t::NEST;
    else if(ordering_str == "RING")
	map.ordering = Healpix::Ordering_t::RING;
    else {
	auto msg = boost::format("unknown ordering for map %1%: \"%2\"")
	    % file_name % ordering_str;
	throw std::runtime_error(msg.str());
    }

    file.getColumn(column, map.pixels, 1, num_of_pixels);
}

void loadConvolutionParametersFromUCDS(SQLite3Connection & ucds,
				       const LfiRadiometer & radiometer,
				       PlanckVelocity & vel);

void saveGainTable(const std::string & file_name,
		   signed short od,
		   const LfiRadiometer & radiometer,
		   const gainTable & gain_table,
		   const std::string & comment = "");

void loadPointingInformation(SQLite3Connection & ucds,
			     int first_od,
			     int last_od,
			     std::vector<Pointing_t> & pointings);
// Wrapper around the previous definition, with the assumption that
// first_od == last_od == od.
void loadPointingInformation(SQLite3Connection & ucds,
			     int od,
			     std::vector<Pointing_t> & pointings);


#endif
