#ifndef IO_HPP
#define IO_HPP

#include "ahf_info.hpp"
#include "healpix_map.hpp"
#include "logging.hpp"
#include "fits_object.hpp"

#include <string>
#include <vector>
#include <boost/format.hpp>

#include <sqlite3.h>

class PlanckVelocity;
class SQLite3Connection;
class gainTable;
class LfiRadiometer;
struct dipoleFit;

////////////////////////////////////////////////////////////////////////////////

/* This function returns its argument. Additionally, assuming that
 * "path" is a valid POSIX path, it creates any directory which does
 * not exist but is referred by "path".
 *
 * Example:
 *
 *     std::cout << ensure_path_exists("/usr/lib/my/test/dir/file.txt");
 *
 * Output:
 *
 *     /usr/lib/my/test/dir/file.txt
 *
 * Moreover, if the directory /usr/lib/my/test/dir/ or any of its
 * parents does not exist, it will be silently created.
 */
const std::string &
ensure_path_exists(const std::string & path);

////////////////////////////////////////////////////////////////////////////////

template<typename T>
void load_map(const std::string & file_name, 
	      int column, 
	      Healpix::Map_t<T> & map)
{
    Logger * log = Logger::get_instance();

    FitsObject file;

    log->info(boost::format("Opening Healpix FITS file %1%") % file_name);
    file.openTable(file_name);

    file.getKey("NSIDE", map.nside);
    const size_t num_of_pixels = 12 * map.nside * map.nside;

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

    log->debug(boost::format("The map has %1% pixels (NSIDE=%2%) and order %3%")
	       % num_of_pixels
	       % map.nside
	       % ordering_str);
    file.getColumn(column, map.pixels, 1, num_of_pixels);
    log->info(boost::format("%1% pixels read from %2%")
	      % num_of_pixels
	      % file_name);
}

void loadConvolutionParametersFromUCDS(SQLite3Connection & ucds,
				       const LfiRadiometer & radiometer,
				       PlanckVelocity & vel);

void loadPointingInformation(SQLite3Connection & ucds,
			     int first_od,
			     int last_od,
			     std::vector<Pointing_t> & pointings);
// Wrapper around the previous definition, with the assumption that
// first_od == last_od == od.
void loadPointingInformation(SQLite3Connection & ucds,
			     int od,
			     std::vector<Pointing_t> & pointings);

void saveGainTable(const std::string & file_name,
		   signed short od,
		   const LfiRadiometer & radiometer,
		   const gainTable & gain_table,
		   const std::string & comment = "");

void save_tod(const std::string & file_name,
	      signed short od,
	      const LfiRadiometer & radiometer,
	      const DifferencedData & datadiff,
	      const std::string & comment = "");

void save_dipole_fit(const std::string & file_name,
		     const LfiRadiometer & radiometer,
		     const dipoleFit & fit,
		     const std::string & comment = "");

#endif
