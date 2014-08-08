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

class Planck_velocity_t;
class Sqlite_connection_t;
class Gain_table_t;
class Lfi_radiometer_t;
struct Dipole_fit_t;

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

void load_convolution_params(Sqlite_connection_t & ucds,
                             const Lfi_radiometer_t & radiometer,
                             Planck_velocity_t & vel);

void load_pointing_information(Sqlite_connection_t & ucds,
                             int first_od,
                             int last_od,
                             std::vector<Pointing_t> & pointings);
// Wrapper around the previous definition, with the assumption that
// first_od == last_od == od.
void load_pointing_information(Sqlite_connection_t & ucds,
                               int od,
                               std::vector<Pointing_t> & pointings);

void save_gain_table(FitsObject & file,
                     const Lfi_radiometer_t & radiometer,
                     const Gain_table_t & gain_table,
                     const std::string & comment = "");

void save_gain_table(const std::string & file_name,
                     const Lfi_radiometer_t & radiometer,
                     const Gain_table_t & gain_table,
                     const std::string & comment = "");

void save_tod(const std::string & file_name,
              signed short od,
              const Lfi_radiometer_t & radiometer,
              const DifferencedData & datadiff,
              const std::string & comment = "");

void load_pointings(const std::string & file_name,
                    PointingData & pointings);

void load_differenced_data(const std::string & file_name,
                           DifferencedData & datadiff);

void load_subsampled_ref_and_hk_data(Sqlite_connection_t & ucds,
                                     const Lfi_radiometer_t & radiometer,
                                     const std::string & hk_sensor_name,
                                     const Range_t<int> & pointing_range,
                                     std::vector<int> & pointing_id,
                                     std::vector<double> & ref_data,
                                     std::vector<double> & hk_data);

#endif
