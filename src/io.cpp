#include "io.hpp"
#include "data_binning.hpp"
#include "dipole_fit.hpp"
#include "planck_velocity.hpp"
#include "fits_object.hpp"
#include "gain_table.hpp"
#include "logging.hpp"
#include "datatypes.hpp"
#include "sqlite3xx.hpp"
#include "squeezer.hpp"

#include <sstream>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

////////////////////////////////////////////////////////////////////////////////

const std::string &
ensure_path_exists(const std::string & path)
{
    boost::filesystem::path p(path);
    boost::filesystem::path dirname = p.parent_path();

    boost::system::error_code err;
    boost::filesystem::create_directories(dirname, err);
    if(err.value() != 0) {
        throw std::runtime_error(err.message());
    }

    return path;
}

////////////////////////////////////////////////////////////////////////////////

void
load_convolution_params(Sqlite_connection_t & ucds,
                        const Lfi_radiometer_t & radiometer,
                        Planck_velocity_t & vel)
{
    Logger * log = Logger::get_instance();

    std::stringstream query;
    query << "\n\tSELECT m100, m010, m001, m200, m110, m101, m020, m011, m002"
          << "\n\t  FROM convolved_dipole"
          << "\n\t WHERE horn = " << radiometer.horn
          << "\n\t   AND radiometer = " << radiometer.radiometer
          << "\n\t";

    Sqlite_statement_t statement(ucds, query.str().c_str());

    int result = statement.step();
    if(result != SQLITE_ROW) {
        auto error_msg =
            boost::format("No match for radiometer %1% in the "
                          "UCDS table \"convolved_dipole\" (error code: %2%)")
            % radiometer.shortName()
            % result;

        log->error(error_msg);
        throw IoError(error_msg.str());
    }

    vel.M100 = statement.column_double(0);
    vel.M010 = statement.column_double(1);
    vel.M001 = statement.column_double(2);
    vel.M200 = statement.column_double(3);
    vel.M110 = statement.column_double(4);
    vel.M101 = statement.column_double(5);
    vel.M020 = statement.column_double(6);
    vel.M011 = statement.column_double(7);
    vel.M002 = statement.column_double(8);
}

////////////////////////////////////////////////////////////////////////////////

void
save_gain_table(FitsObject & gain_file,
                const Lfi_radiometer_t & radiometer,
                const Gain_table_t & gain_table,
                const std::string & comment)
{
    std::vector<fitscolumn> columns {
        { "PID", "", fitsTypeC<int32_t>(), 1 },
        { "GAIN", "", fitsTypeC<double>(), 1 },
        { "OFFSET", "", fitsTypeC<double>(), 1 } };
    gain_file.insertBINtable(columns, radiometer.shortName());
    gain_file.writeColumn(1, gain_table.pointingIds);
    gain_file.writeColumn(2, gain_table.gain);
    gain_file.writeColumn(3, gain_table.offset);

    gain_file.setKey("HORN", radiometer.horn);
    gain_file.setKey("RAD", radiometer.radiometer);

    if(! comment.empty()) {
        gain_file.setComment(comment);
    }
}

////////////////////////////////////////////////////////////////////////////////

void
save_gain_table(const std::string & file_name,
                const Lfi_radiometer_t & radiometer,
                const Gain_table_t & gain_table,
                const std::string & comment)
{
    FitsObject gain_file;
    gain_file.create(file_name, true);

    save_gain_table(gain_file, radiometer, gain_table, comment);
}

////////////////////////////////////////////////////////////////////////////////

void
load_pointing_information(Sqlite_connection_t & db,
                          int first_od,
                          int last_od,
                          std::vector<Pointing_t> & pointings)
{
    Logger * log = Logger::get_instance();

    std::stringstream query;
    query << "\n\tSELECT pointingID, start_pointing, start_time, end_time, "
          << "\n\t       spin_ecl_lon, spin_ecl_lat, od_int"
          << "\n\t  FROM pointings"
          << "\n\t WHERE od_int >= " << first_od
          << "\n\t   AND od_int <= " << last_od
          << "\n\t";

    Sqlite_statement_t statement(db, query.str().c_str());

    int result = statement.step();
    pointings.resize(0);
    while(result == SQLITE_ROW) {
        Pointing_t cur_pointing;

        cur_pointing.id             = statement.column_int   (0);
        cur_pointing.start_pointing = statement.column_double(1);
        cur_pointing.start_time     = statement.column_double(2);
        cur_pointing.end_time       = statement.column_double(3);
        cur_pointing.spin_ecl_lon   = statement.column_double(4);
        cur_pointing.spin_ecl_lat   = statement.column_double(5);
        cur_pointing.od             = statement.column_int   (6);

        pointings.push_back(cur_pointing);
        result = statement.step();
    }
}

////////////////////////////////////////////////////////////////////////////////

void
load_pointing_information(Sqlite_connection_t & ucds,
                          int od,
                          std::vector<Pointing_t> & pointings)
{
    load_pointing_information(ucds, od, od, pointings);
}

////////////////////////////////////////////////////////////////////////////////

void save_tod(const std::string & file_name,
              signed short od,
              const Lfi_radiometer_t & radiometer,
              const DifferencedData & datadiff,
              const std::string & comment)
{
    Logger * log = Logger::get_instance();
    log->info(boost::format("Going to write a TOD into %1%") % file_name);

    FitsObject tod_file;
    tod_file.create(file_name, true);

    std::vector<fitscolumn> columns {
        { "OBT", "", fitsTypeC<double>(), 1 },
        { "SCET", "", fitsTypeC<double>(), 1 },
        { "SKYLOAD", "", fitsTypeC<double>(), 1 },
        { "FLAGS", "", fitsTypeC<int>(), 1 } };
    tod_file.insertBINtable(columns, radiometer.shortName());
    tod_file.writeColumn(1, datadiff.obt_time);
    tod_file.writeColumn(2, datadiff.scet_time);
    tod_file.writeColumn(3, datadiff.sky_load);
    tod_file.writeColumn(4, datadiff.flags);

    tod_file.setKey("HORN", radiometer.horn);
    tod_file.setKey("RAD", radiometer.radiometer);
    tod_file.setKey("OD", od);

    if(! comment.empty()) {
        tod_file.setComment(comment);
    }
}

////////////////////////////////////////////////////////////////////////////////

void
load_pointings(const std::string & file_name,
               PointingData & pointings)
{
    Logger * log = Logger::get_instance();

    if(is_a_squeezer_file(file_name)) {
        decompress_pointings(file_name, pointings);
        return;
    }

    log->info(boost::format("Reading pointings from FITS file %1%")
              % file_name);
    FitsObject fits_file;
    fits_file.openTable(file_name);

    fits_file.getColumn(1, pointings.obt_time);
    fits_file.getColumn(2, pointings.scet_time);
    fits_file.getColumn(3, pointings.theta);
    fits_file.getColumn(4, pointings.phi);
    fits_file.getColumn(5, pointings.psi);

    fits_file.close();
}

////////////////////////////////////////////////////////////////////////////////

void
load_differenced_data(const std::string & file_name,
                      DifferencedData & datadiff)
{
    Logger * log = Logger::get_instance();

    if(is_a_squeezer_file(file_name)) {
        decompress_differenced_data(file_name, datadiff);
        return;
    }

    log->info(boost::format("Reading differenced data from FITS file %1%")
              % file_name);
    FitsObject fits_file;
    fits_file.openTable(file_name);

    fits_file.getColumn(1, datadiff.obt_time);
    fits_file.getColumn(2, datadiff.scet_time);
    fits_file.getColumn(3, datadiff.sky_load);
    fits_file.getColumn(4, datadiff.flags);

    fits_file.close();
}

////////////////////////////////////////////////////////////////////////////////


void
load_subsampled_ref_and_hk_data(Sqlite_connection_t & ucds,
                                const Lfi_radiometer_t & radiometer,
                                const std::string & hk_sensor_name,
                                const Range_t<int> & pointing_range,
                                std::vector<int> & pointing_id,
                                std::vector<double> & ref_data,
                                std::vector<double> & hk_data)
{
    const std::string table_name((boost::format("sci%1%%2%_weighted")
                                  % radiometer.horn
                                  % radiometer.radiometer).str());

    auto query =
        boost::format("SELECT %2%.pointingID as pID, mean_ref, hk.%1% AS temp "
                      "FROM %2% JOIN lfi_hk_temperatures AS hk "
                      "USING (pointingID) "
                      "WHERE pID >= %3% AND pID <= %4% "
                      "AND mean_ref > 0 AND temp > 0")
        % hk_sensor_name
        % table_name
        % pointing_range.start
        % pointing_range.end;

    Sqlite_statement_t statement(ucds, query.str().c_str());

    int result = statement.step();
    pointing_id.clear();
    ref_data.clear();
    while(result == SQLITE_ROW) {
        pointing_id.push_back(statement.column_int(0));
        ref_data.push_back(statement.column_double(1));
        hk_data.push_back(statement.column_double(2));

        result = statement.step();
    }
}
