#include "io.hpp"
#include "dipole_fit.hpp"
#include "planck_velocity.hpp"
#include "fits_object.hpp"
#include "gain_table.hpp"
#include "logging.hpp"
#include "datatypes.hpp"
#include "sqlite3xx.hpp"

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
loadConvolutionParametersFromUCDS(SQLite3Connection & ucds,
				  const Lfi_radiometer_t & radiometer,
				  PlanckVelocity & vel)
{
    Logger * log = Logger::get_instance();

    std::stringstream query;
    query << "\n\tSELECT m100, m010, m001, m200, m110, m101, m020, m011, m002"
	  << "\n\t  FROM convolved_dipole"
	  << "\n\t WHERE horn = " << radiometer.horn
	  << "\n\t   AND radiometer = " << radiometer.radiometer
	  << "\n\t";

    SQLite3Statement statement(ucds, query.str().c_str());

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
saveGainTable(const std::string & file_name,
	      signed short od,
	      const Lfi_radiometer_t & radiometer,
	      const Gain_table_t & gain_table,
	      const std::string & comment)
{
    FitsObject gain_file;

    gain_file.create(file_name, true);

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
    gain_file.setKey("OD", od);

    if(! comment.empty()) {
	gain_file.setComment(comment);
    }
}

////////////////////////////////////////////////////////////////////////////////

void
loadPointingInformation(SQLite3Connection & db,
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

    SQLite3Statement statement(db, query.str().c_str());

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
loadPointingInformation(SQLite3Connection & ucds,
			int od,
			std::vector<Pointing_t> & pointings)
{
    loadPointingInformation(ucds, od, od, pointings);
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
save_dipole_fit(const std::string & file_name,
		const Lfi_radiometer_t & radiometer,
		const Dipole_fit_t & fit,
		const std::string & comment)
{
    Logger * log = Logger::get_instance();
    log->info(boost::format("Going to write a dipoleFit result "
			    "into %1%") % file_name);

    FitsObject dipole_fit_file;
    dipole_fit_file.create(file_name, true);

    std::vector<fitscolumn> columns {
	{ "PIXIDX", "", fitsTypeC<int>(), 1 },
	{ "PIXHITS", "", fitsTypeC<int>(), 1 },
	{ "PIXDATA", "", fitsTypeC<double>(), 1 },
	{ "PIXDIP", "", fitsTypeC<float>(), 1 } };
    dipole_fit_file.insertBINtable(columns, radiometer.shortName());
    dipole_fit_file.writeColumn(1, fit.pixIndex);
    dipole_fit_file.writeColumn(2, fit.pixSumHits);
    dipole_fit_file.writeColumn(3, fit.pixSumData);
    dipole_fit_file.writeColumn(4, fit.pixSumDipole);

    dipole_fit_file.setKey("NSIDE", fit.nSide);
    dipole_fit_file.setKey("PID", fit.pointingID);
    dipole_fit_file.setKey("FLAG", fit.qualityFlag);
    dipole_fit_file.setKey("GAINV", fit.gainv);
    dipole_fit_file.setKey("OFFSET", fit.offset);
    dipole_fit_file.setKey("DIPMIN", fit.minDipole);
    dipole_fit_file.setKey("DIPMAX", fit.maxDipole);

    if(! comment.empty()) {
	dipole_fit_file.setComment(comment);
    }
}

