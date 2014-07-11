#include "io.hpp"
#include "planck_velocity.hpp"
#include "logging.hpp"
#include "datatypes.hpp"
#include "sqlite3xx.hpp"

#include <sstream>
#include <boost/format.hpp>

////////////////////////////////////////////////////////////////////////////////

void
loadConvolutionParametersFromUCDS(SQLite3Connection & ucds,
				  const LfiRadiometer & radiometer,
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
	      const std::string & detectorId,
	      const std::string & instrument)
{
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

        cur_pointing.pointing_id    = statement.column_int   (0);
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
