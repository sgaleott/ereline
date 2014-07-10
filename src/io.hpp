#ifndef IO_HPP
#define IO_HPP

#include "ahf_info.hpp"

#include <string>
#include <vector>

#include <sqlite3.h>

class PlanckVelocity;
class SQLite3Connection;
class LfiRadiometer;

void loadConvolutionParametersFromUCDS(SQLite3Connection & ucds,
				       const LfiRadiometer & radiometer,
				       PlanckVelocity & vel);

void saveGainTable(const std::string & file_name,
		   signed short od,
		   const std::string & detectorId,
		   const std::string & instrument);

void loadPointingInformation(SQLite3Connection & ucds,
			     int first_od,
			     int last_od,
			     std::vector<Pointing_t> & pointings);

#endif
