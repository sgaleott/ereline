#ifndef IO_HPP
#define IO_HPP

#include "ahf_info.hpp"

#include <string>
#include <vector>

#include <sqlite3.h>

class PlanckVelocity;
class SQLite3Connection;
class gainTable;
class LfiRadiometer;

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
