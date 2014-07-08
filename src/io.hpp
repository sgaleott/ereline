#ifndef IO_HPP
#define IO_HPP

#include <sqlite3.h>

class PlanckVelocity;

void loadConvolutionParametersFromUCDS(sqlite3 * ucds,
				       PlanckVelocity & vel);
void saveGainTable(const std::string & file_name,
		   signed short od,
		   const std::string & detectorId,
		   const std::string & instrument);

#endif
