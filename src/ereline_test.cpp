#include <string>
#include <iostream>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "logging.hpp"
#include "configuration.hpp"

#include "io.hpp"
#include "dipole_fit.hpp"
#include "da_capo.hpp"
#include "smooth_gains.hpp"
#include "sqlite3xx.hpp"

int main(int argc, const char ** argv)
{
  if(argc != 2) {
      std::cout << "Usage: dx11_pipeline CONFIGURATION_FILE\n";
      return 1;
  }
  Configuration ereline_config;
  ereline_config.read_from_json(argv[1]);
  Configuration storage_config;
  storage_config.read_from_json(ereline_config.getWithSubst("common.data_storage_layout"));

  Logger * log = Logger::get_instance();
  ereline_config.configure_logging();

  std::string ucds_file_path = storage_config.getWithSubst("ucds_file_path");
  log->info(boost::format("UCDS file path: %1%") % ucds_file_path);
  SQLite3Connection ucds(ucds_file_path.c_str());

  std::vector<Pointing_t> pointings;
  loadPointingInformation(ucds,
			  ereline_config.get<int>("common.first_od"),
			  ereline_config.get<int>("common.last_od"),
			  pointings);

  run_dipole_fit(ereline_config);
  run_da_capo(ereline_config);
  run_smooth_gains(ereline_config);

  return 0;
}
