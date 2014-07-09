#include <string>
#include <iostream>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "logging.hpp"
#include "configuration.hpp"

#include "dipole_fit.hpp"
#include "da_capo.hpp"
#include "smooth_gains.hpp"

int main(int argc, const char ** argv)
{
  if(argc != 2) {
      std::cout << "Usage: dx11_pipeline CONFIGURATION_FILE\n";
      return 1;
  }
  Configuration ereline_config;
  ereline_config.read_from_json(argv[1]);

  Logger * log = Logger::get_instance();
  log->configure(ereline_config);

  run_dipole_fit(ereline_config);
  run_da_capo(ereline_config);
  run_smooth_gains(ereline_config);

  return 0;
}
