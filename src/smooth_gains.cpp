#include "smooth_gains.hpp"
#include "configuration.hpp"
#include "logging.hpp"

void
run_smooth_gains(const Configuration & program_conf,
		 const Configuration & storage_conf)
{
  Logger * log = Logger::get_instance();

  log->info("Starting module smoothGains");
  log->info("Quitting module smoothGains");
}
