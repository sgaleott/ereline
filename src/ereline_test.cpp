#include <string>
#include <iostream>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "logging.hpp"
#include "configuration.hpp"

int main()
{
  std::string hello("Hello, world!");

  Logger * log = Logger::get_instance();

  Configuration conf;
  conf.read_from_json("/home/tomasi/work/planck/ereline/configuration.json");

  log->info(boost::str(boost::format("smooth_gains.run = %1%") % 
		       conf.get<bool>("smooth_gains.run")));
  log->info(boost::str(boost::format("da_capo.input_gains = %1%\n") % 
		       conf.get<std::string>("da_capo.input_gains")));

  return 0;
}
