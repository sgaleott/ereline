#include <string>
#include <iostream>
#include <boost/foreach.hpp>

#include "common/logging.hpp"
#include "common/configuration.hpp"

int main()
{
  std::string hello("Hello, world!");

  Logger * log = Logger::get_instance();

  Configuration conf;
  conf.read_from_json("/home/tomasi/work/planck/ereline/configuration.json");

  std::cout << "smooth_gains.run = " 
	    << conf.get<bool>("smooth_gains.run") 
	    << std::endl;

  std::cout << "da_capo.input_gains = " 
	    << conf.get<std::string>("da_capo.input_gains") 
	    << std::endl;

  return 0;
}
