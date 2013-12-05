#include "configuration.hpp"

#include <boost/property_tree/json_parser.hpp>

////////////////////////////////////////////////////////////////////////

void
init_associations(std::map<std::string, std::string> & ass)
{
    ass["dipole_fit.input_tod"] = "apply_r.output_tod";
    ass["da_capo.input_gains"] = "dipole_fit.output_gains";
    ass["smooth_gains.input_gains"] = "da_capo.output_gains";
}

////////////////////////////////////////////////////////////////////////

Configuration::Configuration()
{
    init_associations(associations);
}

////////////////////////////////////////////////////////////////////////

Configuration::~Configuration()
{
}

////////////////////////////////////////////////////////////////////////

void
Configuration::read_from_json(const std::string & file_name)
{
    read_json(file_name, ptree);
}
