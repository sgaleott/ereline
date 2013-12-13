#include "configuration.hpp"

#include <boost/property_tree/json_parser.hpp>
#include <boost/format.hpp>

////////////////////////////////////////////////////////////////////////

void
init_fallbacks(std::map<std::string, std::string> & ass)
{
    ass["dipole_fit.input_tod"] = "apply_r.output_tod";
    ass["da_capo.input_gains"] = "dipole_fit.output_gains";
    ass["smooth_gains.input_gains"] = "da_capo.output_gains";
}

////////////////////////////////////////////////////////////////////////

Configuration::Configuration()
{
    init_fallbacks(fallbacks);
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

////////////////////////////////////////////////////////////////////////

void 
Configuration::set_variable(const std::string & name,
			    const std::string & value)
{
    variables[name] = value;
}

////////////////////////////////////////////////////////////////////////

std::string
Configuration::substitute_variables(const std::string & str)
{
    std::string result;

    for(auto cur_char = str.begin(); cur_char != str.end(); ++cur_char)
    {
	if(*cur_char == '{') {
	    cur_char++; // Skip the '{'
	    auto cur_pos = cur_char - str.begin();

	    // Read the variable name
	    auto end_pos = str.find('}', cur_pos);
	    if(end_pos == str.npos) {
		auto msg =
		    boost::format("No closing '}' in \"%1%\"")
		    % str;
		throw Configuration_error(boost::str(msg));
	    }

	    std::string var_name = str.substr(cur_pos, end_pos - cur_pos);
	    auto var_match = variables.find(var_name);
	    if(var_match == variables.end()) {
		auto msg =
		    boost::format("Undefined variable %1% in \"%2%\"")
		    % var_name
		    % str;
		throw Configuration_error(boost::str(msg));
	    }

	    result += var_match->second;
	    cur_char += end_pos - cur_pos;
	} else {
	    result += *cur_char;
	}
    }

    return result;
}
