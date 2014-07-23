#include "configuration.hpp"

#include "datatypes.hpp"
#include "logging.hpp"

#include <boost/property_tree/json_parser.hpp>
#include <boost/format.hpp>

////////////////////////////////////////////////////////////////////////

void
init_fallbacks(std::map<std::string, std::string> & ass)
{
    ass["dipole_fit.input_tod"] = "apply_r.output_tod";
    ass["da_capo.solar_dipole"] = "dipole_fit.solar_dipole";
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
Configuration::read_from_json(std::istream & input_stream)
{
    read_json(input_stream, ptree);
    fill_with_standard_variables("");
}

////////////////////////////////////////////////////////////////////////

void
Configuration::read_from_json(const std::string & file_name)
{
    read_json(file_name, ptree);
    fill_with_standard_variables("");
}

////////////////////////////////////////////////////////////////////////

void
Configuration::fill_with_standard_variables(const std::string & start_path)
{
    boost::property_tree::ptree children = ptree.get_child(start_path);

    for(const auto & kv : children)
    {
        if(kv.first.empty())
            break;

        std::string cur_path;
        if(start_path.empty())
            cur_path = kv.first;
        else
            cur_path = start_path + "." + kv.first;

        set_variable(cur_path, ptree.get<std::string>(cur_path));
        fill_with_standard_variables(cur_path);
    }
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
Configuration::substitute_variables(const std::string & str) const
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
                    boost::format("Variable %1% in \"%2%\" is undefined or has the wrong type")
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

////////////////////////////////////////////////////////////////////////////////

void
Configuration::configure_logging() const
{
    Logger * log = Logger::get_instance();

    auto log_file_name = getWithSubst("common.log_file");
    std::ofstream * log_stream = new std::ofstream(log_file_name);
    log->append_stream(log_stream);

    Logger::Log_level log_level;
    switch(get<int>("common.log_level")) {
    case 1: log_level = Logger::Log_level::ERROR; break;
    case 2: log_level = Logger::Log_level::WARNING; break;
    case 3: log_level = Logger::Log_level::INFO; break;
    default: log_level = Logger::Log_level::DEBUG;
    }

    log->set_log_level(log_level);
}

////////////////////////////////////////////////////////////////////////////////

void setup_od_variable(int od, Configuration & conf)
{
    conf.set_variable("od", (boost::format("%d") % od).str());
    conf.set_variable("odNNNN", (boost::format("%04d") % od).str());
    conf.set_variable("odNNNNNN", (boost::format("%06d") % od).str());
    conf.set_variable("odNNNNNNNN", (boost::format("%08d") % od).str());
}

////////////////////////////////////////////////////////////////////////////////

void
setup_variables_for_radiometer(const Lfi_radiometer_t & rad,
                               Configuration & conf)
{
    conf.set_variable("horn",
                      (boost::format("%1%") % rad.horn).str());
    conf.set_variable("arm",
                      rad.armName());
    conf.set_variable("frequency_GHz",
                      (boost::format("%1%") % rad.frequencyInGhz()).str());
}
