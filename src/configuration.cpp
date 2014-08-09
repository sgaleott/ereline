#include "configuration.hpp"

#include "datatypes.hpp"
#include "logging.hpp"

#include <fstream>
#include <boost/format.hpp>
#include "rapidjson/error/en.h"

////////////////////////////////////////////////////////////////////////

void
init_fallbacks(std::map<std::string, std::string> & ass)
{
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
    std::string text;
    while(input_stream.good()) {
        std::string cur_line;
        std::getline(input_stream, cur_line);
        text += cur_line;
    }
    document.Parse(text.c_str());
    if(document.HasParseError()) {
        const size_t error_pos = document.GetErrorOffset();
        size_t start_pos, len;
        if(error_pos > 10)
            start_pos = error_pos - 10;
        else
            start_pos = 0;

        len = text.size() - start_pos;
        if(len > 20)
            len = 20;

        std::string msg =
            (boost::format("Error parsing the JSON parameter file: %1% near «%2%»")
             % rapidjson::GetParseError_En(document.GetParseError())
             % text.substr(start_pos, len)).str();
        throw std::runtime_error(msg);
    }

    fill_with_standard_variables(document, "");
}

////////////////////////////////////////////////////////////////////////

void
Configuration::read_from_json(const std::string & file_name)
{
    std::ifstream is(file_name);
    read_from_json(is);
}

////////////////////////////////////////////////////////////////////////

// Look for a member of the JSON document using a dotted-path syntax
rapidjson::Value::ConstMemberIterator
Configuration::find_member_from_key(const rapidjson::Value & start,
                                    const std::string & key) const
{
    if(! start.IsObject())
        return start.MemberEnd();

    std::string s(key);
    const std::string delimiter(".");

    size_t pos = s.find(delimiter);
    std::string token;
    if(pos == std::string::npos) {
        token = s;
        s = "";
    } else {
        token = s.substr(0, pos);
        s.erase(0, pos + delimiter.length());
    }

    const auto & sub_element = start.FindMember(token.c_str());
    if(s.empty() || sub_element == start.MemberEnd())
        return sub_element;

    if(sub_element->value.IsObject()) {
        auto result = find_member_from_key(sub_element->value, s);
        if(result == sub_element->value.MemberEnd())
            return start.MemberEnd();
        else
            return result;
    }

    return start.MemberEnd();
}

////////////////////////////////////////////////////////////////////////

void
Configuration::fill_with_standard_variables(const rapidjson::Value & start,
                                            const std::string & start_path)
{
    if(! start.IsObject())
        return;

    for (auto itr = start.MemberBegin();
         itr != start.MemberEnd(); 
         ++itr) 
    {
        std::string sub_path;
        if(start_path.empty())
            sub_path = itr->name.GetString();
        else
            sub_path = 
                (boost::format("%s.%s") 
                 % start_path 
                 % itr->name.GetString()).str();

        if(itr->value.IsString()) {
            set_variable(sub_path, itr->value.GetString());
        } else if(itr->value.IsObject()) {
            fill_with_standard_variables(itr->value, sub_path);
        }
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
