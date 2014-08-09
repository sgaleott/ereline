#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

#include <istream>
#include <map>
#include <stdexcept>
#include <string>

#include <boost/format.hpp>
#include "rapidjson/document.h"

////////////////////////////////////////////////////////////////////////

class Configuration_error : public std::runtime_error {
public:
    Configuration_error(const std::string & what)
        : std::runtime_error(what) {}
    virtual ~Configuration_error() throw() {}
};

////////////////////////////////////////////////////////////////////////

template<typename T> inline T 
conf_getter(const rapidjson::Value & value);

template<> inline bool
conf_getter<bool>(const rapidjson::Value & value)
{
    return value.GetBool();
}

template<> inline std::string 
conf_getter<std::string>(const rapidjson::Value & value)
{
    return value.GetString();
}

template<> inline int
conf_getter<int>(const rapidjson::Value & value)
{
    return value.GetInt();
}

template<> inline unsigned
conf_getter<unsigned>(const rapidjson::Value & value)
{
    return value.GetUint();
}

template<> inline double
conf_getter<double>(const rapidjson::Value & value)
{
    return value.GetDouble();
}

////////////////////////////////////////////////////////////////////////

struct Configuration {
    rapidjson::Document document;
    std::map<std::string, std::string> fallbacks;
    std::map<std::string, std::string> variables;

    Configuration();
    virtual ~Configuration();

    void read_from_json(std::istream & input_stream);
    void read_from_json(const std::string & file_name);

    void set_variable(const std::string & name,
                      const std::string & value);
    void fill_with_standard_variables(const rapidjson::Value & start,
                                      const std::string & start_path);

    std::string substitute_variables(const std::string & str) const;

    rapidjson::Value::ConstMemberIterator
    find_member_from_key(const rapidjson::Value & start,
                         const std::string & key) const;

    template<typename T> T get(const std::string & key) const
    {
        const rapidjson::Value::ConstMemberIterator & member =
            find_member_from_key(document, key);

        if(member == document.MemberEnd())
        {
            auto alternative = fallbacks.find(key);
            if(alternative != fallbacks.end())
            {
                // Recursive call
                return get<T>(alternative->second);
            } else {
                auto msg =
                    boost::format("Unable to find %1% in the configuration file")
                    % key;
                throw Configuration_error(boost::str(msg));
            }
        } else {
            return conf_getter<T>(member->value);
        }
    }
    template<typename T> T get(const std::string & key, T default_val) const
    {
        try {
            return get<T>(key);
        }
        catch(const Configuration_error & err) {
            return default_val;
        }
    }

    std::string getWithSubst(const std::string & key) const
    {
        return substitute_variables(get<std::string>(key));
    }
    std::string getWithSubst(const std::string & key,
                             const std::string & default_val) const
    {
        try {
            return getWithSubst(key);
        }
        catch(const Configuration_error & err) {
            return default_val;
        }
    }

    void configure_logging() const;
};

void setup_od_variable(int od, Configuration & conf);

struct Lfi_radiometer_t;
void setup_variables_for_radiometer(const Lfi_radiometer_t & rad,
                                    Configuration & conf);

#endif
