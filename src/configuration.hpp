#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

#include <istream>
#include <map>
#include <stdexcept>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/format.hpp>
#include <boost/optional.hpp>

////////////////////////////////////////////////////////////////////////

class Configuration_error : public std::runtime_error {
public:
    Configuration_error(const std::string & what)
        : std::runtime_error(what) {}
    virtual ~Configuration_error() throw() {}
};

////////////////////////////////////////////////////////////////////////

struct Configuration {
    boost::property_tree::ptree ptree;
    std::map<std::string, std::string> fallbacks;
    std::map<std::string, std::string> variables;

    Configuration();
    virtual ~Configuration();

    void read_from_json(std::istream & input_stream);
    void read_from_json(const std::string & file_name);

    void set_variable(const std::string & name,
                      const std::string & value);
    void fill_with_standard_variables(const std::string & start_path);

    std::string substitute_variables(const std::string & str) const;

    template<typename T> T get(const std::string & key) const
    {
        boost::optional<T> value = ptree.get_optional<T>(key);
        if(! value)
        {
            auto alternative = fallbacks.find(key);
            if(alternative != fallbacks.end())
            {
                // Recursive call
                return get<T>(alternative->second);
            } else {
                auto msg =
                    boost::format("Unable to find the path %1%")
                    % key;
                throw Configuration_error(boost::str(msg));
            }
        } else {
            return value.get();
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
