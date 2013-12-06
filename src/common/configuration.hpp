#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

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
    virtual ~Configuration_error() {}
};

////////////////////////////////////////////////////////////////////////

class Configuration {
private:
    boost::property_tree::ptree ptree;
    std::map<std::string, std::string> fallbacks;
    std::map<std::string, std::string> variables;

public:
    Configuration();
    virtual ~Configuration();

    void read_from_json(const std::string & file_name);

    void set_variable(const std::string & name,
		      const std::string & value);

    std::string substitute_variables(const std::string & str);

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
};

#endif
