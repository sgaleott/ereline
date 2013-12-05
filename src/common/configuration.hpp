#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

#include <string>
#include <map>

#include <boost/property_tree/ptree.hpp>
#include <boost/optional.hpp>

class Configuration {
private:
    boost::property_tree::ptree ptree;
    std::map<std::string, std::string> associations;

public:
    Configuration();
    virtual ~Configuration();

    void read_from_json(const std::string & file_name);

    template<typename T> T get(const std::string & key) const
    {
	boost::optional<T> value = ptree.get_optional<T>(key);
	if(! value)
	{
	    auto alternative = associations.find(key);
	    if(alternative != associations.end())
	    {
		return get<T>(alternative->second);
	    } else {
		throw boost::property_tree::ptree_bad_path("Unable to find the path",
							   boost::property_tree::ptree::path_type(key));
	    }
	} else {
	    return value.get();
	}
    }
};

#endif
