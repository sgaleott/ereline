#ifndef DATATYPES_HPP
#define DATATYPES_HPP

#include <exception>
#include <string>
#include <boost/format.hpp>

////////////////////////////////////////////////////////////////////////////////

class ConfigurationError : public std::exception {
    std::string description;

public:
    ConfigurationError(const std::string & a_description) noexcept
        : description(a_description) { }
    ConfigurationError(const boost::format & a_description) noexcept
	: description(a_description.str()) { }

    const char * what() noexcept { return description.c_str(); }
};

////////////////////////////////////////////////////////////////////////////////

class IoError : public std::exception {
    std::string description;

public:
    IoError(const std::string & a_description) noexcept
	: description(a_description) { }
    IoError(const boost::format & a_description) noexcept
        : description(a_description.str()) { }

    const char * what() noexcept { return description.c_str(); }
};

////////////////////////////////////////////////////////////////////////////////

struct LfiRadiometer {
    int horn;
    int radiometer;

    LfiRadiometer(const std::string & name);
    LfiRadiometer(int a_horn = 18, int a_radiometer = 0);

    void assign(const std::string & name);
    void assign(int a_horn, int a_radiometer);

    std::string shortName() const; // E.g. LFI28M
    std::string fullNameWithDetector(int detector = 0) const; // E.g. LFI28M-00
    std::string armName() const;
    LfiRadiometer twinRadiometer() const;
};

#endif
