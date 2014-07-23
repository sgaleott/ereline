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

struct Lfi_radiometer_t {
    int horn;
    int radiometer;

    Lfi_radiometer_t(const std::string & name);
    Lfi_radiometer_t(int a_horn = 18, int a_radiometer = 0);

    void assign(const std::string & name);
    void assign(int a_horn, int a_radiometer);

    int frequencyInGhz() const;
    std::string shortName() const; // E.g. LFI28M
    std::string fullNameWithDetector(int detector = 0) const; // E.g. LFI28M-00
    std::string armName() const;
    Lfi_radiometer_t twinRadiometer() const;
};

////////////////////////////////////////////////////////////////////////////////

template<typename T>
struct Range_t {
    T start;
    T end;

    // This is useful for unit tests
    bool operator!=(const Range_t<T> & other) const {
	return (start != other.start) || (end != other.end);
    }
};

template<typename T>
std::ostream & operator<<(std::ostream & os, Range_t<T> const & yt)
{
    os << "(" << yt.start << ", " << yt.end << ")";
    return os;
}

////////////////////////////////////////////////////////////////////////////////

struct PointingData {
    std::vector<uint64_t> obt_time;
    std::vector<double> scet_time;
    std::vector<double> theta;
    std::vector<double> phi;
    std::vector<double> psi;

    PointingData() {}
    PointingData(const PointingData & obj, 
		 Range_t<uint64_t> obt_range);
};

////////////////////////////////////////////////////////////////////////////////

struct DifferencedData {
    std::vector<uint64_t> obt_time;
    std::vector<double> scet_time;
    std::vector<double> sky_load;
    std::vector<uint32_t> flags;

    DifferencedData() {}
    DifferencedData(const DifferencedData & obj, 
		    Range_t<uint64_t> obt_range);
};

////////////////////////////////////////////////////////////////////////////////

Range_t<size_t>
find_boundaries_in_obt_times(const std::vector<uint64_t> & source,
			     Range_t<uint64_t> obt_range);

#endif
