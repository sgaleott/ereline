#ifndef SQUEEZER_HPP
#define SQUEEZER_HPP

#include <string>
#include <exception>
#include <boost/format.hpp>

struct PointingData;
struct DifferencedData;

class SqueezerError : public std::exception {
    std::string description;

public:
    SqueezerError(const std::string & a_description) noexcept
        : description(a_description) { }
    SqueezerError(const boost::format & a_description) noexcept
	: description(a_description.str()) { }

    const char * what() noexcept { return description.c_str(); }
};

void decompress_pointings(const std::string & file_name,
			  PointingData & pointings);

void decompress_differenced_data(const std::string & file_name,
				 DifferencedData & datadiff);

#endif
