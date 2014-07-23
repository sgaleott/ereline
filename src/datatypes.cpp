#include "datatypes.hpp"

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

Lfi_radiometer_t::Lfi_radiometer_t(const std::string & name)
{
    assign(name);
}

Lfi_radiometer_t::Lfi_radiometer_t(int a_horn, int a_radiometer)
{
    assign(a_horn, a_radiometer);
}

void
Lfi_radiometer_t::assign(const std::string & name)
{
    size_t horn_str_start_pos;
    size_t arm_char_pos;

    switch(name.size()) {
    case 3: // E.g. "28M"
	horn_str_start_pos = 0;
	arm_char_pos = 2;
	break;
	
    case 6: // E.g. "LFI28M"
    case 8: // E.g. "LFI28M-0"
    case 9: // E.g. "LFI28M-00"
	horn_str_start_pos = 3;
	arm_char_pos = 5;
	break;

    default:
	throw ConfigurationError(boost::format("Invalid radiometer name \"%1%\"")
				 % name);
    }

    try {
	horn = boost::lexical_cast<int>(name.substr(horn_str_start_pos, 2));
    }
    catch(const boost::bad_lexical_cast & exc) {
	throw ConfigurationError(boost::format("No horn in radiometer "
					       "name \"%1%\"")
				 % name);
    }

    if(horn < 18 || horn > 28) {
	throw ConfigurationError(boost::format("Invalid horn number (%2%) in "
					       "radiometer name \"%1%\"")
				 % name % horn);
    }

    if(name[arm_char_pos] == 'M' || name[arm_char_pos] == 'm') {
	radiometer = 0;
    } else if(name[arm_char_pos] == 'S' || name[arm_char_pos] == 's') {
	radiometer = 1;
    } else {
	throw ConfigurationError(boost::format("Invalid arm in radiometer name \"%1%\"")
				 % name);
    }
}

void 
Lfi_radiometer_t::assign(int a_horn, int a_radiometer)
{
    horn = a_horn;
    radiometer = a_radiometer;
}

int Lfi_radiometer_t::frequencyInGhz() const
{
    if(horn <= 23)
	return 70;
    else if(horn <= 26)
	return 44;
    else
	return 30;
}

std::string Lfi_radiometer_t::shortName() const
{
    auto f = boost::format("LFI%1%%2%") % horn % armName();
    return f.str();
}

std::string
Lfi_radiometer_t::fullNameWithDetector(int detector) const
{
    auto f = 
	boost::format("LFI%1%%2%-%3%%4%") 
	% horn 
	% armName() 
	% radiometer 
	% detector;

    return f.str();
}

std::string Lfi_radiometer_t::armName() const
{
    if(radiometer == 0)
	return "M";
    else
	return "S";
}

Lfi_radiometer_t 
Lfi_radiometer_t::twinRadiometer() const
{
    return Lfi_radiometer_t(horn, radiometer == 0 ? 1 : 0);
}

////////////////////////////////////////////////////////////////////////////////

Range_t<size_t>
find_boundaries_in_obt_times(const std::vector<uint64_t> & source,
			     Range_t<uint64_t> obt_range)
{
    const auto first_element = std::lower_bound(source.begin(),
						source.end(),
						obt_range.start);
    const auto last_element = std::upper_bound(source.begin(),
					       source.end(),
					       obt_range.end);

    Range_t<size_t> result { 0, source.size() - 1 };

    if(first_element != source.end())
	result.start = std::distance(source.begin(), first_element);

    if(last_element != source.end())
	result.end = std::distance(source.begin(), last_element);

    return result;
}

////////////////////////////////////////////////////////////////////////////////

PointingData::PointingData(const PointingData & obj, 
			   Range_t<uint64_t> obt_range)
{
    const Range_t<size_t> index(
	find_boundaries_in_obt_times(obj.obt_time, obt_range));

    obt_time.assign(obj.obt_time.begin() + index.start, 
		    obj.obt_time.begin() + index.end);
    scet_time.assign(obj.scet_time.begin() + index.start, 
		     obj.scet_time.begin() + index.end);
    theta.assign(obj.theta.begin() + index.start, 
		 obj.theta.begin() + index.end);
    phi.assign(obj.phi.begin() + index.start, 
	       obj.phi.begin() + index.end);
    psi.assign(obj.psi.begin() + index.start, 
	       obj.psi.begin() + index.end);
}

////////////////////////////////////////////////////////////////////////////////

DifferencedData::DifferencedData(const DifferencedData & obj, 
				 Range_t<uint64_t> obt_range)
{
    const Range_t<size_t> index(
	find_boundaries_in_obt_times(obj.obt_time, obt_range));

    obt_time.assign(obj.obt_time.begin() + index.start, 
		    obj.obt_time.begin() + index.end);
    scet_time.assign(obj.scet_time.begin() + index.start, 
		     obj.scet_time.begin() + index.end);
    sky_load.assign(obj.sky_load.begin() + index.start, 
		    obj.sky_load.begin() + index.end);
    flags.assign(obj.flags.begin() + index.start, 
		 obj.flags.begin() + index.end);
}
