#include "logging.hpp"
#include "configuration.hpp"
#include <iostream>
#include <fstream>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/format.hpp>

bool Logger::exist_instance_flag = false;
Logger * Logger::singleton = nullptr;

////////////////////////////////////////////////////////////////////////

Logger * 
Logger::get_instance()
{
    if(! exist_instance_flag) {
	singleton = new Logger();
	exist_instance_flag = true;
    }

    return singleton;
}

////////////////////////////////////////////////////////////////////////

static std::string
current_date_and_time()
{
    boost::posix_time::ptime now(boost::date_time::second_clock<boost::posix_time::ptime>::universal_time());

    return boost::posix_time::to_iso_extended_string(now);
}

////////////////////////////////////////////////////////////////////////

void 
Logger::log(Log_level level, const std::string & string) const
{
    if(level < log_level || log_stream_list.empty())
	return; // Do not log anything

    std::string error_type;
    switch(level) {
    case Log_level::DEBUG:   error_type = "  DEBUG"; break;
    case Log_level::INFO:    error_type = "   INFO"; break;
    case Log_level::WARNING: error_type = "WARNING"; break;
    case Log_level::ERROR:   error_type = "  ERROR"; break;
    }

    auto msg = (boost::format("%1% %2%: %3%%4%") 
		% current_date_and_time() 
		% error_type 
		% std::string(indent_level, '\t') // Indent
		% string);

    for(auto & log_stream : log_stream_list) {
	if(log_stream.get() != NULL) {
	    (*log_stream) << msg << std::endl;
	    log_stream->flush();
	}
    }
}
