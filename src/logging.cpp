#include "logging.hpp"
#include "configuration.hpp"
#include <ctime>
#include <iostream>
#include <fstream>

#include <boost/format.hpp>

bool Logger::exist_instance_flag = false;
Logger * Logger::singleton = NULL;

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
    time_t t = std::time(NULL);
    struct tm * tmp = std::localtime(&t);
    if (tmp == NULL)
        std::runtime_error("unable to call \"std::localtime\"");

    char datetime_str[64];
    if(std::strftime(datetime_str, sizeof(datetime_str) - 1,
                     "%FT%T", tmp) == 0)
        std::runtime_error("unable to call \"std::strftime\"");

    return std::string(datetime_str);
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

    boost::format msg;
    if(mpi_rank >= 0 && mpi_size >= 0) {
        // Include the MPI rank of this process and the number of
        // processes in the log message
        msg = (boost::format("%1% %2%: [%5%/%6%] %3%%4%")
               % current_date_and_time()
               % error_type
               % std::string(indent_level, '\t') // Indent
               % string
               % (mpi_rank + 1)
               % mpi_size);
    } else {
        msg = (boost::format("%1% %2%: %3%%4%")
               % current_date_and_time()
               % error_type
               % std::string(indent_level, '\t') // Indent
               % string);
    }

    for(auto log_stream : log_stream_list) {
        if(log_stream != NULL) {
            (*log_stream) << msg << std::endl;
            log_stream->flush();
        }
    }
}
