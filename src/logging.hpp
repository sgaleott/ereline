#ifndef LOGGING_HPP
#define LOGGING_HPP

#include <string>
#include <memory>
#include <ostream>
#include <boost/format.hpp>

class Logger {
public:
    enum class Log_level { DEBUG, INFO, WARNING, ERROR };

    static Logger * get_instance();

    void log(Log_level level, const std::string & string) const;

    void increase_indent() { 
	indent_level++; 
    }
    void decrease_indent() { 
	if(indent_level > 0) 
	    indent_level--; 
    }

    void append_stream(std::ostream * new_log_stream) { 
	log_stream_list.push_back(std::auto_ptr<std::ostream>(new_log_stream)); 
    }

    void debug(const std::string & string) const { 
	log(Log_level::DEBUG, string); 
    }
    void debug(const boost::format & string_fmt) const {
	debug(boost::str(string_fmt)); 
    }

    void info(const std::string & string) const {
	log(Log_level::INFO, string); 
    }
    void info(const boost::format & string_fmt) const { 
	info(boost::str(string_fmt)); 
    }

    void warning(const std::string & string) const { 
	log(Log_level::WARNING, string); 
    }
    void warning(const boost::format & string_fmt) const { 
	warning(boost::str(string_fmt)); 
    }

    void error(const std::string & string) const { 
	log(Log_level::ERROR, string); 
    }
    void error(const boost::format & string_fmt) const { 
	error(boost::str(string_fmt)); 
    }

private:
    static bool exist_instance_flag;
    static Logger * singleton;

    std::vector<std::unique_ptr<std::ostream> > log_stream_list;

    int indent_level;
    Log_level log_level;

    Logger() // Private constructor
	{
	    append_stream(&std::cerr);
	    indent_level = 0;
	    log_level = Log_level::INFO;
	}

public:
    void set_log_level(Log_level new_level) { 
	log_level = new_level; 
    }

    Log_level get_log_level() const { 
	return log_level; 
    }
};

#endif
