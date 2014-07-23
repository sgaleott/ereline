#ifndef LOGGING_HPP
#define LOGGING_HPP

#include <string>
#include <memory>
#include <iostream>
#include <boost/format.hpp>

class Logger {
public:
    enum Log_level { DEBUG, INFO, WARNING, ERROR };

    static Logger * get_instance();

    void set_mpi_rank(int rank, int size) {
        mpi_rank = rank;
        mpi_size = size;
    }

    void log(Log_level level, const std::string & string) const;

    void increase_indent() {
        indent_level++;
    }
    void decrease_indent() {
        if(indent_level > 0)
            indent_level--;
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
    int mpi_rank;
    int mpi_size;

    std::vector<std::ostream *> log_stream_list;

    int indent_level;
    Log_level log_level;

    Logger() // Private constructor
        {
            append_stream(&std::cerr);
            indent_level = 0;
            log_level = Log_level::INFO;
            mpi_rank = mpi_size = -1;
        }

public:
    void append_stream(std::ostream * new_log_stream) {
        log_stream_list.push_back(new_log_stream);
    }

    void set_log_level(Log_level new_level) {
        log_level = new_level;
    }

    Log_level get_log_level() const {
        return log_level;
    }
};

#endif
