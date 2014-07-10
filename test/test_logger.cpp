#include "logging.hpp"

#include <sstream>

#define BOOST_TEST_MODULE "Logging"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(logging_test)
{
    Logger * log = Logger::get_instance();
    std::stringstream log_stream;

    log->append_stream(&log_stream);
    log->set_log_level(Logger::Log_level::INFO);

    const std::string debug_str("This is a debug message");
    const std::string info_str("This is just some plain information");
    const std::string warning_str("Warning! Something suspicious has happened");
    const std::string error_str("Error!");

    log->debug(debug_str); // This should be ignored
    log->info(info_str);
    log->warning(warning_str);
    log->error(error_str);

    log_stream.seekg(0); // Move the stream back to start
    std::cout << "The string is:\n\"" << log_stream.str() << "\"\n\n";
    log_stream.seekg(0); // Move the stream back to start

    std::string cur_line;

    // Do not check for the debug message, Logger has ignored it
    const int length_of_log_header = 29;
    std::getline(log_stream, cur_line);
    BOOST_CHECK_EQUAL(info_str, cur_line.substr(length_of_log_header));
    std::getline(log_stream, cur_line);
    BOOST_CHECK_EQUAL(warning_str, cur_line.substr(length_of_log_header));
    std::getline(log_stream, cur_line);
    BOOST_CHECK_EQUAL(error_str, cur_line.substr(length_of_log_header));
}

