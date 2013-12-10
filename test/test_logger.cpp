#include "common/logging.hpp"

#include <sstream>

#define BOOST_TEST_MODULE Logging
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE(logging_test)
{
    Logger * log = Logger::get_instance();
    std::stringstream str;

    log->set_stream(&str);
}
