#include "common/configuration.hpp"

#define BOOST_TEST_MODULE "Configuration"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(variable_substitution_test)
{
    Configuration conf;

    conf.set_variable("test1", "value1");
    conf.set_variable("test2", "value2");

    // No substitution at all
    BOOST_CHECK_EQUAL(conf.substitute_variables("-- test1 is test2 --"),
		      "-- test1 is test2 --");
    // Basic variable substitution
    BOOST_CHECK_EQUAL(conf.substitute_variables("-- {test1} is {test2} --"),
		      "-- value1 is value2 --");
    // Check that the substitution works at the beginning of the string
    BOOST_CHECK_EQUAL(conf.substitute_variables("{test1} --"),
		      "value1 --");
    // Check that the substitution works at the end of the string
    BOOST_CHECK_EQUAL(conf.substitute_variables("{test1} --"),
		      "value1 --");
    // Check that multiple substitution work
    BOOST_CHECK_EQUAL(conf.substitute_variables("{test1}{test1}{test1}"),
		      "value1value1value1");

    // Check that undefined variables are properly detected
    BOOST_CHECK_THROW(conf.substitute_variables("{does_not_exist}"),
		      Configuration_error);

    // Check that unmatched '{}' are detected
    BOOST_CHECK_THROW(conf.substitute_variables("{test1"),
		      Configuration_error);
}
