#include "configuration.hpp"

#include <stdexcept>
#include <sstream>

#define BOOST_TEST_MODULE "Configuration"
#include <boost/test/unit_test.hpp>

////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(config_file_parsing)
{
    Configuration conf;
    conf.read_from_json(TEST_DATA_DIR "configuration_test.js");

    // Check that information retrieval works
    BOOST_CHECK_EQUAL(conf.get<int>("common.od_range.first"),
		      91);

    // Check that the fallback mechanism works
    BOOST_CHECK_EQUAL(conf.get<std::string>("dipole_fit.input_tod"),
		      conf.get<std::string>("apply_r.output_tod"));

    // Check that wrong property names are flagged properly
    BOOST_CHECK_THROW(conf.get<int>("this.does.not.exist"),
		      std::runtime_error);
}
