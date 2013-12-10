check_PROGRAMS = \
	test_configuration \
	test_logger

test_configuration_SOURCES = test_configuration.cpp $(top_srcdir)/src/common/configuration.cpp
test_logger_SOURCES = test_logger.cpp $(top_srcdir)/src/common/logging.cpp

AM_CPPFLAGS = -I$(top_srcdir)/src
AM_CXXFLAGS = $(BOOST_CPPFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

TESTS = $(check_PROGRAMS)
