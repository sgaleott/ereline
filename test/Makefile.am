AUTOMAKE_OPTIONS = subdir-objects

check_PROGRAMS = \
	test_configuration \
	test_logger
EXTRA_DIST = ./test_data

test_configuration_SOURCES = \
	test_configuration.cpp \
	../src/configuration.cpp
test_configuration_LDADD =  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_logger_SOURCES = \
	test_logger.cpp \
	../src/logging.cpp
test_logger_LDADD =  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

AM_CPPFLAGS = -I$(top_srcdir)/src -DTEST_DATA_DIR=\"${srcdir}/test_data/\"
AM_CXXFLAGS = $(BOOST_CPPFLAGS) -DBOOST_TEST_DYN_LINK
AM_LDFLAGS = $(BOOST_LDFLAGS)

TESTS = $(check_PROGRAMS)
