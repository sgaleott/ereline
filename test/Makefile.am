check_PROGRAMS = \
	test_configuration \
	test_logger

test_configuration_SOURCES = \
	test_configuration.cpp \
	$(top_srcdir)/src/common/configuration.cpp
test_configuration_LDADD =  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_logger_SOURCES = \
	test_logger.cpp \
	$(top_srcdir)/src/common/logging.cpp
test_logger_LDADD =  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

AM_CPPFLAGS = -I$(top_srcdir)/src
AM_CXXFLAGS = $(BOOST_CPPFLAGS) -DBOOST_TEST_DYN_LINK
AM_LDFLAGS = $(BOOST_LDFLAGS)

TESTS = $(check_PROGRAMS)
