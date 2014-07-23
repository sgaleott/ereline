AUTOMAKE_OPTIONS = subdir-objects

check_PROGRAMS = \
	test_configuration \
	test_logger \
	test_datatypes \
	test_io \
	test_misc \
	test_project_tod \
	test_squeezer
EXTRA_DIST = ./test_data

test_configuration_SOURCES = \
	test_configuration.cpp \
	../src/logging.cpp \
	../src/datatypes.cpp \
	../src/configuration.cpp
test_configuration_LDADD =  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_logger_SOURCES = \
	test_logger.cpp \
	../src/logging.cpp
test_logger_LDADD =  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_datatypes_SOURCES = \
	test_datatypes.cpp \
	../src/datatypes.cpp
test_datatypes_LDADD =  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_project_tod_SOURCES = \
	test_project_tod.cpp \
	../src/datatypes.cpp \
	../src/logging.cpp \
	../src/project_tod.cpp
test_project_tod_LDADD =  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_misc_SOURCES = \
	test_misc.cpp \
	../src/ahf_info.cpp \
	../src/mpi_processes.cpp
test_misc_LDADD =  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_io_SOURCES = \
	test_io.cpp \
	../src/datatypes.cpp \
	../src/dipole_parameters.cpp \
	../src/fits_object.cpp \
	../src/io.cpp \
	../src/logging.cpp \
	../src/misc.cpp \
	../src/planck_velocity.cpp \
	../src/sqlite3xx.cpp
test_io_LDADD = $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
	$(BOOST_SYSTEM_LIB) \
	$(BOOST_FILESYSTEM_LIB) \
	$(SQLITE3_LDFLAGS)

test_squeezer_SOURCES = \
	test_squeezer.cpp \
	../src/datatypes.cpp \
	../src/logging.cpp \
	../src/squeezer.cpp
test_squeezer_LDADD =  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

AM_CPPFLAGS = -I$(top_srcdir)/src -DTEST_DATA_DIR=\"${srcdir}/test_data/\"
AM_CXXFLAGS = $(BOOST_CPPFLAGS) -DBOOST_TEST_DYN_LINK
AM_LDFLAGS = $(BOOST_LDFLAGS)

TESTS = $(check_PROGRAMS)
