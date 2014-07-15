#include "ahf_info.hpp"
#include "mpi_processes.hpp"

#include <stdexcept>
#include <sstream>

#define BOOST_TEST_MODULE "Miscellaneous"
#include <boost/test/unit_test.hpp>

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(OdList)
{
    const std::vector<Pointing_t> list_of_pointings {
	// Od  Unused  Unused  Unused  Unused  Unused  Od
	{   1,    0.0,    0.1,    1.0,    0.0,    0.0, 91 },
	{   2,    1.0,    1.1,    2.0,    0.0,    0.0, 91 },
	{   3,    2.0,    2.1,    3.0,    0.0,    0.0, 91 },
	{   4,    3.0,    3.1,    4.0,    0.0,    0.0, 91 },
	{   5,    4.0,    4.1,    5.0,    0.0,    0.0, 92 },
	{   6,    5.0,    5.1,    6.0,    0.0,    0.0, 93 },
	{   7,    6.0,    6.1,    7.0,    0.0,    0.0, 93 },
	{   8,    7.0,    7.1,    8.0,    0.0,    0.0, 94 },
    };

    auto od_list = build_od_list(list_of_pointings);
    BOOST_CHECK_EQUAL(od_list.size(), 4);

    BOOST_CHECK_EQUAL(od_list[0].od, 91);
    BOOST_CHECK_EQUAL(od_list[0].first_pointing_id, 1);
    BOOST_CHECK_EQUAL(od_list[0].last_pointing_id, 4);
    BOOST_CHECK_EQUAL(od_list[0].num_of_pointings, 4);

    BOOST_CHECK_EQUAL(od_list[1].od, 92);
    BOOST_CHECK_EQUAL(od_list[1].first_pointing_id, 5);
    BOOST_CHECK_EQUAL(od_list[1].last_pointing_id, 5);
    BOOST_CHECK_EQUAL(od_list[1].num_of_pointings, 1);

    BOOST_CHECK_EQUAL(od_list[2].od, 93);
    BOOST_CHECK_EQUAL(od_list[2].first_pointing_id, 6);
    BOOST_CHECK_EQUAL(od_list[2].last_pointing_id, 7);
    BOOST_CHECK_EQUAL(od_list[2].num_of_pointings, 2);

    BOOST_CHECK_EQUAL(od_list[3].od, 94);
    BOOST_CHECK_EQUAL(od_list[3].first_pointing_id, 8);
    BOOST_CHECK_EQUAL(od_list[3].last_pointing_id, 8);
    BOOST_CHECK_EQUAL(od_list[3].num_of_pointings, 1);
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(MpiProcesses)
{
    const std::vector<Od_t> list_of_ods {
	// Od  firstPid  lastPid  numOfPids
	{   1,        1,     10,         10 },
	{   2,       11,     24,         14 },
	{   3,       25,     33,          9 },
	{   4,       34,     45,         11 },
	{   5,       46,     56,         11 },
	{   6,       57,     66,          8 },
	{   7,       67,     79,         13 },
	{   8,       80,     96,         17 },
	{   9,       97,    108,         12 }
    };

    BOOST_CHECK((splitOdsIntoMpiProcesses(2, list_of_ods) == 
		 std::vector<int>{ 105 }));
    BOOST_CHECK((splitOdsIntoMpiProcesses(4, list_of_ods) == 
		 std::vector<int>{ 55, 50 }));
    BOOST_CHECK((splitOdsIntoMpiProcesses(6, list_of_ods) == 
		 std::vector<int>{ 44, 49, 12 }));
    BOOST_CHECK((splitOdsIntoMpiProcesses(8, list_of_ods) == 
		 std::vector<int>{ 33, 30, 42, 0 }));
}
