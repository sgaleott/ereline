#include "ahf_info.hpp"
#include "misc.hpp"
#include "mpi_processes.hpp"
#include "datatypes.hpp"

#include <stdexcept>
#include <sstream>

#define BOOST_TEST_MODULE "Miscellaneous"
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

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

BOOST_AUTO_TEST_CASE(FindRange)
{
    std::vector<int> vec { 0, 10, 20, 30, 40, 50, 60 };
    int first_idx, last_idx;

    find_range_of_indexes(vec, 15, 35, first_idx, last_idx);
    BOOST_CHECK_EQUAL(first_idx, 2);
    BOOST_CHECK_EQUAL(last_idx, 3);

    find_range_of_indexes(vec, 40, 50, first_idx, last_idx);
    BOOST_CHECK_EQUAL(first_idx, 4);
    BOOST_CHECK_EQUAL(last_idx, 5);

    find_range_of_indexes(vec, 100, 101, first_idx, last_idx);
    BOOST_CHECK_EQUAL(first_idx, -1);
    BOOST_CHECK_EQUAL(last_idx, -1);
}

////////////////////////////////////////////////////////////////////////////////

#include <iostream>
BOOST_AUTO_TEST_CASE(MpiProcesses)
{
    const std::vector<Od_t> list_of_ods {
	// Od  firstPid  lastPid  numOfPids
	{   1,        1,     10,         10 },
	{   2,       11,     24,         14 },
	{   3,       25,     33,          9 },
	{   4,       34,     44,         11 },
	{   5,       45,     55,         11 },
	{   6,       56,     65,          8 },
	{   7,       66,     78,         13 },
	{   8,       79,     95,         17 },
	{   9,       96,    107,         12 }
    };

    std::vector<Data_range_t> result;
    std::vector<Data_range_t> expected;

    result = splitOdsIntoMpiProcesses(2, list_of_ods);
    expected = std::vector<Data_range_t> { { { 1, 9 }, { 1, 107 }, 9, 105 } };
    BOOST_CHECK_EQUAL_COLLECTIONS(result.begin(), result.end(),
				  expected.begin(), expected.end());

    result = splitOdsIntoMpiProcesses(4, list_of_ods);
    expected = std::vector<Data_range_t>{ { { 1, 5 }, { 1, 55 }, 5, 55 },
					  { { 6, 9 }, { 56, 107 }, 4, 50 } };
    BOOST_CHECK_EQUAL_COLLECTIONS(result.begin(), result.end(),
				  expected.begin(), expected.end());

    result = splitOdsIntoMpiProcesses(6, list_of_ods);
    expected = std::vector<Data_range_t>{ { { 1, 4 }, { 1, 44 }, 4, 44 },
					  { { 5, 8 }, { 45, 95 }, 4, 49 },
					  { { 9, 9 }, { 96, 107 }, 1, 12 } };
    BOOST_CHECK_EQUAL_COLLECTIONS(result.begin(), result.end(),
				  expected.begin(), expected.end());

    result = splitOdsIntoMpiProcesses(8, list_of_ods);
    expected = std::vector<Data_range_t>{ { { 1, 3 }, { 1, 33 }, 3, 33 },
					  { { 4, 6 }, { 34, 65 }, 3, 30 },
					  { { 7, 9 }, { 66, 107 }, 3, 42 },
					  { { -1, -1 }, { -1, -11 }, 0, 0 } };
    BOOST_CHECK_EQUAL_COLLECTIONS(result.begin(), result.end(),
				  expected.begin(), expected.end());
}
