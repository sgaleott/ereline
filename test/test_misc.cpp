#include "mpi_processes.hpp"

#include <stdexcept>
#include <sstream>

#define BOOST_TEST_MODULE "Miscellaneous"
#include <boost/test/unit_test.hpp>

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
