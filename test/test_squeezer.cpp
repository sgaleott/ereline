#include "ahf_info.hpp"
#include "datatypes.hpp"
#include "squeezer.hpp"

#include <stdexcept>
#include <sstream>

#define BOOST_TEST_MODULE "Squeezer"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(Pointings)
{
    PointingData pointings;
    decompress_pointings(TEST_DATA_DIR "LFI27M_0118_pointings.sqz", pointings);

    BOOST_CHECK_EQUAL(pointings.obt_time.at(0),  106897764866642);
    BOOST_CHECK_CLOSE(pointings.scet_time.at(0), 1631130448754574.0, 1e-7);
    BOOST_CHECK_CLOSE(pointings.theta.at(0),     1.055340647697449, 1e-7);
    BOOST_CHECK_CLOSE(pointings.phi.at(0),       4.556172370910645, 1e-7);
    BOOST_CHECK_CLOSE(pointings.psi.at(0),       1.022354483604431, 1e-7);

    BOOST_CHECK_EQUAL(pointings.obt_time.at(9),  106897764884786);
    BOOST_CHECK_CLOSE(pointings.scet_time.at(9), 1631130449031430, 1e-7);
    BOOST_CHECK_CLOSE(pointings.theta.at(9),     1.083954110477151, 1e-7);
    BOOST_CHECK_CLOSE(pointings.phi.at(9),       4.550970192741829, 1e-7);
    BOOST_CHECK_CLOSE(pointings.psi.at(9),       1.025531066471899, 1e-7);

    BOOST_CHECK_EQUAL(pointings.obt_time.at(1686385),  106901164618802);
    BOOST_CHECK_CLOSE(pointings.scet_time.at(1686385), 1631182324885342, 1e-7);
    BOOST_CHECK_CLOSE(pointings.theta.at(1686385),     1.0321191264936829, 1e-7);
    BOOST_CHECK_CLOSE(pointings.phi.at(1686385),       1.2142055233289284, 1e-7);
    BOOST_CHECK_CLOSE(pointings.psi.at(1686385),       4.4921789040816051, 1e-7);
 
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(Datadiff)
{
    DifferencedData datadiff;
    decompress_differenced_data(TEST_DATA_DIR "LFI27M_0118_datadiff.sqz", datadiff);

    BOOST_CHECK_EQUAL(datadiff.obt_time.at(0),  106897764866642);
    BOOST_CHECK_CLOSE(datadiff.scet_time.at(0), 1631130448754574, 1e-7);
    BOOST_CHECK_CLOSE(datadiff.sky_load.at(0),  3.798198886215687e-04, 1e-7);
    BOOST_CHECK_EQUAL(datadiff.flags.at(0),     1105);

    BOOST_CHECK_EQUAL(datadiff.obt_time.at(9),  106897764884786);
    BOOST_CHECK_CLOSE(datadiff.scet_time.at(9), 1631130449031430, 1e-7);
    BOOST_CHECK_CLOSE(datadiff.sky_load.at(9),  1.214335439726710e-03, 1e-7);
    BOOST_CHECK_EQUAL(datadiff.flags.at(9),     1105);

    BOOST_CHECK_EQUAL(datadiff.obt_time.at(1686385),  106901164618802);
    BOOST_CHECK_CLOSE(datadiff.scet_time.at(1686385), 1631182324885342, 1e-7);
    BOOST_CHECK_CLOSE(datadiff.sky_load.at(1686385),  4.522699397057295e-03, 1e-7);
    BOOST_CHECK_EQUAL(datadiff.flags.at(1686385),     1089);
}

