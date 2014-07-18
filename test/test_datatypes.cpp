#include "datatypes.hpp"

#define BOOST_TEST_MODULE "Datatypes"
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(lfi_radiometer_type)
{
    LfiRadiometer rad1(18, 1);
    LfiRadiometer rad2("LFI24S");
    LfiRadiometer rad3("LFI26M-00");
    LfiRadiometer rad4("28S");
    LfiRadiometer rad5("LFI21M-0");
    LfiRadiometer rad6;

    BOOST_CHECK_EQUAL(rad1.horn, 18); BOOST_CHECK_EQUAL(rad1.radiometer, 1);
    BOOST_CHECK_EQUAL(rad2.horn, 24); BOOST_CHECK_EQUAL(rad2.radiometer, 1);
    BOOST_CHECK_EQUAL(rad3.horn, 26); BOOST_CHECK_EQUAL(rad3.radiometer, 0);
    BOOST_CHECK_EQUAL(rad4.horn, 28); BOOST_CHECK_EQUAL(rad4.radiometer, 1);
    BOOST_CHECK_EQUAL(rad5.horn, 21); BOOST_CHECK_EQUAL(rad5.radiometer, 0);

    BOOST_CHECK_THROW(rad6.assign("dummy"), ConfigurationError);
    BOOST_CHECK_THROW(rad6.assign("LFI16M"), ConfigurationError);
    BOOST_CHECK_THROW(rad6.assign("LFI26x"), ConfigurationError);
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(pointings)
{
    PointingData pointings;
    pointings.obt_time = std::vector<uint64_t> { 1, 2, 3, 4, 5 };
    pointings.scet_time = std::vector<double> { 11., 12., 13., 14., 15. };
    pointings.theta = std::vector<double> { 101., 102., 103., 104., 105. };
    pointings.phi = std::vector<double> { 201., 202., 203., 204., 205. };
    pointings.psi = std::vector<double> { 301., 302., 303., 304., 305. };

    const PointingData subset(pointings, 2, 4);

    PointingData expected;
    expected.obt_time = std::vector<uint64_t> { 2, 3, 4 };
    expected.scet_time = std::vector<double> { 12., 13., 14. };
    expected.theta = std::vector<double> { 102., 103., 104. };
    expected.phi = std::vector<double> { 202., 203., 204. };
    expected.psi = std::vector<double> { 302., 303., 304. };

    BOOST_CHECK_EQUAL_COLLECTIONS(subset.obt_time.begin(), subset.obt_time.end(),
				  expected.obt_time.begin(), expected.obt_time.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(subset.scet_time.begin(), subset.scet_time.end(),
				  expected.scet_time.begin(), expected.scet_time.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(subset.theta.begin(), subset.theta.end(),
				  expected.theta.begin(), expected.theta.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(subset.phi.begin(), subset.phi.end(),
				  expected.phi.begin(), expected.phi.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(subset.psi.begin(), subset.psi.end(),
				  expected.psi.begin(), expected.psi.end());
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(datadiff)
{
    DifferencedData datadiff;
    datadiff.obt_time = std::vector<uint64_t> { 1, 2, 3, 4, 5 };
    datadiff.scet_time = std::vector<double> { 11., 12., 13., 14., 15. };
    datadiff.sky_load = std::vector<double> { 101., 102., 103., 104., 105. };
    datadiff.flags = std::vector<uint32_t> { 201, 202, 203, 204, 205 };

    const DifferencedData subset(datadiff, 2, 4);

    DifferencedData expected;
    expected.obt_time = std::vector<uint64_t> { 2, 3, 4 };
    expected.scet_time = std::vector<double> { 12., 13., 14. };
    expected.sky_load = std::vector<double> { 102., 103., 104. };
    expected.flags = std::vector<uint32_t> { 202, 203, 204 };

    BOOST_CHECK_EQUAL_COLLECTIONS(subset.obt_time.begin(), subset.obt_time.end(),
				  expected.obt_time.begin(), expected.obt_time.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(subset.scet_time.begin(), subset.scet_time.end(),
				  expected.scet_time.begin(), expected.scet_time.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(subset.sky_load.begin(), subset.sky_load.end(),
				  expected.sky_load.begin(), expected.sky_load.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(subset.flags.begin(), subset.flags.end(),
				  expected.flags.begin(), expected.flags.end());
}
