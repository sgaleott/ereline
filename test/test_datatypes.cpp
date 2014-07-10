#include "datatypes.hpp"

#define BOOST_TEST_MODULE "Datatypes"
#include <boost/test/unit_test.hpp>

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
