#include "io.hpp"
#include "ahf_info.hpp"
#include "datatypes.hpp"
#include "planck_velocity.hpp"
#include "sqlite3xx.hpp"

#include <stdexcept>
#include <sstream>

#define BOOST_TEST_MODULE "I/O"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

////////////////////////////////////////////////////////////////////////////////

#include "logging.hpp"
BOOST_AUTO_TEST_CASE(load_convolved_dipole)
{
    Logger * log = Logger::get_instance();
    log->set_log_level(Logger::Log_level::DEBUG);

    SQLite3Connection ucds(TEST_DATA_DIR "/ucds.db");
    std::vector<Pointing_t> pointings;
    PlanckVelocity vel;
    loadConvolutionParametersFromUCDS(ucds, LfiRadiometer("LFI24M"), vel);

    BOOST_CHECK_CLOSE(vel.M100,  0.00119669520223024,  1e-7);
    BOOST_CHECK_CLOSE(vel.M010, -1.35903112318943e-06, 1e-7);
    BOOST_CHECK_CLOSE(vel.M001,  0.998945455192864,    1e-7);
    BOOST_CHECK_CLOSE(vel.M200,  0.00106781574504692,  1e-7);
    BOOST_CHECK_CLOSE(vel.M110, -2.31094686909285e-06, 1e-7);
    BOOST_CHECK_CLOSE(vel.M101,  0.000231187225785091, 1e-7);
    BOOST_CHECK_CLOSE(vel.M020,  0.000136825401708578, 1e-7);
    BOOST_CHECK_CLOSE(vel.M011,  8.95456618142487e-07, 1e-7);
    BOOST_CHECK_CLOSE(vel.M002,  0.998795358853245,    1e-7);
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(load_pointings)
{
    SQLite3Connection ucds(TEST_DATA_DIR "/ucds.db");
    std::vector<Pointing_t> pointings;
    loadPointingInformation(ucds, 943, 945, pointings);

    BOOST_CHECK_EQUAL(pointings.size(), 5);

    BOOST_CHECK_EQUAL(pointings[0].pointing_id, 25619);
    BOOST_CHECK_EQUAL(pointings[0].start_pointing, 111568338130683);
    BOOST_CHECK_EQUAL(pointings[0].start_time, 111568338294523);
    BOOST_CHECK_EQUAL(pointings[0].end_time, 111570160752379);
    BOOST_CHECK_CLOSE(pointings[0].spin_ecl_lon, 77.7838, 1e-7);
    BOOST_CHECK_CLOSE(pointings[0].spin_ecl_lat, 0.0024, 1e-7);
    BOOST_CHECK_EQUAL(pointings[0].od, 943);

    BOOST_CHECK_EQUAL(pointings[1].pointing_id, 25620);
    BOOST_CHECK_EQUAL(pointings[1].start_pointing, 111570160752379);
    BOOST_CHECK_EQUAL(pointings[1].start_time, 111570278770430);
    BOOST_CHECK_EQUAL(pointings[1].end_time, 111570278770430);
    BOOST_CHECK_CLOSE(pointings[1].spin_ecl_lon, 0, 1e-7);
    BOOST_CHECK_CLOSE(pointings[1].spin_ecl_lat, 0, 1e-7);
    BOOST_CHECK_EQUAL(pointings[1].od, 943);

    BOOST_CHECK_EQUAL(pointings[2].pointing_id, 25621);
    BOOST_CHECK_EQUAL(pointings[2].start_pointing, 111570278770430);
    BOOST_CHECK_EQUAL(pointings[2].start_time, 111570390734587);
    BOOST_CHECK_EQUAL(pointings[2].end_time, 111573961135874);
    BOOST_CHECK_CLOSE(pointings[2].spin_ecl_lon, 82.0504, 1e-7);
    BOOST_CHECK_CLOSE(pointings[2].spin_ecl_lat, 0.0048, 1e-7);
    BOOST_CHECK_EQUAL(pointings[2].od, 943);

    BOOST_CHECK_EQUAL(pointings[3].pointing_id, 25622);
    BOOST_CHECK_EQUAL(pointings[3].start_pointing, 111573961135874);
    BOOST_CHECK_EQUAL(pointings[3].start_time, 111573961283330);
    BOOST_CHECK_EQUAL(pointings[3].end_time, 111579584059131);
    BOOST_CHECK_CLOSE(pointings[3].spin_ecl_lon, 82.051, 1e-7);
    BOOST_CHECK_CLOSE(pointings[3].spin_ecl_lat, 0.0108, 1e-7);
    BOOST_CHECK_EQUAL(pointings[3].od, 944);

    BOOST_CHECK_EQUAL(pointings[4].pointing_id, 25623);
    BOOST_CHECK_EQUAL(pointings[4].start_pointing, 111579584059131);
    BOOST_CHECK_EQUAL(pointings[4].start_time, 111579584272123);
    BOOST_CHECK_EQUAL(pointings[4].end_time, 111585207047931);
    BOOST_CHECK_CLOSE(pointings[4].spin_ecl_lon, 82.0518, 1e-7);
    BOOST_CHECK_CLOSE(pointings[4].spin_ecl_lat, 0.011, 1e-7);
    BOOST_CHECK_EQUAL(pointings[4].od, 945);
}
