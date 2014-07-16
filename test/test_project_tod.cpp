#include "project_tod.hpp"

#include <vector>
#include <stdexcept>
#include <sstream>

#define BOOST_TEST_MODULE "Project TOD"
#include <boost/test/unit_test.hpp>

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(ProjectTOD)
{
    PointingData pnt_data {
	std::vector<uint64_t> { 0, 1, 2, 3, 4, 5, 6, 7, 8 },
        std::vector<double> { 0., 1., 2., 3., 4., 5., 6., 7., 8., },
        std::vector<double> { 0., 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08 },
	std::vector<double> { 0., 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04 },
        std::vector<double> { 0., 0., 0., 0., 0., 0., 0., 0., 0. }
    };

    DifferencedData diff_data {
	std::vector<uint64_t> { 0, 1, 2, 3, 4, 5, 6, 7, 8 },
        std::vector<double> { 0., 1., 2., 3., 4., 5., 6., 7., 8., },
        std::vector<double> { 10., 20., 30., 40., 50., 60., 70., 80., 90. },
        std::vector<uint32_t> { 1, 1, 1, 2, 1, 1, 1, 1 }
    };

    Projected_tod_t proj_tod;
    proj_tod.pointing_id = 1;
    proj_tod.nside = 64;
    Range_t<double> extrema;
    project_tod(2, Range_t<uint64_t> { 1, 7 }, pnt_data, diff_data, extrema, proj_tod);

    BOOST_CHECK_EQUAL(proj_tod.pixel_idx.size(), 5);
    BOOST_CHECK_EQUAL(proj_tod.sum.size(), 5);
    BOOST_CHECK_EQUAL(proj_tod.hits.size(), 5);

    BOOST_CHECK_EQUAL(proj_tod.pixel_idx.at(0), 4078);
    BOOST_CHECK_EQUAL(proj_tod.pixel_idx.at(1), 4079);
    BOOST_CHECK_EQUAL(proj_tod.pixel_idx.at(2), 4090);
    BOOST_CHECK_EQUAL(proj_tod.pixel_idx.at(3), 4094);
    BOOST_CHECK_EQUAL(proj_tod.pixel_idx.at(4), 4095);

    BOOST_CHECK_EQUAL(proj_tod.sum.at(0), 80.0);
    BOOST_CHECK_EQUAL(proj_tod.sum.at(1), 70.0);
    BOOST_CHECK_EQUAL(proj_tod.sum.at(2), 110.0);
    BOOST_CHECK_EQUAL(proj_tod.sum.at(3), 30.0);
    BOOST_CHECK_EQUAL(proj_tod.sum.at(4), 20.0);

    BOOST_CHECK_EQUAL(proj_tod.hits.at(0), 1);
    BOOST_CHECK_EQUAL(proj_tod.hits.at(1), 1);
    BOOST_CHECK_EQUAL(proj_tod.hits.at(2), 2);
    BOOST_CHECK_EQUAL(proj_tod.hits.at(3), 1);
    BOOST_CHECK_EQUAL(proj_tod.hits.at(4), 1);
}
