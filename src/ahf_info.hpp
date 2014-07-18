#ifndef AHF_INFO_HPP
#define AHF_INFO_HPP

#include <vector>
#include "datatypes.hpp"

struct Pointing_t {
    int id;
    double start_pointing;
    double start_time;
    double end_time;
    double spin_ecl_lon;
    double spin_ecl_lat;
    int od;
};

struct Od_t {
    int od;
    int first_pointing_id;
    int last_pointing_id;
    int num_of_pointings;
};

std::vector<Od_t> build_od_list(const std::vector<Pointing_t> & pointings);
void get_pid_iterators_for_range(const std::vector<Pointing_t> & list_of_pointings,
				 const Range_t<int> & pid_range,
				 std::vector<Pointing_t>::const_iterator & first,
				 std::vector<Pointing_t>::const_iterator & last);

#endif
