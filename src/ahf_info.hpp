#ifndef AHF_INFO_HPP
#define AHF_INFO_HPP

struct Pointing_t {
    int pointing_id;
    double start_pointing;
    double start_time;
    double end_time;
    double spin_ecl_lon;
    double spin_ecl_lat;
    int od;
};

#endif
