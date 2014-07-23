#ifndef DIPOLE_FIT_RESULTS_HPP
#define DIPOLE_FIT_RESULTS_HPP

struct Dipole_fit_results_t {
    Gain_table_t gain_table;
    std::vector<Dipole_fit_t> list_of_fits;
    Healpix::Map_t<float> mask;
};

#endif
