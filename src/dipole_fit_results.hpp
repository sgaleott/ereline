#ifndef DIPOLE_FIT_RESULTS_HPP
#define DIPOLE_FIT_RESULTS_HPP

#include "dipole_fit.hpp"
#include "gain_table.hpp"
#include "healpix_map.hpp"

#include <vector>

struct Dipole_fit_results_t {
    Gain_table_t gain_table;
    std::vector<Dipole_fit_t> list_of_fits;
    Healpix::Map_t<float> mask;

    // This is needed to figure out which gains in "gain_table"
    // belong to the M/S arms.
    std::vector<int> pids_per_process;
};

#endif
