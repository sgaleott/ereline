#ifndef DIPOLE_FIT_RESULTS_HPP
#define DIPOLE_FIT_RESULTS_HPP

#include "data_binning_results.hpp"
#include "dipole_fit.hpp"
#include "gain_table.hpp"
#include "healpix_map.hpp"

#include <vector>

struct Dipole_fit_results_t {
    std::vector<Dipole_fit_t> dipole_fits;
    Healpix::Map_t<float> mask;

    // This is needed to figure out which gains in "gain_table"
    // belong to the M/S arms.
    std::vector<int> pids_per_process;
};

struct Lfi_radiometer_t;

void save_dipole_fit_results(const std::string & file_name,
                             const Lfi_radiometer_t & radiometer,
                             const Dipole_fit_results_t & results,
                             const std::string & comment = "");

#endif
