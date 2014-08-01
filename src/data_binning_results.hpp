#ifndef DATA_BINNING_RESULTS_HPP
#define DATA_BINNING_RESULTS_HPP

#include <vector>
#include "healpix_map.hpp"
#include "data_binning.hpp"

struct Data_binning_results_t {
    std::vector<Binned_data_t> binned_pids;
    Healpix::Map_t<float> mask;

    // This is needed to figure out which gains in "gain_table"
    // belong to the M/S arms.
    std::vector<int> pids_per_process;
};

#endif
