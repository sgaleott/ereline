#ifndef DATA_BINNING_RESULTS_HPP
#define DATA_BINNING_RESULTS_HPP

#include <vector>
#include "healpix_map.hpp"
#include "data_binning.hpp"

class Configuration;

struct Data_binning_results_t {
    Lfi_radiometer_t radiometer;
    std::vector<Binned_data_t> binned_pids;
    Healpix::Map_t<float> mask;

    // These are needed to figure out which gains in "gain_table"
    // belong to the M/S arms.
    int mpi_size;
    std::vector<int> pids_per_process;

    void save_to_disk(Configuration & program_conf) const;
    void load_from_disk(const Lfi_radiometer_t &a_radiometer,
                        Configuration & program_conf);
};

#endif
