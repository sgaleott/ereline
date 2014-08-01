#ifndef DATA_BINNING_HPP
#define DATA_BINNING_HPP

#include <vector>

#include "ahf_info.hpp"
#include "gain_table.hpp"
#include "healpix_map.hpp"
#include "planck_velocity.hpp"

/* This class contains data for one pointing period binned on a
 * Healpix map. It stores the pixels efficiently (as most of the sky
 * is not covered, keeping the full Healpix map in memory would be a
 * waste of space). */
struct Binned_data_t
{
    int quality_flag;
    int nside;
    int pointing_id;

    double max_dipole;
    double min_dipole;

    std::vector<int> pix_index;
    std::vector<int> pix_num_of_hits;
    std::vector<double> pix_data_sum;
    std::vector<float> pix_model_mean;

    Binned_data_t(int a_quality_flag,
                  int a_nside,
                  int a_pointing_id) :
        quality_flag(a_quality_flag),
        nside(a_nside),
        pointing_id(a_pointing_id),
        max_dipole(),
        min_dipole() { }

    double getDipolePeakToPeak() const { return max_dipole - min_dipole; }
};

struct Lfi_radiometer_t;
class Configuration;
struct Sqlite_connection_t;
struct Data_binning_results_t;

void run_data_binning(Sqlite_connection_t & ucds,
                      const Lfi_radiometer_t & rad,
                      Configuration & program_conf,
                      Configuration & storage_conf,
                      const std::vector<Pointing_t> & list_of_pointings,
                      Data_binning_results_t & result);

#endif
