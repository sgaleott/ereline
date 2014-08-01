#ifndef _DIPOLE_FIT_
#define _DIPOLE_FIT_

#include <vector>

#include "ahf_info.hpp"
#include "data_binning.hpp"
#include "gain_table.hpp"
#include "healpix_map.hpp"
#include "planck_velocity.hpp"

struct Lfi_radiometer_t;

/* class calib iter*/
struct Dipole_fit_t
{
    Binned_data_t binned_data;
    double gainv;
    double offset;

    Dipole_fit_t(const Binned_data_t & a_binned_data) :
        binned_data(a_binned_data),
        gainv(),
        offset() {}
};

class Configuration;
struct Sqlite_connection_t;
struct Data_binning_results_t;
struct Dipole_fit_results_t;

void run_dipole_fit(Sqlite_connection_t & ucds,
                    const Lfi_radiometer_t & rad,
                    Configuration & program_conf,
                    Configuration & storage_conf,
                    const std::vector<Pointing_t> & list_of_pointings,
                    const Data_binning_results_t & binned_data,
                    Dipole_fit_results_t & result);

#endif
