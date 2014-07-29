#ifndef SMOOTH_GAINS
#define SMOOTH_GAINS

#include <vector>

class Configuration;
struct Pointing_t;
struct Sqlite_connection_t;
struct Lfi_radiometer_t;
struct Dipole_fit_results_t;
struct Da_capo_results_t;
struct Smooth_gains_results_t;

void run_smooth_gains(Sqlite_connection_t & ucds,
                      const Configuration & program_conf,
                      const Configuration & storage_conf,
                      const Lfi_radiometer_t & rad,
                      const std::vector<Pointing_t> & list_of_pointings,
                      const Dipole_fit_results_t & fit_results,
                      const Da_capo_results_t & da_capo_results,
                      Smooth_gains_results_t & smooth_results);

#endif
