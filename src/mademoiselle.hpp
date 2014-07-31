#ifndef MADEMOISELLE_HPP
#define MADEMOISELLE_HPP

#include <vector>
#include "dipole_fit.hpp"

class Configuration;
struct Dipole_fit_results_t;
struct Da_capo_results_t;
struct Lfi_radiometer_t;

void run_mademoiselle(const Configuration & program_conf,
                      const Configuration & storage_conf,
                      const Lfi_radiometer_t & user_rad,
                      const std::vector<Pointing_t> & list_of_pointings,
                      Dipole_fit_results_t & dipole_fit_results,
                      Da_capo_results_t & da_capo_results);

#endif
