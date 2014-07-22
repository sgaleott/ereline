#ifndef DIPOLE_PARAMETERS_HPP
#define DIPOLE_PARAMETERS_HPP

struct Dipole_parameters_t {
    // Axis of the dipole (length = 1)
    double axis[3];
    // Velocity vector (equal to axis * solar_speed)
    double solar_velocity[3]; 
    // Speed of the Sun wrt CMB rest frame [m/s]
    double solar_speed; 
    // Thermodynamic temperature of the monopole [K_CMB]
    double monopole; 

    Dipole_parameters_t(double ecl_dir_theta, 
			double ecl_dir_phi, 
			double a_solar_speed, 
			double a_monopole);
};

#include <iostream>
std::ostream & operator<<(std::ostream & os, Dipole_parameters_t const & yt);

struct Configuration;
Dipole_parameters_t read_dipole_fit_params(const Configuration & conf);

#endif
