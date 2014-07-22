#include "dipole_parameters.hpp"
#include "configuration.hpp"
#include <cmath>

////////////////////////////////////////////////////////////////////////////////

Dipole_parameters_t::Dipole_parameters_t(double ecl_dir_theta, 
					 double ecl_dir_phi, 
					 double a_solar_speed, 
					 double a_monopole)
    : solar_speed(a_solar_speed),
      monopole(a_monopole)
{
    axis[0] = std::sin(ecl_dir_theta) * std::cos(ecl_dir_phi);
    axis[1] = std::sin(ecl_dir_theta) * std::sin(ecl_dir_phi);
    axis[2] = std::cos(ecl_dir_theta);

    for(int i = 0; i < 3; ++i)
	solar_velocity[i] = axis[i] * solar_speed;
}

////////////////////////////////////////////////////////////////////////////////

Dipole_parameters_t
read_dipole_fit_params(const Configuration & conf)
{
    const double ecl_theta = 
	conf.get<double>("dipole_fit.solar_dipole.theta_ecl_rad");
    const double ecl_phi = 
	conf.get<double>("dipole_fit.solar_dipole.phi_ecl_rad");
    const double speed = 
	conf.get<double>("dipole_fit.solar_dipole.speed_m_s");
    const double monopole = 
	conf.get<double>("dipole_fit.solar_dipole.monopole_K_CMB");

    return Dipole_parameters_t(ecl_theta, ecl_phi, speed, monopole);
}
