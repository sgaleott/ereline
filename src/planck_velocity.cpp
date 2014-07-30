#include "planck_velocity.hpp"
#include "configuration.hpp"
#include "dipole_parameters.hpp"
#include "logging.hpp"
#include "rotmatrix.hpp"
#include "fits_object.hpp"

#include <cmath>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////

// Open object and read velocities
Planck_velocity_t::Planck_velocity_t (const std::string & file_name,
                                      const Dipole_parameters_t & a_dipole_params)
    : dipole_params(a_dipole_params)
{
    Logger * log = Logger::get_instance();
    std::stringstream s;
    s << a_dipole_params;
    log->debug("Planck_velocity_t initialized with the following dipole: "
               + s.str());

    FitsObject velocityFile;
    velocityFile.openTable(file_name);

    int numOfRows;
    velocityFile.getKey("NAXIS2", numOfRows);

    velocityFile.getColumn("SCET", scet, 1, numOfRows);
    velocityFile.getColumn("XVEL", xvel, 1, numOfRows);
    velocityFile.getColumn("YVEL", yvel, 1, numOfRows);
    velocityFile.getColumn("ZVEL", zvel, 1, numOfRows);
}

// Return velocity at given SCET
std::vector<double>
Planck_velocity_t::getVelocity (double scetTime) const
{
  std::vector<double> vSat(3,0.0);

  size_t idx = searchArr(scet, scetTime);
  if (idx < (scet.size()-1))
    {
      vSat[0] = xvel[idx]+((xvel[idx+1]-xvel[idx])/(scet[idx+1]-scet[idx]))*(scetTime-scet[idx]);
      vSat[1] = yvel[idx]+((yvel[idx+1]-yvel[idx])/(scet[idx+1]-scet[idx]))*(scetTime-scet[idx]);
      vSat[2] = zvel[idx]+((zvel[idx+1]-zvel[idx])/(scet[idx+1]-scet[idx]))*(scetTime-scet[idx]);
    }
  else
    {
      vSat[0] = xvel[idx-1]+((xvel[idx]-xvel[idx-1])/(scet[idx]-scet[idx-1]))*(scetTime-scet[idx-1]);
      vSat[1] = yvel[idx-1]+((yvel[idx]-yvel[idx-1])/(scet[idx]-scet[idx-1]))*(scetTime-scet[idx-1]);
      vSat[2] = zvel[idx-1]+((zvel[idx]-zvel[idx-1])/(scet[idx]-scet[idx-1]))*(scetTime-scet[idx-1]);
    }

  return vSat;
}

std::vector<double>
Planck_velocity_t::getAbsoluteVelocity(double scetTime) const
{
  std::vector<double> vSat = getVelocity(scetTime);

  // Add the solar system velocity
  for(int i = 0; i < 3; ++i)
      vSat[i] += dipole_params.solar_velocity[i];

  return vSat;
}

double
Planck_velocity_t::dipole(const std::vector<double> & velocity,
                          double theta,
                          double phi) const
{
  double velocityLength = sqrt(velocity[0]*velocity[0] +
                               velocity[1]*velocity[1] +
                               velocity[2]*velocity[2]);

  double velocity_versor[3];
  velocity_versor[0] = velocity[0]/velocityLength;
  velocity_versor[1] = velocity[1]/velocityLength;
  velocity_versor[2] = velocity[2]/velocityLength;

  double beta = velocityLength/SPEED_OF_LIGHT;
  double gamma = 1./sqrt(1.-beta*beta);

  double detDir[3];
  angToCart(theta, phi, detDir);

  double cosDir = (velocity_versor[0]*detDir[0] +
                   velocity_versor[1]*detDir[1] +
                   velocity_versor[2]*detDir[2]);

  return (1./(gamma*(1. - beta*cosDir )) -1.) * dipole_params.monopole;
}

double
Planck_velocity_t::convolvedDipole (const std::vector<double> & velocity,
                                    double theta,
                                    double phi,
                                    double psi) const
{
  double xv = velocity[0]/SPEED_OF_LIGHT;
  double yv = velocity[1]/SPEED_OF_LIGHT;
  double zv = velocity[2]/SPEED_OF_LIGHT;

  // Rotate velocity
  double x1 = cos(phi)*xv+sin(phi)*yv;
  double y1 = -sin(phi)*xv+cos(phi)*yv;
  double z1 = zv;

  double x2 = cos(theta)*x1-sin(theta)*z1;
  double y2 = y1;
  double z2 = sin(theta)*x1+cos(theta)*z1;

  double x3 = cos(psi)*x2+sin(psi)*y2;
  double y3 = -sin(psi)*x2+cos(psi)*y2;
  double z3 = z2;

  // Compute dipole amplitude
  double scalarProduct = (x3*M100+y3*M010+z3*M001);

  double relativisticCorrection = (M200*x3+M110*y3+M101*z3)*x3+
    (M110*x3+M020*y3+M011*z3)*y3+
    (M101*x3+M011*y3+M002*z3)*z3;

  double totalDipole = scalarProduct + relativisticCorrection;
  return totalDipole * dipole_params.monopole;
}

double
Planck_velocity_t::getConvolvedDipole(double scetTime, double theta, double phi, double psi) const
{
  // Compute velocity versor
  const std::vector<double> vSatAbsolute = getAbsoluteVelocity(scetTime);
  return convolvedDipole (vSatAbsolute, theta, phi, psi);
}

std::vector<double>
Planck_velocity_t::getConvolvedDipole(const std::vector<double> & scetTime,
                                      const std::vector<double> & theta,
                                      const std::vector<double> & phi,
                                      const std::vector<double> & psi) const
{
  std::vector<double> local_dipole(scetTime.size());
  for (size_t idx = 0; idx < local_dipole.size(); ++idx)
    {
      std::vector<double> vSatAbsolute = getAbsoluteVelocity(scetTime[idx]);
      local_dipole[idx] = convolvedDipole(vSatAbsolute, theta[idx], phi[idx], psi[idx]);
    }

  return local_dipole;
}

void
Planck_velocity_t::use_pencil_beam()
{
    M100 = M010 = M200 = M110 = M101 = M020 = M011 = 0.0;
    M001 = M002 = 1.0;
}
