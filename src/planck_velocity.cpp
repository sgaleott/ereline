#include <cmath>

#include "planck_velocity.hpp"
#include "rotmatrix.hpp"
#include "fits_object.hpp"

// Open object and read velocities
PlanckVelocity::PlanckVelocity (const std::string & file_name)
{
    FitsObject velocityFile;
    velocityFile.openTable(file_name);
    
    int numOfRows;
    velocityFile.getKey("NAXIS2", numOfRows);

    velocityFile.getColumn("SCET", scet, 1, numOfRows);
    velocityFile.getColumn("XVEL", xvel, 1, numOfRows);
    velocityFile.getColumn("YVEL", yvel, 1, numOfRows);
    velocityFile.getColumn("ZVEL", zvel, 1, numOfRows);

    // Initialize common variables
    std::vector<double> SOLSYSDIR_V=angToCart(SOLSYSDIR_ECL_THETA, SOLSYSDIR_ECL_PHI);
    vSolSys.push_back(SOLSYSDIR_V[0]*SOLSYSSPEED);
    vSolSys.push_back(SOLSYSDIR_V[1]*SOLSYSSPEED);
    vSolSys.push_back(SOLSYSDIR_V[2]*SOLSYSSPEED);
}

// Initialize V SolSys only
PlanckVelocity::PlanckVelocity ()
{
  // Initialize common variables
  std::vector<double> SOLSYSDIR_V=angToCart(SOLSYSDIR_ECL_THETA, SOLSYSDIR_ECL_PHI);
  vSolSys.push_back(SOLSYSDIR_V[0]*SOLSYSSPEED);
  vSolSys.push_back(SOLSYSDIR_V[1]*SOLSYSSPEED);
  vSolSys.push_back(SOLSYSDIR_V[2]*SOLSYSSPEED);
}

// Return velocity at given SCET
std::vector<double> 
PlanckVelocity::getVelocity (double scetTime)
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
PlanckVelocity::getAbsoluteVelocity(double scetTime)
{
  std::vector<double> vSat = getVelocity(scetTime);

  // Add the solar system velocity
  vSat[0] += vSolSys[0];
  vSat[1] += vSolSys[1];
  vSat[2] += vSolSys[2];
  
  return vSat;
}

double 
PlanckVelocity::dipole(std::vector<double> velocity, double theta, double phi)
{
  double velocityLength = sqrt(velocity[0]*velocity[0] + 
			       velocity[1]*velocity[1] + 
			       velocity[2]*velocity[2]);

  std::vector<double> velocity_versor(3,0.0);
  velocity_versor[0] = velocity[0]/velocityLength;
  velocity_versor[1] = velocity[1]/velocityLength;
  velocity_versor[2] = velocity[2]/velocityLength;
  
  double beta = velocityLength/SPEED_OF_LIGHT;
  double gamma = 1./sqrt(1.-beta*beta);
  
  std::vector<double> detDir = angToCart(theta, phi);
    
  double cosDir = (velocity_versor[0]*detDir[0] +
		   velocity_versor[1]*detDir[1] +
		   velocity_versor[2]*detDir[2]);
  
  return (1./(gamma*(1. - beta*cosDir )) -1.)*TCMB;
  //return cosDir*TCMB*beta;
}

double 
PlanckVelocity::getTotalDipole(double scetTime, double theta, double phi)
{  
  std::vector<double> vSatAbsolute = getAbsoluteVelocity(scetTime);

  return dipole(vSatAbsolute, theta, phi);
}

std::vector<double> 
PlanckVelocity::getTotalDipole(const std::vector<double> & scetTime, const std::vector<double> & theta, 
			       const std::vector<double> & phi)
{
  std::vector<double> local_dipole;
  for (size_t idx=0; idx<scetTime.size(); ++idx)
    {
      std::vector<double> vSatAbsolute = getAbsoluteVelocity(scetTime[idx]);
      local_dipole.push_back(dipole(vSatAbsolute,theta[idx],phi[idx]));
    }
  return local_dipole;
}

double 
PlanckVelocity::getSolarDipole(double theta, double phi)
{  
  return dipole(vSolSys, theta, phi);
}

std::vector<double> 
PlanckVelocity::getSolarDipole(const std::vector<double> & theta, const std::vector<double> & phi)
{
  std::vector<double> local_dipole;
  for (size_t idx=0; idx<theta.size(); ++idx)
    {
      local_dipole.push_back(dipole(vSolSys, theta[idx], phi[idx]));
    }

  return local_dipole;
}

double 
PlanckVelocity::getOrbitalDipole(double scetTime, double theta, double phi)
{  
  std::vector<double> vSat = getVelocity(scetTime);

  return dipole(vSat, theta, phi);
}

double 
PlanckVelocity::getTotalDipoleTAnt(double scetTime, double theta, double phi, double hnydk)
{
  double totalDip = getTotalDipole(scetTime, theta, phi);
  totalDip=totalDip+TCMB;

  double tcmbAnt = hnydk/(exp(hnydk/TCMB)-1);
  double totalDipAnt = hnydk/(exp(hnydk/totalDip)-1);

  return totalDipAnt-tcmbAnt;
}

double 
PlanckVelocity::getSolarDipoleTAnt(double theta, double phi, double hnydk)
{
  double solarDip = getSolarDipole(theta, phi);
  solarDip=solarDip+TCMB;

  double tcmbAnt = hnydk/(exp(hnydk/TCMB)-1);
  double solarDipAnt = hnydk/(exp(hnydk/solarDip)-1);

  return solarDipAnt-tcmbAnt;
}

double 
PlanckVelocity::getOrbitalDipoleTAnt(double scetTime, double theta, double phi, double hnydk)
{
  double orbitalDip = getOrbitalDipole(scetTime, theta, phi);
  orbitalDip=orbitalDip+TCMB;

  double tcmbAnt = hnydk/(exp(hnydk/TCMB)-1);
  double orbitalDipAnt = hnydk/(exp(hnydk/orbitalDip)-1);

  return orbitalDipAnt-tcmbAnt;
}

double 
PlanckVelocity::convolvedDipole (std::vector<double> velocity, double theta, double phi, double psi)
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

  double totalDipole = scalarProduct+relativisticCorrection;
  return totalDipole*TCMB;
}

double 
PlanckVelocity::getConvolvedDipole(double scetTime, double theta, double phi, double psi)
{  
  // Compute velocity versor
  std::vector<double> vSatAbsolute = getAbsoluteVelocity(scetTime);
  return convolvedDipole (vSatAbsolute,theta,phi,psi);
}

std::vector<double> 
PlanckVelocity::getConvolvedDipole(const std::vector<double> & scetTime, const std::vector<double> & theta, 
				   const std::vector<double> & phi, const std::vector<double> & psi)
{
  std::vector<double> local_dipole;
  for (size_t idx=0; idx<scetTime.size(); ++idx)
    {
      std::vector<double> vSatAbsolute = getAbsoluteVelocity(scetTime[idx]);
      local_dipole.push_back(convolvedDipole(vSatAbsolute,theta[idx],phi[idx],psi[idx]));
    }
  return local_dipole;
}

double 
PlanckVelocity::getSolarConvolvedDipole(double theta, double phi, double psi)
{  
  return convolvedDipole (vSolSys,theta,phi,psi);
}

std::vector<double> 
PlanckVelocity::getSolarConvolvedDipole(const std::vector<double> & theta, const std::vector<double> & phi, const std::vector<double> & psi)
{
  std::vector<double> local_dipole;
  for (size_t idx=0; idx<theta.size(); ++idx)
    {
      local_dipole.push_back(convolvedDipole(vSolSys,theta[idx],phi[idx],psi[idx]));
    }
  return local_dipole;
}

double 
PlanckVelocity::getOrbitalConvolvedDipole(double scetTime, double theta, double phi, double psi)
{  
  // Compute velocity versor
  std::vector<double> vSatAbsolute = getVelocity(scetTime);
  return convolvedDipole (vSatAbsolute,theta,phi,psi);
}

std::vector<double>
PlanckVelocity::getOrbitalConvolvedDipole(const std::vector<double> & scetTime, const std::vector<double> & theta, 
					  const std::vector<double> & phi, const std::vector<double> & psi)
{  
  std::vector<double> local_dipole;
  for (size_t idx=0; idx<scetTime.size(); ++idx)
    {
      std::vector<double> vSatAbsolute = getVelocity(scetTime[idx]);
      local_dipole.push_back(convolvedDipole(vSatAbsolute,theta[idx],phi[idx],psi[idx]));
    }
  return local_dipole;
}

double 
PlanckVelocity::getConvolvedDipoleTAnt(double scetTime, double theta, double phi, double psi, double hnydk)
{
  double totalDip = getConvolvedDipole(scetTime, theta, phi, psi);
  totalDip=totalDip+TCMB;

  double tcmbAnt = hnydk/(exp(hnydk/TCMB)-1);
  double totalDipAnt = hnydk/(exp(hnydk/totalDip)-1);

  return totalDipAnt-tcmbAnt;
}

double 
PlanckVelocity::getSolarConvolvedDipoleTAnt(double theta, double phi, double psi, double hnydk)
{
  double totalDip = getSolarConvolvedDipole(theta, phi, psi);
  totalDip=totalDip+TCMB;

  double tcmbAnt = hnydk/(exp(hnydk/TCMB)-1);
  double totalDipAnt = hnydk/(exp(hnydk/totalDip)-1);

  return totalDipAnt-tcmbAnt;
}

double 
PlanckVelocity::getOrbitalConvolvedDipoleTAnt(double scetTime, double theta, double phi, double psi, double hnydk)
{
  double totalDip = getOrbitalConvolvedDipole(scetTime, theta, phi, psi);
  totalDip=totalDip+TCMB;

  double tcmbAnt = hnydk/(exp(hnydk/TCMB)-1);
  double totalDipAnt = hnydk/(exp(hnydk/totalDip)-1);

  return totalDipAnt-tcmbAnt;
}

double 
PlanckVelocity::constrainedDipole (std::vector<double> velocity, double theta, double phi, double psi)
{
  double xv = velocity[0]/SOLSYSSPEED;
  double yv = velocity[1]/SOLSYSSPEED;
  double zv = velocity[2]/SOLSYSSPEED;

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

  return x3*M100+y3*M010+z3*M001;
}

std::vector<double> 
PlanckVelocity::getSolarConvolvedConstraint(const std::vector<double> & theta, const std::vector<double> & phi, const std::vector<double> & psi)
{
  std::vector<double> local_dipole;
  for (size_t idx=0; idx<theta.size(); ++idx)
    {
      local_dipole.push_back(constrainedDipole(vSolSys,theta[idx],phi[idx],psi[idx]));
    }
  return local_dipole;
}
