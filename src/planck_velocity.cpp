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
PlanckVelocity::getVelocity (double scetTime) const
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
PlanckVelocity::getAbsoluteVelocity(double scetTime) const
{
  std::vector<double> vSat = getVelocity(scetTime);

  // Add the solar system velocity
  vSat[0] += vSolSys[0];
  vSat[1] += vSolSys[1];
  vSat[2] += vSolSys[2];
  
  return vSat;
}

double 
PlanckVelocity::dipole(const std::vector<double> & velocity, 
		       double theta, 
		       double phi) const
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
PlanckVelocity::convolvedDipole (const std::vector<double> & velocity, 
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

  double totalDipole = scalarProduct+relativisticCorrection;
  return totalDipole*TCMB;
}

double 
PlanckVelocity::getConvolvedDipole(double scetTime, double theta, double phi, double psi) const
{  
  // Compute velocity versor
  const std::vector<double> vSatAbsolute = getAbsoluteVelocity(scetTime);
  return convolvedDipole (vSatAbsolute,theta,phi,psi);
}

std::vector<double> 
PlanckVelocity::getConvolvedDipole(const std::vector<double> & scetTime, 
				   const std::vector<double> & theta, 
				   const std::vector<double> & phi, 
				   const std::vector<double> & psi) const
{
  std::vector<double> local_dipole;
  for (size_t idx=0; idx<scetTime.size(); ++idx)
    {
      std::vector<double> vSatAbsolute = getAbsoluteVelocity(scetTime[idx]);
      local_dipole.push_back(convolvedDipole(vSatAbsolute,theta[idx],phi[idx],psi[idx]));
    }
  return local_dipole;
}
