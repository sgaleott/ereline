/*
 *  Class for holding information velocity of the Planck satellite
 *
 *  Copyright (C) 2011 LFI-DPC
 *  \author Samuele Galeotta
 */

#ifndef PLANCK_VELOCITY_H
#define PLANCK_VELOCITY_H

#include <string>
#include <vector>

struct PlanckVelocity
{
  double M100;
  double M010;
  double M001;
  double M200;
  double M110;
  double M101;
  double M020;
  double M011;
  double M002;

  std::vector<double> scet;
  std::vector<double> xvel;
  std::vector<double> yvel;
  std::vector<double> zvel;

  std::vector<double> vSolSys;

  double dipole(std::vector<double> velocity, double theta, double phi);
  double convolvedDipole(std::vector<double> velocity, double theta, double phi, double psi);

  PlanckVelocity (const std::string & file_name);
  PlanckVelocity ();
  
  std::vector<double> getVelocity (double scetTime);  
  std::vector<double> getAbsoluteVelocity(double scetTime);

  double getTotalDipole(double scetTime, double theta, double phi);
  std::vector<double> getTotalDipole(const std::vector<double> & scetTime, const std::vector<double> & theta, 
				const std::vector<double> & phi);
  double getSolarDipole(double theta, double phi);
  std::vector<double> getSolarDipole(const std::vector<double> & theta, const std::vector<double> & phi);
  double getOrbitalDipole(double scetTime, double theta, double phi);

  double getTotalDipoleTAnt(double scetTime, double theta, double phi, double hnydk);
  double getSolarDipoleTAnt(double theta, double phi, double hnydk);
  double getOrbitalDipoleTAnt(double scetTime, double theta, double phi, double hnydk);

  double getConvolvedDipole(double scetTime, double theta, double phi, double psi);
  std::vector<double> getConvolvedDipole(const std::vector<double> & scetTime, const std::vector<double> & theta, 
				    const std::vector<double> & phi, const std::vector<double> & psi);
  double getSolarConvolvedDipole(double theta, double phi, double psi);
  std::vector<double> getSolarConvolvedDipole(const std::vector<double> & theta, const std::vector<double> & phi, const std::vector<double> & psi);
  double getOrbitalConvolvedDipole(double scetTime, double theta, double phi, double psi);
  std::vector<double> getOrbitalConvolvedDipole(const std::vector<double> & scetTime, const std::vector<double> & theta, 
					   const std::vector<double> & phi, const std::vector<double> & psi);

  double getConvolvedDipoleTAnt(double scetTime, double theta, double phi, double psi, double hnydk);
  double getSolarConvolvedDipoleTAnt(double theta, double phi, double psi, double hnydk);
  double getOrbitalConvolvedDipoleTAnt(double scetTime, double theta, double phi, double psi, double hnydk);

  double constrainedDipole (std::vector<double> velocity, double theta, double phi, double psi);
  std::vector<double> getSolarConvolvedConstraint(const std::vector<double> & theta, const std::vector<double> & phi, const std::vector<double> & psi);
};

#endif
