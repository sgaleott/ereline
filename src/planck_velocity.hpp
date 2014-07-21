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

  double dipole(const std::vector<double> & velocity, 
		double theta, 
		double phi) const;
  double convolvedDipole(const std::vector<double> & velocity, 
			 double theta, 
			 double phi, 
			 double psi) const;

  PlanckVelocity (const std::string & file_name);
  PlanckVelocity ();

  std::vector<double> getVelocity (double scetTime) const;
  std::vector<double> getAbsoluteVelocity(double scetTime) const;

  double getConvolvedDipole(double scetTime,
			    double theta,
			    double phi, double psi) const;
  std::vector<double> getConvolvedDipole(const std::vector<double> & scetTime,
					 const std::vector<double> & theta,
					 const std::vector<double> & phi,
					 const std::vector<double> & psi) const;
};

#endif
