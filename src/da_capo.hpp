#ifndef _DACAPO_HPP_
#define _DACAPO_HPP_

#include <vector>
#include <algorithm>

#include "dipole_fit.hpp"

struct Dipole_parameters_t;

//Class for holding a set of 4 baseline vectors,
//two of which represent the usual destriping baselines, and the remaining two the gain.
class basevec
{
public:
  std::vector<double> base;
  std::vector<double> gain;
  int nbase;

public:

  basevec(int npp)
  {
    nbase=npp;

    base = std::vector<double> (npp, 0.);
    gain = std::vector<double> (npp, 0.);
  }

  void SetValues(double g, double b, int index)
  {
    base[index]=b;
    gain[index]=g;
  }

  //Initialize
  void SetZero()
  {
    base.assign (base.size(), 0.);
    gain.assign (gain.size(), 0.);
  }

  void Fill(double g, double b)
  {
    base.assign (base.size(), b);
    gain.assign (gain.size(), g);
  }

  //Scale by constant
  void Scale(const double scal)
  {
    for (int i=0; i<nbase; i++)
      {
        base[i] *= scal;
        gain[i] *= scal;
      }
  }

  //Add another vec
  void Add(const basevec other)
  {
    for (int i=0; i<nbase; i++)
      {
        base[i] += other.base[i];
        gain[i] += other.gain[i];
      }
  }

  //Dot product between two vecs
  double Dotprod(const basevec other) const
  {
    double dnorm=0;
    for (int i=0; i<nbase; i++)
      {
        dnorm += base[i]*other.base[i];
        dnorm += gain[i]*other.gain[i];
      }
    return dnorm;
  }
};

class daCapo
{
 protected:

  int nPixelMap;
  int nside;

  int sizeMPI;
  int rankMPI;

  double rzinit;

  std::vector<double> constraintMap;
  std::vector<double> localMap;
  std::vector<double> fullMap;
  std::vector<double> dipolenorm;
  std::vector<double> ccFull;
  std::vector< std::vector<double> > preconditioner;

  std::vector<int> pixelIndexLocal;
  std::vector< std::vector<int> > pixelIndexLocalMap;
  std::vector<int> pixelIndexFull;
  std::vector<int> pixelIndexFullMap;
  std::vector<int> sendcnt;

  void initializeConstraint(bool constraint,
                            const Dipole_parameters_t & solar_dipole);
  void initializeConstraint(const std::vector<double> & constraint);
  void initializeLocmap(const std::vector<Dipole_fit_t> & binnedData);
  void initializeFullmap();
  void applyMask(const std::vector<Dipole_fit_t> & binnedData,
                 const std::vector<float> & mask);

 public:
  daCapo(const std::vector<Dipole_fit_t> & locallyBinnedData,
         const std::vector<float> & mask,
         bool constraint,
         const Dipole_parameters_t & dipole_params);
  daCapo (const std::vector<Dipole_fit_t> & locallyBinnedData,
          const std::vector<float> & mask,
          const std::vector<double> & constraint);

  void constructCCmatrix(const std::vector<Dipole_fit_t> & binnedData);
  void updateDipolenorm();
  void buildPreconditioner(const std::vector<Dipole_fit_t> & binnedData);
  void toiToLocmap(const std::vector<Dipole_fit_t> & binnedData);
  void locToFullmap();
  void ccMultiply();
  void applyConstraint();
  void fullToLocmap();
  void applyCC();
  void subtractMapFromTod(const std::vector<Dipole_fit_t> & binnedData, basevec &p);
  void baseToLocmap(const std::vector<Dipole_fit_t> & binnedData, const basevec &p);
  void subtractMapFromBase(const std::vector<Dipole_fit_t> & binnedData, basevec &p);
  void applyPreconditioner(const basevec &r, basevec &z);
  void updateSignal(std::vector<Dipole_fit_t> & binnedData);

  double iterativeCalibration(std::vector<Dipole_fit_t> & binnedData,
                              bool firstLoop);
};

class Configuration;
struct Dipole_fit_results_t;
struct Da_capo_results_t;
struct Lfi_radiometer_t;

void run_da_capo(const Configuration & program_conf,
                 const Configuration & storage_conf,
                 const Lfi_radiometer_t & user_rad,
                 const std::vector<Pointing_t> & list_of_pointings,
                 Dipole_fit_results_t & dipole_fit_results,
                 Da_capo_results_t & da_capo_results);

#endif
