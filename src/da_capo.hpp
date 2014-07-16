#ifndef _DACAPO_HPP_
#define _DACAPO_HPP_

#include <vector>
#include <algorithm>

#include "dipole_fit.hpp"

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
  int nSide;

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

  void initializeConstraint(bool constraint);
  void initializeConstraint(std::vector<double> & constraint);
  void initializeLocmap(const std::vector<dipoleFit> & binnedData);
  void initializeFullmap();
  void applyMask(const std::vector<dipoleFit> & binnedData, 
		 const std::vector<float> & mask);

 public:
  daCapo(std::vector<dipoleFit> & binnedData, std::vector<float> & mask, bool constraint);
  daCapo (std::vector<dipoleFit> & binnedData, std::vector<float> & mask, std::vector<double> & constraint);

  void constructCCmatrix(const std::vector<dipoleFit> & binnedData);
  void updateDipolenorm();
  void buildPreconditioner(const std::vector<dipoleFit> & binnedData);
  void toiToLocmap(const std::vector<dipoleFit> & binnedData);
  void locToFullmap();
  void ccMultiply();
  void applyConstraint();
  void fullToLocmap();
  void applyCC();
  void subtractMapFromTod(const std::vector<dipoleFit> & binnedData, basevec &p);
  void baseToLocmap(const std::vector<dipoleFit> & binnedData, const basevec &p);
  void subtractMapFromBase(const std::vector<dipoleFit> & binnedData, basevec &p);
  void applyPreconditioner(const basevec &r, basevec &z);
  void updateSignal(std::vector<dipoleFit> & binnedData);

  double iterativeCalibration(std::vector<dipoleFit> & binnedData, 
			      bool firstLoop);
};

class Configuration;

void run_da_capo(const Configuration & program_conf,
		 const Configuration & storage_conf);

#endif
