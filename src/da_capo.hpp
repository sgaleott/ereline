#ifndef _DACAPO_HPP_
#define _DACAPO_HPP_

#include <vector>
#include <algorithm>

#include "dipole_fit.hpp"

using namespace std;

//Class for holding a set of 4 baseline vectors,
//two of which represent the usual destriping baselines, and the remaining two the gain.
class basevec
{
public:
  vector<double> base;
  vector<double> gain;
  int nbase;
  
public:
  
  basevec(int npp)
  {
    nbase=npp;

    base = vector<double> (npp, 0.);
    gain = vector<double> (npp, 0.);
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

  vector<double> constraintMap;
  vector<double> localMap;
  vector<double> fullMap;
  vector<double> dipolenorm;
  vector<double> ccFull;
  vector< vector<double> > preconditioner;

  vector<int> pixelIndexLocal;
  vector< vector<int> > pixelIndexLocalMap;
  vector<int> pixelIndexFull;
  vector<int> pixelIndexFullMap;
  vector<int> sendcnt;

  void initializeConstraint(bool constraint);
  void initializeConstraint(vector<double> & constraint);
  void initializeLocmap(vector<dipoleFit> & binnedData);
  void initializeFullmap();
  void applyMask(vector<dipoleFit> & binnedData, vector<float> & mask);

 public:
  daCapo(vector<dipoleFit> & binnedData, vector<float> & mask, bool constraint);
  daCapo (vector<dipoleFit> & binnedData, vector<float> & mask, vector<double> & constraint);

  void constructCCmatrix(vector<dipoleFit> & binnedData);
  void updateDipolenorm();
  void buildPreconditioner(vector<dipoleFit> & binnedData);
  void toiToLocmap(vector<dipoleFit> & binnedData);
  void locToFullmap();
  void ccMultiply();
  void applyConstraint();
  void fullToLocmap();
  void applyCC();
  void subtractMapFromTod(vector<dipoleFit> & binnedData, basevec &p);
  void baseToLocmap(vector<dipoleFit> & binnedData, const basevec &p);
  void subtractMapFromBase(vector<dipoleFit> & binnedData, basevec &p);
  void applyPreconditioner(const basevec &r, basevec &z);
  void updateSignal(vector<dipoleFit> & binnedData);

  double iterativeCalibration(vector<dipoleFit> & binnedData, bool firstLoop);
};

#endif
