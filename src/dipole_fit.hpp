#ifndef _DIPOLE_FIT_
#define _DIPOLE_FIT_

#include <vector>

#include "planck_velocity.hpp"

/* class calib iter*/
class dipoleFit
{
 private:
  int qualityFlag;
  int nSide;
  int pointingID;

  double gainv;
  double offset;

  double maxDipole;
  double minDipole;

  std::vector<int> pixIndex;
  std::vector<int> pixSumHits;
  std::vector<double> pixSumData;
  std::vector<float> pixSumDipole;
  std::vector<float> inputMap;

 public:
  
  dipoleFit(const int a_qualityFlag, const int a_nSide, const int a_pointingID);

  ~dipoleFit(){};

  bool binData(const std::vector<double> & data, const std::vector<int>& flag,
	       const std::vector<double> & theta, const std::vector<double> & phi, const std::vector<double> & dipole,
	       const std::vector<size_t> & pidRange, const std::vector<double> & sidelobes);
  bool fitData(const std::vector<float> & maskMap);
  bool fit(const std::vector<double> & data, const std::vector<int> & flag,
	   const std::vector<double> & theta, const std::vector<double> & phi, const std::vector<double> & dipole,
	   const std::vector<size_t> & pidRange, const std::vector<float> & maskMap, const std::vector<double> & sidelobes);

  void setGainV(double a_gainv);
  void setOffset(double a_offset);
  void setPixSumDipole(const std::vector<float> inpArr);

  double getGainV() const;
  double getGain() const;
  double getOffset() const;

  int getPointingID() const;
  int getNSide() const;

  const std::vector<int> & getPixIndex() const;
  const std::vector<double> & getPixSumData() const;
  const std::vector<int> & getPixSumHits() const;
  const std::vector<float> & getPixSumDipole() const;

  double getDipoleVariance() const;
  double getMaxDipole() const;
  double getMinDipole() const;

  void unload();
};

class Configuration;

#endif

  


