#ifndef _DIPOLE_FIT_
#define _DIPOLE_FIT_

#include <vector>

#include "planck_velocity.hpp"

using namespace std;

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

  vector<int> pixIndex;
  vector<int> pixSumHits;
  vector<double> pixSumData;
  vector<float> pixSumDipole;
  vector<float> inputMap;

 public:
  
  dipoleFit(const int a_qualityFlag, const int a_nSide, const int a_pointingID);

  ~dipoleFit(){};

  bool binData(const vector<double> & data, const vector<int>& flag,
	       const vector<double> & theta, const vector<double> & phi, const vector<double> & dipole,
	       const vector<size_t> & pidRange, const vector<double> & sidelobes);
  bool fitData(const vector<float> & maskMap);
  bool fit(const vector<double> & data, const vector<int> & flag,
	   const vector<double> & theta, const vector<double> & phi, const vector<double> & dipole,
	   const vector<size_t> & pidRange, const vector<float> & maskMap, const vector<double> & sidelobes);

  void setGainV(double a_gainv);
  void setOffset(double a_offset);
  void setPixSumDipole(vector<float> inpArr);

  double getGainV();
  double getGain();
  double getOffset();

  int getPointingID();
  int getNSide();

  vector<int> & getPixIndex();
  vector<double> & getPixSumData();
  vector<int> & getPixSumHits();
  vector<float> & getPixSumDipole();

  double getDipoleVariance();
  double getMaxDipole();
  double getMinDipole();

  void unload();
};


#endif

  


