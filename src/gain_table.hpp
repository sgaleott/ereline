/*
 * The class implements the reading and storing
 * of the toi.attitude.HighFrequency objects.
 *
 * Copyright (C) 2011 LFI-DPC
 *
 * Author: Samuele Galeotta
 *
 */

#ifndef GAIN_TABLE_HPP
#define GAIN_TABLE_HPP

#include <string>
#include <vector>

using namespace std;

class gainTable
{
private:
  vector<int> pointingIds;
  vector<double> gain;
  vector<double> offset;
  vector<int> windowVector;

  vector<int> createWindowsVector(int windowLen1, int windowLen2,
				  double minRangeDipole, double maxRangeDipole, vector<double> & dipole);
  vector<double> variableWMAlocal(vector<int> & windowLen, vector<double> & dipole, vector<double> & sectorGain);
  vector<double> variableWMA(vector<int> & windowLen, vector<double> & dipole);
public:
  gainTable();

  vector<int> & getPids();
  vector<double> & getGains();
  vector<double> & getOffsets();
  vector<int> & getWindowVector();
  double getGainValue (int index);
  double getOffsetValue (int index);
  int getPidIndex (int pid);

  void setPidValue (int a_pid);
  void setGainValue (double a_gain);
  void setOffsetValue (double a_offset);
  void setWindowVector (vector<int> & a_windowVector);

  void mergeResults();

  vector<double> wma(int windowLen, vector<double> weights, vector<double> raw);
  vector<double> offsetSmoothing(int windowLenMinima, int windowLenMaxima, 
				 double minRangeDipole, double maxRangeDipole, vector<double> & dipole);
  vector<double> gainSmoothing(int windowLenMinima, int windowLenMaxima, 
			       int windowLenSlowSmoothing, double percentSlowVariations,
			       double minRangeDipole, double maxRangeDipole,
			       vector<double> & dipole);
  vector<double> zeroing(int windowLen, double percent, vector<double> & dipole);

  void selectRadiometerGains(int detectorIdIdx, size_t detectorIdsSize, 
			     vector<int> nIdsRange);
  void selectDiodeGains(int detectorIdIdx, size_t detectorIdsSize, 
			vector<int> nIdsRange, vector<double> & outGain,
			vector<double> & outOffset);
};

#endif
