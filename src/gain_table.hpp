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

#include <cstdlib>
#include <vector>

class gainTable
{
private:
  std::vector<int> pointingIds;
  std::vector<double> gain;
  std::vector<double> offset;
  std::vector<int> windowVector;

  std::vector<int> createWindowsVector(int windowLen1, 
				       int windowLen2,
				       double minRangeDipole, 
				       double maxRangeDipole, 
				       const std::vector<double> & dipole);
  std::vector<double> variableWMAlocal(const std::vector<int> & windowLen, 
				       const std::vector<double> & dipole, 
				       const std::vector<double> & sectorGain);
  std::vector<double> variableWMA(const std::vector<int> & windowLen, 
				  const std::vector<double> & dipole);
public:
  gainTable();

  std::vector<int> & getPids();
  std::vector<double> & getGains();
  std::vector<double> & getOffsets();
  std::vector<int> & getWindowVector();
  double getGainValue (int index);
  double getOffsetValue (int index);
  int getPidIndex (int pid);

  void setPidValue (int a_pid);
  void setGainValue (double a_gain);
  void setOffsetValue (double a_offset);
  void setWindowVector (std::vector<int> & a_windowVector);

  void mergeResults();

  std::vector<double> wma(int windowLen, 
			  const std::vector<double> weights, 
			  const std::vector<double> raw);
  std::vector<double> offsetSmoothing(int windowLenMinima, int windowLenMaxima, 
				 double minRangeDipole, double maxRangeDipole, std::vector<double> & dipole);
  std::vector<double> gainSmoothing(int windowLenMinima, int windowLenMaxima, 
			       int windowLenSlowSmoothing, double percentSlowVariations,
			       double minRangeDipole, double maxRangeDipole,
			       std::vector<double> & dipole);
  std::vector<double> zeroing(int windowLen, double percent, std::vector<double> & dipole);

  void selectRadiometerGains(int detectorIdIdx, 
			     size_t detectorIdsSize, 
			     const std::vector<int> & nIdsRange);
  void selectDiodeGains(int detectorIdIdx, size_t detectorIdsSize, 
			const std::vector<int> nIdsRange, 
			std::vector<double> & outGain,
			std::vector<double> & outOffset);
};

#endif
