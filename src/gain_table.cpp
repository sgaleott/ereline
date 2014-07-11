#include "gain_table.hpp"
#include "misc.hpp"

#include <algorithm>
#include <iomanip>
#include <numeric>

#include <mpi.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

gainTable::gainTable()
{
}

vector<int> & 
gainTable::getPids()
{
  return pointingIds;
}

vector<double> & 
gainTable::getGains()
{
  return gain;
}

vector<double> & 
gainTable::getOffsets()
{
  return offset;
}

vector<int> & 
gainTable::getWindowVector()
{
  return windowVector;
}

int
gainTable::getPidIndex (int pid)
{
  unsigned int res = static_cast<unsigned int>(lower_bound(pointingIds.begin(), pointingIds.end(), pid) - pointingIds.begin());
  return res;
}

double 
gainTable::getGainValue (int index)
{
  return gain[index];
}

double 
gainTable::getOffsetValue (int index)
{
  return offset[index];
}

void 
gainTable::setPidValue (int a_pid)
{
  pointingIds.push_back(a_pid);
}

void 
gainTable::setGainValue (double a_gain)
{
  gain.push_back(a_gain);
}

void 
gainTable::setOffsetValue (double a_offset)
{
  offset.push_back(a_offset);
}

void 
gainTable::setWindowVector (vector<int> & a_windowVector)
{
  windowVector = a_windowVector;
}

void
gainTable::mergeResults()
{
    int rankMPI = MPI::COMM_WORLD.Get_rank();
    int sizeMPI = MPI::COMM_WORLD.Get_size();

    // Collect the number of pointing periods processed by each MPI
    // process into the "lengths" vector, and the total number of
    // periods into "overallLength".
    int length = static_cast<int>(pointingIds.size());
    int overallLength;
    MPI::COMM_WORLD.Allreduce(&length, &overallLength, 1, MPI::INT, MPI::SUM);

    std::vector<int> lengths(sizeMPI);
    MPI::COMM_WORLD.Gather(&length, 1, MPI::INT, 
			   lengths.data(), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(lengths.data(), lengths.size(), MPI::INT, 0);

    // Now use MPI's "gatherv" function to concatenate pIDs, gains, and
    // offsets into "overallPointings, overallGains, and overallOffsets.
    std::vector<int> displacement(sizeMPI);
    for(int i = 0; i < sizeMPI; ++i) {
	if(i == 0)
	    displacement[i] = 0;
	else
	    displacement[i] = displacement.at(i - 1) + lengths.at(i - 1);
    }

    std::vector<int> overallPointings(overallLength);
    std::vector<double> overallGains(overallLength);
    std::vector<double> overallOffsets(overallLength);

    MPI::COMM_WORLD.Gatherv(pointingIds.data(), pointingIds.size(), MPI::INT,
			    overallPointings.data(), lengths.data(),
			    displacement.data(), MPI::INT, 0);

    MPI::COMM_WORLD.Gatherv(gain.data(), gain.size(), MPI::DOUBLE,
			    overallGains.data(), lengths.data(),
			    displacement.data(), MPI::DOUBLE, 0);

    MPI::COMM_WORLD.Gatherv(offset.data(), offset.size(), MPI::DOUBLE,
			    overallOffsets.data(), lengths.data(), 
			    displacement.data(), MPI::DOUBLE, 0);

    if(rankMPI == 0) {
	pointingIds = overallPointings;
	gain = overallGains;
	offset = overallOffsets;
    } else {
	pointingIds.resize(overallLength);
	gain.resize(overallLength);
	offset.resize(overallLength);
    }

    // So far overallPointings, overallGains, and overallOffsets have
    // been set up in the root process only. Broadcast them to every
    // other process.
    MPI::COMM_WORLD.Bcast(pointingIds.data(), pointingIds.size(), MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(gain.data(), gain.size(), MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(offset.data(), offset.size(), MPI::DOUBLE, 0);
}

vector<double> 
gainTable::wma(int windowLen, 
	       const std::vector<double> weights, 
	       const std::vector<double> raw)
{
  // Moving window smoothing
  vector<double> windowVec;
  vector<double> paddedRaw;
  vector<double> paddedWeights;
  for (size_t idx=windowLen/2; idx>0; --idx)
    {
      paddedRaw.push_back(raw[idx]);
      paddedWeights.push_back(weights[idx]);
    }
  for (size_t idx=0; idx<raw.size(); ++idx)
    {
      paddedRaw.push_back(raw[idx]);
      paddedWeights.push_back(weights[idx]);
    }
  for (int idx=0; idx<=windowLen/2; ++idx)
    {
      paddedRaw.push_back(raw[raw.size()-idx-2]);
      paddedWeights.push_back(weights[weights.size()-idx-2]);
    }

  for (size_t extIdx=0; extIdx<raw.size(); extIdx++) 
    { 
      windowVec.push_back(0);
      double wAcc = 0;
      for (int intIdx=0; intIdx<=windowLen; intIdx++) 
	{ 
	    if ((paddedRaw[extIdx+intIdx]!=0)&&(!std::isinf(paddedRaw[extIdx+intIdx])))
	    {
	      windowVec[extIdx] += paddedRaw[extIdx+intIdx]*paddedWeights[extIdx+intIdx];
	      wAcc += paddedWeights[extIdx+intIdx];
	    }
	}
      windowVec[extIdx] /= wAcc;
    }

  return windowVec;
}

vector<int>
gainTable::createWindowsVector(int windowLen1, 
			       int windowLen2,
			       double minRangeDipole, 
			       double maxRangeDipole, 
			       const std::vector<double> & dipole)
{
  vector<int> result;
  size_t count = 0;
  size_t beginCount = 0;
  bool inTheMiddle = false;
  for (size_t idx=0; idx<dipole.size(); ++idx)
    {
      if (dipole[idx] == 0)
	{
	  if (result.size() == 0)
	    {
	      beginCount++;
	      continue;
	    }

	  if (inTheMiddle)
	    {
	      count++;
	      continue;
	    }
	  result.push_back(result[idx-1]);
	  continue;
	}

      if (dipole[idx] < minRangeDipole)
	{
	  if (result.size() == 0)
	    {
	      for (size_t intIdx=0; intIdx<beginCount; intIdx++)
		result.push_back(windowLen1);
	      
	      beginCount=0;
	    }
	      
	  if (inTheMiddle)
	    {
	      for (size_t boundIdx=0; boundIdx<count; boundIdx++)
		{
		  int step = static_cast<int>((windowLen1-result.back())/(count-boundIdx));
		  result.push_back(result.back()+step);
		}

	      inTheMiddle=false;
	      count=0;
	    }
	  result.push_back(windowLen1);
	  continue;
	}

      if (dipole[idx] > maxRangeDipole)
	{
	  if (result.size() == 0)
	    {
	      for (size_t intIdx=0; intIdx<beginCount; intIdx++)
		result.push_back(windowLen1);
	      
	      beginCount=0;
	    }

	  if (inTheMiddle)
	    {
	      for (size_t boundIdx=0; boundIdx<count; boundIdx++)
		{
		  int step = static_cast<int>((result.back()-windowLen2)/(count-boundIdx));
		  result.push_back(result.back()-step);
		}

	      inTheMiddle=false;
	      count=0;
	    }
	  result.push_back(windowLen2);
	  continue;
	}

      count++;
      inTheMiddle=true;
    }

  if (inTheMiddle)
    {
      if (result.back() == windowLen2)
	{
	  for (size_t boundIdx=0; boundIdx<count; boundIdx++)
	    {
	      int step = static_cast<int>((windowLen1-result.back())/(count-boundIdx));
	      result.push_back(result.back()+step);
	    }
	}
      else
	{
	  for (size_t boundIdx=0; boundIdx<count; boundIdx++)
	    {
	      int step = static_cast<int>((result.back()-windowLen2)/(count-boundIdx));
	      result.push_back(result.back()-step);
	    }
	}
    }
  
  return result;
}

vector<double> 
gainTable::variableWMA(const vector<int> & windowLen, 
		       const vector<double> & dipole)
{  
  // Moving window smoothing
  vector<double> windowVec;
  vector<double> paddedRaw;
  vector<double> paddedWeights;
  int maxWindow = (*std::max_element(windowLen.begin(), windowLen.end()))/2;
  for (size_t idx=maxWindow; idx>0; --idx)
    {
      paddedRaw.push_back(gain[idx]);
      paddedWeights.push_back(dipole[idx]);
    }
  for (size_t idx=0; idx<gain.size(); ++idx)
    {
      paddedRaw.push_back(gain[idx]);
      paddedWeights.push_back(dipole[idx]);
    }
  for (int idx=0; idx<=maxWindow; ++idx)
    {
      paddedRaw.push_back(gain[gain.size()-idx-2]);
      paddedWeights.push_back(dipole[dipole.size()-idx-2]);
    }

  for (size_t extIdx=0; extIdx<gain.size(); extIdx++) 
    {
      double value = 0;
      double wAcc = 0;
      int expanded = 0;
      while (((value==0)||(wAcc==0))&&(static_cast<size_t>(windowLen[extIdx]/2)+expanded<paddedRaw.size()))
	{
	  for (int intIdx=-windowLen[extIdx]/2-expanded; intIdx<=windowLen[extIdx]/2+expanded; intIdx++) 
	    { 
	      int idx = static_cast<int>(extIdx)+maxWindow+intIdx;
	      if ((paddedRaw[idx]==0)||(paddedWeights[idx]==0)||(std::isinf(paddedRaw[idx])))
		continue;
	      
	      value += paddedRaw[idx]*paddedWeights[idx];
	      wAcc += paddedWeights[idx];
	    }
	  expanded++;
	}
      windowVec.push_back(value/wAcc);
    }

  return windowVec;
}

vector<double> 
gainTable::variableWMAlocal(const vector<int> & windowLen, 
			    const vector<double> & dipole, 
			    const vector<double> & sectorGain)
{  
  // Moving window smoothing
  vector<double> windowVec;
  vector<double> paddedRaw;
  vector<double> paddedWeights;

  int maxWindow = (*max_element(windowLen.begin(), windowLen.end()))/2;
  for (size_t idx=maxWindow; idx>0; --idx)
    {
      paddedRaw.push_back(sectorGain[idx]);
      paddedWeights.push_back(dipole[idx]);
    }
  for (size_t idx=0; idx<sectorGain.size(); ++idx)
    {
      paddedRaw.push_back(sectorGain[idx]);
      paddedWeights.push_back(dipole[idx]);
    }
  for (int idx=0; idx<=maxWindow; ++idx)
    {
      paddedRaw.push_back(sectorGain[sectorGain.size()-idx-2]);
      paddedWeights.push_back(dipole[dipole.size()-idx-2]);
    }

  for (size_t extIdx=0; extIdx<sectorGain.size(); extIdx++) 
    { 
      double value = 0;
      double wAcc = 0;
      int expanded = 0;
      while (((value==0)||(wAcc==0))&&(static_cast<size_t>(windowLen[extIdx]/2)+expanded<paddedRaw.size()))
	{
	  for (int intIdx=-windowLen[extIdx]/2-expanded; intIdx<=windowLen[extIdx]/2+expanded; intIdx++) 
	    { 
	      int idx = static_cast<int>(extIdx)+maxWindow+intIdx;
	      if ((paddedRaw[idx]==0)||(paddedWeights[idx]==0)||(std::isinf(paddedRaw[idx])))
		continue;
	      
	      value += paddedRaw[idx]*paddedWeights[idx];
	      wAcc += paddedWeights[idx];
	    }
	  expanded++;
	}
      windowVec.push_back(value/wAcc);
    }

  return windowVec;
}

vector<double>
gainTable::gainSmoothing(int windowLenMinima, int windowLenMaxima, 
			 int windowLenSlowSmoothing, double percentSlowVariations,
			 double minRangeDipole, double maxRangeDipole,
			 vector<double> & dipole)
{
  vector<int> dipoleWindowVector = createWindowsVector(windowLenMinima, windowLenMaxima, minRangeDipole, maxRangeDipole, dipole);

  // Rough Smoothing of the Raw Gains
  vector<double> windowedSlowSmoothing = wma(windowLenSlowSmoothing, dipole, gain);

  // Padding
  vector<double> paddedRaw;
  for (size_t idx=windowLenSlowSmoothing/2; idx>0; --idx)
    {
      paddedRaw.push_back(windowedSlowSmoothing[idx]);
    }
  for (size_t idx=0; idx<windowedSlowSmoothing.size(); ++idx)
    {
      paddedRaw.push_back(windowedSlowSmoothing[idx]);
    }
  for (int idx=0; idx<=windowLenSlowSmoothing/2; ++idx)
    {
      paddedRaw.push_back(windowedSlowSmoothing[windowedSlowSmoothing.size()-idx-2]);
    }
  
  // Compute variance every windowLenSlowSmoothing samples and multiply by dipole
  vector<double> locVariance;
  for (size_t extIdx=0; extIdx<windowedSlowSmoothing.size(); extIdx++)
    { 
      vector<double> localGain;
      for (int intIdx=-windowLenSlowSmoothing/2; intIdx<=windowLenSlowSmoothing/2; intIdx++)
	  if ((paddedRaw[extIdx+windowLenSlowSmoothing/2+intIdx]!=0)&&(!std::isinf(paddedRaw[extIdx+windowLenSlowSmoothing/2+intIdx])))
	  {
	    localGain.push_back(paddedRaw[extIdx+windowLenSlowSmoothing/2+intIdx]);
	  }

      locVariance.push_back(computeVariance(localGain)*dipole[extIdx]);
    }
  
  // Sort Variance and find the percentile percentSlowVariations
  vector<double> sortedVariance = locVariance;
  sort(sortedVariance.begin(),sortedVariance.end());
  double percentileSlowVariations = sortedVariance[static_cast<int>(percentSlowVariations*static_cast<double>(sortedVariance.size()))];

  // Find sudden jumps using variances
  vector<size_t> startIdx;
  vector<size_t> endIdx;
  startIdx.push_back(0);
  for (size_t idx=0; idx<locVariance.size(); idx++)
    { 
      if (locVariance[idx]>percentileSlowVariations)
	{
	  vector<double> upVariance;
	  size_t rangeStartIdx = idx;
	  while ((locVariance[idx]>percentileSlowVariations)||(locVariance[idx]==0)||(idx==locVariance.size()))
	    {
	      upVariance.push_back(locVariance[idx]);
	      idx++;
	    }

	  size_t position = distance (upVariance.begin(), max_element(upVariance.begin(), upVariance.end()));
	  if ((position+rangeStartIdx-startIdx.back()) < static_cast<size_t>(dipoleWindowVector[idx]))
	    continue;

	  if ((gain.size()-(position+rangeStartIdx)) < static_cast<size_t>(dipoleWindowVector[idx]))
	    continue;

	  endIdx.push_back(position+rangeStartIdx);
	  startIdx.push_back(position+rangeStartIdx+1);
	}
    }
  endIdx.push_back(gain.size()-1);

  // Smooth gains for each interval between jumps
  vector<double> windowedG = vector<double>(dipoleWindowVector.size(), 0.);
  for (size_t idx=0; idx<startIdx.size(); idx++)
    {
      size_t len = endIdx[idx]-startIdx[idx]+1;

      vector<int> sectorWinVector = vector<int>(len, 0);
      vector<double> sectorDipole = vector<double>(len, 0.);
      vector<double> sectorRaw = vector<double>(len, 0.);

      // Fill sector windows, dipole and raw gains vectors
      size_t secIdx = 0;
      while ((secIdx < static_cast<size_t>(dipoleWindowVector[startIdx[idx]+secIdx]/2))&&(secIdx < len))
	{
	  sectorWinVector[secIdx]=dipoleWindowVector[startIdx[idx]+secIdx]/2+static_cast<int>(secIdx);
	  sectorDipole[secIdx]=dipole[startIdx[idx]+secIdx];
	  sectorRaw[secIdx]=gain[startIdx[idx]+secIdx];
	  secIdx++;
	}
      while (secIdx < len-static_cast<size_t>(dipoleWindowVector[startIdx[idx]+secIdx]/2))
	{
	  sectorWinVector[secIdx]=dipoleWindowVector[startIdx[idx]+secIdx];
	  sectorDipole[secIdx]=dipole[startIdx[idx]+secIdx];
	  sectorRaw[secIdx]=gain[startIdx[idx]+secIdx];
	  secIdx++;
	}
      while (secIdx < len)
	{
	  sectorWinVector[secIdx]=dipoleWindowVector[startIdx[idx]+secIdx]/2+static_cast<int>(len-secIdx);
	  sectorDipole[secIdx]=dipole[startIdx[idx]+secIdx];
	  sectorRaw[secIdx]=gain[startIdx[idx]+secIdx];
	  secIdx++;
	}

      // Interpolate raw gains for zeros
      vector<int> p2i;
      vector<int> pC;
      vector<double> g2i;
      for (size_t lIdx=0; lIdx<sectorRaw.size(); lIdx++)
	{ 
	  pC.push_back(pointingIds[startIdx[idx]+lIdx]);
	  if ((sectorRaw[lIdx] != 0)&&(!std::isinf(sectorRaw[lIdx])))
	    {
	      p2i.push_back(pointingIds[startIdx[idx]+lIdx]);
	      g2i.push_back(sectorRaw[lIdx]);
	    }
	}  

      double data[sectorRaw.size()];
      unsigned int locIdx=0;
      for (unsigned int lIdx=0; lIdx<pC.size(); ++lIdx)
	{
	  if (locIdx==p2i.size())
	    {
	      data[lIdx]=g2i[locIdx-1];
	      continue;
	    }
	  
	  if (pC[lIdx]==p2i[locIdx])
	    {
	      data[lIdx]=g2i[locIdx];
	      ++locIdx;
	    }
	  else
	    {
	      if (locIdx>0)
		data[lIdx]=g2i[locIdx-1]+((g2i[locIdx]-g2i[locIdx-1])/(p2i[locIdx]-p2i[locIdx-1])*(pC[lIdx]-p2i[locIdx-1]));
	      else
		data[lIdx]=g2i[locIdx];
	    }
	}

      // Compute FFT of the RAW gains
      gsl_fft_real_wavetable * real;
      gsl_fft_halfcomplex_wavetable * hc;
      gsl_fft_real_workspace * work;
      work = gsl_fft_real_workspace_alloc (sectorRaw.size());
      real = gsl_fft_real_wavetable_alloc (sectorRaw.size());
      gsl_fft_real_transform (data, 1, sectorRaw.size(), 
			      real, work);
      gsl_fft_real_wavetable_free (real);

      // Low filtering (5 percent)
      size_t tenPercent = static_cast<size_t>(0.05*static_cast<double>(sectorRaw.size()));
      for (size_t lIdx=tenPercent; lIdx<sectorRaw.size(); lIdx++)
	{
	  data[lIdx] = 0;
	}

      // Inverse FFT to find filtered gains
      hc = gsl_fft_halfcomplex_wavetable_alloc (sectorRaw.size());
      gsl_fft_halfcomplex_inverse (data, 1, sectorRaw.size(), 
				   hc, work);
      gsl_fft_halfcomplex_wavetable_free (hc);      
      gsl_fft_real_workspace_free (work);

      vector<double> sectorGain;
      for(size_t sectorIdx=0; sectorIdx<len; ++sectorIdx)
	sectorGain.push_back(data[sectorIdx]);
      
      // Smooth gains using filtered gains, variable windows and dipole variance
      vector<double> sectorWMA = variableWMAlocal(sectorWinVector, sectorDipole, sectorGain);

      // Store smoothed gains
      for(size_t sectorIdx=0; sectorIdx<len; ++sectorIdx)
	windowedG[startIdx[idx]+sectorIdx]=sectorWMA[sectorIdx];
    }

  return windowedG;
}

vector<double> 
gainTable::zeroing(int windowLen, double percent, vector<double> & dipole)
{
  vector<double> paddedRaw;
  for (size_t idx=windowLen/2; idx>0; --idx)
    {
      paddedRaw.push_back(gain[idx]);
    }
  for (size_t idx=0; idx<gain.size(); ++idx)
    {
      paddedRaw.push_back(gain[idx]);
    }
  for (int idx=0; idx<=windowLen/2; ++idx)
    {
      paddedRaw.push_back(gain[gain.size()-idx-2]);
    }

  vector<double> locVariance;
  for (size_t extIdx=0; extIdx<gain.size(); extIdx++)
    { 
      vector<double> localGain;
      for (int intIdx=-windowLen/2; intIdx<=windowLen/2; intIdx++)
	localGain.push_back(paddedRaw[extIdx+windowLen/2+intIdx]);
      
      locVariance.push_back(computeVariance(localGain));
    }

  vector<double> sortedVariance = locVariance;
  sort(sortedVariance.begin(),sortedVariance.end());
  double mulFactor = sortedVariance[static_cast<int>(percent*static_cast<double>(sortedVariance.size()))]*windowLen;

  for (size_t idx=0; idx<gain.size(); idx++)
    { 
      int locWindowLen = static_cast<int>(mulFactor/locVariance[idx]);
      if (locWindowLen > windowLen)
	windowVector.push_back(windowLen);
      else
	windowVector.push_back(locWindowLen);
    }

  // Weighted Moving Average
  vector<double> windowedG = variableWMA(windowVector, dipole);

  // Subtract windowed
  vector<double> retVec;
  for (size_t idx=0; idx<windowedG.size(); ++idx)
    retVec.push_back(gain[idx]-windowedG[idx]);

  return retVec;
}

vector<double>
gainTable::offsetSmoothing(int windowLenMinima, 
			   int windowLenMaxima, 
			   double minRangeDipole, 
			   double maxRangeDipole, 
			   vector<double> & dipole)
{
  vector<int> dipoleWindowVector = createWindowsVector(windowLenMinima, windowLenMaxima, minRangeDipole, maxRangeDipole, dipole);
  vector<double> retVec = variableWMAlocal(dipoleWindowVector, dipole, offset);
  return retVec;
}

void
gainTable::selectRadiometerGains(int armNumber, 
				 size_t detectorIdsSize, 
				 const vector<int> & nIdsRange)
{
    // Select values
    vector<int> tmpPointingIds;
    vector<double> tmpGain;
    vector<double> tmpOffset;

    size_t startPoint = 0;
    for (auto chunkSize : nIdsRange)
    {
	size_t chunkOffset = armNumber * chunkSize;
	size_t chunkStart = startPoint + chunkOffset;

	tmpPointingIds.insert(tmpPointingIds.end(),
			      pointingIds.begin() + chunkStart,
			      pointingIds.begin() + chunkStart + chunkSize);
	tmpGain.insert(tmpGain.end(),
		       gain.begin() + chunkStart,
		       gain.begin() + chunkStart + chunkSize);
	tmpOffset.insert(tmpOffset.end(),
			 offset.begin() + chunkStart,
			 offset.begin() + chunkStart + chunkSize);

	startPoint += detectorIdsSize * chunkSize;
    }

  // Swap vectors
  pointingIds.swap(tmpPointingIds);
  gain.swap(tmpGain);
  offset.swap(tmpOffset);
}

void
gainTable::selectDiodeGains(int detectorIdIdx, size_t detectorIdsSize, 
			    const std::vector<int> nIdsRange, 
			    std::vector<double> & outGain,
			    std::vector<double> & outOffset)
{
  // Select values
  vector<int> tmpPointingIds;
  vector<double> tmpGain;
  vector<double> tmpOffset;
  for (size_t idx=0; idx<nIdsRange.size(); ++idx)
    {
      int offsetIdx=detectorIdIdx*nIdsRange[idx];
      int startPoint=0;
      for (size_t intIdx=0; intIdx<idx; ++intIdx)
	startPoint += nIdsRange[intIdx]*static_cast<int>(detectorIdsSize);

      for (int intIdx=0; intIdx<nIdsRange[idx]; ++intIdx)
	{
	  tmpPointingIds.push_back(pointingIds[startPoint+intIdx+offsetIdx]);
	  tmpGain.push_back(gain[startPoint+intIdx+offsetIdx]);
	  tmpOffset.push_back(offset[startPoint+intIdx+offsetIdx]);	  
	}
    }

  // Swap vectors
  pointingIds.swap(tmpPointingIds);
  gain.swap(tmpGain);
  offset.swap(tmpOffset);

  int otherIdIdx = detectorIdIdx;
  if (detectorIdIdx % 2 == 0)
    otherIdIdx += 1;
  else
    otherIdIdx -= 1;

  for (size_t idx=0; idx<nIdsRange.size(); ++idx)
    {
      int offsetIdx=otherIdIdx*nIdsRange[idx];
      int startPoint=0;
      for (size_t intIdx=0; intIdx<idx; ++intIdx)
	startPoint += nIdsRange[intIdx]*static_cast<int>(detectorIdsSize);

      for (int intIdx=0; intIdx<nIdsRange[idx]; ++intIdx)
	{
	  outGain.push_back(tmpGain[startPoint+intIdx+offsetIdx]);
	  outOffset.push_back(tmpOffset[startPoint+intIdx+offsetIdx]);	  
	}
    }
}
