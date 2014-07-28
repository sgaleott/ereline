#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_fit.h>
#include <mpi.h>
#include <algorithm>
#include <cmath>

#include "deltavv.hpp"
#include "misc.hpp"

////////////////////////////////////////////////////////////////////////////////

Delta_vv_t::Delta_vv_t()
{
}

////////////////////////////////////////////////////////////////////////////////

void
Delta_vv_t::hybridFit(const std::vector<double> & subData,
                      const std::vector<double> & subSensor)
{
    // Remove wrong data
    const double stdDevData = computeVariance(subData);
    const double meanSubData = computeMean(subData);
    const double stdDevSensor = computeVariance(subSensor);
    const double meanSubSensor = computeMean(subSensor);

    // Mask missing pid
    std::vector<double> locData;
    std::vector<double> locGains;
    std::vector<double> locDipole;
    std::vector<double> locSensor;
    for (unsigned int idx=0; idx<pid.size(); ++idx) {
        if ((gain[idx]==0)||(std::isinf(gain[idx]))
            ||(std::isnan(subSensor[idx]))||(subSensor[idx]==0)
            ||(subSensor[idx]>(meanSubSensor+5*stdDevSensor))||(subSensor[idx]<(meanSubSensor-5*stdDevSensor))
            ||(std::isnan(subData[idx]))||(subData[idx]==0)
            ||(subData[idx]>(meanSubData+5*stdDevData))||(subData[idx]<(meanSubData-5*stdDevData)))
            continue;

        locData.push_back(subData[idx]);
        locGains.push_back(gain[idx]);
        locDipole.push_back(dipole[idx]);
        locSensor.push_back(subSensor[idx]);
    }
    double meanSensor = computeMean(locSensor);

    // Multi-parameters fitting
    // Prepare data
    gsl_matrix *data, *cov;
    gsl_vector *loadVolt, *weights, *coeff;

    data = gsl_matrix_alloc (locData.size(), 2);
    loadVolt = gsl_vector_alloc (locData.size());
    weights = gsl_vector_alloc (locData.size());

    coeff = gsl_vector_alloc (2);
    cov = gsl_matrix_alloc (2, 2);

    for(size_t idx=0; idx<locData.size(); idx++) {
        gsl_matrix_set (data, idx, 0, 1.0);
        gsl_matrix_set (data, idx, 1, locSensor[idx]-meanSensor);

        gsl_vector_set (loadVolt, idx, locData[idx]*locGains[idx]);
        gsl_vector_set (weights, idx, abs(locDipole[idx]));
    }

    // Fit
    double chisq;
    gsl_multifit_linear_workspace * work =
        gsl_multifit_linear_alloc (locData.size(), 2);
    gsl_multifit_wlinear (data, weights, loadVolt, coeff, cov,
                          &chisq, work);
    gsl_multifit_linear_free (work);

    // Save results
    std::vector<double> returnGains;
    std::vector<int> returnPids;
    for (unsigned int idx=0; idx<pid.size(); ++idx) {
        if ((gain[idx]==0)||(std::isinf(gain[idx]))
            ||(std::isnan(subSensor[idx]))||(subSensor[idx]==0)
            ||(subSensor[idx]>(meanSubSensor+5*stdDevSensor))||(subSensor[idx]<(meanSubSensor-5*stdDevSensor))
            ||(std::isnan(subData[idx]))||(subData[idx]==0)
            ||(subData[idx]>(meanSubData+5*stdDevData))||(subData[idx]<(meanSubData-5*stdDevData)))
            continue;

        returnGains.push_back((gsl_vector_get(coeff,0)+gsl_vector_get(coeff,1)*(subSensor[idx]-meanSensor))/subData[idx]);
        returnPids.push_back(pid[idx]);
    }

    // Memory free
    gsl_matrix_free (data);
    gsl_vector_free (loadVolt);
    gsl_vector_free (weights);
    gsl_vector_free (coeff);
    gsl_matrix_free (cov);

    gain.swap(returnGains);
    pid.swap(returnPids);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double>
Delta_vv_t::eval(const std::vector<int> & pointingIDs)
{
    std::vector<double> interpGain;
    unsigned int locIdx=0;
    for (unsigned int idx=0; idx<pointingIDs.size(); ++idx) {
        if (locIdx==pid.size())
            break;

        if (pointingIDs[idx]==pid[locIdx]) {
            interpGain.push_back(gain[locIdx]);
            ++locIdx;
        } else {
            if (locIdx>0)
                interpGain.push_back(gain[locIdx-1]+((gain[locIdx]-gain[locIdx-1])/(pid[locIdx]-pid[locIdx-1])*(pointingIDs[idx]-pid[locIdx-1])));
            else
                interpGain.push_back(gain[locIdx]+((gain[locIdx+1]-gain[locIdx])/(pid[locIdx+1]-pid[locIdx])*(pointingIDs[idx]-pid[locIdx])));
        }
    }

  return interpGain;
}

////////////////////////////////////////////////////////////////////////////////

void
Delta_vv_t::setPidValue (int a_pid)
{
    pid.push_back(a_pid);
}

////////////////////////////////////////////////////////////////////////////////

void
Delta_vv_t::setGainValue (double a_gain)
{
    gain.push_back(a_gain);
}

////////////////////////////////////////////////////////////////////////////////

void
Delta_vv_t::setDipoleValue (double a_dipole)
{
    dipole.push_back(a_dipole);
}

////////////////////////////////////////////////////////////////////////////////

void
Delta_vv_t::mergeResults()
{
    // Get MPI infos
    int rankMPI = MPI::COMM_WORLD.Get_rank();
    int sizeMPI = MPI::COMM_WORLD.Get_size();

    // Set gain and pid arrays
    int length = static_cast<int>(pid.size());
    int * lengths = new int[sizeMPI];
    std::fill (lengths,lengths+sizeMPI,0);
    for (int jj=0; jj<sizeMPI; ++jj) {
        if (jj==rankMPI)
            lengths[jj] = length;
        else
            lengths[jj] = 0.0;
    }

    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, lengths, sizeMPI, MPI::INT, MPI::SUM);
    MPI::COMM_WORLD.Barrier();

    MPI::COMM_WORLD.Allreduce(&length, &length, 1, MPI::INT, MPI::SUM);
    MPI::COMM_WORLD.Barrier();

  // Retrieve data from other processes
    int startPid=0;
    for (int jj=0; jj<rankMPI; ++jj)
        startPid+=lengths[jj];

    int * pids = new int[length];
    std::fill (pids, pids+length, 0);
    double * gainArr = new double[length];
    std::fill (gainArr, gainArr+length, 0.0);
    double * dipoleArr = new double[length];
    std::fill (dipoleArr, dipoleArr+length, 0);

    for (unsigned int jj = 0; jj < pid.size(); ++jj) {
      pids[startPid+jj] = pid[jj];
      gainArr[startPid+jj] = gain[jj];
      dipoleArr[startPid+jj] = dipole[jj];
    }

    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, pids, length, MPI::INT, MPI::SUM);
    MPI::COMM_WORLD.Barrier();

    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, gainArr, length, MPI::DOUBLE, MPI::SUM);
    MPI::COMM_WORLD.Barrier();

    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, dipoleArr, length, MPI::DOUBLE, MPI::SUM);
    MPI::COMM_WORLD.Barrier();

    std::vector<int>().swap(pid);
    std::vector<double>().swap(gain);
    std::vector<double>().swap(dipole);

    pid.assign(pids, pids+length);
    gain.assign(gainArr, gainArr+length);
    dipole.assign(dipoleArr, dipoleArr+length);

    delete [] lengths;
    delete [] gainArr;
    delete [] pids;
    delete [] dipoleArr;
}

////////////////////////////////////////////////////////////////////////////////

void
Delta_vv_t::selectRadiometerGains(int detectorIdIdx,
                                  size_t detectorIdsSize,
                                  const std::vector<int> & nIdsRange)
{
    // Select values
    std::vector<int> tmpPid;
    std::vector<double> tmpGain;
    std::vector<double> tmpDipole;

    for (size_t idx = 0; idx < nIdsRange.size(); ++idx) {
        int offsetIdx = detectorIdIdx * nIdsRange[idx];
        int startPoint = 0;
        for (size_t intIdx = 0; intIdx < idx; ++intIdx)
            startPoint += nIdsRange[intIdx] * detectorIdsSize;

        for (int intIdx=0; intIdx<nIdsRange[idx]; ++intIdx) {
            tmpPid.push_back(pid[startPoint+intIdx+offsetIdx]);
            tmpGain.push_back(gain[startPoint+intIdx+offsetIdx]);
            tmpDipole.push_back(dipole[startPoint+intIdx+offsetIdx]);
        }
    }

    // Swap vectors
    pid.swap(tmpPid);
    gain.swap(tmpGain);
    dipole.swap(tmpDipole);
}
