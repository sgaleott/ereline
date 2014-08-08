#include <mpi.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_fit.h>
#include <mpi.h>
#include <algorithm>
#include <cmath>

#include "deltavv.hpp"
#include "logging.hpp"
#include "misc.hpp"
#include "mpi_processes.hpp"

////////////////////////////////////////////////////////////////////////////////

Delta_vv_t::Delta_vv_t()
{
}

////////////////////////////////////////////////////////////////////////////////

void
Delta_vv_t::hybridFit(const std::vector<double> & subData,
                      const std::vector<double> & subSensor)
{
    Logger * log = Logger::get_instance();

    log->debug(boost::format("Entering function %1%")
              % __PRETTY_FUNCTION__);
    log->increase_indent();

    // These constants are used to flag "bad" data
    const double stdDevData = computeVariance(subData);
    const double meanSubData = computeMean(subData);
    const double stdDevSensor = computeVariance(subSensor);
    const double meanSubSensor = computeMean(subSensor);

    log->debug(boost::format("stdDevData = %1%, meanSubData = %2%, "
                             "stdDevSensor = %3%, meanSubSensor = %4%")
               % stdDevData
               % meanSubData
               % stdDevSensor
               % meanSubSensor);

    // Mask missing pid
    log->debug("Going to mask those pIDs for which gains are too «strange»");
    log->increase_indent();

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
        {
            log->debug(boost::format("Skipping pid %1%")
                       % pid[idx]);
            continue;
        }

        locData.push_back(subData[idx]);
        locGains.push_back(gain[idx]);
        locDipole.push_back(dipole[idx]);
        locSensor.push_back(subSensor[idx]);
    }
    log->decrease_indent();
    log->debug(boost::format("Done, %1% elements kept for the fit")
               % locData.size());

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

    for(size_t idx = 0; idx < locData.size(); idx++) {
        gsl_matrix_set (data, idx, 0, 1.0);
        gsl_matrix_set (data, idx, 1, locSensor[idx] - meanSensor);

        gsl_vector_set (loadVolt, idx, locData[idx] * locGains[idx]);
        gsl_vector_set (weights, idx, std::fabs(locDipole[idx]));
    }

    // Fit
    double chisq;
    gsl_multifit_linear_workspace * work =
        gsl_multifit_linear_alloc (locData.size(), 2);
    gsl_multifit_wlinear (data, weights, loadVolt, coeff, cov,
                          &chisq, work);
    gsl_multifit_linear_free (work);

    log->debug(boost::format("Interpolation coefficients: c_0 = %1%, "
                             "c_1 = %2%, \u03c7\u00b2 = %3%")
               % gsl_vector_get(coeff, 0)
               % gsl_vector_get(coeff, 1)
               % chisq);

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

    log->decrease_indent();
    log->debug(boost::format("Quitting function %1%")
              % __PRETTY_FUNCTION__);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double>
Delta_vv_t::eval(const std::vector<int> & pointingIDs) const
{
    std::vector<double> interpGain(pointingIDs.size());
    size_t locIdx = 0;
    for (size_t idx = 0; idx < pointingIDs.size(); ++idx) {
        if (pointingIDs[idx]==pid[locIdx]) {
            interpGain[idx] = gain[locIdx];
            ++locIdx;
        } else {
            if (locIdx>0)
                interpGain[idx] = gain[locIdx-1]+((gain[locIdx]-gain[locIdx-1])/(pid[locIdx]-pid[locIdx-1])*(pointingIDs[idx]-pid[locIdx-1]));
            else
                interpGain[idx] = gain[locIdx]+((gain[locIdx+1]-gain[locIdx])/(pid[locIdx+1]-pid[locIdx])*(pointingIDs[idx]-pid[locIdx]));
        }
    }

  return interpGain;
}

////////////////////////////////////////////////////////////////////////////////

void
Delta_vv_t::appendPidValue (int a_pid)
{
    pid.push_back(a_pid);
}

////////////////////////////////////////////////////////////////////////////////

void
Delta_vv_t::appendGainValue (double a_gain)
{
    gain.push_back(a_gain);
}

////////////////////////////////////////////////////////////////////////////////

void
Delta_vv_t::appendDipoleValue (double a_dipole)
{
    dipole.push_back(a_dipole);
}

////////////////////////////////////////////////////////////////////////////////

void
Delta_vv_t::mergeResults()
{
    merge_tables(pid, gain, dipole);
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
