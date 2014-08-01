#include <mpi.h>
#include "mademoiselle.hpp"
#include "configuration.hpp"
#include "da_capo_results.hpp"
#include "dipole_fit_results.hpp"
#include "io.hpp"
#include "logging.hpp"

#include <boost/format.hpp>
#include <numeric>
#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

////////////////////////////////////////////////////////////////////////////////

template <class T> T
productSum2D(const std::vector<std::vector<T>> & vecIn1,
             const std::vector<std::vector<T>> & vecIn2)
{
    T sum = 0.;
    for(size_t i = 0; i < vecIn1.size(); i++) {
        for(size_t j = 0; j < vecIn1[i].size(); j++) {
            sum += vecIn1[i][j] * vecIn2[i][j];
        }
    }
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

inline static std::string
gain_table_file_path(const Configuration & program_conf,
                     const Lfi_radiometer_t & radiometer)
{
    return (boost::format("%s/mademoiselle/%s_mademoiselle_gains.fits")
        % program_conf.getWithSubst("common.base_output_dir")
        % radiometer.shortName()).str();
}

////////////////////////////////////////////////////////////////////////////////

void
ffBuild(unsigned point,
        const std::vector<int> & hits,
        const std::vector<float> & dipole,
        std::vector<std::vector<double>>& ff)
{
    /* build ff std::vector SQRTFF  */
    // calculate sum of hits & dipole
    double sum_h = 0.0;//accumulate(hits.begin(),hits.end(),0);
    double sum_d = 0.0; //accumulate(dipole.begin(),dipole.end(),0.0);
    double sum_4 = 0.0;
    for (size_t i=0; i<dipole.size(); ++i) {
        sum_h += hits[i];
        sum_d += dipole[i];
        sum_4 += dipole[i]*dipole[i]/hits[i];
    }

    // allocate matrix space
    gsl_matrix* components   = gsl_matrix_alloc(2,2);
    gsl_matrix* eigenvectors = gsl_matrix_alloc(2,2);
    gsl_vector* eigenvalues  = gsl_vector_alloc(2);
    gsl_eigen_symmv_workspace* work = gsl_eigen_symmv_alloc(2);

    // set the matrix
    gsl_matrix_set(components, 0, 0, sum_h);
    gsl_matrix_set(components, 1, 0, sum_d);
    gsl_matrix_set(components, 0, 1, sum_d);
    gsl_matrix_set(components, 1, 1, sum_4);
    // solve the matrix
    gsl_eigen_symmv (components, eigenvalues, eigenvectors, work);
    gsl_eigen_symmv_free(work);
    gsl_matrix_free(components);

    // compute ff
    double d0 = sqrt(gsl_vector_get(eigenvalues,0));
    double d1 = sqrt(gsl_vector_get(eigenvalues,1));
    double v00 = gsl_matrix_get (eigenvectors, 0, 0);
    double v01 = gsl_matrix_get (eigenvectors, 0, 1);
    double v10 = gsl_matrix_get (eigenvectors, 1, 0);
    double v11 = gsl_matrix_get (eigenvectors, 1, 1);

    if ((d0 > 1.0e-10) && (d1 > 1.0e-10)) {
        ff[0][point] = v00*v00/d0 + v01*v01/d1;
        ff[1][point] = v00*v10/d0 + v01*v11/d1;
        ff[2][point] = v10*v10/d0 + v11*v11/d1;
    }
    /* END build ff std::vector SQRTFF  */
}

////////////////////////////////////////////////////////////////////////////////

/*
 * Mademoiselle destriping - single detector -
 */
double
mademoiselle(const std::vector<float> & maskMap,
             std::vector<Dipole_fit_t> & calibratedPID)
{
    Logger * log = Logger::get_instance();
    int mpi_rank = MPI::COMM_WORLD.Get_rank();

    // check not empty calibrated data
    if (calibratedPID.empty()) {
        const std::string msg("No binned data fed to Mademoiselle, aborting");
        log->error(msg);
        throw std::runtime_error(msg);
    }

    // get the nside of the binned tod map
    const int nside = calibratedPID[0].binned_data.nside;
    const int nPixelMap = 12 * nside * nside;

    // build initial map && ffstd::vector for each PID
    int *hMap = new int[nPixelMap];
    memset(hMap,0,nPixelMap*sizeof(int));
    float *dMap = new float[nPixelMap];
    memset(dMap,0.0,nPixelMap*sizeof(float));

    std::vector< std::vector<double> > ff (3,std::vector<double>(calibratedPID.size(),0.));
    for (size_t point = 0; point < calibratedPID.size(); point++) {
        std::vector<float> dipole = calibratedPID[point].binned_data.pix_model_mean;
        std::vector<int> hits      = calibratedPID[point].binned_data.pix_num_of_hits;
        std::vector<int> iPix      = calibratedPID[point].binned_data.pix_index;

        std::vector<double> data   = calibratedPID[point].binned_data.pix_data_sum;
        const double gain   = 1.0 / calibratedPID[point].gainv;
        const double gainV  = calibratedPID[point].gainv;
        const double offset = calibratedPID[point].offset;
        if (gain != 0.0) {
            for (size_t sample = 0; sample < hits.size(); sample++) {
                // compute TOD
                const double tod =
                    (data[sample] - gainV*dipole[sample] - hits[sample]*offset)*gain;

                // build map
                dMap[iPix[sample]] += tod;
                hMap[iPix[sample]] += hits[sample];
            }

            // calculate ff values for each point
            ffBuild(point, hits, dipole, ff);
        }
    }

    // collect map from nodes
    float *ddMap = new float[nPixelMap];
    MPI::COMM_WORLD.Allreduce(dMap, ddMap, nPixelMap, MPI_FLOAT, MPI_SUM);
    MPI::COMM_WORLD.Barrier();
    memcpy (dMap, ddMap, nPixelMap*sizeof(float));
    delete[] ddMap;

    int *hhMap = new int[nPixelMap];
    MPI::COMM_WORLD.Allreduce(hMap, hhMap, nPixelMap, MPI_INT, MPI_SUM);
    MPI::COMM_WORLD.Barrier();
    memcpy (hMap, hhMap, nPixelMap*sizeof(int));
    delete[] hhMap;

    // normalize map
    for (int i=0; i<nPixelMap; i++)
        if (hMap[i] > 0)
            dMap[i] /= hMap[i];

    /* build initial bb array */
    std::vector< std::vector<double> > bb(2,std::vector<double>(calibratedPID.size(),0.0));
    for (size_t point=0; point<calibratedPID.size(); point++) {
        double p1 = 0.;
        double p2 = 0.;

        auto dipole = calibratedPID[point].binned_data.pix_model_mean;
        auto hits   = calibratedPID[point].binned_data.pix_num_of_hits;
        auto iPix   = calibratedPID[point].binned_data.pix_index;

        auto data = calibratedPID[point].binned_data.pix_data_sum;
        double gain   = 1.0 / calibratedPID[point].gainv;
        double gainV  = calibratedPID[point].gainv;
        double offset = calibratedPID[point].offset;
        if (gain != 0.0) {
            for (size_t sample=0; sample < hits.size(); sample++) {
                if (maskMap[iPix[sample]] != 0) {
                    // compute TOD
                    double ddtod = (data[sample] - gainV*dipole[sample] - hits[sample]*offset)*gain;
                    double dtod = ddtod - (dMap[iPix[sample]]*hits[sample]);
                    p1 += dtod;
                    p2 += dtod*(dipole[sample]/hits[sample]);
                }
            }
        }

        bb[0][point] = ff[0][point] * p1 +ff[1][point] * p2;
        bb[1][point] = ff[1][point] * p1 +ff[2][point] * p2;
    }

    /*  ITERATE */
    std::vector<std::vector<double>> aa(2,std::vector<double>(calibratedPID.size(),0.0)); // sono le baselines
    std::vector<std::vector<double>> r = bb;
    std::vector<std::vector<double>> p = bb;

    // free space
    bb.erase(bb.begin(),bb.end());

    // get rz init collecting from nodes
    double rz = productSum2D<double>(r,r);
    MPI::COMM_WORLD.Allreduce(&rz, &rz, 1, MPI_DOUBLE, MPI_SUM);
    MPI::COMM_WORLD.Barrier();

    double rzinit = rz;
    log->info(boost::format("rzinit = %.12f") % rzinit);

    /* iteration */
    for (int step = 0; step < 800; step++) {
        // reset map
        delete[] dMap;
        dMap = new float[nPixelMap];
        memset(dMap, 0.0, nPixelMap * sizeof(float));

        // calculate new map
        for (size_t point = 0; point < calibratedPID.size(); point++) {
            auto hits   = calibratedPID[point].binned_data.pix_num_of_hits;
            auto dipole = calibratedPID[point].binned_data.pix_model_mean;
            auto iPix   = calibratedPID[point].binned_data.pix_index;

            if (calibratedPID[point].gainv != 0.0) {
                for (size_t sample=0; sample < hits.size(); sample++) {
                    if (maskMap[iPix[sample]] == 0)
                        continue;

                    double p0 = p[0][point];
                    double p1 = p[1][point];

                    dMap[iPix[sample]] += ( (p0*ff[0][point] + p1*ff[1][point])*hits[sample] +
                        (p0*ff[1][point] + p1*ff[2][point])*dipole[sample]);
                }
            }
        }

        // collect map from nodes
        ddMap = new float[nPixelMap];
        MPI::COMM_WORLD.Allreduce(dMap, ddMap, nPixelMap, MPI_FLOAT, MPI_SUM);
        MPI::COMM_WORLD.Barrier();
        memcpy (dMap, ddMap, nPixelMap*sizeof(float));
        delete [] ddMap;

        // normalize map
        for (int i=0; i<nPixelMap; i++) {
            if (hMap[i] > 0)
                dMap[i] /= hMap[i];
        }

        // BOH part
        std::vector< std::vector<double> > ap(2,std::vector<double>(calibratedPID.size(),0.0));
        for (size_t point = 0; point < calibratedPID.size(); point++) {
            double pw1 = 0.0;
            double pw2 = 0.0;

            auto hits   = calibratedPID[point].binned_data.pix_num_of_hits;
            auto dipole = calibratedPID[point].binned_data.pix_model_mean;
            auto iPix   = calibratedPID[point].binned_data.pix_index;

            if (calibratedPID[point].gainv != 0.0) {
                for (size_t sample = 0; sample < hits.size(); sample++) {
                    pw1 += hits[sample]*dMap[iPix[sample]];
                    pw2 += dipole[sample]*dMap[iPix[sample]];
                }
            }

            ap[0][point] = p[0][point] - pw1*ff[0][point] - pw2*ff[1][point];
            ap[1][point] = p[1][point] - pw1*ff[1][point] - pw2*ff[2][point];
        }

        // get global pap from nodes
        double pap = productSum2D<double>(p,ap);
        MPI::COMM_WORLD.Allreduce(&pap, &pap, 1, MPI_DOUBLE, MPI_SUM);
        MPI::COMM_WORLD.Barrier();

        // update aa and r
        for (size_t i = 0; i < calibratedPID.size(); i++) {
            aa[0][i] += rz/pap*p[0][i];
            aa[1][i] += rz/pap*p[1][i];

            r[0][i] -= rz/pap*ap[0][i];
            r[1][i] -= rz/pap*ap[1][i];
        }

        // compute new rz global
        double rzo = rz;
        rz = productSum2D<double>(r,r);
        MPI::COMM_WORLD.Allreduce(&rz, &rz, 1, MPI_DOUBLE, MPI_SUM);
        MPI::COMM_WORLD.Barrier();

        // exit condition
        if ((rz/rzinit < 1.e-12) && (step > 4)) {
            if (mpi_rank == 0)
                log->info(boost::format("rz/rzInit = %.12e at step %d")
                          % (rz / rzinit) % step);
            break;
        }

        // update p
        for (size_t i = 0; i < calibratedPID.size(); i++) {
            p[0][i] = r[0][i] + (rz/rzo)*p[0][i];
            p[1][i] = r[1][i] + (rz/rzo)*p[1][i];
        }
        MPI::COMM_WORLD.Barrier();

    }// end iteration

    for (size_t i = 0; i < calibratedPID.size(); i++) {
        double pw1 = aa[0][i];
        double pw2 = aa[1][i];
        aa[0][i] = pw1*ff[0][i] + pw2*ff[1][i];
        aa[1][i] = pw1*ff[1][i] + pw2*ff[2][i];
    }
    /* END ITERATE */

    // set the new gain
    for (size_t point = 0; point < calibratedPID.size(); point++) {
        double offset = calibratedPID[point].offset;
        double gainV  = calibratedPID[point].gainv;

        calibratedPID[point].offset = offset + aa[0][point] * gainV;
        double gainNew = gainV + aa[1][point] * gainV;
        if (gainNew > 1.e-20)
            calibratedPID[point].gainv = gainNew;
        else
            calibratedPID[point].gainv = 0.0;
    }

    delete [] hMap;
    delete [] dMap;
    return rzinit;
}

////////////////////////////////////////////////////////////////////////////////

void
run_mademoiselle(const Configuration & program_conf,
                 const Configuration & storage_conf,
                 const Lfi_radiometer_t & user_rad,
                 const std::vector<Pointing_t> & list_of_pointings,
                 Dipole_fit_results_t & dipole_fit_results,
                 Da_capo_results_t & da_capo_results)
{
    Logger * log = Logger::get_instance();
    const int mpi_size = MPI::COMM_WORLD.Get_size();
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();

    log->info("Starting module Mademoiselle");
    log->increase_indent();

    for (size_t iteration=0; iteration<60; ++iteration) {
        log->debug(boost::format("CG iteration %1%") % iteration);

        const double init = mademoiselle(dipole_fit_results.mask.pixels,
                                         dipole_fit_results.dipole_fits);
        if (init < 1.e-10) {
            if (mpi_rank == 0)
                log->info(boost::format("Break at CG iteration %1%")
                          % iteration);
            break;
        }
    }
    MPI::COMM_WORLD.Barrier();

    {
        std::string pointings("[");
        for(auto cur_fit = dipole_fit_results.dipole_fits.begin();
            cur_fit != dipole_fit_results.dipole_fits.end();
            ++cur_fit)
        {
            if(cur_fit != dipole_fit_results.dipole_fits.begin())
                pointings += ", ";
            pointings += 
                (boost::format("%1%") 
                 % cur_fit->binned_data.pointing_id).str();
        }

        log->info("The pointings contained in the list of dipole fits are "
                  + pointings);
    }

    Gain_table_t & gain_table = da_capo_results.gain_table;
    gain_table.pointingIds.resize(list_of_pointings.size());
    gain_table.gain.assign(list_of_pointings.size(), 0.0);
    gain_table.offset.assign(list_of_pointings.size(), 0.0);

    for(size_t idx = 0; idx < gain_table.pointingIds.size(); ++idx) {
        gain_table.pointingIds[idx] = list_of_pointings[idx].id;
    }

    for(size_t fit_idx = 0;
        fit_idx < dipole_fit_results.dipole_fits.size();
        ++fit_idx)
    {
        auto const & fit = dipole_fit_results.dipole_fits[fit_idx];
        auto const item =
            std::lower_bound(gain_table.pointingIds.begin(),
                             gain_table.pointingIds.end(),
                             fit.binned_data.pointing_id);
        if(item == gain_table.pointingIds.end() ||
           *item != fit.binned_data.pointing_id)
        {
            log->error(boost::format("Mademoiselle fitted pointing ID"
                                     "%1%, but this was not in the"
                                     "list of pIDs to process")
                       % fit.binned_data.pointing_id);
            continue;
        }

        const int pid_idx = std::distance(gain_table.pointingIds.begin(), item);
        gain_table.gain.at(pid_idx) = 1.0 / fit.gainv;
        gain_table.offset.at(pid_idx) = fit.offset;
    }

    // To understand the following lines, have a look at the implementation
    // of "run_dipole_fit" (in dipole_fit.cpp).

    Lfi_radiometer_t real_radiometer;
    if(mpi_rank % 2 == 0)
        real_radiometer = user_rad;
    else
        real_radiometer = user_rad.twinRadiometer();

    gain_table.mergeResults();
    gain_table.selectRadiometerGains(mpi_rank % 2, 2,
                                     dipole_fit_results.pids_per_process);

    if(mpi_rank == 0 || mpi_rank == 1) {
        const std::string gain_file_path(gain_table_file_path(program_conf,
                                                              real_radiometer));
        log->info(boost::format("Saving dipoleFit gains for "
                                "radiometer %1% into %2%")
                  % real_radiometer.shortName() % gain_file_path);
        save_gain_table(ensure_path_exists(gain_file_path),
                        real_radiometer, gain_table);
    }

    log->decrease_indent();
    log->info("Quitting module Mademoiselle");
}
