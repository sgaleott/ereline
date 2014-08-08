#include <mpi.h>

#include "da_capo.hpp"
#include "da_capo_results.hpp"
#include "dipole_fit_results.hpp"
#include "io.hpp"
#include "misc.hpp"
#include "configuration.hpp"
#include "logging.hpp"

#include <algorithm>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>

#define MPI_BARRIER { \
    Logger * private_log = Logger::get_instance(); \
    private_log->debug(boost::format("Entering MPI::COMM_WORLD.Barrier() at %1%, line %2% (%3%)") \
               % __FILE__ % __LINE__ % __PRETTY_FUNCTION__); \
    MPI::COMM_WORLD.Barrier(); \
    private_log->debug(boost::format("Exiting MPI::COMM_WORLD.Barrier() at %1%, line %2% (%3%)") \
               % __FILE__ % __LINE__ % __PRETTY_FUNCTION__); \
}

extern "C" {
#include "chealpix.h"
}

/**
 * Take binned data and dipolefit gains as input.
 * It also take the number of streams (single diode/radiometer, horn using all streams)
 **/
daCapo::daCapo(const std::vector<Dipole_fit_t> & locallyBinnedData,
               const std::vector<float> & mask,
               bool constraint,
               const Dipole_parameters_t & dipole_params)
{
  Logger * log = Logger::get_instance();
  log->debug(boost::format("Creating daCapo object using daCapo::daCapo("
               "locallyBinnedData, mask, %1%, dipole_params)")
         % constraint);

  sizeMPI = MPI::COMM_WORLD.Get_size();
  rankMPI = MPI::COMM_WORLD.Get_rank();

  nside = locallyBinnedData.at(0).binned_data.nside;
  nPixelMap = 12 * nside * nside;
  log->debug(boost::format("DaCapo will create maps with NSIDE = %1%") % nside);
  rzinit=0;

  initializeLocmap(locallyBinnedData);
  initializeFullmap();
  applyMask(locallyBinnedData, mask);
  initializeConstraint(constraint, dipole_params);

  localMap = std::vector<double>(pixelIndexLocal.size()+1, 0.);
  fullMap = std::vector<double>(pixelIndexFull.size(), 0.);
  ccFull = std::vector<double>(pixelIndexFull.size(), 0.);
  dipolenorm = std::vector<double>(3, 0.);
  preconditioner = std::vector< std::vector<double> >(locallyBinnedData.size(), std::vector<double>(3, 0.));
}

daCapo::daCapo (const std::vector<Dipole_fit_t> & locallyBinnedData,
                const std::vector<float> & mask,
                const std::vector<double> & constraint)
{
  Logger * log = Logger::get_instance();
  log->debug("Creating daCapo object using daCapo::daCapo("
         "locallyBinnedData, mask, constraint)");

  sizeMPI = MPI::COMM_WORLD.Get_size();
  rankMPI = MPI::COMM_WORLD.Get_rank();

  nside = locallyBinnedData.at(0).binned_data.nside;
  nPixelMap = 12 * nside * nside;
  rzinit=0;

  initializeLocmap(locallyBinnedData);
  initializeFullmap();
  applyMask(locallyBinnedData, mask);
  initializeConstraint(constraint);

  localMap = std::vector<double>(pixelIndexLocal.size()+1, 0.);
  fullMap = std::vector<double>(pixelIndexFull.size(), 0.);
  ccFull = std::vector<double>(pixelIndexFull.size(), 0.);
  dipolenorm = std::vector<double>(3, 0.);
  preconditioner = std::vector< std::vector<double> >(locallyBinnedData.size(), std::vector<double>(3, 0.));
}

/**
 * Initialize a Solar Dipole Map
 **/
void
daCapo::initializeConstraint(bool constraint,
                             const Dipole_parameters_t & solar_dipole)
{
  Logger * log = Logger::get_instance();

  if (!constraint)
    return;

  log->debug(boost::format("Running daCapo::initializeConstraint"
                           "(%1%, solar_dipole)")
             % constraint);

  // Build Dipole Constraint Map
  for (size_t idx=0; idx<pixelIndexFull.size(); idx++)
    {
      double theta, phi;
      pix2ang_nest(nside, pixelIndexFull[idx], &theta, &phi);
      double cartesianPixel[3];
      angToCart(theta, phi, cartesianPixel);
      constraintMap.push_back(solar_dipole.axis[0] * cartesianPixel[0] +
                              solar_dipole.axis[1] * cartesianPixel[1] +
                              solar_dipole.axis[2] * cartesianPixel[2]);
    }
}

/**
 * Initialize a Solar Dipole Map
 **/
void
daCapo::initializeConstraint(const std::vector<double> & constraint)
{
  Logger * log = Logger::get_instance();
  log->debug("Running daCapo::initializeConstraint(constraint)");

  // Build Dipole Constraint Map
  for (size_t idx=0; idx<pixelIndexFull.size(); idx++)
    {
      constraintMap.push_back(constraint[pixelIndexFull[idx]]);
    }
}

/**
 * Construct a "local map" as sequence of pixels hit by the local toi.
 * A Local map consists of the pixels hit by the toi of the current process.
 * At return, array iloc, of same size as ipix, contains an index pointing from toi to the local map.
 * Output: ipix_loc, iloc
 **/
void
daCapo::initializeLocmap(const std::vector<Dipole_fit_t> & binnedData)
{
  Logger * log = Logger::get_instance();
  log->debug(boost::format("Running daCapo::initializeLocmap(binnedData), "
                           "binnedData has %1% elements")
             % binnedData.size());

  std::vector<int> tmpPixels;
  for (size_t ipp=0; ipp<binnedData.size(); ipp++)
    {
      std::vector<int> pix_indexPid = binnedData[ipp].binned_data.pix_index;
      for (size_t i=0; i<pix_indexPid.size(); i++)
        {
          tmpPixels.push_back(pix_indexPid[i]);
        }
    }

  //Compress. This now holds a list of pixels hit by the local toi.
  pixelIndexLocal = sortAndCount(tmpPixels);

  //Build a search table from pixel number to local map, to speed up the next step
  std::vector<int> search_table(nPixelMap, 0);
  for (size_t i=0; i<pixelIndexLocal.size(); i++)
    search_table[pixelIndexLocal[i]]=static_cast<int>(i);

  //Construct a pointer from toi to the local map
  for (size_t ipp=0; ipp<binnedData.size(); ipp++)
    {
      std::vector<int> pix_indexPid = binnedData[ipp].binned_data.pix_index;
      std::vector<int> pixMap;
      for (size_t i=0; i<pix_indexPid.size(); i++)
        {
          pixMap.push_back(search_table[pix_indexPid[i]]);
        }
      pixelIndexLocalMap.push_back(pixMap);
    }
}

/**
 * Construct the full_map structure, and auxiliary tables that control the transfer between
 * local and full maps. A full map consists of all pixels covered by the data.
 * Output:
 * ipix_full: list of pixels covered by the "full map" of this process
 * sendcnt[irank]:numer of adjacent pixels of the local map projected to the full map of process IRANK
 * ifull[ipix]: a pointer from a locmap pixel to a fullmap pixel
 **/
void
daCapo::initializeFullmap()
{
  Logger * log = Logger::get_instance();
  log->debug("Running daCapo::initializeFullmap");

  int ipixMax=pixelIndexLocal[pixelIndexLocal.size()-1];
  MPI::COMM_WORLD.Allreduce(&ipixMax, &ipixMax, 1, MPI_INT, MPI_MAX);
  MPI_BARRIER;
  int pixPerProc = 1+ipixMax/sizeMPI;

  sendcnt = std::vector<int>(sizeMPI, 0);
  for (size_t i=0; i<pixelIndexLocal.size(); ++i)
    sendcnt[pixelIndexLocal[i]/pixPerProc]++;

  std::vector<int> dummy(sizeMPI, 0);
  MPI::COMM_WORLD.Alltoall (sendcnt.data(), static_cast<int>(sendcnt.size()/sizeMPI), MPI_INT,
                            dummy.data(), static_cast<int>(sendcnt.size()/sizeMPI), MPI_INT);
  MPI_BARRIER;

  std::vector<int> disin(sizeMPI, 0);
  std::vector<int> disout(sizeMPI, 0);
  for (int i=1; i<sizeMPI; ++i)
    {
      disin [i]=disin [i-1]+sendcnt[i-1];
      disout[i]=disout[i-1]+dummy[i-1];
    }

  std::vector<int> tmpPixels(disout[sizeMPI-1]+dummy[sizeMPI-1], 0);
  MPI::COMM_WORLD.Alltoallv(pixelIndexLocal.data(), sendcnt.data(), disin.data(), MPI_INT,
                            tmpPixels.data(), dummy.data(), disout.data(), MPI_INT);

  // Compress. This now holds a list of pixels hit by the local toi.
  pixelIndexFull = sortAndCount(tmpPixels);

  // Construct an index table pointing from local map to full map
  int offset=0;
  pixelIndexFullMap = std::vector<int>(pixelIndexLocal.size(), 0);
  for (int irank=0; irank<sizeMPI; irank++) //Loop over fullmaps
    {
      int npix_tmp = static_cast<int>(pixelIndexFull.size());

      MPI::COMM_WORLD.Bcast(&npix_tmp, 1, MPI_INT, irank);
      MPI::COMM_WORLD.Barrier();

      std::vector<int> ipix_tmp(npix_tmp);
      if (rankMPI==irank)
        for (int i=0; i<npix_tmp; i++)
          ipix_tmp[i]=pixelIndexFull[i];

      MPI::COMM_WORLD.Bcast(ipix_tmp.data(), npix_tmp, MPI_INT, irank);
      MPI_BARRIER;

      int m=0;
      for (int i=0; i<sendcnt[irank]; i++)
        {
          while (pixelIndexLocal[offset+i]>ipix_tmp[m])
            {
              m++;
            }
          pixelIndexFullMap[offset+i]=m;
        }
      offset +=sendcnt[irank];
    }
}

/**
 * Apply an input mask.
 * Masked samples are repointed to a dummy pixel.
 **/
void
daCapo::applyMask(const std::vector<Dipole_fit_t> & binnedData,
                  const std::vector<float> & mask)
{
  Logger * log = Logger::get_instance();
  log->debug(boost::format("Running daCapo::applyMask(binnedData, mask), "
                           "binnedData has %1% elements")
             % binnedData.size());

  int dummy_pixel = static_cast<int>(pixelIndexLocal.size());
  for (size_t ipp=0; ipp<binnedData.size(); ipp++)
    {
      std::vector<int> pix_indexPid = binnedData[ipp].binned_data.pix_index;
      for (size_t i=0; i<pix_indexPid.size(); i++)
        {
          if (mask[pix_indexPid[i]]==0) pixelIndexLocalMap[ipp][i]=dummy_pixel;
        }
    }
}

/**
 * Construct and invert the hit map.
 **/
void
daCapo::constructCCmatrix(const std::vector<Dipole_fit_t> & binnedData)
{
  Logger * log = Logger::get_instance();
  log->debug(boost::format("Running daCapo::constructCCmatrix(binnedData), "
                           "binnedData has %1% elements")
             % binnedData.size());

  std::vector<double> ccLoc(pixelIndexFullMap.size()+1, 0.);

  for (size_t ipp=0; ipp<binnedData.size(); ipp++)
    {
      std::vector<int> hits = binnedData[ipp].binned_data.pix_num_of_hits;
      double gain = binnedData[ipp].gainv;
      for (size_t i=0; i<hits.size(); i++)
        {
          ccLoc[pixelIndexLocalMap[ipp][i]] += gain*gain*hits[i];
        }
    }

  // Local CC to Full CC
  int offset=0;
  ccFull.assign(ccFull.size(), 0.);
  for (int irank=0; irank<sizeMPI; irank++) //Loop over fullmaps
    {
      int nbf = static_cast<int>(pixelIndexFull.size());

      MPI::COMM_WORLD.Bcast(&nbf, 1, MPI_INT, irank);
      MPI_BARRIER;

      std::vector<double> cbuf(nbf, 0.);
      for (int i=0; i<sendcnt[irank]; i++)
        cbuf[pixelIndexFullMap[offset+i]] = ccLoc[offset+i];

      if(! cbuf.empty()) {
        MPI::COMM_WORLD.Reduce (cbuf.data(), ccFull.data(), static_cast<int>(cbuf.size()),
                                MPI_DOUBLE, MPI_SUM, irank);
      } else {
        std::string msg;
        if(binnedData.empty()) {
            msg = (boost::format("Empty \"cbuf\" in daCapo::constructCCmatrix, "
                                 "irank=%1%, sendcnt[irank]=%2%")
                   % irank
                   % sendcnt.at(irank)).str();
        } else {
            msg = (boost::format("Empty \"cbuf\" in daCapo::constructCCmatrix, "
                                 "pID range is [%1%, %2%], irank=%3%, "
                                 "sendcnt[irank]=%4%")
                   % binnedData.front().binned_data.pointing_id
                   % binnedData.back().binned_data.pointing_id
                   % irank
                   % sendcnt.at(irank)).str();
        }
        log->warning(msg);
      }
      MPI_BARRIER;

      offset += sendcnt.at(irank);
    }

  for (size_t ip=0; ip<pixelIndexFull.size(); ip++)
    ccFull[ip] = ccFull[ip]>1.e-30 ? 1/ccFull[ip] : 0;
}

/**
 * Update dipolenorm with the current cc matrix.
 **/
void
daCapo::updateDipolenorm()
{
  Logger * log = Logger::get_instance();
  log->debug("Running daCapo::updateDipolenorm()");

  if (constraintMap.size()==0)
    return;

  dipolenorm.assign(3, 0.);
  for (size_t ip=0; ip<constraintMap.size(); ip++)
    {
      dipolenorm[0] += ccFull[ip];
      dipolenorm[1] += ccFull[ip]*constraintMap[ip];
      dipolenorm[2] += ccFull[ip]*constraintMap[ip]*constraintMap[ip];
    }

  MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, dipolenorm.data(), static_cast<int>(dipolenorm.size()), MPI_DOUBLE, MPI_SUM);
  MPI_BARRIER;

  invert2_eig(dipolenorm);
}

/**
 * Construct the preconditioner as the 2x2 block diagonal of the full matrix.
 * The same matrix gives an error estimate for the gain&baseline solution.
 * Output: prec
 * The 2x2 block is a symmetric matrix and is stored as a 3 element vector:
 * 0: base x base
 * 1: cross base x gain
 * 2: gain x gain
 **/
void
daCapo::buildPreconditioner(const std::vector<Dipole_fit_t> & binnedData)
{
  Logger * log = Logger::get_instance();
  log->debug(boost::format("Running daCapo::buildPreconditioner(binnedData), "
                           "binnedData has %1% elements")
             % binnedData.size());

  for (size_t ipp=0; ipp<binnedData.size(); ipp++)
    {
      std::vector<int> hits = binnedData[ipp].binned_data.pix_num_of_hits;
      std::vector<float> signal = binnedData[ipp].binned_data.pix_model_mean;

       std::vector<double> ippPrec(3, 0.);
       for (size_t i=0; i<hits.size(); i++)
         {
           ippPrec[0] += static_cast<float>(hits[i]);
           ippPrec[1] += static_cast<float>(hits[i])*signal[i];
           ippPrec[2] += static_cast<float>(hits[i])*signal[i]*signal[i];
         }
       invert2_eig(ippPrec);
       preconditioner[ipp]=ippPrec;
    }
}

/**
 * Coadd toi on an existing locmap.
 * Locmap must be preallocated and initialized.
 * Output: locmap
 **/
void
daCapo::toiToLocmap(const std::vector<Dipole_fit_t> & binnedData)
{
  Logger * log = Logger::get_instance();
  log->debug(boost::format("Running daCapo::toiToLocmap(binnedData), "
                           "binnedData has %1% elements")
             % binnedData.size());

  for (size_t ipp=0; ipp<pixelIndexLocalMap.size(); ipp++)
    {
      std::vector<double> signal = binnedData[ipp].binned_data.pix_data_sum;

      const double gain = binnedData[ipp].gainv;

      for (size_t i=0; i<signal.size(); i++)
        {
          localMap[pixelIndexLocalMap[ipp][i]] += gain*signal[i];
        }
    }
}

/**
 * This routine collects local maps (locmap) from processes and coadds them into a full map (map_full).
 * Ouput array map_full must be preallocated.
 * Output: map_full
 **/
void
daCapo::locToFullmap()
{
  Logger * log = Logger::get_instance();
  log->debug("Running daCapo::locToFullmap()");

  fullMap.assign(fullMap.size(), 0.);

  int offset = 0;
  for (int irank=0; irank<sizeMPI; irank++) // loop over receiving process
    {
      int npix_buff = static_cast<int>(pixelIndexFull.size());
      MPI::COMM_WORLD.Bcast(&npix_buff, 1, MPI_INT, irank);
      MPI_BARRIER;

      std::vector<double> mbuff(npix_buff, 0.);
      for (int i=0; i<sendcnt[irank]; i++)
        mbuff[pixelIndexFullMap[offset+i]]=localMap[offset+i];

      if(! mbuff.empty()) {
        MPI::COMM_WORLD.Reduce (mbuff.data(), fullMap.data(), static_cast<int>(mbuff.size()),
                                MPI_DOUBLE, MPI_SUM, irank);
      }
      MPI_BARRIER;

      offset +=sendcnt[irank];
  }
  MPI_BARRIER;
}

/**
 * Multiply a full sky map by the pre-inverted cc matrix
 * Output: map
 **/
void
daCapo::ccMultiply()
{
  Logger * log = Logger::get_instance();
  log->debug("Running daCapo::ccMultiply()");

  for (size_t ip=0; ip<fullMap.size(); ip++)
    {
      if (ccFull[ip]>0)
        fullMap[ip] *= ccFull[ip];
      else
        fullMap[ip]=HEALPIX_UNDEF;
    }
  MPI_BARRIER;
}

/**
 * Apply a constraint to force the dipole of the sky map to zero.
 **/
void
daCapo::applyConstraint()
{
  Logger * log = Logger::get_instance();
  log->debug("Running daCapo::applyConstraint()");

  if (constraintMap.size()==0)
    return;

  double ddotmap1=0;
  double ddotmap2=0;
  for (size_t ip=0; ip<constraintMap.size(); ip++)
    if (ccFull[ip]>0)
      {
        ddotmap1 += fullMap[ip];
        ddotmap2 += fullMap[ip]*constraintMap[ip];
      }

  MPI::COMM_WORLD.Allreduce(&ddotmap1, &ddotmap1, 1, MPI_DOUBLE, MPI_SUM);
  MPI_BARRIER;
  MPI::COMM_WORLD.Allreduce(&ddotmap2, &ddotmap2, 1, MPI_DOUBLE, MPI_SUM);
  MPI_BARRIER;

  double d0 = ddotmap1;
  double d1 = ddotmap2;
  ddotmap1 = d0*dipolenorm[0] +d1*dipolenorm[1];
  ddotmap2 = d0*dipolenorm[1] +d1*dipolenorm[2];

  for (size_t ip=0; ip<constraintMap.size(); ip++)
    {
      double dm = ddotmap1 + constraintMap[ip]*ddotmap2;
      fullMap[ip] -= ccFull[ip]*dm;
    }
  MPI_BARRIER;
}

/**
 * Scan a locmap from a full map.
 * Locmap must be preallocated.
 * Output: locmap
 **/
void
daCapo::fullToLocmap()
{
  Logger * log = Logger::get_instance();
  log->debug("Running daCapo::fullToLocmap()");

  localMap.assign(localMap.size(), 0.);

  int offset=0;
  for (int irank=0; irank<sizeMPI; irank++)
    {
      int npix_buff=static_cast<int>(fullMap.size());
      MPI::COMM_WORLD.Bcast(&npix_buff, 1, MPI_INT, irank);
      MPI_BARRIER;

      std::vector<double> mbuff(npix_buff, 0.);

      if (rankMPI==irank)
        for (size_t i=0; i<fullMap.size(); i++)
          mbuff[i]=fullMap[i];

      MPI::COMM_WORLD.Bcast(mbuff.data(), npix_buff, MPI_DOUBLE, irank);
      MPI_BARRIER;

      for (int i=0; i<sendcnt[irank]; i++)
        localMap[offset+i]=mbuff[pixelIndexFullMap[offset+i]];

      offset +=sendcnt[irank];
    }
  MPI_BARRIER;
}

/*
Map binning operation, including constraint
 */
void
daCapo::applyCC()
{
  Logger * log = Logger::get_instance();
  log->debug("Running daCapo::applyCC()");
  log->increase_indent();

  locToFullmap();
  ccMultiply();
  applyConstraint();
  fullToLocmap();

  log->decrease_indent();
  log->debug("daCapo::applyCC() completed");
}

/**
 * Pick baselines from a locmap, subtract them from input toi, and sum into baselines
 * Output: p
 */
void
daCapo::subtractMapFromTod(const std::vector<Dipole_fit_t> & binnedData,
                           basevec &p)
{
  Logger * log = Logger::get_instance();
  log->debug(boost::format("Running daCapo::subtractMapFromTod(binnedData, p), "
                           "binnedData has %1% elements")
             % binnedData.size());

  int dummy_pixel=static_cast<int>(localMap.size())-1;

  p.SetZero();
  for (size_t ipp=0; ipp<pixelIndexLocalMap.size(); ipp++)
    {
      size_t nn = pixelIndexLocalMap[ipp].size();
      std::vector<double> tmp1(nn,0.);

      std::vector<double> signal = binnedData[ipp].binned_data.pix_data_sum;
      std::vector<float> dipole = binnedData[ipp].binned_data.pix_model_mean;
      std::vector<int> hits = binnedData[ipp].binned_data.pix_num_of_hits;
      const double gain = binnedData[ipp].gainv;

      for (size_t i=0; i<nn; i++)
        if (pixelIndexLocalMap[ipp][i]<dummy_pixel)
          {
            tmp1[i]=signal[i]-gain*hits[i]*localMap[pixelIndexLocalMap[ipp][i]];
          }
      for (size_t i=0; i<nn; i++)
        {
          p.base[ipp] += tmp1[i];
          p.gain[ipp] += tmp1[i]*dipole[i];
        }
    }
}

/**
 * Coadd baselines onto a locmap.
 * Locmap must be preinitialized.
 * Output: locmap
 **/
void
daCapo::baseToLocmap(const std::vector<Dipole_fit_t> & binnedData,
                     const basevec &p)
{
  Logger * log = Logger::get_instance();
  log->debug(boost::format("Running daCapo::baseToLocmap(binnedData, p), "
                           "binnedData has %1% elements")
             % binnedData.size());

  for (size_t ipp=0; ipp<pixelIndexLocalMap.size(); ipp++)
    {
      std::vector<int> hits = binnedData[ipp].binned_data.pix_num_of_hits;
      std::vector<float> dipole = binnedData[ipp].binned_data.pix_model_mean;

      const double gain = binnedData[ipp].gainv;

      for (size_t i=0; i<pixelIndexLocalMap[ipp].size(); i++)
        {
          localMap[pixelIndexLocalMap[ipp][i]] +=
            gain*hits[i]*(p.base[ipp]+p.gain[ipp]*dipole[i]);
        }
    }
}

/**
 * Pick baselines from a locmap and subtract them from input baselines
 * Output: p
 **/
void
daCapo::subtractMapFromBase(const std::vector<Dipole_fit_t> & binnedData,
                            basevec &p)
{
  Logger * log = Logger::get_instance();
  log->debug(boost::format("Running daCapo::subtractMapFromBase(binnedData, p), "
                           "binnedData has %1% elements")
             % binnedData.size());

  int dummy_pixel = static_cast<int>(localMap.size())-1;

  basevec r=p;
  p.SetZero();

  for (size_t ipp=0; ipp<binnedData.size(); ipp++)
    {
      std::vector<int> hits = binnedData[ipp].binned_data.pix_num_of_hits;
      std::vector<float> dipole = binnedData[ipp].binned_data.pix_model_mean;

      const double gain = binnedData[ipp].gainv;

      size_t nn=pixelIndexLocalMap[ipp].size();
      std::vector<double> tmp1(nn, 0.);

      for (size_t i=0; i<nn; i++)
        if (pixelIndexLocalMap[ipp][i]<dummy_pixel)
          {
            tmp1[i] = r.base[ipp]+r.gain[ipp]*dipole[i]
              -gain*localMap[pixelIndexLocalMap[ipp][i]];
          }

      for (size_t i=0; i<nn; i++)
        {
          p.base[ipp] += tmp1[i]*static_cast<float>(hits[i]);
          p.gain[ipp] += tmp1[i]*static_cast<float>(hits[i])*static_cast<float>(dipole[i]);
        }
    }
}

/**
 * Apply the preconditioner on the r vector.
 * Output: z
 **/
void
daCapo::applyPreconditioner(const basevec &r, basevec &z)
{
  Logger * log = Logger::get_instance();
  log->debug("Running daCapo::applyPreconditioner(r, z)");

  for (size_t i=0; i<preconditioner.size(); i++)
    {
      z.base[i] = preconditioner[i][0]*r.base[i] +preconditioner[i][1]*r.gain[i];
      z.gain[i] = preconditioner[i][1]*r.base[i] +preconditioner[i][2]*r.gain[i];
    }
}

/**
 * Update signal with the current map
 **/
void
daCapo::updateSignal(std::vector<Dipole_fit_t> & binnedData)
{
  Logger * log = Logger::get_instance();
  log->debug(boost::format("Running daCapo::updateSignal(binnedData), "
                           "binnedData has %1% elements")
             % binnedData.size());

  for (size_t ipp=0; ipp<pixelIndexLocalMap.size(); ipp++)
    {
      std::vector<float> dipole = binnedData[ipp].binned_data.pix_model_mean;
      std::vector<float> localDipole;
      for (size_t i=0; i<pixelIndexLocalMap[ipp].size(); i++)
        {
          localDipole.push_back(dipole[i] + static_cast<float>(localMap[pixelIndexLocalMap[ipp][i]]));
        }
      binnedData[ipp].binned_data.pix_model_mean = localDipole;
    }
}

/**
 * Run Iterative calibration
 **/
double
daCapo::iterativeCalibration(std::vector<Dipole_fit_t> & binnedData,
                             bool firstLoop)
{
  Logger * log = Logger::get_instance();
  log->debug(boost::format("Running daCapo::iterativeCalibration(binnedData, %2%), "
                           "binnedData has %1% elements")
             % binnedData.size() % firstLoop);
  log->increase_indent();

  constructCCmatrix(binnedData);
  updateDipolenorm();
  buildPreconditioner(binnedData);
  localMap.assign(localMap.size(), 0.);
  toiToLocmap(binnedData);

  applyCC();

  basevec yb(static_cast<int>(binnedData.size()));
  subtractMapFromTod(binnedData, yb);

  basevec aa(static_cast<int>(binnedData.size()));
  for (size_t point=0; point<binnedData.size(); point++)
    {
      double gain = binnedData[point].gainv;
      double offset = binnedData[point].offset;

      aa.SetValues(gain, offset, static_cast<int>(point));
    }

  localMap.assign(localMap.size(), 0.);
  baseToLocmap(binnedData, aa);

  applyCC();

  basevec r=aa;
  subtractMapFromBase(binnedData, r);

  r.Scale(-1);
  r.Add(yb);

  basevec z=r;
  applyPreconditioner(r, z);

  double rzo;
  double rz=r.Dotprod(z);

  MPI::COMM_WORLD.Allreduce(&rz, &rz, 1, MPI_DOUBLE, MPI_SUM);
  MPI_BARRIER;

  if (firstLoop == true)
    rzinit=rz;
  basevec p=z;

  for (int istep=0; istep<1000; istep++)
    {
      localMap.assign(localMap.size(), 0.);
      baseToLocmap(binnedData, p);

      applyCC();

      basevec ap=p;
      subtractMapFromBase(binnedData, ap);

      double pap=p.Dotprod(ap);
      MPI::COMM_WORLD.Allreduce(&pap, &pap, 1, MPI_DOUBLE, MPI_SUM);
      MPI_BARRIER;

      double alpha=rz/pap;
      ap.Scale(-alpha);
      r.Add(ap);

      applyPreconditioner(r, z);

      rzo = rz;
      rz = r.Dotprod(z);
      MPI::COMM_WORLD.Allreduce(&rz, &rz, 1, MPI_DOUBLE, MPI_SUM);
      MPI_BARRIER;

      basevec phelp=p;
      phelp.Scale(alpha);
      aa.Add(phelp);

      if (rz/rzinit<1e-12)
        {
          if (rankMPI == 0)
              std::cout << "   rz/rzinit " << std::setw(12) << rz/rzinit
                        << "   at step " << istep << std::endl << std::flush;
          break;
        }

      double beta=rz/rzo;
      p.Scale(beta);
      p.Add(z);
    }
  MPI_BARRIER;

  localMap.assign(localMap.size(), 0.);
  toiToLocmap(binnedData);

  p = aa;
  p.Scale(-1);
  baseToLocmap(binnedData, p);
  applyCC();

  updateSignal(binnedData);

  double dmax=0;
  for (size_t i=0; i<aa.gain.size(); i++)
    {
      double gainDiff = aa.gain[i]-binnedData[i].gainv;
      dmax = std::max(dmax, std::abs(gainDiff));
    }

  MPI::COMM_WORLD.Allreduce(&dmax, &dmax, 1, MPI_DOUBLE, MPI_MAX);
  MPI_BARRIER;

  // Set the new gain
  for (size_t point=0; point<binnedData.size(); point++)
    {
      binnedData[point].offset = aa.base[point];
      binnedData[point].gainv = aa.gain[point];
    }

  log->decrease_indent();
  log->debug(boost::format("daCapo::iterativeCalibration completed "
                           "with result %1%")
             % dmax);

  return dmax;
}

////////////////////////////////////////////////////////////////////////////////

inline static std::string
gain_table_file_path(const Configuration & program_conf,
                     const Lfi_radiometer_t & radiometer)
{
    return (boost::format("%s/da_capo/%s_da_capo_gains.fits")
        % program_conf.getWithSubst("common.base_output_dir")
        % radiometer.shortName()).str();
}

////////////////////////////////////////////////////////////////////////////////

/* Note that this function will manipulate the values in
 * "dipole_fit_results". */
void
run_da_capo(const Configuration & program_conf,
            const Configuration & storage_conf,
            const Lfi_radiometer_t & user_rad,
            const std::vector<Pointing_t> & list_of_pointings,
            Dipole_fit_results_t & dipole_fit_results,
            Da_capo_results_t & da_capo_results)
{
    Logger * log = Logger::get_instance();
    const int mpi_size = MPI::COMM_WORLD.Get_size();
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();

    log->info("Starting module daCapo");
    log->increase_indent();

    daCapo calibrator(dipole_fit_results.dipole_fits,
                      dipole_fit_results.mask.pixels,
                      program_conf.get<bool>("da_capo.constraint"),
                      read_dipole_fit_params(program_conf));
    bool firstLoop = true;
    for (size_t iteration=0; iteration<10000; ++iteration)
    {
        if (mpi_rank == 0)
            log->debug(boost::format("CG iteration %d") % iteration);

        double init =
            calibrator.iterativeCalibration(dipole_fit_results.dipole_fits,
                                            firstLoop);
        firstLoop = false;
        if (init < 1.e-6)
        {
            if (mpi_rank == 0)
                log->debug(boost::format("Break at CG iteration %1%, "
                                         "rzInit = %2%")
                           % iteration % init);
            break;
        }
    }
    MPI_BARRIER;

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
            log->error(boost::format("Da Capo fitted pointing ID"
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
    da_capo_results.radiometer = real_radiometer;

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
    log->info("Quitting module daCapo");
}
