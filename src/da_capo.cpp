#include "da_capo.hpp"
#include "misc.hpp"
#include "configuration.hpp"
#include "logging.hpp"

#include <mpi.h>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>

extern "C" {
#include "chealpix.h"
}

/**
 * Take binned data and dipolefit gains as input.
 * It also take the number of streams (single diode/radiometer, horn using all streams)
 **/
daCapo::daCapo(std::vector<dipoleFit> & binnedData, 
	       std::vector<float> & mask, 
	       bool constraint)
{
  sizeMPI = MPI::COMM_WORLD.Get_size();
  rankMPI = MPI::COMM_WORLD.Get_rank();

  nSide = binnedData[0].getNSide();
  nPixelMap = 12 * nSide * nSide;
  rzinit=0;

  initializeLocmap(binnedData);
  initializeFullmap();
  applyMask(binnedData, mask);
  initializeConstraint(constraint);

  localMap = std::vector<double>(pixelIndexLocal.size()+1, 0.);
  fullMap = std::vector<double>(pixelIndexFull.size(), 0.);
  ccFull = std::vector<double>(pixelIndexFull.size(), 0.);
  dipolenorm = std::vector<double>(3, 0.);
  preconditioner = std::vector< std::vector<double> >(binnedData.size(), std::vector<double>(3, 0.));
}

daCapo::daCapo (std::vector<dipoleFit> & binnedData, std::vector<float> & mask, std::vector<double> & constraint)
{
  sizeMPI = MPI::COMM_WORLD.Get_size();
  rankMPI = MPI::COMM_WORLD.Get_rank();

  nSide = binnedData[0].getNSide();
  nPixelMap = 12 * nSide * nSide;
  rzinit=0;

  initializeLocmap(binnedData);
  initializeFullmap();
  applyMask(binnedData, mask);
  initializeConstraint(constraint);

  localMap = std::vector<double>(pixelIndexLocal.size()+1, 0.);
  fullMap = std::vector<double>(pixelIndexFull.size(), 0.);
  ccFull = std::vector<double>(pixelIndexFull.size(), 0.);
  dipolenorm = std::vector<double>(3, 0.);
  preconditioner = std::vector< std::vector<double> >(binnedData.size(), std::vector<double>(3, 0.));
}

/**
 * Initialize a Solar Dipole Map
 **/
void
daCapo::initializeConstraint(bool constraint)
{
  if (!constraint)
    return;

  // Build Dipole Constraint Map 
  std::vector<double> SOLSYSDIR_V=angToCart(SOLSYSDIR_ECL_THETA, SOLSYSDIR_ECL_PHI);
  for (size_t idx=0; idx<pixelIndexFull.size(); idx++) 
    {
      double theta, phi;
      pix2ang_nest(nSide, pixelIndexFull[idx], &theta, &phi);
      std::vector<double> cartesianPixel=angToCart(theta, phi);
      constraintMap.push_back(SOLSYSDIR_V[0]*cartesianPixel[0]+SOLSYSDIR_V[1]*cartesianPixel[1]+SOLSYSDIR_V[2]*cartesianPixel[2]);
    }
}

/**
 * Initialize a Solar Dipole Map
 **/
void
daCapo::initializeConstraint(std::vector<double> & constraint)
{
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
daCapo::initializeLocmap(const std::vector<dipoleFit> & binnedData)
{
  std::vector<int> tmpPixels;
  for (size_t ipp=0; ipp<binnedData.size(); ipp++)
    {
      std::vector<int> pixIndexPid = binnedData[ipp].getPixIndex();
      for (size_t i=0; i<pixIndexPid.size(); i++)
	{
	  tmpPixels.push_back(pixIndexPid[i]);
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
      std::vector<int> pixIndexPid = binnedData[ipp].getPixIndex();
      std::vector<int> pixMap;
      for (size_t i=0; i<pixIndexPid.size(); i++)
	{
	  pixMap.push_back(search_table[pixIndexPid[i]]);
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
  int ipixMax=pixelIndexLocal[pixelIndexLocal.size()-1];
  MPI::COMM_WORLD.Allreduce(&ipixMax, &ipixMax, 1, MPI_INT, MPI_MAX);
  MPI::COMM_WORLD.Barrier();
  int pixPerProc = 1+ipixMax/sizeMPI;

  sendcnt = std::vector<int>(sizeMPI, 0);
  for (size_t i=0; i<pixelIndexLocal.size(); ++i)
    sendcnt[pixelIndexLocal[i]/pixPerProc]++;

  std::vector<int> dummy(sizeMPI, 0);
  MPI::COMM_WORLD.Alltoall (sendcnt.data(), static_cast<int>(sendcnt.size()/sizeMPI), MPI_INT,
			    dummy.data(), static_cast<int>(sendcnt.size()/sizeMPI), MPI_INT);
  MPI::COMM_WORLD.Barrier();

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
      MPI::COMM_WORLD.Barrier();

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
daCapo::applyMask(const std::vector<dipoleFit> & binnedData, 
		  const std::vector<float> & mask)
{
  int dummy_pixel = static_cast<int>(pixelIndexLocal.size());
  for (size_t ipp=0; ipp<binnedData.size(); ipp++)
    {
      std::vector<int> pixIndexPid = binnedData[ipp].getPixIndex();
      for (size_t i=0; i<pixIndexPid.size(); i++)
	{
	  if (mask[pixIndexPid[i]]==0) pixelIndexLocalMap[ipp][i]=dummy_pixel;
	}
    }
}

/**
 * Construct and invert the hit map.
 **/
void 
daCapo::constructCCmatrix(const std::vector<dipoleFit> & binnedData)
{
  std::vector<double> ccLoc(pixelIndexFullMap.size()+1, 0.);
  
  for (size_t ipp=0; ipp<binnedData.size(); ipp++)
    {
      std::vector<int> hits = binnedData[ipp].getPixSumHits();
      double gain = binnedData[ipp].getGainV();
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
      MPI::COMM_WORLD.Barrier();

      std::vector<double> cbuf(nbf, 0.);
      for (int i=0; i<sendcnt[irank]; i++)
	cbuf[pixelIndexFullMap[offset+i]] = ccLoc[offset+i];
      
      MPI::COMM_WORLD.Reduce (cbuf.data(), ccFull.data(), static_cast<int>(cbuf.size()),
			      MPI_DOUBLE, MPI_SUM, irank);
      MPI::COMM_WORLD.Barrier();

      offset += sendcnt[irank];
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
  MPI::COMM_WORLD.Barrier();

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
daCapo::buildPreconditioner(const std::vector<dipoleFit> & binnedData)
{
  for (size_t ipp=0; ipp<binnedData.size(); ipp++)
    {
      std::vector<int> hits = binnedData[ipp].getPixSumHits();
      std::vector<float> signal = binnedData[ipp].getPixSumDipole();

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
daCapo::toiToLocmap(const std::vector<dipoleFit> & binnedData)
{
  for (size_t ipp=0; ipp<pixelIndexLocalMap.size(); ipp++)
    {
      std::vector<double> signal = binnedData[ipp].getPixSumData();
      
      double gain = binnedData[ipp].getGainV();
      
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
  fullMap.assign(fullMap.size(), 0.);

  int offset = 0;
  for (int irank=0; irank<sizeMPI; irank++) // loop over receiving process
    {
      int npix_buff = static_cast<int>(pixelIndexFull.size());
      MPI::COMM_WORLD.Bcast(&npix_buff, 1, MPI_INT, irank);
      MPI::COMM_WORLD.Barrier();

      std::vector<double> mbuff(npix_buff, 0.);
      for (int i=0; i<sendcnt[irank]; i++)
	mbuff[pixelIndexFullMap[offset+i]]=localMap[offset+i];

      MPI::COMM_WORLD.Reduce (mbuff.data(), fullMap.data(), static_cast<int>(mbuff.size()),
			      MPI_DOUBLE, MPI_SUM, irank);
      MPI::COMM_WORLD.Barrier();
  
      offset +=sendcnt[irank];
  }
  MPI::COMM_WORLD.Barrier();
}

/**
 * Multiply a full sky map by the pre-inverted cc matrix
 * Output: map
 **/
void 
daCapo::ccMultiply()
{
  for (size_t ip=0; ip<fullMap.size(); ip++)
    {
      if (ccFull[ip]>0)
	fullMap[ip] *= ccFull[ip];
      else
	fullMap[ip]=HEALPIX_UNDEF;
    }
  MPI::COMM_WORLD.Barrier();
}

/**
 * Apply a constraint to force the dipole of the sky map to zero.
 **/
void 
daCapo::applyConstraint()
{
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
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Allreduce(&ddotmap2, &ddotmap2, 1, MPI_DOUBLE, MPI_SUM);
  MPI::COMM_WORLD.Barrier();

  double d0 = ddotmap1;
  double d1 = ddotmap2;
  ddotmap1 = d0*dipolenorm[0] +d1*dipolenorm[1];
  ddotmap2 = d0*dipolenorm[1] +d1*dipolenorm[2];

  for (size_t ip=0; ip<constraintMap.size(); ip++)
    {
      double dm = ddotmap1 + constraintMap[ip]*ddotmap2;
      fullMap[ip] -= ccFull[ip]*dm;
    }
  MPI::COMM_WORLD.Barrier();
}

/**
 * Scan a locmap from a full map.
 * Locmap must be preallocated.
 * Output: locmap
 **/
void 
daCapo::fullToLocmap()
{
  localMap.assign(localMap.size(), 0.);

  int offset=0;
  for (int irank=0; irank<sizeMPI; irank++)
    {
      int npix_buff=static_cast<int>(fullMap.size());
      MPI::COMM_WORLD.Bcast(&npix_buff, 1, MPI_INT, irank);
      MPI::COMM_WORLD.Barrier();

      std::vector<double> mbuff(npix_buff, 0.);

      if (rankMPI==irank)
	for (size_t i=0; i<fullMap.size(); i++)
	  mbuff[i]=fullMap[i];

      MPI::COMM_WORLD.Bcast(mbuff.data(), npix_buff, MPI_DOUBLE, irank);
      MPI::COMM_WORLD.Barrier();

      for (int i=0; i<sendcnt[irank]; i++)
	localMap[offset+i]=mbuff[pixelIndexFullMap[offset+i]];

      offset +=sendcnt[irank];
    }
  MPI::COMM_WORLD.Barrier();
}

/*
Map binning operation, including constraint
 */
void 
daCapo::applyCC()
{
  locToFullmap();
  ccMultiply();
  applyConstraint();
  fullToLocmap();
}

/**
 * Pick baselines from a locmap, subtract them from input toi, and sum into baselines
 * Output: p
 */
void 
daCapo::subtractMapFromTod(const std::vector<dipoleFit> & binnedData,
			   basevec &p)
{
  int dummy_pixel=static_cast<int>(localMap.size())-1;

  p.SetZero();
  for (size_t ipp=0; ipp<pixelIndexLocalMap.size(); ipp++)
    {
      size_t nn = pixelIndexLocalMap[ipp].size();
      std::vector<double> tmp1(nn,0.);

      std::vector<double> signal = binnedData[ipp].getPixSumData();
      std::vector<float> dipole = binnedData[ipp].getPixSumDipole();
      std::vector<int> hits = binnedData[ipp].getPixSumHits();
      double gain = binnedData[ipp].getGainV();

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
daCapo::baseToLocmap(const std::vector<dipoleFit> & binnedData, 
		     const basevec &p)
{
  for (size_t ipp=0; ipp<pixelIndexLocalMap.size(); ipp++)
    {
      std::vector<int> hits = binnedData[ipp].getPixSumHits();
      std::vector<float> dipole = binnedData[ipp].getPixSumDipole();

      double gain = binnedData[ipp].getGainV();
      
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
daCapo::subtractMapFromBase(const std::vector<dipoleFit> & binnedData, 
			    basevec &p)
{
  int dummy_pixel = static_cast<int>(localMap.size())-1;

  basevec r=p;
  p.SetZero();

  for (size_t ipp=0; ipp<binnedData.size(); ipp++)
    {
      std::vector<int> hits = binnedData[ipp].getPixSumHits();
      std::vector<float> dipole = binnedData[ipp].getPixSumDipole();

      double gain = binnedData[ipp].getGainV();

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
daCapo::updateSignal(std::vector<dipoleFit> & binnedData)
{
  for (size_t ipp=0; ipp<pixelIndexLocalMap.size(); ipp++)
    {
      std::vector<float> dipole = binnedData[ipp].getPixSumDipole();
      std::vector<float> localDipole;
      for (size_t i=0; i<pixelIndexLocalMap[ipp].size(); i++)
	{
	  localDipole.push_back(dipole[i] + static_cast<float>(localMap[pixelIndexLocalMap[ipp][i]]));
	}
      binnedData[ipp].setPixSumDipole(localDipole);
    }
}

/**
 * Run Iterative calibration
 **/
double
daCapo::iterativeCalibration (std::vector<dipoleFit> & binnedData, 
			      bool firstLoop)
{
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
      double gain = binnedData[point].getGainV();
      double offset = binnedData[point].getOffset();
      
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
  MPI::COMM_WORLD.Barrier();
  
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
      MPI::COMM_WORLD.Barrier();
      
      double alpha=rz/pap;
      ap.Scale(-alpha);
      r.Add(ap);
      
      applyPreconditioner(r, z);
      
      rzo = rz;
      rz = r.Dotprod(z);
      MPI::COMM_WORLD.Allreduce(&rz, &rz, 1, MPI_DOUBLE, MPI_SUM);
      MPI::COMM_WORLD.Barrier();
      
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
  MPI::COMM_WORLD.Barrier();
  
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
      double gainDiff = aa.gain[i]-binnedData[i].getGainV();
      dmax = std::max(dmax, std::abs(gainDiff));
    }
  
  MPI::COMM_WORLD.Allreduce(&dmax, &dmax, 1, MPI_DOUBLE, MPI_MAX);
  MPI::COMM_WORLD.Barrier();

  // Set the new gain
  for (size_t point=0; point<binnedData.size(); point++) 
    {
      binnedData[point].setOffset(aa.base[point]);  
      binnedData[point].setGainV(aa.gain[point]);
    }
  
  return dmax;
}

void
run_da_capo(const Configuration & program_conf,
	    const Configuration & storage_conf)
{
  Logger * log = Logger::get_instance();

  log->info("Starting module daCapo");
  log->info("Quitting module daCapo");
}
