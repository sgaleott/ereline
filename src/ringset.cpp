/*
 * The Class is fixed on the Ecliptic to Galactic Conversion
 * in order to be compliant with Galactic Sidelobes beam convolution
 *
 * The Class is able to produce white noise given a specific instrumentDb
 * and the radiometer
 */

#include "ringset.hpp"
#include "logging.hpp"

#include <fitsio.h>
#include <math.h>
#include <cstdio>

const double PI = 3.141592653589793238462643383279502884197;
const double HALFPI = 1.570796326794896619231321691639751442099;
const double DEG_TO_RAD = PI / 180.0;

////////////////////////////////////////////////////////////////////////

ringset::ringset ()
{
    loaded = false;
    ok = true;
}

////////////////////////////////////////////////////////////////////////

/*
 * Initializiation reading the ringset.LS_ringset object
 * from dB. Initialize the pointing conversion from ecliptic
 * to galactic and the order of the interpolation.
 */
ringset::ringset (std::string fileName, int order, bool feedback_flag)
{
    loaded = false;
    ok = true;
    init (fileName, order, feedback_flag);
}

////////////////////////////////////////////////////////////////////////

void
ringset::init (std::string fileName, int order, bool feedback_flag)
{
    // Initialize Rotmatrix Ecliptic to Galactic
    ecl2gal = rotmatrix(-0.054882486, -0.993821033, -0.096476249,
			 0.494116468, -0.110993846,  0.862281440,
			-0.867661702, -0.000346354,  0.497154957);

    Logger * log = Logger::get_instance();
    log->info(boost::format("Reading ringsets from FITS file %1%")
	      % fileName);

    fitsfile * fptr;
    int status = 0;
    fits_open_file(&fptr, fileName.c_str(), READONLY, &status);
    if(status != 0) {
	char msg[30];
	fits_get_errstatus(status, msg);
	throw std::runtime_error((boost::format("Error reading %1%: %2%") 
				  % fileName % msg).str());
    }

    fits_read_key(fptr, TINT, "beam_mmax", &beammMax, NULL, &status);
    fits_read_key(fptr, TINT, "nphi", &nPhi, NULL, &status);

    double dPhi;
    fits_read_key(fptr, TDOUBLE, "dphi", &dPhi, NULL, &status);
    invDeltaPhi = 1./(DEG_TO_RAD*dPhi);

    double phi0;
    fits_read_key(fptr, TDOUBLE, "phi0", &phi0, NULL, &status);
    phiOffset = (DEG_TO_RAD*phi0)/(DEG_TO_RAD*dPhi);

    fits_read_key(fptr, TINT, "ntheta", &nTheta, NULL, &status);

    double dTheta;
    fits_read_key(fptr, TDOUBLE, "dtheta", &dTheta, NULL, &status);
    invDeltaTheta = 1./(DEG_TO_RAD*dTheta);

    double theta0;
    fits_read_key(fptr, TDOUBLE, "theta0", &theta0, NULL, &status);
    thetaOffset = (DEG_TO_RAD*theta0)/(DEG_TO_RAD*dTheta);

    // Read data
    fits_movabs_hdu(fptr, 3, NULL, &status);
    long numRows;
    fits_get_num_rows(fptr, &numRows, &status);

    nPsi=0;
    std::vector<int> tmpPresent(numRows);
    fits_read_col(fptr, TINT, 1, 1, 1, numRows, NULL, tmpPresent.data(), NULL, &status);
    std::vector<int> tmpVector(beammMax+1, -1);
    psiVector.swap(tmpVector);
    for (size_t idx=0; idx<tmpPresent.size(); ++idx)
    {
	psiVector[tmpPresent[idx]] = nPsi;
	(tmpPresent[idx] == 0) ? nPsi+=1 : nPsi+=2;
    }

    fits_movabs_hdu(fptr, 2, NULL, &status);
    fits_get_num_rows(fptr, &numRows, &status);

    std::vector<std::vector<std::vector<float> > > tmpMat(nPhi, std::vector< std::vector<float> >(nTheta, std::vector<float>(nPsi, 0)));
    sky.swap(tmpMat);
    for (int psiIdx=0; psiIdx<nPsi; ++psiIdx)
    {
	if(feedback_flag && psiIdx % 10 == 0) {
	    const int bar_length = 40;
	    char bar[bar_length + 1];
	    int percent = (psiIdx * bar_length) / (nPsi - 1);
	    int char_idx;

	    bar[bar_length] = 0;
	    for(char_idx = 0; char_idx < bar_length; ++char_idx) {
		bar[char_idx] = char_idx < percent ? '*' : '.';
	    }
	    std::fprintf(stderr,
			 "|%s| %d/%d, %d%%\r", bar, psiIdx + 1, nPsi, (psiIdx * 100) / (nPsi - 1));
	}

	std::vector<float> tmpSky(nTheta * nPhi);
	fits_read_col(fptr, TFLOAT, 1,
		      1 + psiIdx * nTheta * nPhi, 1, nTheta * nPhi,
		      NULL, tmpSky.data(), NULL, &status);

	for (int thetaIdx = 0; thetaIdx < nTheta; ++thetaIdx) {
	    for (int phiIdx = 0; phiIdx < nPhi; ++phiIdx) {
		sky[phiIdx][thetaIdx][psiIdx] = 
		    (psiIdx == 0) 
		    ? tmpSky[thetaIdx * nPhi + phiIdx] 
		    : 2 * tmpSky[thetaIdx * nPhi + phiIdx];
	    }
	}
    }

    if(feedback_flag) {
	std::fputc('\n', stderr);
    }

    fits_close_file(fptr, &status);
    if(status != 0) {
	fits_report_error(stderr, status);
	loaded = false;
	ok = false;
	return;
    }

    initializeWeights(order);
    log->info(boost::format("Ringsets loaded from %1% successfully")
	      % fileName);

    loaded = true;
    ok = true;
}

////////////////////////////////////////////////////////////////////////

void
ringset::initializeWeights(int order)
{
    nPoints = order+1;
    iOffset = order/2;

    std::vector<double> tmpWgt(nPoints, 1);
    baseWgt.swap(tmpWgt);
    for (int extIdx = 0; extIdx < nPoints; ++extIdx) {
	for (int intIdx=0; intIdx < nPoints; ++intIdx) {
	    if (intIdx != extIdx) {
		baseWgt[extIdx] *= extIdx - intIdx;
	    }
	}
    }
  
    for (int idx=0; idx < nPoints; ++idx)
	baseWgt[idx] = 1. / baseWgt[idx];
}

////////////////////////////////////////////////////////////////////////

inline std::vector<double>
ringset::weightN (double x)
{
    std::vector<double> wgt = baseWgt;

    double mul1 = x;
    double mul2 = x - nPoints + 1;

    for (int idx=1; idx<nPoints; ++idx)
    {
	wgt[idx] *= mul1;
	wgt[nPoints - idx - 1] *= mul2;
	mul1 *= x - idx;
	mul2 *= x - nPoints + idx + 1;
    }
    return wgt;
}

////////////////////////////////////////////////////////////////////////

inline std::vector<double>
ringset::interpolN (double theta, double phi)
{
    double frac = theta*invDeltaTheta - thetaOffset;
    int iTheta0 = int (frac) - iOffset;

    if (iTheta0 > (nTheta-nPoints)) 
	iTheta0 = nTheta-nPoints;
    if (iTheta0 < 0) 
	iTheta0 = 0;
    frac -= iTheta0;
    std::vector<double> wgt1 = weightN (frac);

    frac = phi * invDeltaPhi - phiOffset;
    if(frac >= 0) {
	frac = (frac < double(nPhi)) ? frac : remainder(frac, double(nPhi));
    } else {
	frac = remainder(frac, double(nPhi)) + double(nPhi);
    }

    int iPhi0 = int (frac) - iOffset;
    frac -= iPhi0;
    if (iPhi0 >= nPhi) 
	iPhi0 -= nPhi;
    if (iPhi0 < 0) 
	iPhi0 += nPhi;
    std::vector<double> wgt2 = weightN (frac);

    std::vector<double> result(nPsi, 0);
    int iPhi = iPhi0;
    for (int i = 0; i < nPoints; ++i)
    {
	for (int j = 0; j < nPoints; ++j)
	{
	    double weight = wgt2[i] * wgt1[j];
	    const float * ref = &sky[iPhi][iTheta0 + j][0];
	    for (int k=0; k < nPsi; ++k)
		result[k] += weight * ref[k];
	}
	if (++iPhi >= nPhi) 
	    iPhi -= nPhi;
    }
    return result;
}

inline double 
ringset::interpolPsi(double omega, const std::vector<double> & kArr)
{
    double result = (psiVector[0] >= 0) ? kArr[0] : 0;
    if (nPsi <= 1) 
	return result;

    double cosang = 1;
    double sinang = 0;
    double sinomg = sin(omega);
    double cosomg = cos(omega);

    for (int k = 1; k <= beammMax; ++k)
    {
	const double tmp = sinang * cosomg + cosang * sinomg;
	cosang = cosang * cosomg - sinang * sinomg;
	sinang = tmp;
	if (psiVector[k] >= 0)
	    result += cosang * kArr[psiVector[k]] - sinang * kArr[psiVector[k] + 1];
    }
    return result;
}

////////////////////////////////////////////////////////////////////////

std::vector<double>
ringset::getIntensities (const std::vector<double> & theta,
			 const std::vector<double> & phi, 
			 const std::vector<double> & psi)
{ 
    // Init rotation trigonometry functions
    double cosBeta  = ecl2gal.entry[2][2];
    double sinBeta  = sqrt(1 - cosBeta * cosBeta);

    double cosAlpha = -ecl2gal.entry[2][1] / sinBeta;
    double sinAlpha = ecl2gal.entry[2][0] / sinBeta;
    double alpha = atan2(sinAlpha, cosAlpha) + 2 * PI;

    double cosGamma = ecl2gal.entry[1][2] / sinBeta;
    double sinGamma = -ecl2gal.entry[0][2] / sinBeta;
    double gamma = atan2(sinGamma, cosGamma);

    std::vector<double> result;
    for (size_t idx = 0; idx < theta.size(); ++idx)
    {
	// Convert ecliptic to galactic
	double z = -sinBeta * sin(theta[idx]) * sin(phi[idx] - alpha) 
	    + cosBeta * cos(theta[idx]);

	double yh = cosBeta * sin(theta[idx]) * sin(phi[idx] - alpha) 
	    + sinBeta*cos(theta[idx]);
	double xh = sin(theta[idx]) * cos(phi[idx] - alpha);

	std::vector<double> tmpArr = interpolN (acos(z), atan2(yh,xh) + gamma);

	double ys = -cos(phi[idx] - alpha) * sinBeta;
	double xs = sin(theta[idx]) * cosBeta 
	    + cos(theta[idx]) * sin(phi[idx] - alpha) * sinBeta;
	double omegaWg = psi[idx] + atan2(ys,xs) - HALFPI;

	// * 2 Multiplication because of LevelS SCR 234
	result.push_back(interpolPsi (omegaWg, tmpArr) * 2);
    }

    return result;
}

