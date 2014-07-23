/*
 *  Copyright (C) 2003, 2004 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef RINGSET_HPP
#define RINGSET_HPP

#include "rotmatrix.hpp"

#include <string>
#include <vector>

class ringset
{
public:
    rotmatrix ecl2gal;

    double invDeltaPhi;
    double invDeltaTheta;
    double phiOffset;
    double thetaOffset;

    int beammMax;
    int nPsi;
    int nPhi;
    int nTheta;
    int nPoints;
    int iOffset;

    std::vector<int> psiVector;
    std::vector<std::vector<std::vector<float> > > sky;
    std::vector<double> baseWgt;

    inline std::vector<double> weightN (double x) const;
    inline std::vector<double> interpolN (double theta, double phi) const;
    inline double interpolPsi(double omega, const std::vector<double> & kArr) const;
    void init (std::string fileName, int order, bool feedback_flag);

    void initializeWeights(int order);

    bool loaded;
    bool ok;

    ringset ();
    ringset (std::string fileName, int order, bool feedback_flag);

    std::vector<double>
    getIntensities (const std::vector<double> & theta,
                    const std::vector<double> & phi,
                    const std::vector<double> & psi) const;
    void getIntensities (const std::vector<double> & theta,
                         const std::vector<double> & phi,
                         const std::vector<double> & psi,
                         std::vector<double> & I,
                         std::vector<double> & Q,
                         std::vector<double> & U) const;
};

#endif
