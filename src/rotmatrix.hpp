/*
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef ROTMATRIX_HPP_
#define ROTMATRIX_HPP_

#include <vector>
#include <ostream>

struct rotmatrix
{
    double entry[3][3];

    rotmatrix ();
    rotmatrix (double a00, double a01, double a02,
               double a10, double a11, double a12,
               double a20, double a21, double a22);
    void setRotmatrix(double a00, double a01, double a02,
                      double a10, double a11, double a12,
                      double a20, double a21, double a22);
    void makeAxisRotationTransform (const std::vector<double> &axis,
                                    double angle);
    std::vector<double> transform (const std::vector<double> &vec);
};

rotmatrix operator* (const rotmatrix &a, const rotmatrix &b);
std::ostream &operator<< (std::ostream &os, const rotmatrix &mat);


#endif
