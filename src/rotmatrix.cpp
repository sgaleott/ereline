/*
 *  Class for rotation transforms in 3D space
 *
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <cmath>

#include "rotmatrix.hpp"

rotmatrix::rotmatrix ()
{
}

rotmatrix::rotmatrix (double a00, double a01, double a02,
		      double a10, double a11, double a12,
		      double a20, double a21, double a22)
{
    entry[0][0]=a00; entry[0][1]=a01; entry[0][2]=a02;
    entry[1][0]=a10; entry[1][1]=a11; entry[1][2]=a12;
    entry[2][0]=a20; entry[2][1]=a21; entry[2][2]=a22;
}

void rotmatrix::setRotmatrix(double a00, double a01, double a02,
			     double a10, double a11, double a12,
			     double a20, double a21, double a22)
{
    entry[0][0]=a00; entry[0][1]=a01; entry[0][2]=a02;
    entry[1][0]=a10; entry[1][1]=a11; entry[1][2]=a12;
    entry[2][0]=a20; entry[2][1]=a21; entry[2][2]=a22;
}

void
rotmatrix::makeAxisRotationTransform (const std::vector<double> &axis, 
				      double angle)
{
    double sa=sin(angle), ca=cos(angle);
    double ica=1-ca;
    entry[0][0] = axis[0]*axis[0]*ica + ca;
    entry[1][1] = axis[1]*axis[1]*ica + ca;
    entry[2][2] = axis[2]*axis[2]*ica + ca;
    double t1 = axis[0]*axis[1]*ica, t2 = axis[2]*sa;
    entry[1][0] = t1 + t2;
    entry[0][1] = t1 - t2;
    t1 = axis[0]*axis[2]*ica; t2 = axis[1]*sa;
    entry[2][0] = t1 - t2;
    entry[0][2] = t1 + t2;
    t1 = axis[1]*axis[2]*ica; t2 = axis[0]*sa;
    entry[1][2] = t1 - t2;
    entry[2][1] = t1 + t2;
}

std::vector<double>
rotmatrix::transform (const std::vector<double> &vec)
{
    std::vector<double> outVec;
    outVec.push_back(vec[0]*entry[0][0] + vec[1]*entry[0][1] + vec[2]*entry[0][2]);
    outVec.push_back(vec[0]*entry[1][0] + vec[1]*entry[1][1] + vec[2]*entry[1][2]);
    outVec.push_back(vec[0]*entry[2][0] + vec[1]*entry[2][1] + vec[2]*entry[2][2]);

    return outVec;
}

rotmatrix
operator* (const rotmatrix &a, const rotmatrix &b)
{
    rotmatrix res;
    for (int i=0; i<3; ++i)
	for (int j=0; j<3; ++j)
	    res.entry[i][j] = a.entry[i][0] * b.entry[0][j]
		+ a.entry[i][1] * b.entry[1][j]
		+ a.entry[i][2] * b.entry[2][j];
    return res;
}

std::ostream &operator<< (std::ostream &os, const rotmatrix &mat)
{
    for (int i=0;i<3;++i) {
	os << '[' << mat.entry[i][0] << ','
	   << mat.entry[i][1] << ','
	   << mat.entry[i][2] << ']' << std::endl;
    }
    return os;
}
