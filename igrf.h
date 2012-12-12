#ifndef _IGRF_H
#define _IGRF_H

#include "vector.h"

//Read spherical harmonic coefficients from model into array
int getshc(const char *file,int iflag,long int strec,int nmax_of_gh,int gh);
//Extrapolate model
int extrapsh(double date,double dte1,int nmax1,int nmax2);
//Interpolate between models
int interpsh(double date,double dte1,int nmax1,double dte2,int nmax2);
//Calculates field components from models
int shval3(double flat,double flon,double elev,int nmax,VEC *dest);

#endif

