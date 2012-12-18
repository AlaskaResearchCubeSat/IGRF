#ifndef _IGRF_H
#define _IGRF_H

#include "vector.h"

//Read spherical harmonic coefficients from model into array
int getshc(const char *file,int iflag,long int strec,int nmax_of_gh,int gh);
//Extrapolate model
int extrapsh(float date);
//Interpolate between models
int interpsh(float date,float dte1,int nmax1,float dte2,int nmax2);
//Calculates field components from models
int shval3(float flat,float flon,float elev,int nmax,VEC *dest);

#endif

