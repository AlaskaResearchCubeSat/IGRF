#include <stdio.h>

#include "igrf.h"


int main(int argc,char **argv){
    VEC field;
    int nmax;
    //extrapolate model to desired date
    nmax=extrapsh(2013);
    //calculate magnetic field
    shval3(64.9261111,-147.4958333,6371.2,nmax,&field);
    //print
    printf("x = %f\ny = %f\nz = %f\n",field.c.x,field.c.y,field.c.z);
    return 0;
}
