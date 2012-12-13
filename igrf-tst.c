#include <stdio.h>

#include "igrf.h"

#define PI 3.14159265358979323846
#define RAD2DEG (180.0/PI)

int main(int argc,char **argv){
    VEC field;
    int nmax;
    //extrapolate model to desired date
    nmax=extrapsh(2013);
    //calculate magnetic field
    shval3(64.9261111/RAD2DEG,-147.4958333/RAD2DEG,6371.2,nmax,&field);
    //print
    printf("x = %f\ny = %f\nz = %f\n",field.c.x,field.c.y,field.c.z);
    return 0;
}
