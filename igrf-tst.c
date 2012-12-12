#include <stdio.h>

#include "igrf.h"


int main(int argc,char **argv){
    VEC field;
    int nmax;
    if(getshc("IGRFy10",1,0,13,1)){
        printf("Error loading model\n");
        return 1;
    }
    if(getshc("IGRFy10",0,0,13,2)){
        printf("Error loading model\n");
        return 1;
    }
    nmax=extrapsh(2013,2010,13,8);
    shval3(64.9261111,-147.4958333,6371.2,nmax,&field);
    printf("x = %f\ny = %f\nz = %f\n",field.c.x,field.c.y,field.c.z);
    return 0;
}
