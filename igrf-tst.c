#include <stdio.h>

#include "igrf.h"


int main(int argc,char **argv){
    VEC field;
    if(getshc("IGRFy10",1,0,13,1)){
        printf("Error loading model\n");
        return 1;
    }
    if(getshc("IGRFy10",0,0,13,2)){
        printf("Error loading model\n");
        return 1;
    }
    extrapsh(2013,2010,13,8,3);
    shval3(2,64.9261111,-147.4958333,6371.2,13,&field);
    printf("x = %f\ny = %f\nz = %f\n",field.c.x,field.c.y,field.c.z);
    return 0;
}
