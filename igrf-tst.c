#include <stdio.h>

#include "igrf.h"


int main(int argc,char **argv){
    if(getshc("IGRFy10",1,0,13,1)){
        printf("Error loading model\n");
        return 1;
    }
    if(getshc("IGRFy10",0,0,13,2)){
        printf("Error loading model\n");
        return 1;
    }
    extrapsh(2013,2010,13,8,3);
    extrapsh(2014,2010,13,8,4);
    shval3(2,64.9261111,-147.4958333,6371.2,13,3,0,0,0,0);
    printf("x = %f\ny = %f\nz = %f\n",x,y,z);
    return 0;
}
