#include <msp430.h>

//use majority function so the timer
//can be read while it is running
short readTA(void){
  int a=TAR,b=TAR,c=TAR;
  return (a&b)|(a&c)|(b&c);
}

//setup timer A to run off 32.768kHz xtal
void init_timerA(void){
    //setup timer A for ADC usage
  //TACTL=TASSEL_1|ID_3|TACLR;
  TACTL=TASSEL_1|ID_0|TACLR;
}

//start timer A in continuous mode
void start_timerA(void){
//start timer A
  TACTL|=MC_2;
}
