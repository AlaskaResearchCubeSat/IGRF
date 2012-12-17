#include <msp430.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "uart.h"
#include "terminal.h"
#include "igrf.h"


#define PI 3.14159265358979323846
#define RAD2DEG (180.0/PI)

int tst(char **argv,unsigned short argc);
int resetCmd(char **argv,unsigned short argc);

//table of commands with help
CMD_SPEC cmd_tbl[]={{"help"," [command]",helpCmd},
                   {"reset","\r\n\t""Reset the msp430.",resetCmd},
                   {"tst","\r\n\t""Test IGRF conversion",tst},
                   //end of list
                   {NULL,NULL,NULL}};

int tst(char **argv,unsigned short argc){
    VEC field;
    int nmax;
    P7OUT|=BIT0|BIT1;
    //extrapolate model to desired date
    nmax=extrapsh(2013);
    P7OUT&=~BIT1;
    P7OUT|=BIT2;
    //calculate magnetic field
    shval3(64.9261111/RAD2DEG,-147.4958333/RAD2DEG,6371.2,nmax,&field);
    //shval3(64.9261111/RAD2DEG,-147.4958333/RAD2DEG,6771.2,nmax,&field);
    P7OUT&=~(BIT0|BIT2);
    //print
    printf("x = %f\r\ny = %f\r\nz = %f\r\n",field.c.x,field.c.y,field.c.z);
    return 0;
}

int resetCmd(char **argv,unsigned short argc){
  if(argc!=0){
    printf("Error : %s takes no arguments.\r\n");
  }
  //reset by abusing WDT, this causes PUC
  WDTCTL=0;
}

void initCLK(void){
  //set XT1 load caps, do this first so XT1 starts up sooner
  BCSCTL3=XCAP_0;
  //stop watchdog
  WDTCTL = WDTPW|WDTHOLD;
  //setup clocks

  //set DCO to 16MHz from calibrated values
  DCOCTL=0;
  BCSCTL1=CALBC1_16MHZ;
  DCOCTL=CALDCO_16MHZ;
}


void main(void){
  //character from port
  int c=0;
  //buffer for command
  char cmd[30];
  //buffer for last command
  char last[30];
  //command string index
  unsigned int cIdx=0;
  initCLK();
  //initialize UART for USB
  initUART();

  //Turn off LED's
  P7OUT=0;
  //set LED's as outputs
  P7DIR=0xFF;
  
  //enable interrupts
  __enable_interrupt();

  printf("\rIGRF test program ready\r\n>");

    //initialize buffers
  cmd[0]=0;
  last[0]=0;
  for(;;){
    //get character
    c=Getc();
    //wait for UART events
    while(c==EOF){
      LPM0;
      c=Getc();
    }
    //process received character
    switch(c){
      case EOF:
      continue;
      case '\r':
      case '\n':
        //return key run command
        if(cIdx==0){
          //if nothing entered, do last command
          printf(last);       //print command
          strcpy(cmd,last);   //copy into command buffer
        }else{
          //run command from buffer
          cmd[cIdx]=0;    //terminate command string
          cIdx=0;         //reset the command index
          //save this as last command
          strcpy(last,cmd);
        }
        //send carriage return and new line
        printf("\r\n");
        //run command
        doCmd(cmd);
        //print prompt char
        putchar('>');
      continue;
      case '\x7f':
      case '\b':
        //backspace
        if(cIdx==0)continue;
        //backup and write over char
        printf("\b \b");
        //decrement command index
        cIdx--;
      continue;
      case '\t':
      //ignore tab character
      continue;
    }
    if(!iscntrl(c)){
      //echo character
      putchar(c);
      //put character in command buffer
      cmd[cIdx++]=c;
    }
  }
}
