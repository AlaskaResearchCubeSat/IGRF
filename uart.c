
#include <msp430.h>
#include <string.h>
#include <stdio.h>
#include "uart.h"

void initUART(void){
 //init structures
  memset(&TxBuf,0,sizeof(TxBuf));
  memset(&RxBuf,0,sizeof(RxBuf));
  TxBuf.done=0;
   //setup UCA1 for USB-UART operation
  UCA1CTL1=UCSWRST;
  UCA1CTL0=0;
 UCA1CTL1|=UCSSEL_1;
 //UCA1CTL1|=UCSSEL_2;
 //set baud rate to 9600
  UCA1BR0=3;
  UCA1BR1=0;
  UCA1MCTL=UCBRS_3;
  //set baud rate to 38400
  //UCA1BR0=26;
  //UCA1BR1=0;
  //UCA1MCTL=UCBRF_1|UCOS16;
  //set baud rate to 57600
  //UCA1BR0=17;
  //UCA1BR1=0;
  //UCA1MCTL=UCBRF_6|UCBRS_0|UCOS16;
  //setup pins
  P3SEL|=BIT6|BIT7;
  //take UCA1 out of reset mode
  UCA1CTL1&=~UCSWRST;
  //enable interrupts
  UC1IE|=UCA1TXIE|UCA1RXIE;
}

//queue byte to get transmitted
int TxChar(unsigned char c){
	unsigned int t;
	int res=c;
        unsigned short state;
	//disable interrupt
	state=__disable_interrupt();
	//check if transmitting
	if(TxBuf.done){
		//bypass queue for first byte if not transmitting
		UCA1TXBUF=c;
		//clear done flag
		TxBuf.done=0;
	//queue byte
	}else{
		//get next input index
		t=(TxBuf.in+1)%TX_SIZE;
		//check if next input index equals output index
		if(t!=TxBuf.out){
			//not equal, room in buffer
			//put char in buffer
			TxBuf.buf[TxBuf.in]=c;
			//set next input index
			TxBuf.in=t;
		}else{
			//buffer full
			res=EOF;
		}
	}
	//enable interrupt
	__set_interrupt(state);
	//return result
	return res;
}

//get byte from buffer return if buffer is empty
int Getc(void){
	int c;
	unsigned int t;
        unsigned short state;
	//check for bytes
	if(RxBuf.in==RxBuf.out){
		//no bytes return EOF
		return EOF;
	}else{
		//disable interrupt while byte is retrieved
		state=__disable_interrupt();
		//get output index
		t=RxBuf.out;
		//get char from buffer
		c=(unsigned char)(RxBuf.buf[t++]);
		//set new output index
		RxBuf.out=(t)%RX_SIZE;
		//re enable interrupt
		__set_interrupt(state);
	}
	//return result
	return c;
}

