
#include <msp430.h>
#include "uart.h"

//buffers for USB UART
struct Tx TxBuf;
struct Rx RxBuf;

//UART TX ISR called to transmit UART data
void USB_TX(void) __interrupt[USCIAB1TX_VECTOR]{
  unsigned char flags=UC1IFG&(UC1IE);
  //process UART TXIFG
  if(flags&UCA1TXIFG){
    unsigned short t=TxBuf.out;
    //check if their are more chars
    if(TxBuf.in!=t){
      //more chars TX next
      UCA1TXBUF=TxBuf.buf[t++];
      TxBuf.out=(t)%TX_SIZE;
    }else{
      //buffer empty disable TX
      TxBuf.done=1;
      UC1IFG&=~UCA1TXIFG;
    }
    //more room in buffer, exit LPM
    LPM4_EXIT;
  }
}

// receive UART ISR
void USB_rx(void) __interrupt[USCIAB1RX_VECTOR]{
  unsigned char flags=UC1IFG&(UC1IE);
  //process UART RXIFG
  if(flags&UCA1RXIFG){
    unsigned int t;
    t=(RxBuf.in+1)%RX_SIZE;
    //check if there is room
    if(t!=RxBuf.out){
      //write char in buffer
      RxBuf.buf[RxBuf.in]=UCA1RXBUF;
      //advance index
      RxBuf.in=t;
      //new char ready, exit LPM
      LPM4_EXIT;
    }
    //if no room char is lost
  }
}
