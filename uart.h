#ifndef __UART_H
  #define __UART_H

#define RX_SIZE   (8)
#define TX_SIZE   (1024)

  //TX buffer structure
  struct Tx{
          unsigned int in;
          volatile unsigned int out;
          char done;
          char buf[TX_SIZE];
  };
  
  //RX buffer structure
  struct Rx{
          volatile unsigned int in;
          unsigned int out;
          char buf[RX_SIZE];
  };

  #define rxFlush()  (BT_RxBuf.in=BT_RxBuf.out)

  int TxChar(unsigned char c);
  int Getc(void);
  void initUART(void);

  extern struct Tx TxBuf;
  extern struct Rx RxBuf;

#endif
  