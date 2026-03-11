/*SHIFTTB.C-DEMONSTRATES -90 DEGREES SHIFT IN PHASE*/ 
#define Max 8192
#define Points 256
#define Pi  3.1415926 
#include "scdat"        /*look-up table for acos, asin*/

main(void) 
  { 
  int count, Y1; 
  double yo, Y, IN_data; 
  volatile int *IO_INPUT = (volatile int*) 0x804000; 
  volatile int *IO_OUTPUT = (volatile int*) 0x804001; 
  for (count=0; count < Points; count++) 
    { 
    yo = IN_data;                /*store signal value         */ 
    IN_data = *IO_INPUT;         /*input signal in            */ 
    Y1 = ((IN_data/Max)+1)*100;  /*step up array ( 0-> 200 )  */ 
    if (yo >= IN_data)           /*is signal falling or rising*/ 
       *IO_OUTPUT = yc[Y1]*Max;  /*signal falling             */ 
    else 
       *IO_OUTPUT = ys[Y1]*Max;  /*signal rising              */ 
    } 
  }
