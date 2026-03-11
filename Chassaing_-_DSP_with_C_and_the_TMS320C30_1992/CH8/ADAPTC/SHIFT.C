/*SHIFT.C-DEMONSTRATES 90 DEGREES DELAY IN PHASE */
#include <math.h>
#define Pi 3.1415926

main()
{
  int count;
  double Y,X,yo;
  volatile int *IO_INPUT = (volatile int*) 0x804000;
  volatile int *IO_OUTPUT = (volatile int*) 0x804001;
  for (count=0; count < 256; count++)
    {
    Y= *IO_INPUT;                
    Y= Y/8192;              /*Y must be between 1 and -1 */ 
    if (yo>=Y)              /*is signal falling or rising*/
      X=acos(Y);            /*signal is falling          */
    else
      X=asin(Y)-(Pi/2);     /*signal is rising           */
    X=X-(Pi/2);             /*shift by 90 degrees        */
    *IO_OUTPUT=8192*cos(X); /*shifted output Y value     */
    yo=Y;                   /*store Y value              */
    }
}  

