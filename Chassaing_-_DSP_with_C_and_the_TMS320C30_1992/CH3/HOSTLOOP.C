/*HOSTLOOP.C-LOOP PROGRAM WITH AMPLITUDE CONTROL*/
#include "aiccom.c"                           /* AIC communications routines*/
int AICSEC[4] = {0x1428,0x1,0x4A96,0x67};       /* AIC setup data           */
volatile int *host = (volatile int *) 0x804000; /* address of host          */

void main(void)
{
  int data_IN, data_OUT, ampt = 1, i;            /* declare variables       */
  AICSET();                                      /* initialize AIC          */
  do
  {
    i = *host & 0xFF;                            /* i = 8 bit word from host*/
    if (i > 0 && i < 11) ampt = i;               /* i is attenuation value  */
    data_IN = UPDATE_SAMPLE(data_OUT);           /* input output sample     */
    data_OUT = data_IN / ampt;                   /* scale input to output   */
  }
  while (1);                                     /* endless loop            */
}  
  

