/*SINE4AIC.C-SINE PROGRAM WITH 4 POINTS USING INTERRUPTS                */
#include "aiccom.c"                       /*AIC comm routines           */     
#define VEC_ADDR (volatile int *) 0x00;   /*addr of vectors             */
int AICSEC[4] = {0x1428,0x1,0x4A96,0x67}; /*config data for AIC SP0     */
int data_out, loop = 0;                   /*declare global variables    */
int sin_table[4] = {0,1000,0,-1000};      /*values for 4-point sinewave */

void c_int05()                            /*TINT0 interrupt routine     */
{
  PBASE[0x48] = sin_table[loop] << 2;     /*output value from sine table*/
  if (loop < 3) ++loop;                   /*increment loop counter < 3  */
  else loop = 0;                          /*reset loop counter          */
}                                                                       

main()
{
  volatile int *INTVEC = VEC_ADDR;        /*pointer to vectors          */
  INTVEC[5] = (volatile int) c_int05;     /*install interrupt 5 Handler */
  AICSET_I();                             /*function to configure AIC   */ 
  for (;;);                               /*wait for interrupt          */
}  
      
      
