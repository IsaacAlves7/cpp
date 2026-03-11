/*LOOPCI.C-LOOP PROGRAM TO TEST AIC USING INTERRUPTS                 */
#include "aiccom01.c"                      /*AIC comm routines       */        
#define VEC_ADDR (volatile int *) 0x00     /*addr of vectors         */
int AICSEC[4] = {0x1428,0x1,0x4A96,0x67};  /*config data for SP0     */
int AICSEC1[4];                            /*SP1 not used            */
int data_in, data_out;                     /*declare global variables*/

void c_int05()                             /*TINT0 interrupt routine */
{
  data_in = UPDATE_SAMPLE(SP0, data_out);  /*update sample to SP0 AIC*/
  data_out = data_in;                      /*loop input to output    */
}      

main()
{
  volatile int *INTVEC = VEC_ADDR;      /*pointer to vectors         */
  INTVEC[5] = (volatile int) c_int05;   /*install interrupt 5 Handler*/
  AICSET_I(SP0);                        /*configure SP0 of AIC       */ 
  for (;;);                             /*wait for interrupt         */
}  
  
   
