/*LOOPC.C-LOOP PROGRAM USING PORT 0 WITH POLLING */
#include "aiccom.c" /*AIC Communication routines */              
int AICSEC[4] = {0x1428,0x1,0x4A96,0x67}; /*Config data for AIC*/
main()
{
  int  data_in, data_out;             /*Initialize variables  */
  AICSET();                           /*Function to config AIC*/ 
  while (1)                           /*Create endless loop   */
  {
   data_in = UPDATE_SAMPLE(data_out); /*Call function to update sample*/
   data_out = data_in;                /*Loop input to output          */
  }
}  
  
    
      
