/*LOOPC1.C-LOOP PROGRAM TO TEST AIC PORT 1                          */
#include "aiccom01.c"                      /*AIC COMM ROUTINES      */ 
int AICSEC[4];                             /*PORT 0 NOT USED FOR AIC*/
int AICSEC1[4] = {0x1428,0x1,0x4A96,0x67}; /*CONFIG DATA FOR AIC    */
main()
{
  int data_in, data_out;       /*INITIALIZE VARABLES      */
  AICSET(SP1);                 /*FUNCTION TO CONFIGURE AIC*/ 
  while (1)                    /*CREATE ENDLESS LOOP      */
  {
    data_in = UPDATE_SAMPLE(SP1, data_out); /*CALL FUNCT TO UPDATE SAMPLE*/
    data_out = data_in;                     /*LOOP INPUT TO OUTPUT       */
  }
}  
  
    
      
