/*SINECR.C-REAL-TIME SINE GENERATION BY RECURSIVE EQUATION */
#include "aiccom.c"                /*AIC comm routines     */
#include <math.h>                  /*math library function */
#define SAMPLE_FREQ 10000          /*sample frequency      */
#define SINE_FREQ 1000             /*desired frequency     */  
#define PI 3.14159                 /*constant PI           */  
int AICSEC[4] = {0x1428,0x1,0x4A96,0x67}; /*AIC config data*/

void sinewave(float A,float B,float C)
{                        
 float Y[3] = {0.0,0.0,0.0};  /*Y[N] array         */
 float X[3] = {0.0,0.0,1.0};  /*X[N] array         */
 int N = 2, result;           /*declare variables  */
 while(1)
  {                              
   TWAIT;
   Y[N] = A*Y[N-1] + B*Y[N-2] + C*X[N-1];  /*determine Y[N]*/ 
   result = (int)(Y[N]*1000); /*output Y[N] scaled by 1000 */
   PBASE[0x48] = result << 2; /*output to AIC              */
   Y[N-2] = Y[N-1];           /*shift Y's back in array    */
   Y[N-1] = Y[N];                                 
   X[N-2] = X[N-1];           /*shift X's back in array    */
   X[N-1] = X[N];
   X[N] = 0.0;                /*set future X's to 0        */
  } 
}                               

main()
{
 float Fs, Fosc, w, T, A, B, C;  /*declare variables       */
 AICSET();                       /*initialize AIC          */
 Fs = SAMPLE_FREQ;               /*get sampling frequency  */
 Fosc = SINE_FREQ;               /*get oscillator frequency*/
 T = 1/Fs;                       /*determine sample period */
 w = 2*PI*Fosc;                  /*determine angular freq  */
 A = 2 * cos((w * T));           /*determine coefficient A */
 B = -1.0;                       /*coeff B is constant     */
 C = sin((w * T));               /*determine coefficient B */
   sinewave(A,B,C);              /*call sinewave function  */
}                                                      









