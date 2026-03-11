/*SINEC.C-SINE GENERATION BY RECURSIVE EQUATION   */
#include <math.h>          /*math library function*/
#define SAMPLE_FREQ 10000  /*sample frequency     */
#define SINE_FREQ 1000     /*desired frequency    */  
#define PI 3.14159         /*constant PI          */  
volatile int *IO_OUTPUT=(volatile int *)0x804002; 

void sinewave(float A,float B,float C)
{                        
 float Y[3] = {0.0,0.0,0.0};  /*Y[N] array         */
 float X[3] = {0.0,0.0,1.0};  /*X[N] array         */
 int N = 2, i;                /*declare variables  */
 for (i = 0; i < 100; i++)
  {                              
   Y[N] = A*Y[N-1] + B*Y[N-2] + C*X[N-1]; /*determine Y[N]*/ 
   *IO_OUTPUT = Y[N]*1000;   /*output Y[N] scaled by 1000 */
   Y[N-2] = Y[N-1];          /*shift Y's back in array    */
   Y[N-1] = Y[N];                                
   X[N-2] = X[N-1];          /*shift X's back in array    */
   X[N-1] = X[N];
   X[N] = 0.0;               /*set future X's to 0        */
  } 
}                               

main()
{
 float Fs, Fosc, w, T, A, B, C;  /*declare variables       */
 Fs = SAMPLE_FREQ;               /*get sampling frequency  */
 Fosc = SINE_FREQ;               /*get oscillator frequency*/
 T = 1/Fs;                       /*determine sample period */
 w = 2*PI*Fosc;                  /*determine angular freq  */
 A = 2 * cos((w * T));           /*determine coefficient A */
 B = -1.0;                       /*coeff B is constant     */
 C = sin((w * T));               /*determine coefficient B */
   sinewave(A,B,C);              /*call sinewave function  */
}                                                      









