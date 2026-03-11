/*ADAPTC.C-ADAPTATION USING LMS IN C WITHOUT THE TI SIMULATOR*/
#include <stdio.h>
#include <math.h>
#define beta 0.05                      //step factor for coeff update
#define N  20                          //order of filter (N+1)
#define NS  40                         //number of samples
#define FS  8000                       //sampling frequency
#define Pi  3.1415926                     
#define DESIRED 2*cos(2*Pi*T*1000/FS)  //desired signal
#define NOISE sin(2*Pi*T*1000/FS)      //noise signal

main()
{
  long I, T;
  double D, Y, E; 
  double W[N+1] = {0.0};
  double X[N+1] = {0.0};
  FILE *desired, *Y_out, *error, *noise; 
  desired = fopen ("DESIRED", "w++"); //open file for desired samples
  Y_out = fopen ("Y_OUT", "w++");     //open file for Y output samples
  error = fopen ("ERROR", "w++");     //open file for error samples
  noise = fopen ("NOISE", "w++");     //open file for noise samples
  for (T = 0; T < NS; T++)           //start adaptive alogrithm
  {    
    X[0] = NOISE;                     //new noise sample
    D = DESIRED;                      //desired signal
    Y = 0;                            //set output of filter to zero
    for (I = 0; I <= N; I++)
      Y += (W[I] * X[I]);             //calculate filter output
    E = D - Y;                        //calculate error signal
    for (I = N; I >= 0; I--)          //convolve the coeff with noise samples
    {
      W[I] = W[I] + (2*beta*E*X[I]);  //update filter coefficients
      if (I != 0)                     
        X[I] = X[I-1];                //move data sample 
    }   
    fprintf (desired, "\n%10g    %10f", (float) T/FS, D); 
    fprintf (Y_out, "\n%10g    %10f", (float) T/FS, Y);
    fprintf (error, "\n%10g    %10f", (float) T/FS, E);
    fprintf (noise, "\n%10g    %10f", (float) T/FS, X[0]);
  }   
  fclose (desired);
  fclose (Y_out);
  fclose (error);
  fclose (noise);
}


