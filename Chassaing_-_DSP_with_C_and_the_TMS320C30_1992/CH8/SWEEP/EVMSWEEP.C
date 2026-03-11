/*EVMSWEEP.C-GENERATES SINUSOID WITH FREQUENCY SWEEP */                      
#include <math.h>
#include "aiccom.c"
#define PI 3.14159                          
#define SAMPLE_FREQ 10000
volatile int *host = (volatile int *) 0x804000;
volatile int *dma = (volatile int *) 0x808000;
void change_freq(float);                                            
int average_peak(int *, int);
int AICSEC[4] = {0x1428,0x01,0x4A96,0x67};
float A,B,C;

void sinewave()
    {                        
    float Y[3]={0.0,0.0,0.0};           /* Y[N] array                 */
    float X[3]={0.0,0.0,1.0};           /* X[N] array                 */
    float current_freq;
                         
    int inval[200];    
    int peak;
    int counter=0;
    int N=2;                            /* used for array indexing    */

    current_freq=100.0;
    change_freq(current_freq);
    
    while(1)
        {                         
        TWAIT;

        /* calculate and output sine wave value */

        Y[N]=A*Y[N-1] + B*Y[N-2] + C*X[N-1];
        PBASE[0x48]=-4 & (int)(Y[N]*2600);

        Y[N-2]=Y[N-1];                  /* shift Y's back in array    */
        Y[N-1]=Y[N];                                 
        X[N-2]=X[N-1];                  /* shift X's back in array    */
        X[N-1]=X[N];
        X[N]=0.0;                       /* set future X's to 0        */
        inval[counter++]=(PBASE[0x4C]<<16)>>18;
        if (counter==200)                    
            {        
            peak=average_peak(inval, counter);
            if (peak==0)
                peak=1;
            while(*host & 0x0000FFFF);
            *host=peak;
            while(*host!=(peak & 0x0000FFFF));
            *host=0x0;
            peak=0;
            counter=0;
            Y[0]=Y[1]=Y[2]=0.0;
            X[0]=X[1]=0.0;
            X[2]=1.0;
            if (current_freq<2000.0)
                current_freq+=100.0;
            else 
                current_freq=100.0;
            change_freq(current_freq);
            }
        } /* end of loop */
    }                               

int average_peak(int *values, int count)
    {                                     
    long int sum;
    int avg_peak;
    int i;
            
    sum=0;
    for (i=0;i<count;i++)
        {
        if (values[i]<0)
            values[i]=0-values[i];
        sum+=values[i];
        }
    avg_peak=sum/count;
    return(avg_peak);
    }

void change_freq(float newfreq)
    {                         
    float Fs;
    float T;
    float w;    

    Fs=SAMPLE_FREQ;        
    T=1/Fs;                         /* determine sample period    */
    w=2*PI*newfreq;                    /* determine angular freq.    */
    A=2 * cos((w * T));            /* determine coefficients     */
    B=-1.0;
    C=sin((w * T));
    return;
    }                       
                                  

main()
    {
    AICSET();
    INIT_EVM();
    *host=0x0;
    sinewave();
    }              

















