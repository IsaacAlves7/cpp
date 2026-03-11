/*FFTGEN.C-GENERATES TWIDDLE CONSTANTS */
#include <math.h>
#include <stdio.h>
#define N 256  /*to generate 256 complex points*/

main()
  {
  FILE *fptr;
  double sinval[N];
  double cosval[N];
  double arg;
  int i;
  fptr=fopen("twid256.h","w");
  arg=2*3.141592654/512;
  for(i=0;i<N;i++)
    {
    cosval[i]=(float)cos((i*arg));
    sinval[i]=-(float)sin((i*arg));
    }
  fprintf(fptr,"struct\n");
  fprintf(fptr,"    {\n");
  fprintf(fptr,"    double real;\n");
  fprintf(fptr,"    double imag;\n");
  fprintf(fptr,"    }");
  fprintf(fptr," w[]={");
  for(i=0;i<N;i++)
    {
    fprintf(fptr,"%8.5f,%8.5f,\n",cosval[i],sinval[i]);
    fprintf(fptr,"      ");
    }
  fclose(fptr);
  }


