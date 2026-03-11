/*SINEGEN.C-GENERATES SINE VALUES FOR REAL FFT*/
#include <math.h>
#include <stdio.h>
#define N 512
#define pi 3.141592654

main()
{
  FILE *stream;
  int n;
  float result;
  stream = fopen("twid512.asm", "w+");
  fprintf(stream, "\n%s", "          .global    _sine");
  fprintf(stream, "\n%s", "          .data");
  fprintf(stream, "\n%s%7f", "_sine     .float     ", 0.0000000);
  for (n = 1; n < N/2; n++)
  {
    result = sin(n*2*pi/N);
    fprintf(stream, "\n%s%7f", "          .float     ", result);
  }
  fclose(stream);
}
