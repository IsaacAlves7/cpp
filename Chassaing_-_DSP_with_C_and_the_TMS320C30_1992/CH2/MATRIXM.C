/*MATRIXM.C MATRIX MULTIPLICATION.CALLS ASSEMBLY FUNCTION*/
volatile int *IO_OUTPUT=(volatile int *) 0x804001; 
extern void matrixmf (float *, float *, int *);
main()
{
 float A[3][3] = {{1,2,3},
                  {4,5,6},
                  {7,8,9}};
 float B[3] =   {1,2,3};
  matrixmf ((float *) A, (float *)B, (int *) IO_OUTPUT);
}
