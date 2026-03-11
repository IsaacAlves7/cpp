/*MATRIXC.C-MATRIX MULTIPLICATION (3x3)(3x1)=(3x1)*/
main()
{
 volatile int *IO_OUTPUT = (volatile int *) 0x804001;
 float A[3][3] = {{1,2,3},   
                  {4,5,6},  
                  {7,8,9}};  
 float  B[3] = {1,2,3};
 float result;                   
 int i, j;
 for (i = 0; i < 3; i++)
  {
   result = 0; 
   for (j = 0; j < 3; j++)
    {
     result += A[i][j] * B[j];   
    }
   *IO_OUTPUT = (int)result;  /*result=14,32,50*/
  }
}

