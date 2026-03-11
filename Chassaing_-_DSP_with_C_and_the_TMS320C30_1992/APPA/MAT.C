/*MAT.C - CONCEPTS OF ARRAYS AND LOOPS*/

main()
{
  float A[3][3] = {{1.0,2.0,3.0},
                   {4.0,5.0,6.0},
                   {7.0,8.0,9.0}};
  float B[3] =     {1.0,2.0,3.0};
  int i, j;
  float result;
  for (i = 0; i < 3; i++)           /*first for loop*/
  {
    result = 0;
    for (j = 0; j < 3; j++)         /*second for loop*/
      result += A[i][j]*B[j];
    printf("%f\n", result);
  }
} 
