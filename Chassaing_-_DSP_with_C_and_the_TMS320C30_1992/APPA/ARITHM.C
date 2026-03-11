/*ARITHM.C - USE OF VARIABLES AND ARITHMETIC OPERATIONS*/
#include <stdio.h>

main()
{
  float x = 5.0;   /*declare and initialize x to 5.0*/
  float y = 2.0;   /*declare and initialize y to 2.0*/
  float z;         /*declare z*/
  int A,B;
  A = 5;
  B = 2;
  z = x/y;
  printf ("x/y = %f\n", z);                  /*output z*/
  printf ("A/B = %i\n", A/B);                /*output int part of A / B*/
  printf ("A mod B = %i\n", A%B);            /*output modulus of A / B */
  printf ("A/B = %f\n", (float)A/(float)B);  /*floating output of A / B*/
}
